#!/usr/bin/env nextflow
/*
    Barrett Lab Pool-Seq Alignment/Analysis Pipeline
    Authors:
    - Wing-Zheng Ho <wingzhg@gmail.com>
*/
nextflow.enable.dsl=2

/*
    Params
*/
params.publish_dir = './results'

// Define the genome
params.genome = null
genome = file(params.genome)

// Read sample sheet
params.sample_sheet = null


sample_sheet = Channel.fromPath(params.sample_sheet, checkIfExists: true)
                .ifEmpty { exit 1, "sample sheet not found" }
                .splitCsv(header:true, sep: ",")
                .map { row -> row.fq1 = row.seq_folder ? row.fq1 = row.seq_folder + "/" + row.fq1 : row.fq1; row }
                .map { row -> row.fq2 = row.seq_folder ? row.fq2 = row.seq_folder + "/" + row.fq2 : row.fq2; row }
                .map { row -> row.id = row.pool ? row.id = row.pool + "_" + row.id : row.id; row }
                .map { row -> [row, file(row.fq1), file(row.fq2)] }
                .set { sample_info_ch }



// Includes

// Procceses


/* 
    Genome Prep
*/

process genome_prep_bwa {
    publishDir "${params.publish_dir}/01_genome/bwa", mode: 'copy'
    input:
        path(genome)
    output:
        path("${genome}.*"), emit: bwa_idx

    script:
    """
    bwa index -a bwtsw ${genome}
    samtools faidx ${genome}

    """
}


/* 
    FastQC
*/

process fastqc {
    tag {row.id}
    publishDir "${params.publish_dir}/02_QC/FastQC", mode:'copy'
    input:
    tuple val(row), path(fq1), path(fq2)

    output:
    path("*.html")
    path("*.zip"), emit: fastqc_f

    script:
    """
    fastqc -t 8 "${fq1}" "${fq2}"
    """
}


/* 
    MultiQC
*/


process MULTIQC {
    publishDir "${params.publish_dir}/02_QC/MultiQC", mode:'copy'

    input:
    file("*")

    output:
    path('multiqc_report.html')

    script:
    """
    multiqc .
    """
}


/* 
    fastp!
*/
process fastp {

    tag {row.id}
    publishDir "${params.publish_dir}/03_trim", mode:'copy'

    input:
        tuple val(row), path(fq1), path(fq2)

    output:
        tuple val(row), file("trimmed_*_R1.fastq.gz"), file("trimmed_*_R2.fastq.gz"), emit: trimmed_reads

    script:
    """
    fastp -i ${fq1} ${fq2} -o trimmed_${fq1} -O trimmed_${fq2} --unpaired1 unpaired${fq1} --unpaired2 unpaired${fq2} --failed_out ${row.id}_failed.fastq.gz --detect_adapter_for_pe -q 20 -l 50 -g
    """
}



/* 
    Alignment
*/

process alignment {

    tag {row.id}
    publishDir "${params.publish_dir}/04_align", mode: 'copy'
    input:
        tuple val(row), file(reads1), file(reads2)
        path(bwa_idx)
        
    output:
        tuple val(row), file("${row.id}.bam"), file("${row.id}.bam.bai")

	script:
		// Construct read group
        def idxbase = bwa_idx[0].baseName
		RG = ["@RG",
			  "ID:${row.id}",
			  "SM:${row.pool}",
			  "LB:${row.lb}",
			  "PL:${row.pl}"].join("\\t")

    """
        bwa mem -t 8 -R '${RG}' ${idxbase} ${trimmed_R1} ${trimmed_R2} | \\
        sambamba view --nthreads=8 --show-progress --sam-input --format=bam --with-header /dev/stdin | \\
        sambamba sort --nthreads=8 --show-progress --tmpdir=. --out=${row.id}.bam /dev/stdin
        sambamba index --nthreads=8 ${row.id}.bam
        if [[ ! \$(samtools view ${row.id}.bam | head -n 10) ]]; then
            exit 1;
        fi
    """
}


/* 
    Merge ID Bams by Pool
*/
process merge_bam {

    tag { row.pool }
    publishDir "${params.publish_dir}/05_merge", mode: 'copy'
    input:
        tuple val(pool), val(row), path(bam), path(bai), val(n_count)

    output:
        tuple val(pool), val(row), file("${row.pool}.bam"), file("${row.pool}.bam.bai")

    script:
        if (n_count == 1)
            """
                mv ${bam} ${pool}.bam
                mv ${bai} ${pool}.bam.bai
            """
        else
            """
                sambamba merge --nthreads=8 --show-progress ${pool}.bam ${bam}
                sambamba index --nthreads=8 ${pool}.bam
            """
}

/* 
    Remove PCR Duplicate
*/

process mark_dups {

    tag { "${pool}" }

    publishDir "${params.publish_dir}/06_dedup", mode: 'copy', pattern: '*.bam*'

    input:
        tuple val(pool), val(row), path("${pool}.in.bam"), path("${pool}.in.bam.bai")
    output:
        tuple val(row), file("*.bam"), file("*.bam.bai"), emit: "bams"
        tuple val(row), path("${pool}.bam"), path("${pool}.bam.bai"), emit: "pool_sheet"
        path "${pool}.duplicates.txt", emit: "markdups"
        tuple path("${pool}.bam"), path("${pool}.bam.bai"), emit: "npr"

    """
        java -Xmx20g -Xms1g -jar $EBROOTPICARD/picard.jar  MarkDuplicates -I ${pool}.in.bam -O ${pool}.bam -M ${pool}.duplicates.txt --VALIDATION_STRINGENCY STRICT --REMOVE_DUPLICATES false --TAGGING_POLICY All --REMOVE_SEQUENCING_DUPLICATES TRUE --SORTING_COLLECTION_SIZE_RATIO 0.1

        sambamba index --nthreads=8 ${pool}.bam
    """
}

process samtools_mpileup {
    
    //tag { pool }
    publishDir "${params.publish_dir}/07_mpileup", mode: 'copy'

    input:
    //tuple val(row), path("*.bam"), path("${row.pool}.bam.bai")
    file(npr)
    path(genome)
    path(bwa_idx)

    output:
    path("*.mpileup"), emit: "mpileup"
    file ("*.txt")


    script:
    def idxbase = bwa_idx[0].baseName
    //def all_bam = npr
    """
    find . -name "*.bam" > all_bams.txt
    samtools mpileup -d 1000 -Q 20 -B -f ${idxbase} -b all_bams.txt -o all_bams.mpileup
    """
}

process mpileup_to_sync {
    publishDir "${params.publish_dir}/08_sync", mode: 'copy'
    input:
    file(mpileup)
    path(bwa_idx)

    output:
    path("*.sync"), emit: "sync"
    script:
    """
    grenedalf sync --pileup-path *.mpileup --pileup-min-base-qual 20 --threads 8 --reference-genome-fai-file *.fai
    """
}

process fst {
    publishDir "${params.publish_dir}/09_fst", mode: 'copy'
    input:
    file(sync)

    output:
    path("*.sync"), emit: "sync"
    script:
    """
    grenedalf sync --pileup-path *.mpileup --pileup-min-base-qual 20 --threads 8 --reference-genome-fai-file *.fai
    """
}



workflow {
        genome_prep_bwa(params.genome)
        fastqc(sample_info_ch)
        MULTIQC(fastqc.out.fastqc_f.collect())
        fastp(sample_info_ch)
        alignment(fastp.out.trimmed_reads, genome_prep_bwa.out.bwa_idx)

        merge_in = alignment.out.map { row -> [row[0].pool] + row }
                    .groupTuple()
                    .map { pool, row, bam, bai -> [pool, row.sort()[0], bam, bai, row.size()] }
        merge_in | merge_bam | mark_dups

        samtools_mpileup(mark_dups.out.npr.collect(), params.genome, genome_prep_bwa.out.bwa_idx)
        mpileup_to_sync(samtools_mpileup.out.mpileup, genome_prep_bwa.out.bwa_idx)
        fst(mpileup_to_sync.out.sync)
}

