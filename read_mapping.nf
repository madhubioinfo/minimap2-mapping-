#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads     = "${params.reads ?: './reads'}"
params.reference = "${params.reference ?: './reference.fa'}"

workflow {

    reads_ch = Channel.fromPath("${params.reads}/*.fastq")
                      .map { file -> tuple(file.baseName, file) }

    MAP_READS(reads_ch, params.reference)
    SORT_BAM(MAP_READS.out)
    COVERAGE_CSV(SORT_BAM.out)
}

process MAP_READS {
    tag "$sample_id"
    input:
        tuple val(sample_id), path(fastq)
        path reference
    output:
        tuple val(sample_id), path("${sample_id}.bam")
    publishDir "results", mode: 'copy'

    script:
    """
    minimap2 -ax map-ont $reference $fastq | \
        samtools view -bS - > ${sample_id}.bam
    """
}

process SORT_BAM {
    tag "$sample_id"
    input:
        tuple val(sample_id), path(bam)
    output:
        tuple val(sample_id), path("${sample_id}_sorted.bam")
    publishDir "results", mode: 'copy'

    script:
    """
    samtools sort $bam -o ${sample_id}_sorted.bam
    """
}

process COVERAGE_CSV {
    tag "$sample_id"
    input:
        tuple val(sample_id), path(sorted_bam)
    output:
        path("${sample_id}_coverage.csv")
    publishDir "results", mode: 'copy'

    script:
    """
    samtools coverage $sorted_bam | \
        cut -f 1,4 | awk '\$2 > 0' | sort -rnk 2,2 | sed 's/\\t/,/g' > ${sample_id}_coverage.csv
    echo "gene_name,number_of_reads_mapped" | cat - ${sample_id}_coverage.csv > tmp && mv tmp ${sample_id}_coverage.csv
    """
}
