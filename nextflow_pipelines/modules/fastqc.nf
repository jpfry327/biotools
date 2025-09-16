#!/usr/bin/env nextflow

process FASTQC {
    tag "$sample"
    label 'process_single'

    container "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/71/7180d4f5ff01f33209c7b39b31277d416bffe8d221149a25894805afc9755835/data"


    input:
    tuple val(sample), path(read1), path(read2)

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc ${read1} ${read2} \\
        --outdir . \\
        --threads ${task.cpus} \\
        --quiet
    """

    stub:
    """
    touch ${read1.baseName}_fastqc.zip
    touch ${read1.baseName}_fastqc.html
    touch ${read2.baseName}_fastqc.zip
    touch ${read2.baseName}_fastqc.html
    """
}