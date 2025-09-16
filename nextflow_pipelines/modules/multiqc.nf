#!/usr/bin/env nextflow

process MULTIQC {

    container "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/71/7180d4f5ff01f33209c7b39b31277d416bffe8d221149a25894805afc9755835/data"
    executor 'local'

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path '*'
    val(output_name)

    output:
    path "${output_name}.html", emit: report
    path "${output_name}_data", emit: data

    script:
    """
    multiqc . -n ${output_name}.html
    """
}