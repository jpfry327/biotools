/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EXAMPLE_MODULE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Simple example module that demonstrates basic module structure
*/

process EXAMPLE_MODULE {
    tag "$meta.id"
    label 'process_single'
    
    publishDir "${params.outdir}/example_results", mode: 'copy'
    
    input:
    tuple val(meta), path(files)
    
    output:
    tuple val(meta), path("${meta.id}_output.txt"), emit: results
    path("versions.yml"),                           emit: versions
    
    script:
    """
    echo "Hello, this is the main workflow processing sample: ${meta.id}" > ${meta.id}_output.txt
    echo "Input files: ${files}" >> ${meta.id}_output.txt
    echo "Processed at: \$(date)" >> ${meta.id}_output.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        example_module: 1.0.0
    END_VERSIONS
    """
    
    stub:
    """
    touch ${meta.id}_output.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        example_module: 1.0.0
    END_VERSIONS
    """
}