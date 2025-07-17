
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def helpMessage() {
    log.info"""
    ========================================
     ${workflow.manifest.name} v${workflow.manifest.version}
    ========================================
    
    Usage:
      nextflow run main.nf --input samplesheet.csv --outdir results
    
    Required arguments:
      --input               Path to input samplesheet (CSV format)
      --outdir              Output directory
    
    Optional arguments:
      --help                Show this help message and exit
    
    """.stripIndent()
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow PIPELINE_INITIALISATION {

    take:
    help         // boolean: Display help and exit
    input        // string: Path to input file
    outdir       // string: Output directory

    main:

    // Display help message and exit if requested
    if (help) {
        helpMessage()
        exit 0
    }

    // Validate required parameters
    if (!input) {
        log.error "ERROR: --input parameter is required"
        helpMessage()
        exit 1
    }

    // Check input file exists
    if (!file(input).exists()) {
        log.error "ERROR: Input file does not exist: ${input}"
        exit 1
    }

    // Create output directory if it doesn't exist
    file(outdir).mkdirs()

    // Log pipeline start
    log.info """
    ========================================
     Pipeline Started Successfully
    ========================================
     Input: ${input}
     Output: ${outdir}
     Nextflow: ${workflow.nextflow.version}
     Pipeline: ${workflow.manifest.name} ${workflow.manifest.version}
    ========================================
    """.stripIndent()

}


workflow PIPELINE_COMPLETION {
    take:
    outdir
    results      // channel: Final results from main workflow

    main:
    // Access params.outdir directly since it's a global parameter
    ch_final_results = results
    
    workflow.onComplete {
            log.info """
            ========================================
            Pipeline completed successfully!
            ========================================
            Output directory: ${outdir} 
            Completed at: ${new Date()}
            ========================================
            """.stripIndent()
    }

    workflow.onError { error ->
        log.error """
            ========================================
            Pipeline failed!
            ========================================
            Error: ${error.message}
            Exit status: ${workflow.exitStatus}
            Failed at: ${new Date()}
            ========================================
            """.stripIndent()
    }
    // ... rest of your workflow.onError block
}