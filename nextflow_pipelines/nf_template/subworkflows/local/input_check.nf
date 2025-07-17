/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUT_CHECK SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Validate input samplesheet and create sample channels
*/

workflow INPUT_CHECK {
    take:
    samplesheet // file: path to samplesheet

    main:
    // Read and parse samplesheet
    ch_samples = Channel
        .fromPath(samplesheet)
        .splitCsv(header: true)
        .map { row ->
            validateSampleRow(row)
        }
        .map { row ->
            createSampleMeta(row)
        }

    emit:
    samples = ch_samples
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def validateSampleRow(LinkedHashMap row) {
    // Check required columns exist
    def required_cols = ['sample_id']
    required_cols.each { col ->
        if (!row.containsKey(col) || !row[col]) {
            error("Missing required column '${col}' in samplesheet row: ${row}")
        }
    }
    
    return row
}

def createSampleMeta(LinkedHashMap row) {
    // Create meta map and file list
    def meta = [:]
    meta.id = row.sample_id
    
    // For this simple example, create a dummy file list
    def files = []
    if (row.containsKey('file_path') && row.file_path) {
        files = [file(row.file_path)]
    } else {
        // Create empty file for demonstration
        files = []
    }
    
    return [meta, files]
}