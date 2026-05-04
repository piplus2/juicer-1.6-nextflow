include { samplesheetToList } from 'plugin/nf-schema'

workflow PIPELINE_INITIALISATION {

    take:
    input
    outdir
    nextflow_args
    version

    main:

    ch_versions = channel.empty()

    // Validate the input here

    ch_samplesheet = channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
    .map {
        meta, fastq_1, fastq_2 ->
        if (!fastq_2) {
            return [ meta.id, meta + [ single_end:true ], [fastq_1]]
            }
        else {
            return [ meta.id, meta + [ single_end:false ], [fastq_1, fastq_2]]
            }
    }
    .groupTuple()
    .map { samplesheet -> validateInputSamplesheet(samplesheet)
    }
    .map {
        meta, fastqs ->
            return [ meta, fastqs.flatten() ]
    }

    emit:
    samplesheet = ch_samplesheet

}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    // Emit a warning if `single_end` is true
    if (metas[0].single_end == true) {
        log.warn "Sample ${metas[0].id} is detected as single-end reads (fastq_1 only)."
        metas[0].single_end = false
    }

    return [ metas[0], fastqs ]
}