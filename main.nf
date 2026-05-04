include { JUICER } from './workflows/juicer.nf'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_pipeline/main.nf'

workflow {
    main:

        // TODO: complete the pipeline initialisation
        PIPELINE_INITIALISATION(
            params.input,
            params.outdir,
            args,
            params.version
        )

        // TODO: add workflow for preparing the genome

        JUICER()
}
