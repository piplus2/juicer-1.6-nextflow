nextflow.enable.dsl = 2

include { NFCORE_JUICER } from './workflow/main.nf'

workflow {
    main:
        NFCORE_JUICER()
}
