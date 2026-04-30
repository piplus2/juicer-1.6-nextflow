nextflow.enable.dsl = 2

include { JUICER } from './workflow/main.nf'

workflow {
    main:
        JUICER()
}
