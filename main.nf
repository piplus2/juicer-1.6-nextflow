include { NFCORE_JUICER } from './workflow/main.nf'

workflow {
    // TODO: make a PREPARE_GENOME process that prepares the genome for juicer, and then pass the prepared genome to NFCORE_JUICER

    NFCORE_JUICER()
}
