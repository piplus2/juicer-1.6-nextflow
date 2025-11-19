// Process: Arrowhead (contact domain calling)
process ARROWHEAD {
    tag "${sample}"
    label "highmemory"
    label "juicertools_1_22"

    publishDir "${params.outdir}/${sample}/aligned", mode: 'copy'

    input:
    tuple val(sample), path(inter_filt_hic)

    output:
    path contact_domains, type: 'dir'

    script:
    contact_domains = "${inter_filt_hic.simpleName}_contact_domains"
    """
    export LC_ALL=en_US.UTF-8
    export _JAVA_OPTIONS="-Xmx${params.java_mem}"

    mkdir -p ${contact_domains}

    juicer_tools arrowhead \\
        --threads 0 \\
        --ignore-sparsity \\
        ${inter_filt_hic} \\
        ${contact_domains}
    """
}
