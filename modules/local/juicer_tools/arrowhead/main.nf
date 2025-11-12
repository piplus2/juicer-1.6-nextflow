// Process: Arrowhead (contact domain calling)
process ARROWHEAD {
    tag "${sample}"
    label "gpu"

    publishDir "${params.outdir}/${sample}/aligned", mode: 'copy'

    input:
    tuple val(sample), path(inter_30_hic)

    output:
    path "${inter_30_hic.simpleName}_contact_domains", type: 'dir'

    script:
    def contact_domains = "${inter_30_hic.simpleName}_contact_domains"
    """
    export _JAVA_OPTIONS=-Xmx${params.java_mem}
    export LC_ALL=en_US.UTF-8

    mkdir -p ${contact_domains}

    juicer_tools arrowhead \\
        --threads 0 \\
        --ignore-sparsity \\
        ${inter_30_hic} \\
        ${contact_domains}
    """
}
