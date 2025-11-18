// Process: Arrowhead (contact domain calling)
process ARROWHEAD {
    tag "${sample}"
    label "highcpu"
    label "juicertools"

    env [
        'LC_ALL': 'en_US.UTF-8',
        '_JAVA_OPTIONS': "-Xmx${params.java_mem}",
    ]

    publishDir "${params.outdir}/${sample}/aligned", mode: 'copy'

    input:
    tuple val(sample), path(inter_30_hic)

    output:
    path "${inter_30_hic.simpleName}_contact_domains", type: 'dir'

    script:
    def contact_domains = "${inter_30_hic.simpleName}_contact_domains"
    """
    mkdir -p ${contact_domains}

    juicer_tools arrowhead \\
        --threads 0 \\
        --ignore-sparsity \\
        ${inter_30_hic} \\
        ${contact_domains}
    """
}
