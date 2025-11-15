// Process: Arrowhead (contact domain calling)
process ARROWHEAD {
    tag "${sample}"
    label "gpu"

    container 'docker://pinglese6022/juicer_tools:1.22.01'

    publishDir "${params.outdir}/${sample}/aligned", mode: 'copy'

    input:
    tuple val(sample), path(inter_30_hic)

    output:
    path "${inter_30_hic.simpleName}_contact_domains", type: 'dir'

    script:
    def contact_domains = "${inter_30_hic.simpleName}_contact_domains"
    """
    export LC_ALL=en_US.UTF-8

    mkdir -p ${contact_domains}

    /usr/local/bin/juicer_tools -Xmx${params.java_mem} arrowhead \\
        --threads 0 \\
        --ignore-sparsity \\
        ${inter_30_hic} \\
        ${contact_domains}
    """
}
