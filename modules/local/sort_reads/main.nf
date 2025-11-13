process SORT {
    tag "${sample}-${name}"
    label "mediumcpu"

    publishDir "${params.outdir}/${sample}/splits", mode: 'copy'

    input:
    tuple val(sample), val(name), path(frag_txt)

    output:
    tuple val(sample), path(output_sort_txt)

    script:
    def tmpdir = "HIC_tmp"
    output_sort_txt = "${name}${params.ext}_sort.txt"
    """
    mkdir -p ${tmpdir}

    sort -T ${tmpdir} -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n ${frag_txt} > ${output_sort_txt}

    rm -rf ${tmpdir}
    """
}
