process MERGE_SORT {
    tag { sample }
    label "mediumcpu"

    publishDir "${params.outdir}/${sample}/aligned", mode: 'copy'

    input:
    tuple val(sample), path(sort_files)

    output:
    tuple val(sample), path(output_merged_sort_txt)

    script:
    def tmpdir = "HIC_tmp"
    output_merged_sort_txt = "${sample}${params.ext}_merged_sort.txt"
    """
    mkdir -p ${tmpdir}

    sort -T ${tmpdir} -m \\
        -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n \\
        ${sort_files.join(' ')} > ${output_merged_sort_txt}

    rm -rf ${tmpdir}
    """
}
