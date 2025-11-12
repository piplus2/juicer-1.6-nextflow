process MERGE_SORT {
    tag "${sample}"
    label "mediumcpu"

    publishDir "${params.outdir}/${sample}/aligned", mode: 'copy'

    input:
    tuple val(sample), path(sort_files)

    output:
    tuple val(sample), path("${sample}${params.ext}_merged_sort.txt")

    script:
    def tmpdir = "HIC_tmp"
    def output_merged_sort_txt = "${sample}${params.ext}_merged_sort.txt"
    """
    mkdir -p ${tmpdir}

    sort -T ${tmpdir} -m \\
        -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n \\
        ${sort_files.join(' ')} > ${output_merged_sort_txt}

    rm -rf ${tmpdir}
    """
}

process REMOVE_DUPLICATES {
    tag "${sample}"
    label "mediumcpu"

    publishDir "${params.outdir}/${sample}/aligned", mode: 'copy'

    input:
    tuple val(sample), path(merged_sort)

    output:
    tuple val(sample), path("merged_nodups.txt"), path("dups.txt"), path("opt_dups.txt")

    script:
    def merged_nodups = "merged_nodups.txt"
    def dups = "dups.txt"
    def optdups = "optdups.txt"
    def opt_dups = "opt_dups.txt"
    """
    mkdir -p aligned

    touch ${dups}
    touch ${optdups}
    touch ${merged_nodups}

    if [ ${params.justexact} -eq 1 ]; then
        dups.awk \\
            -v name="" \\
            -v nowobble=1 \\
            ${merged_sort}
    else
        dups.awk \\
            -v name="" \\
            ${merged_sort}
    fi
    mv ${optdups} ${opt_dups}
    """
}
