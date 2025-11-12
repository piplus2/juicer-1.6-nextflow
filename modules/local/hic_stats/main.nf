process HIC_STATS {
    tag "${sample}"
    label 'mediumcpu'

    publishDir "${params.output_dir}/${sample}/aligned", mode: 'copy'

    input:
    tuple val(sample), path(inter_30_txt), path(merged_nodups)

    output:
    tuple val(sample), path(inter_30_txt), path("inter_30_hists.m")

    script:
    """
    cp ${inter_30_txt} "${inter_30_txt.simpleName}_local"

    statistics.pl \\
        -s ${params.site_file} \\
        -l ${params.ligation} \\
        -o ${inter_30_txt.simpleName}_local \\
        -q 30 \\
        ${merged_nodups}

    mv "${inter_30_txt.simpleName}_local_hists.m" inter_30_hists.m
    mv "${inter_30_txt.simpleName}_local" ${inter_30_txt}
    """
}
