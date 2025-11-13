process HIC_STATS {
    tag "${sample}"
    label 'mediumcpu'

    publishDir "${params.output_dir}/${sample}/aligned", mode: 'copy'

    input:
    tuple val(sample), path(inter_30_txt), path(merged_nodups)

    output:
    tuple val(sample), path(inter_30_txt), path(inter_30_hists_m)

    script:
    inter_30_hists_m = "inter_30_hists.m"
    def inter_30_local_hists_m = "${inter_30_txt.simpleName}_local_hists.m"
    def inter_30_local_txt = "${inter_30_txt.simpleName}_local"
    """
    cp ${inter_30_txt} "${inter_30_local_txt}"

    statistics.pl \\
        -s ${params.site_file} \\
        -l ${params.ligation} \\
        -o ${inter_30_local_txt} \\
        -q 30 \\
        ${merged_nodups}

    mv "${inter_30_local_hists_m}" ${inter_30_hists_m}
    mv "${inter_30_local_txt}" ${inter_30_txt}
    """
}
