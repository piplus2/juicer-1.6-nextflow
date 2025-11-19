process HIC_STATS {
    tag "${sample}"
    label 'mediumcpu'

    publishDir "${params.outdir}/${sample}/aligned", mode: 'copy'

    input:
    tuple val(sample), path(inter_filt_txt), path(merged_nodups)

    output:
    tuple val(sample), path(inter_filt_txt), path(inter_filt_hists_m)

    script:
    inter_filt_hists_m = "inter_${params.mapq}_hists.m"
    def inter_filt_local_hists_m = "${inter_filt_txt.simpleName}_local_hists.m"
    def inter_filt_local_txt = "${inter_filt_txt.simpleName}_local"
    """
    cp ${inter_filt_txt} "${inter_filt_local_txt}"

    statistics.pl \\
        -s ${params.site_file} \\
        -l ${params.ligation} \\
        -o ${inter_filt_local_txt} \\
        -q ${params.mapq} \\
        ${merged_nodups}

    mv "${inter_filt_local_hists_m}" ${inter_filt_hists_m}
    mv "${inter_filt_local_txt}" ${inter_filt_txt}
    """
}
