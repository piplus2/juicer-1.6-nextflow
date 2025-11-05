// Process: Generate .hic files (normal and inter_30.hic)
process GEN_HIC_FILES {
    tag "${sample}"
    label "hpc"

    publishDir "${params.outdir}/${sample}/aligned", mode: 'copy'

    input:
    tuple val(sample), path(inter), path(inter_30), path(inter_hists), path(merged_nodups)

    output:
    tuple val(sample), path("inter.hic"), path("inter_30.hic"), path("inter_30_hists.m")

    script:
    def genome = params.genome_id

    def inter_30_hists = "inter_30_hists.m"
    def inter_hic = "inter.hic"
    def inter_30_hic = "inter_30.hic"
    """
    export _JAVA_OPTIONS=-Xmx${params.java_mem}
    export LC_ALL=en_US.UTF-8

    pre_args=()
    if [[ "${params.nofrag}" != "1" && -n "${params.site_file}" ]]; then
        pre_args+=(-f "${params.site_file}")
    fi
    pre_args+=(-s "${inter}" -g "${inter_hists}" -q 1)
    if [[ -n "${params.resolutions}" ]]; then
        pre_args+=(-r "${params.resolutions}")
    fi
    pre_args+=("${merged_nodups}" "${inter_hic}" "${genome}")

    juicer_tools pre "${pre_args[@]}"

    stat_args=(-l "${params.ligation}" -o "${inter_30}" -q 30 "${merged_nodups}")
    if [[ -n "${params.site_file}" ]]; then
        stat_args=(-s "${params.site_file}" "${stat_args[@]}")
    fi

    statistics.pl "${stat_args[@]}"

    pre_30_args=()
    if [[ "${params.nofrag}" != "1" && -n "${params.site_file}" ]]; then
        pre_30_args+=(-f "${params.site_file}")
    fi
    pre_30_args+=(-s "${inter_30}" -g "${inter_30_hists}" -q 30)
    if [[ -n "${params.resolutions}" ]]; then
        pre_30_args+=(-r "${params.resolutions}")
    fi
    pre_30_args+=("${merged_nodups}" "${inter_30_hic}" "${genome}")

    juicer_tools pre "${pre_30_args[@]}"
    """
}
