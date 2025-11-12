process JUICER_TOOLS_PRE_Q30 {
    tag "${sample}"
    label 'highcpu'

    publishDir "${params.output_dir}/${sample}/aligned", mode: 'copy', pattern: "*.hic"

    input:
    tuple val(sample), path(merged_no_dups), path(inter_30_txt), path(inter_30_hists_m)

    output:
    tuple val(sample), path("${inter_30_txt.baseName}.hic")

    script:
    def site_opt = params.nofrag.toString() == "1" ? "" : "-f ${params.site_file}"
    def resolutions = params.resolutions == null || params.resolutions.toString() == "" ? "" : "-r ${params.resolutions}"
    def output_file = "${inter_30_txt.baseName}.hic"
    """
    export _JAVA_OPTIONS=-Xmx${params.java_mem}
    export LC_ALL=en_US.UTF-8

    juicer_tools pre \\
        ${site_opt} \\
        -s ${inter_30_txt} \\
        -g ${inter_30_hists_m} \\
        -q 30 \\
        ${resolutions} \\
        ${merged_no_dups} \\
        ${output_file} \\
        ${params.genomeID}
    """
}
