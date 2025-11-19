process JUICER_TOOLS_PRE_QFILT {
    tag "${sample}"
    label 'highcpu'
    label "juicertools"

    publishDir "${params.outdir}/${sample}/aligned", mode: 'copy', pattern: "*.hic"

    input:
    tuple val(sample), path(merged_no_dups), path(inter_filt_txt), path(inter_filt_hists_m)
    path site_file

    output:
    tuple val(sample), path(output_file)

    script:
    def site_opt = params.nofrag.toString() == "1" ? "" : "-f ${site_file}"
    def resolutions = params.resolutions == null || params.resolutions.toString() == "" ? "" : "-r ${params.resolutions}"
    output_file = "inter_${params.mapq}.hic"
    """
    export LC_ALL=en_US.UTF-8
    export _JAVA_OPTIONS="-Xmx${params.java_mem}"

    juicer_tools pre \\
        ${site_opt} \\
        -s ${inter_filt_txt} \\
        -g ${inter_filt_hists_m} \\
        -q ${params.mapq} \\
        ${resolutions} \\
        ${merged_no_dups} \\
        ${output_file} \\
        ${params.genome_id}
    """
}
