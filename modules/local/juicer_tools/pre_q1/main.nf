process JUICER_TOOLS_PRE_Q1 {
    tag "${sample}"
    label 'highcpu'
    label "juicertools_1_22"

    publishDir "${params.outdir}/${sample}/aligned", mode: 'copy', pattern: "*.hic"

    input:
    tuple val(sample), path(input_file), path(inter), path(inter_hists)
    path site_file

    output:
    tuple val(sample), path(output_file)

    script:
    def site_opt = params.nofrag.toString() == "1" ? "" : "-f ${site_file}"
    def resolutions = params.resolutions == null || params.resolutions.toString() == "" ? "" : "-r ${params.resolutions}"
    output_file = "inter.hic"
    """
    export LC_ALL=en_US.UTF-8
    export _JAVA_OPTIONS="-Xmx${params.java_mem}"

    juicer_tools pre ${site_opt} -s ${inter} -g ${inter_hists} -q 1 ${resolutions} ${input_file} ${output_file} ${params.genome_id}
    """
}
