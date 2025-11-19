process JUICER_TOOLS_PRE_Q1 {
    tag "${sample}"
    label 'highcpu'
    label "juicertools"

    publishDir "${params.outdir}/${sample}/aligned", mode: 'copy', pattern: "*.hic"

    input:
    tuple val(sample), path(input_file), path(inter), path(inter_hists)
    path site_file

    output:
    tuple val(sample), path("${input_file.baseName}.hic")

    script:
    def site_opt = params.nofrag.toString() == "1" ? "" : "-f ${site_file}"
    def resolutions = params.resolutions == null || params.resolutions.toString() == "" ? "" : "-r ${params.resolutions}"
    def output_file = "${input_file.baseName}.hic"
    """
    export LC_ALL=en_US.UTF-8
    export _JAVA_OPTIONS="-Xmx${params.java_mem}"

    juicer_tools pre ${site_opt} -s ${inter} -g ${inter_hists} -q 1 ${resolutions} ${input_file} ${output_file} ${params.genome_id}
    """
}
