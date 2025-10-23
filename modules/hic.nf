// Process: Generate .hic files (normal and inter_30.hic)
process GEN_HIC_FILES {
    tag "${sample}"
    label "hpc"

    publishDir "${params.output_dir}/${sample}/aligned", mode: 'copy'

    input:
    tuple val(sample), path(inter), path(inter_30), path(inter_hists), path(merged_nodups)

    output:
    tuple val(sample), path("inter.hic"), path("inter_30.hic"), path("inter_30_hists.m")

    script:
    def genome = params.genomeID

    def inter_30_hists = "inter_30_hists.m"
    def inter_hic = "inter.hic"
    def inter_30_hic = "inter_30.hic"
    def site_opt = params.nofrag.toString() == "1" ? "" : "-f ${params.site_file}"
    def resolutions = params.resolutions == null || params.resolutions.toString() == "" ? "" : "-r ${params.resolutions}"
    """
    export _JAVA_OPTIONS=-Xmx${params.java_mem}
    export LC_ALL=en_US.UTF-8

    juicer_tools pre \\
        ${site_opt} \\
        -s ${inter} \\
        -g ${inter_hists} \\
        -q 1 \\
        ${resolutions} \\
        ${merged_nodups} \\
        ${inter_hic} \\
        ${genome}

    statistics.pl \\
        -s ${params.site_file} \\
        -l ${params.ligation} \\
        -o ${inter_30} \\
        -q 30 \\
        ${merged_nodups}

    juicer_tools pre \\
        ${site_opt} \\
        -s ${inter_30} \\
        -g ${inter_30_hists} \\
        -q 30 \\
        ${resolutions} \\
        ${merged_nodups} \\
        ${inter_30_hic} \\
        ${genome}
    """
}
