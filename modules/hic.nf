// Process: Generate .hic files (normal and inter_30.hic)
process GEN_HIC_FILES {
    tag "${sample}"
    label "hpc"
    cache 'deep'

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


process JUICER_TOOLS_PRE_Q1 {
    tag "${sample}"
    label 'hpc'
    cache 'deep'

    publishDir "${params.output_dir}/${sample}/aligned", mode: 'copy', pattern: "*.hic"

    input:
    tuple val(sample), path(input_file), path(inter), path(inter_hists)

    output:
    tuple val(sample), path("${input_file.baseName}.hic")

    script:
    def site_opt = params.nofrag.toString() == "1" ? "" : "-f ${params.site_file}"
    def resolutions = params.resolutions == null || params.resolutions.toString() == "" ? "" : "-r ${params.resolutions}"
    def output_file = "${input_file.baseName}.hic"
    """
    export _JAVA_OPTIONS=-Xmx${params.java_mem}
    export LC_ALL=en_US.UTF-8

    juicer_tools pre \\
        ${site_opt} \\
        -s ${inter} \\
        -g ${inter_hists} \\
        -q 1 \\
        ${resolutions} \\
        ${input_file} \\
        ${output_file} \\
        ${params.genomeID}
    """
}


process JUICER_TOOLS_PRE_Q30 {
    tag "${sample}"
    label 'hpc'
    cache 'deep'

    publishDir "${params.output_dir}/${sample}/aligned", mode: 'copy', pattern: "*.hic"

    input:
    tuple val(sample), path(merged_no_dups), path(inter_30), path(inter_hists)

    output:
    tuple val(sample), path("${inter_30.baseName}.hic")

    script:
    def site_opt = params.nofrag.toString() == "1" ? "" : "-f ${params.site_file}"
    def resolutions = params.resolutions == null || params.resolutions.toString() == "" ? "" : "-r ${params.resolutions}"
    def output_file = "${inter_30.baseName}.hic"
    """
    export _JAVA_OPTIONS=-Xmx${params.java_mem}
    export LC_ALL=en_US.UTF-8

    juicer_tools pre \\
        ${site_opt} \\
        -s ${inter_30} \\
        -g ${inter_hists} \\
        -q 30 \\
        ${resolutions} \\
        ${merged_no_dups} \\
        ${output_file} \\
        ${params.genomeID}
    """
}


process HIC_STATS {
    tag "${sample}"
    label 'hpc'
    cache 'deep'

    publishDir "${params.output_dir}/${sample}/aligned", mode: 'copy'

    input:
    tuple val(sample), path(inter_30), path(merged_nodups)

    output:
    tuple val(sample), path(inter_30), path("inter_30_hists.m")

    script:
    """
    cp ${inter_30} "${inter_30.simpleName}_local"

    statistics.pl \\
        -s ${params.site_file} \\
        -l ${params.ligation} \\
        -o ${inter_30.simpleName}_local \\
        -q 30 \\
        ${merged_nodups}

    mv "${inter_30.simpleName}_local.hists.m" inter_30_hists.m
    mv "${inter_30.simpleName}_local" ${inter_30}
    """
}


workflow hic {
    take:
    inter_ch

    main:
    //     tuple val(sample), path(inter_txt), path(inter_30_txt), path(inter_hists), path(merged_nodups)

    pre_q1_ch = inter_ch.map { sample, inter_txt, _inter_30_txt, inter_hists_m, merged_no_dups ->
        tuple(sample, merged_no_dups, inter_txt, inter_hists_m)
    }

    // output : tuple val(sample), path(inter_hic)
    JUICER_TOOLS_PRE_Q1(pre_q1_ch)

    hic_stats_input = inter_ch.map { sample, _inter_txt, inter_30_txt, _inter_hists_m, merged_nodups ->
        tuple(sample, inter_30_txt, merged_nodups)
    }
    // output: tuple val(sample), path(inter_30_hic), path(inter_30_hists)
    hic_stats_out = HIC_STATS(hic_stats_input)
    // output : tuple val(sample), path(inter_30_hic)
    pre_q30_input = hic_stats_out.map { sample, inter_30, inter_30_hists ->
        tuple(sample, inter_30, inter_30_hists)
    }
    // output : tuple val(sample), path(inter_30_hic)
    inter_30_hic_ch = JUICER_TOOLS_PRE_Q30(pre_q30_input)

    emit:
    hic_out_ch = inter_30_hic_ch
}
