process STATS {
    tag "${sample}"
    label "hpc"

    publishDir "${params.outdir}/${sample}/aligned", mode: 'copy'

    input:
    tuple val(sample), path(header), path(res_txts), path(abnorm_sams), path(unmapped_sams), path(merged_nodups), path(dups_txt), path(opt_dups_txt)

    output:
    tuple val(sample), path("inter.txt"), path("inter_30.txt"), path("inter_hists.m"), path("collisions.txt"), path("abnormal.sam"), path("unmapped.sam")

    script:
    def inter = "inter.txt"
    def inter_30 = "inter_30.txt"
    def collisions = "collisions.txt"
    def abnorm_sam = "abnormal.sam"
    def unmapped_sam = "unmapped.sam"
    """
    export _JAVA_OPTIONS=-Xmx${params.java_mem}
    export LC_ALL=en_US.UTF-8

    tail -n1 ${header} | awk '{printf "%-1000s\\n", \$0}' > ${inter}
    cat ${res_txts.join(' ')} | stats_sub.awk >> ${inter}
    juicer_tools LibraryComplexity "./" ${inter} >> ${inter}

    cp ${inter} ${inter_30}

    stat_cmd="statistics.pl"
    if [[ -n "${params.site_file}" ]]; then
        stat_cmd+=" -s \"${params.site_file}\""
    fi
    stat_cmd+=" -l \"${params.ligation}\" -o \"${inter}\" -q 1 \"${merged_nodups}\""

    eval "\${stat_cmd}"

    cat ${abnorm_sams.join(' ')} > ${abnorm_sam}
    cat ${unmapped_sams.join(' ')} > ${unmapped_sam}

    collisions.awk ${abnorm_sam} > ${collisions}

    collisions_dedup_rearrange_cols.awk -v fname=${collisions} \\
        ${collisions} \\
        | sort -k3,3n -k4,4n -k10,10n -k11,11n -k17,17n -k18,18n -k24,24n -k25,25n -k31,31n -k32,32n \\
        | collisions_dups.awk -v name=.
    """
}
