process STATS {
    tag "${sample}"
    label "highcpu"

    container 'docker://pinglese6022/juicer_tools:1.22.01'
    // TODO: add the mount option for docker when nextflow supports it
    containerOptions "-B ${baseDir}/bin:/workflow/bin"

    publishDir "${params.outdir}/${sample}/aligned", mode: 'copy'

    input:
    tuple val(sample), path(header), path(res_txts), path(abnorm_sams), path(unmapped_sams), path(merged_nodups), path(dups_txt), path(opt_dups_txt)
    path site_file

    output:
    tuple val(sample), path(inter_txt), path(inter_30_txt), path(inter_30_hists_m), path(collisions_txt), path(abnormal_sam), path(unmapped_sam)

    script:
    inter_txt = "inter.txt"
    inter_30_hists_m = "inter_30_hists.m"
    inter_30_txt = "inter_30.txt"
    collisions_txt = "collisions.txt"
    abnormal_sam = "abnormal.sam"
    unmapped_sam = "unmapped.sam"
    """
    export LC_ALL=en_US.UTF-8

    tail -n1 ${header} | awk '{printf "%-1000s\\n", \$0}' > ${inter_txt}
    cat ${res_txts.join(' ')} | /workflow/bin/stats_sub.awk >> ${inter_txt}
    /usr/local/bin/juicer_tools -Xmx${params.java_mem} LibraryComplexity "./" ${inter_txt} >> ${inter_txt}

    cp ${inter_txt} ${inter_30_txt}

    if [[ -n ${site_file} ]]; then
        perl /workflow/bin/statistics.pl -s ${site_file} -l ${params.ligation} -o ${inter_30_txt} -q 30 ${merged_nodups}
    else
        perl /workflow/bin/statistics.pl -l ${params.ligation} -o ${inter_30_txt} -q 30 ${merged_nodups}
    fi

    cat ${abnorm_sams.join(' ')} > ${abnormal_sam}
    cat ${unmapped_sams.join(' ')} > ${unmapped_sam}

    /workflow/bin/collisions.awk ${abnormal_sam} > ${collisions_txt}

    /workflow/bin/collisions_dedup_rearrange_cols.awk -v fname=${collisions_txt} \\
        ${collisions_txt} \\
        | sort -k3,3n -k4,4n -k10,10n -k11,11n -k17,17n -k18,18n -k24,24n -k25,25n -k31,31n -k32,32n \\
        | /workflow/bin/collisions_dups.awk -v name=.
    """
}
