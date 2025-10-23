process COUNT_LIGATIONS {
    tag "${sample}-${name}"
    label "hpc"

    publishDir "${params.output_dir}/${sample}/splits", mode: 'copy'

    input:
    tuple val(sample), val(name), path(read1), path(read2)

    output:
    path ("${name}${params.ext}_initial_norm.txt.res.txt"), emit: init_norm_res
    path ("${name}${params.ext}_linecount.txt"), emit: linecount

    script:
    def output_norm_txt = "${name}${params.ext}_initial_norm.txt.res.txt"
    def output_linecount_txt = "${name}${params.ext}_linecount.txt"
    """
    export LC_ALL=C
    export LC_COLLATE=C

    usegzip=0
    if [[ "${read1}" == *.gz ]]; then
        usegzip=1
    fi

    if [[ -z "${params.ligation}" || "${params.ligation}" == "XXXX" ]]; then
        echo "Skipping ligation match counting due to unset or 'XXXX' ligation motif"
        num1=0
    else
        if [[ \$usegzip -eq 1 ]]; then
            num1=\$(paste <(gunzip -c ${read1}) <(gunzip -c ${read2}) | awk '!((NR+2)%4)' | grep -cE ${params.ligation})
        else
            num1=\$(paste ${read1} ${read2} | awk '!((NR+2)%4)' | grep -cE "${params.ligation}")
        fi
    fi

    if [[ \$usegzip -eq 1 ]]; then
        num2=\$(gunzip -c ${read1} | wc -l | awk '{print \$1}')
    else
        num2=\$(wc -l ${read1} | awk '{print \$1}')
    fi
    echo -ne "\${num1} " > ${output_norm_txt}
    echo "\${num2}" > ${output_linecount_txt}
    """
}


process BWA_ALIGN {
    tag "${sample}-${name}"
    label "hpc"

    container "docker://pinglese6022/bwa-mem2:2.3"

    publishDir "${params.output_dir}/${sample}/splits", mode: 'copy'

    input:
    tuple val(sample), val(name), path(read1), path(read2)

    output:
    path ("${name}${params.ext}.sam"), emit: aligned_sam

    script:
    def output_aligned = "${name}${params.ext}.sam"
    """
    bwa-mem2 mem -SP5M -t ${params.threads} ${params.refSeq} \\
        ${read1} ${read2} > ${output_aligned}
    """
}


process FRAGMENT {
    tag "${sample}-${name}"
    label "hpc"

    publishDir "${params.output_dir}/${sample}/splits", mode: 'copy'

    input:
    tuple val(sample), val(name), path(norm_txt)

    output:
    path "${name}${params.ext}.frag.txt", emit: frag_txt

    script:
    def output_frag_txt = "${name}${params.ext}.frag.txt"
    """
    if [ "${params.site}" != "none" ] && [ -e "${params.site_file}" ]; then
        fragment.pl \\
            ${norm_txt} \\
            ${output_frag_txt} \\
            ${params.site_file}
    elif [ "${params.site}" == "none" ] || [ "${params.nofrag}" -eq 1 ]; then
        awk '{printf("%s %s %s %d %s %s %s %d", \$1, \$2, \$3, 0, \$4, \$5, \$6, 1); for (i=7; i<=NF; i++) {printf(" %s",\$i);}printf("\n");}' ${norm_txt} > ${output_frag_txt}
    else
        echo "***! No ${norm_txt} file created"
        exit 1
    fi
    """
}


process SORT {
    tag "${sample}-${name}"
    label "hpc"

    publishDir "${params.output_dir}/${sample}/splits", mode: 'copy'

    input:
    val sample
    val name
    path frag_txt

    output:
    path "${name}${params.ext}_sort.txt", emit: sorted_frag_txt

    script:
    def tmpdir = "HIC_tmp"
    def output_sort_txt = "${name}${params.ext}_sort.txt"
    """
    mkdir -p ${tmpdir}

    sort -T ${tmpdir} -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n ${frag_txt} > ${output_sort_txt}

    rm -rf ${tmpdir}
    """
}


process CHIMERIC {
    tag "${sample}-${name}"
    label "hpc"

    publishDir "${params.output_dir}/${sample}/splits", mode: 'copy'

    input:
    val sample
    val name
    path aligned_sam
    path norm_res_txt_initial

    output:
    path "${name}${params.ext}_norm.txt", emit: norm_txt
    path "${name}${params.ext}_abnorm.sam", emit: abnorm_sam
    path "${name}${params.ext}_unmapped.sam", emit: unmapped_sam
    path "${name}${params.ext}_norm.sam", emit: norm_sam
    path "${name}${params.ext}_norm.txt.res.txt", emit: norm_res_txt

    script:
    def output_norm = "${name}${params.ext}_norm.txt"
    def output_abnorm = "${name}${params.ext}_abnorm.sam"
    def output_unmapped = "${name}${params.ext}_unmapped.sam"
    def output_norm_sam = "${name}${params.ext}_norm.sam"
    """
    touch ${output_abnorm} ${output_unmapped}

    cp ${norm_res_txt_initial} ${name}${params.ext}_norm.txt.res.txt

    chimeric_blacklist.awk -v "fname1"=${output_norm} \\
        -v "fname2"=${output_abnorm} \\
        -v "fname3"=${output_unmapped} \\
        -v "fname_sam"=${output_norm_sam} \\
        ${aligned_sam}
    """
}


workflow process_fragments {
    take:
    reads

    main:
    sample = reads.map { it -> it[0] }
    name = reads.map { it -> it[1] }

    init_norm_res = COUNT_LIGATIONS(reads).init_norm_res
    aligned = BWA_ALIGN(reads).aligned_sam
    CHIMERIC(sample, name, aligned, init_norm_res)

    // Create the final chimeric output channel: (sample, name, Path1, Path2, Path3, Path4)
    def chimeric_output_channel = CHIMERIC.out.norm_txt
        .join(CHIMERIC.out.abnorm_sam)
        .join(CHIMERIC.out.unmapped_sam)
        .join(CHIMERIC.out.norm_res_txt)
        .join(CHIMERIC.out.norm_sam)
        .map { norm_txt, abnorm_sam, unmapped_sam, norm_res_txt, norm_sam ->
            tuple(sample, name, norm_txt, abnorm_sam, unmapped_sam, norm_res_txt, norm_sam)
        }

    frag_txt = FRAGMENT(chimeric_output_channel.map { it -> it[0..2] }).frag_txt
    SORT(sample, name, frag_txt)
    sorted_fragments_channel = SORT.out.sorted_frag_txt.map { it ->
        tuple(sample, it)
    }

    emit:
    chimeric_output  = chimeric_output_channel
    sorted_fragments = sorted_fragments_channel
}
