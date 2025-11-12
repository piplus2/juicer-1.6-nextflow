process COUNT_LIGATIONS {
    tag "${sample}-${name}"
    label "smallcpu"

    publishDir "${params.outdir}/${sample}/splits", mode: 'copy'

    input:
    tuple val(sample), val(name), path(read1), path(read2)

    output:
    tuple val(sample), val(name), path("${name}${params.ext}_initial_norm.txt.res.txt"), path("${name}${params.ext}_linecount.txt")

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
    label "highcpu"
    label "bwa"

    publishDir "${params.outdir}/${sample}/splits", mode: 'copy'

    input:
    tuple val(sample), val(name), path(read1), path(read2)

    output:
    tuple val(sample), val(name), path("${name}${params.ext}.sam")

    script:
    def output_aligned = "${name}${params.ext}.sam"
    """
    bwa-mem2 mem -SP5M -t ${task.cpus} ${params.reference} \\
        ${read1} ${read2} > ${output_aligned}
    """
}


process FRAGMENT {
    tag "${sample}-${name}"
    label "smallcpu"

    publishDir "${params.outdir}/${sample}/splits", mode: 'copy'

    input:
    tuple val(sample), val(name), path(norm_txt)

    output:
    tuple val(sample), val(name), path("${name}${params.ext}.frag.txt")

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
    label "mediumcpu"

    publishDir "${params.outdir}/${sample}/splits", mode: 'copy'

    input:
    tuple val(sample), val(name), path(frag_txt)

    output:
    tuple val(sample), path("${name}${params.ext}_sort.txt")

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
    label "mediumcpu"

    publishDir "${params.outdir}/${sample}/splits", mode: 'copy'

    input:
    tuple val(sample), val(name), path(aligned_sam), path(norm_res_txt_initial)

    output:
    tuple val(sample), val(name), path("${name}${params.ext}_norm.txt"), path("${name}${params.ext}_abnorm.sam"), path("${name}${params.ext}_unmapped.sam"), path("${name}${params.ext}_norm.sam"), path("${name}${params.ext}_norm.txt.res.txt")

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

    // output = (sample, name, init_norm_res, linecount)
    init_norm_res = COUNT_LIGATIONS(reads)

    // output = (sample, name, aligned_sam)
    aligned_sams = BWA_ALIGN(reads)

    // Prepare CHIMERIC inputs
    chimeric_input_ch = init_norm_res
        .map { sample, name, norm_res_txt, _linecount ->
            tuple(sample, name, norm_res_txt)
        }
        .join(aligned_sams, by: [0, 1])

    chimeric_input_ch = chimeric_input_ch.map { sample, name, norm_txt, aligned_sam ->
        tuple(sample, name, aligned_sam, norm_txt)
    }

    // output = (sample, name, norm_txt, abnorm_sam, unmapped_sam, norm_sam, norm_res_txt)
    chimeric = CHIMERIC(chimeric_input_ch)

    fragment_input_ch = chimeric.map { sample, name, norm_txt, _abnorm_sam, _unmapped_sam, _norm_sam, _norm_res_txt ->
        tuple(sample, name, norm_txt)
    }

    // output = (sample, sort_txt)
    sorted_files = SORT(FRAGMENT(fragment_input_ch))

    emit:
    chimeric_output  = chimeric
    sorted_fragments = sorted_files
}
