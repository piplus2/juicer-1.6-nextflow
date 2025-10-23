process MAKE_HEADERFILE {
    tag "${sample}"
    label "hpc"

    publishDir "${params.output_dir}/${sample}/aligned", mode: 'copy'

    input:
    val sample

    output:
    tuple val(sample), path("header")

    script:
    def header = "header"
    """
    date > ${header}
    echo -ne 'Experiment description: ${params.about};' >> ${header}
    echo -ne "Juicer version: ${params.juicer_version};" >> ${header}
    bwa 2>&1 | awk '\$1=="Version:"{printf(" BWA %s; ", \$2)}' >> ${header}
    echo -ne "${params.threads} threads; " >> ${header}
    java -version 2>&1 | awk 'NR==1{printf("%s; ", \$0);}' >> ${header}
    juicer_tools -V 2>&1 | awk '\$1=="Juicer" && \$2=="Tools"{printf("%s; ", \$0);}' >> ${header}
    echo "" >> ${header}
    """
}
