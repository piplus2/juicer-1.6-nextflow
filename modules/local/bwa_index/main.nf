process PREPARE_BWA_INDEX {
    tag "${reference_fasta.simpleName}"
    label "highcpu"
    label "bwa"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bwa-mem2%3A2.3--he70b90d_0'
        : 'pinglese6022/bwa-mem2:2.3'}"

    input:
    path reference_fasta

    output:
    path index_dir, type: 'dir'

    script:
    index_dir = "bwa_index_files"
    """
    bwa-mem2 index -p bwa_index ${reference_fasta}
    mkdir -p ${index_dir}
    mv bwa_index.* ${index_dir}/
    """
}
