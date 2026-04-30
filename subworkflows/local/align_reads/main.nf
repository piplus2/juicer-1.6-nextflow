include { BLACKLIST_CHIMERIC } from '../../../modules/local/chimeric'
include { SORT               } from '../../../modules/local/sort_reads'
include { CONVERT_FRAGMENTS  } from '../../../modules/local/fragments'
include { COUNT_LIGATIONS    } from '../../../modules/local/count_ligations'
include { BWA_INDEX          } from '../../../modules/nf-core/bwa/index'
include { BWA_MEM            } from '../../../modules/nf-core/bwa/mem/main.nf'


workflow PROCESS_FRAGMENTS {
    take:
    reads     // Expecting: [sample, name, r1, r2]
    reference // Expecting a reference fasta file

    main:

    // TODO: this must go into the future PROCESS_GENOME module, but for now we need it here to prepare the BWA index
    ch_index = BWA_INDEX([[id: reference.baseName], reference]).index
    ch_fasta = [[id: reference.baseName], reference]

    // BWA_MEM expects [ [id: sample_name ], [r1, r2] ]
    ch_bwa_input = reads.map { sample, name, r1, r2 ->
        def meta = [id: "${sample}-${name}", sample: sample, name: name]
        [meta, [r1, r2]]
    }

    BWA_MEM(ch_bwa_input, ch_index, ch_fasta, false)

    // The aligned files are in the SAM format although they have a .bam extension
    aligned_sams = BWA_MEM.out.sam.map { meta, sam ->
        tuple(meta.sample, meta.name, sam)
    }

    // Prepare CHIMERIC inputs
    init_norm_res = COUNT_LIGATIONS(reads)
    chimeric_input_ch = init_norm_res
        .map { sample, name, norm_res_txt, _linecount ->
            tuple(sample, name, norm_res_txt)
        }
        .join(aligned_sams, by: [0, 1])

    chimeric_input_ch = chimeric_input_ch.map { sample, name, norm_txt, aligned_sam ->
        tuple(sample, name, aligned_sam, norm_txt)
    }

    // output = (sample, name, norm_txt, abnorm_sam, unmapped_sam, norm_sam, norm_res_txt)
    chimeric = BLACKLIST_CHIMERIC(chimeric_input_ch)

    fragment_input_ch = chimeric.map { sample, name, norm_txt, _abnorm_sam, _unmapped_sam, _norm_sam, _norm_res_txt ->
        tuple(sample, name, norm_txt)
    }

    // output = (sample, sort_txt)
    sorted_files = SORT(CONVERT_FRAGMENTS(fragment_input_ch))

    emit:
    chimeric_output  = chimeric
    sorted_fragments = sorted_files
    aligned_sams     = aligned_sams
}
