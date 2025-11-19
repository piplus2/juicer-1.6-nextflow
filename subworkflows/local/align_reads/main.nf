include { BLACKLIST_CHIMERIC } from '../../../modules/local/chimeric'
include { SORT               } from '../../../modules/local/sort_reads'
include { CONVERT_FRAGMENTS  } from '../../../modules/local/fragments'
include { BWA_ALIGN          } from '../../../modules/local/bwa'
include { COUNT_LIGATIONS    } from '../../../modules/local/count_ligations'
include { PREPARE_BWA_INDEX  } from '../../../modules/local/bwa_index'


workflow process_fragments {
    take:
    reads

    main:

    // output = (sample, name, init_norm_res, linecount)
    init_norm_res = COUNT_LIGATIONS(reads)

    index_dir = PREPARE_BWA_INDEX(params.reference)

    // output = (sample, name, aligned_sam)
    aligned_sams = BWA_ALIGN(reads.combine(index_dir))

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
    chimeric = BLACKLIST_CHIMERIC(chimeric_input_ch)

    fragment_input_ch = chimeric.map { sample, name, norm_txt, _abnorm_sam, _unmapped_sam, _norm_sam, _norm_res_txt ->
        tuple(sample, name, norm_txt)
    }

    // output = (sample, sort_txt)
    sorted_files = SORT(CONVERT_FRAGMENTS(fragment_input_ch))

    emit:
    chimeric_output  = chimeric
    sorted_fragments = sorted_files
}
