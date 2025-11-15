include { JUICER_TOOLS_PRE_Q1  } from '../../../modules/local/juicer_tools/pre_q1'
include { HIC_STATS            } from '../../../modules/local/hic_stats'
include { JUICER_TOOLS_PRE_Q30 } from '../../../modules/local/juicer_tools/pre_q30'


workflow hic {
    take:
    inter_ch // tuple val(sample), path(inter_txt), path(inter_30_txt), path(inter_hists_m), path(merged_nodups)

    main:

    pre_q1_ch = inter_ch.map { sample, inter_txt, _inter_30_txt, inter_hists_m, merged_no_dups ->
        tuple(sample, merged_no_dups, inter_txt, inter_hists_m)
    }

    // output : tuple val(sample), path(inter_hic)
    site_file = file(params.site_file)
    JUICER_TOOLS_PRE_Q1(pre_q1_ch, site_file)

    hic_stats_input = inter_ch.map { sample, _inter_txt, inter_30_txt, _inter_hists_m, merged_nodups ->
        tuple(sample, inter_30_txt, merged_nodups)
    }
    // output: tuple val(sample), path(inter_30_txt), path(inter_30_hists_m)
    hic_stats_out = HIC_STATS(hic_stats_input)

    pre_q30_ch = inter_ch
        .map { sample, _inter_txt, _inter_30_txt, _inter_hists_m, merged_nodups ->
            tuple(sample, merged_nodups)
        }
        .join(
            hic_stats_out,
            by: 0
        )

    // Input: tuple val(sample), path(merged_no_dups), path(inter_30_txt), path(inter_30_hists_m)
    // Output: tuple val(sample), path(inter_30_hic)
    inter_30_hic_ch = JUICER_TOOLS_PRE_Q30(pre_q30_ch, site_file)

    emit:
    hic_out_ch = inter_30_hic_ch
}
