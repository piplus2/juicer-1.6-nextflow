#!/usr/bin/awk -f
##########
#The MIT License (MIT)
#
# Copyright (c) 2015 Aiden Lab
# ... (License content)
##########
#
# SAM Filter for Chimeric and Normal Hi-C Reads (Juicer Component)
# ----------------------------------------------------------------
# Requires gawk for bitwise operations (and()).
# Takes a SAM file, assumed to be sorted by read name, and classifies alignments.
#
# Output file fname1: Stores "Normal" pairs/chimeras in custom Juicer format
#                     (EXCLUDING quality string).
# Output file fname2: Stores "Abnormal" chimeric reads in original SAM format.
# Output file fname3: Stores fully unmapped read groups in original SAM format.
# Output file fname_sam: Stores Normal/Normal Chimeric pairs in original SAM format.
# Output file fname1.res.txt: Summary statistics of read classifications.
# ----------------------------------------------------------------

# ====================================================================
# UTILITY FUNCTIONS
# ====================================================================

# returns absolute value
function abs(value) {
  return (value<0?-value:value);
}

# returns minimum of two values
function min(value1, value2) {
  return (value1<value2?value1:value2);
}

# Determines canonical pair order for fname1 output
function less_than(strand1, chrom1, pos1, strand2, chrom2, pos2) {
  if (chrom1 < chrom2) return 1;
  if (chrom1 > chrom2) return 0;
  if (strand1 < strand2) return 1;
  if (strand1 > strand2) return 0;
  if (pos1 < pos2) return 1;
  if (pos1 > pos2) return 0;
  return 1;
}

# ====================================================================
# HELPER FUNCTION: PARSE AND ADJUST POSITION
#
# Parses SAM fields and adjusts alignment position based on CIGAR clipping.
# ====================================================================
function parse_and_adjust(ind, sam_record,   sam_fields, clip_parts, match_position, clip_length, aligned_length, current_cigar_segment)
{
    split(sam_record, sam_fields);

    # Parse primary fields (SAM fields: $1, $2, $3, $4, $5, $6, $10, $11)
    split(sam_fields[1], read_name_parts, "/");
    read_mate_index[ind] = read_name_parts[2];
    full_read_name[ind] = sam_fields[1];
    strand_flag[ind] = and(sam_fields[2], 16); # Bit 16 (0x10): Reverse strand
    chrom[ind] = sam_fields[3];
    adjusted_pos[ind] = sam_fields[4];
    mapq[ind] = sam_fields[5];
    cigar[ind] = sam_fields[6];
    sequence[ind] = sam_fields[10];
    quality_string[ind] = sam_fields[11];      # Captured, but NOT used in fname1 print
    is_mapped[ind] = and(sam_fields[2], 4) == 0; # Bit 4 (0x4): Read unmapped

    # Position Adjustment for Clipping

    # Forward strand (strand_flag == 0)
    if (strand_flag[ind] == 0) {
        # Leading Soft Clip (S) or Hard Clip (H): Adjust position backward
        if (sam_fields[6] ~/^[0-9]+S/) {
            split(sam_fields[6], clip_parts, "S");
            clip_length = clip_parts[1] + 0;
            adjusted_pos[ind] = adjusted_pos[ind] - clip_length;
        }
        else if (sam_fields[6] ~/^[0-9]+H/) {
            split(sam_fields[6], clip_parts, "H");
            clip_length = clip_parts[1] + 0;
            adjusted_pos[ind] = adjusted_pos[ind] - clip_length;
        }
        # Ensure position is at least 1
        if (adjusted_pos[ind] <= 0) {
            adjusted_pos[ind] = 1;
        }
    }
    # Reverse strand (strand_flag == 16): Adjust position forward to the end
    else if (strand_flag[ind] == 16) {
        aligned_length = 0;
        current_cigar_segment = sam_fields[6];
        # Sum all aligned operations (M, D, N, X, =) for alignment length
        match_position = match(current_cigar_segment, /[0-9]+[M|D|N|X|=]/);
        while (match_position > 0) {
            aligned_length += substr(current_cigar_segment, match_position, RLENGTH-1) + 0;
            current_cigar_segment = substr(current_cigar_segment, match_position + RLENGTH);
            match_position = match(current_cigar_segment, /[0-9]+[M|D|N|X|=]/);
        }

        # New position = POS + aligned_length - 1
        adjusted_pos[ind] = adjusted_pos[ind] + aligned_length - 1;

        # Add trailing Soft/Hard clips to the end position
        if (sam_fields[6] ~ /[0-9]+S$/) {
            match_position = match(sam_fields[6],/[0-9]+S$/);
            clip_length = substr(sam_fields[6], match_position, RLENGTH-1) + 0;
            adjusted_pos[ind] = adjusted_pos[ind] + clip_length;
        }
        else if (sam_fields[6] ~ /[0-9]+H$/) {
            match_position = match(sam_fields[6],/[0-9]+H$/);
            clip_length = substr(sam_fields[6], match_position, RLENGTH-1) + 0;
            adjusted_pos[ind] = adjusted_pos[ind] + clip_length;
        }
    }
}


# ====================================================================
# AWK BLOCKS
# ====================================================================

BEGIN{
  OFS="\t";
  total_read_groups = -1;
  count_unmapped = 0;
  count_normal_pairs = 0;
  count_normal_chimeras = 0;
  count_abnormal_chimeras = 0;
}

# Main processing block (called for every line/SAM record)
{
  # Print SAM header lines (@HD, @SQ, @PG, etc.) to all SAM output files
  # and skip processing them as alignment records.
  if ($1 ~ /^@/) {
      print $0 > fname_sam;
      next;
  }

  # Grouping records by base read name (QNAME, field 1)
  split($1, read_name_parts, "/");
  current_read_group_id = read_name_parts[1];

  if (current_read_group_id == prev_read_group_id) {
    sam_record_count++;
  }
  else {
    # Process previous read group
    total_read_groups++;

    # --------------------------------------------------------
    # CLASSIFY AND PROCESS PREVIOUS READ GROUP
    # --------------------------------------------------------

    # CHIMERIC READS: 3 or 4 alignments
    if (sam_record_count == 3 || sam_record_count == 4) {

      # 1. Parse fields and calculate adjusted positions
      for (i=1; i <= sam_record_count; i++) {
        parse_and_adjust(i, read_group_records[i]);
      }

      # COUNT == 4 logic (A/B, A/B)
      if (sam_record_count == 4) {
        inter_chr_penalty = 10000000;
	    pair_distance[12] = abs(chrom[1]-chrom[2]) * inter_chr_penalty + abs(adjusted_pos[1]-adjusted_pos[2]);
	    pair_distance[13] = abs(chrom[1]-chrom[3]) * inter_chr_penalty + abs(adjusted_pos[1]-adjusted_pos[3]);
	    pair_distance[14] = abs(chrom[1]-chrom[4]) * inter_chr_penalty + abs(adjusted_pos[1]-adjusted_pos[4]);
	    pair_distance[23] = abs(chrom[2]-chrom[3]) * inter_chr_penalty + abs(adjusted_pos[2]-adjusted_pos[3]);
	    pair_distance[24] = abs(chrom[2]-chrom[4]) * inter_chr_penalty + abs(adjusted_pos[2]-adjusted_pos[4]);
	    pair_distance[34] = abs(chrom[3]-chrom[4]) * inter_chr_penalty + abs(adjusted_pos[3]-adjusted_pos[4]);

	    mate1_index = 0; mate2_index = 0;

	    if ((pair_distance[13] < 1000 && pair_distance[24] < 1000) || (pair_distance[14] < 1000 && pair_distance[23] < 1000)) {
	      mate1_index = 1; mate2_index = 2;
	    }
	    else if ((pair_distance[12] < 1000 && pair_distance[34] < 1000)) {
	      mate1_index = 1; mate2_index = 3;
	    }

	    if (mate1_index != 0) {
	      if (is_mapped[mate1_index] && is_mapped[mate2_index]) {
	        count_normal_chimeras++;

            # Output original SAM records to fname_sam (1-to-1 match)
            for (i in read_group_records) { print read_group_records[i] > fname_sam; }

            # Output in Juicer format to fname1 (NOW WITHOUT QUALITY)
	        if (less_than(strand_flag[mate1_index], chrom[mate1_index], adjusted_pos[mate1_index], strand_flag[mate2_index], chrom[mate2_index], adjusted_pos[mate2_index])) {
	          print strand_flag[mate1_index], chrom[mate1_index], adjusted_pos[mate1_index], strand_flag[mate2_index], chrom[mate2_index], adjusted_pos[mate2_index], mapq[mate1_index], cigar[mate1_index], sequence[mate1_index], mapq[mate2_index], cigar[mate2_index], sequence[mate2_index], full_read_name[mate1_index], full_read_name[mate2_index] > fname1;
	        } else {
	          print strand_flag[mate2_index], chrom[mate2_index], adjusted_pos[mate2_index], strand_flag[mate1_index], chrom[mate1_index], adjusted_pos[mate1_index], mapq[mate2_index], cigar[mate2_index], sequence[mate2_index], mapq[mate1_index], cigar[mate1_index], sequence[mate1_index], full_read_name[mate2_index], full_read_name[mate1_index] > fname1;
	        }
	      } else {
	        # Unmapped
	        for (i in read_group_records) { print read_group_records[i] > fname3; }
	        count_unmapped++;
	      }
	    } else {
	      # Abnormal chimera (too far apart)
	      count_abnormal_chimeras++;
	      for (i in read_group_records) { print read_group_records[i] > fname2; }
	    }
      }
      # COUNT == 3 logic (A/B...B)
      else {
        inter_chr_penalty = 10000000;
	    pair_distance[12] = abs(chrom[1]-chrom[2]) * inter_chr_penalty + abs(adjusted_pos[1]-adjusted_pos[2]);
	    pair_distance[23] = abs(chrom[2]-chrom[3]) * inter_chr_penalty + abs(adjusted_pos[2]-adjusted_pos[3]);
	    pair_distance[13] = abs(chrom[1]-chrom[3]) * inter_chr_penalty + abs(adjusted_pos[1]-adjusted_pos[3]);

	    if (min(pair_distance[12], min(pair_distance[23], pair_distance[13])) < 1000) {
	      if (read_mate_index[1] == read_mate_index[2]) { mate2_index = 3; mate1_index = pair_distance[13] > pair_distance[23] ? 1:2; }
	      else if (read_mate_index[1] == read_mate_index[3]) { mate2_index = 2; mate1_index = pair_distance[12] > pair_distance[23] ? 1:3; }
	      else if (read_mate_index[2] == read_mate_index[3]) { mate2_index = 1; mate1_index = pair_distance[12] > pair_distance[13] ? 2:3; }
	      else { printf("Error: Read mate indices strange for count=3 group.\n") > "/dev/stderr"; exit 1; }

	      if (is_mapped[mate1_index] && is_mapped[mate2_index]) {
	        count_normal_chimeras++;

            # Output original SAM records to fname_sam (1-to-1 match)
            for (i in read_group_records) { print read_group_records[i] > fname_sam; }

	        # Output in Juicer format to fname1 (NOW WITHOUT QUALITY)
	        if (less_than(strand_flag[mate1_index], chrom[mate1_index], adjusted_pos[mate1_index], strand_flag[mate2_index], chrom[mate2_index], adjusted_pos[mate2_index])) {
	          print strand_flag[mate1_index], chrom[mate1_index], adjusted_pos[mate1_index], strand_flag[mate2_index], chrom[mate2_index], adjusted_pos[mate2_index], mapq[mate1_index], cigar[mate1_index], sequence[mate1_index], mapq[mate2_index], cigar[mate2_index], sequence[mate2_index], full_read_name[mate1_index], full_read_name[mate2_index] > fname1;
	        } else {
	          print strand_flag[mate2_index], chrom[mate2_index], adjusted_pos[mate2_index], strand_flag[mate1_index], chrom[mate1_index], adjusted_pos[mate1_index], mapq[mate2_index], cigar[mate2_index], sequence[mate2_index], mapq[mate1_index], cigar[mate1_index], sequence[mate1_index], full_read_name[mate2_index], full_read_name[mate1_index] > fname1;
	        }
	      } else {
	        # Unmapped
	        for (i in read_group_records) { print read_group_records[i] > fname3; }
	        count_unmapped++;
	      }
	    } else {
	      # Abnormal chimera (too far apart)
	      count_abnormal_chimeras++;
	      for (i in read_group_records) { print read_group_records[i] > fname2; }
	    }
      }
    }

    # ABNORMAL/UNEXPECTED READS: > 4 or 1 alignment
    else if (sam_record_count > 4 || sam_record_count == 1) {
      count_abnormal_chimeras++;
      for (i in read_group_records) { print read_group_records[i] > fname2; }
    }

    # NORMAL PAIRED-END READ: 2 alignments
    else if (sam_record_count == 2) {
      # 1. Parse fields (use 0 and 1 for array indexing)
      mate_index = 0;
      for (i in read_group_records) {
        parse_and_adjust(mate_index, read_group_records[i]);
        mate_index++;
      }

      # 2. Process: Both mates must be mapped
      if (is_mapped[0] && is_mapped[1]) {
	    count_normal_pairs++;

        # Output original SAM records to fname_sam (1-to-1 match)
        for (i in read_group_records) { print read_group_records[i] > fname_sam; }

	    # Output in Juicer format to fname1 (NOW WITHOUT QUALITY)
	    if (less_than(strand_flag[0], chrom[0], adjusted_pos[0], strand_flag[1], chrom[1], adjusted_pos[1])) {
	      print strand_flag[0], chrom[0], adjusted_pos[0], strand_flag[1], chrom[1], adjusted_pos[1], mapq[0], cigar[0], sequence[0], mapq[1], cigar[1], sequence[1], full_read_name[0], full_read_name[1] > fname1;
	    } else {
	      print strand_flag[1], chrom[1], adjusted_pos[1], strand_flag[0], chrom[0], adjusted_pos[0], mapq[1], cigar[1], sequence[1], mapq[0], cigar[0], sequence[0], full_read_name[1], full_read_name[0] > fname1;
	    }
      } else {
	    # Unmapped
	    for (i in read_group_records) { print read_group_records[i] > fname3; }
	    count_unmapped++;
      }
    }

    # Reset variables and store current record for new group
    # Note: Clearing arrays used in parse_and_adjust is good practice
    delete read_group_records;
    delete read_mate_index;
    delete full_read_name;
    delete strand_flag;
    delete chrom;
    delete adjusted_pos;
    delete mapq;
    delete cigar;
    delete sequence;
    delete quality_string;
    delete is_mapped;

    sam_record_count = 1;
    prev_read_group_id = current_read_group_id;
  }

  # Store current SAM line
  read_group_records[sam_record_count] = $0;
}

# END Block: Process the very last read group
END{
  if (sam_record_count > 0) {
    total_read_groups++;
  }

  # --------------------------------------------------------
  # CLASSIFY AND PROCESS LAST READ GROUP
  # --------------------------------------------------------

  # CHIMERIC READS: 3 or 4 alignments
  if (sam_record_count == 3 || sam_record_count == 4) {
    # 1. Parse fields and calculate adjusted positions
    for (i=1; i <= sam_record_count; i++) {
        parse_and_adjust(i, read_group_records[i]);
    }

    # COUNT == 4 logic (A/B, A/B)
    if (sam_record_count == 4) {
      inter_chr_penalty = 10000000;
      pair_distance[12] = abs(chrom[1]-chrom[2]) * inter_chr_penalty + abs(adjusted_pos[1]-adjusted_pos[2]);
      pair_distance[13] = abs(chrom[1]-chrom[3]) * inter_chr_penalty + abs(adjusted_pos[1]-adjusted_pos[3]);
      pair_distance[14] = abs(chrom[1]-chrom[4]) * inter_chr_penalty + abs(adjusted_pos[1]-adjusted_pos[4]);
      pair_distance[23] = abs(chrom[2]-chrom[3]) * inter_chr_penalty + abs(adjusted_pos[2]-adjusted_pos[3]);
      pair_distance[24] = abs(chrom[2]-chrom[4]) * inter_chr_penalty + abs(adjusted_pos[2]-adjusted_pos[4]);
      pair_distance[34] = abs(chrom[3]-chrom[4]) * inter_chr_penalty + abs(adjusted_pos[3]-adjusted_pos[4]);

      mate1_index = 0; mate2_index = 0;

      if ((pair_distance[13] < 1000 && pair_distance[24] < 1000) || (pair_distance[14] < 1000 && pair_distance[23] < 1000)) {
        mate1_index = 1; mate2_index = 2;
      }
      else if ((pair_distance[12] < 1000 && pair_distance[34] < 1000)) {
        mate1_index = 1; mate2_index = 3;
      }

      if (mate1_index != 0) {
        if (is_mapped[mate1_index] && is_mapped[mate2_index]) {
          count_normal_chimeras++;
          # Output original SAM records to fname_sam (1-to-1 match)
          for (i in read_group_records) { print read_group_records[i] > fname_sam; }

          # Output in Juicer format to fname1 (NOW WITHOUT QUALITY)
          if (less_than(strand_flag[mate1_index], chrom[mate1_index], adjusted_pos[mate1_index], strand_flag[mate2_index], chrom[mate2_index], adjusted_pos[mate2_index])) {
            print strand_flag[mate1_index], chrom[mate1_index], adjusted_pos[mate1_index], strand_flag[mate2_index], chrom[mate2_index], adjusted_pos[mate2_index], mapq[mate1_index], cigar[mate1_index], sequence[mate1_index], mapq[mate2_index], cigar[mate2_index], sequence[mate2_index], full_read_name[mate1_index], full_read_name[mate2_index] > fname1;
          } else {
            print strand_flag[mate2_index], chrom[mate2_index], adjusted_pos[mate2_index], strand_flag[mate1_index], chrom[mate1_index], adjusted_pos[mate1_index], mapq[mate2_index], cigar[mate2_index], sequence[mate2_index], mapq[mate1_index], cigar[mate1_index], sequence[mate1_index], full_read_name[mate2_index], full_read_name[mate1_index] > fname1;
          }
        } else {
          for (i in read_group_records) { print read_group_records[i] > fname3; }
          count_unmapped++;
        }
      }
      else {
        count_abnormal_chimeras++;
        for (i in read_group_records) { print read_group_records[i] > fname2; }
      }
    }
    # COUNT == 3 logic (A/B...B)
    else {
      inter_chr_penalty = 10000000;
      pair_distance[12] = abs(chrom[1]-chrom[2]) * inter_chr_penalty + abs(adjusted_pos[1]-adjusted_pos[2]);
      pair_distance[23] = abs(chrom[2]-chrom[3]) * inter_chr_penalty + abs(adjusted_pos[2]-adjusted_pos[3]);
      pair_distance[13] = abs(chrom[1]-chrom[3]) * inter_chr_penalty + abs(adjusted_pos[1]-adjusted_pos[3]);

      if (min(pair_distance[12], min(pair_distance[23], pair_distance[13])) < 1000) {
        if (read_mate_index[1] == read_mate_index[2]) { mate2_index = 3; mate1_index = pair_distance[13] > pair_distance[23] ? 1:2; }
        else if (read_mate_index[1] == read_mate_index[3]) { mate2_index = 2; mate1_index = pair_distance[12] > pair_distance[23] ? 1:3; }
        else if (read_mate_index[2] == read_mate_index[3]) { mate2_index = 1; mate1_index = pair_distance[12] > pair_distance[13] ? 2:3; }
        else { printf("Error: Read mate indices strange in END block for count=3 group.\n") > "/dev/stderr"; exit 1; }

        if (is_mapped[mate1_index] && is_mapped[mate2_index]) {
          count_normal_chimeras++;
          # Output original SAM records to fname_sam (1-to-1 match)
          for (i in read_group_records) { print read_group_records[i] > fname_sam; }

          # Output in Juicer format to fname1 (NOW WITHOUT QUALITY)
          if (less_than(strand_flag[mate1_index], chrom[mate1_index], adjusted_pos[mate1_index], strand_flag[mate2_index], chrom[mate2_index], adjusted_pos[mate2_index])) {
            print strand_flag[mate1_index], chrom[mate1_index], adjusted_pos[mate1_index], strand_flag[mate2_index], chrom[mate2_index], adjusted_pos[mate2_index], mapq[mate1_index], cigar[mate1_index], sequence[mate1_index], mapq[mate2_index], cigar[mate2_index], sequence[mate2_index], full_read_name[mate1_index], full_read_name[mate2_index] > fname1;
          } else {
            print strand_flag[mate2_index], chrom[mate2_index], adjusted_pos[mate2_index], strand_flag[mate1_index], chrom[mate1_index], adjusted_pos[mate1_index], mapq[mate2_index], cigar[mate2_index], sequence[mate2_index], mapq[mate1_index], cigar[mate1_index], sequence[mate1_index], full_read_name[mate2_index], full_read_name[mate1_index] > fname1;
          }
        } else {
          for (i in read_group_records) { print read_group_records[i] > fname3; }
          count_unmapped++;
        }
      }
      else {
        count_abnormal_chimeras++;
        for (i in read_group_records) { print read_group_records[i] > fname2; }
      }
    }
  }

  # ABNORMAL/UNEXPECTED READS: > 4 or 1 alignment
  else if (sam_record_count > 4 || sam_record_count == 1) {
    count_abnormal_chimeras++;
    for (i in read_group_records) { print read_group_records[i] > fname2; }
  }

  # NORMAL PAIRED-END READ: 2 alignments
  else if (sam_record_count == 2) {
    # 1. Parse fields
    mate_index = 0;
    for (i in read_group_records) {
        parse_and_adjust(mate_index, read_group_records[i]);
        mate_index++;
    }

    # 2. Process: Both mates must be mapped
    if (is_mapped[0] && is_mapped[1]) {
	  count_normal_pairs++;
      # Output original SAM records to fname_sam (1-to-1 match)
      for (i in read_group_records) { print read_group_records[i] > fname_sam; }

      # Output in Juicer format to fname1 (NOW WITHOUT QUALITY)
	  if (less_than(strand_flag[0], chrom[0], adjusted_pos[0], strand_flag[1], chrom[1], adjusted_pos[1])) {
	    print strand_flag[0], chrom[0], adjusted_pos[0], strand_flag[1], chrom[1], adjusted_pos[1], mapq[0], cigar[0], sequence[0], mapq[1], cigar[1], sequence[1], full_read_name[0], full_read_name[1] > fname1;
	  } else {
	    print strand_flag[1], chrom[1], adjusted_pos[1], strand_flag[0], chrom[0], adjusted_pos[0], mapq[1], cigar[1], sequence[1], mapq[0], cigar[0], sequence[0], full_read_name[1], full_read_name[0] > fname1;
	  }
    } else {
	  for (i in read_group_records) { print read_group_records[i] > fname3; }
	  count_unmapped++;
    }
  }

  # Write final statistics to a results file
  results_file = fname1 ".res.txt";
  printf("%d %d %d %d %d\n", total_read_groups, count_unmapped, count_normal_pairs, count_normal_chimeras, count_abnormal_chimeras) >> results_file;
}
