#!/usr/bin/gawk -f
#
# Filters a SAM file (ARGV[2]) by performing an exact match against the
# full QNAMEs found in the deduplicated list (ARGV[1]).
#
# Usage: gawk -f filter_sam_by_readname.awk merged_nodups.txt normal.sam > final_nodups.sam
#

BEGIN {
    # Set FS to one or more spaces for the *first* file (merged_nodups.txt)
    FS = "[[:space:]]+";
}

# --- Pass 1: Read Allowed Read Names from merged_nodups.txt (ARGV[1]) ---
FNR == NR {
    # The deduplicated file contains the FULL, unique QNAMEs in fields 15 and 16.
    # We store the full string for an exact match.
    allowed_reads[$15] = 1;
    allowed_reads[$16] = 1;

    next;
}

# --- Pass 2: Filter the SAM file (normal.sam, ARGV[2]) ---
{
    # IMPORTANT: Change FS to tab for the SAM file (ARGV[2]).
    if (FNR == 1) {
        FS = "\t";
    }

    # 1. Handle SAM Header
    if ($1 ~ /^@/) {
        print $0;
        next;
    }

    # 2. Check: Use the full SAM QNAME ($1) for direct lookup.
    # This guarantees that only records whose QNAME exactly matches one of the
    # QNAMEs stored in the 'allowed_reads' hash are outputted.
    if ($1 in allowed_reads) {
        print $0; # Print the full, original SAM record
    }
}