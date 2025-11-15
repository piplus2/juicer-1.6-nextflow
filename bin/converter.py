#!/usr/bin/env python3
# join_normal_with_sam.py
# Rebuild "normal" pairs as original SAM records (preserving QUAL + tags).
# Inputs:
#   --norm: output of chimeric_blacklist (fname1)
#   --sam : original SAM (same input that produced --norm)
# Output:
#   --out: SAM containing only the two original records per normal pair, in normal-file order.

import argparse
import gzip
import sys

# -------------------------------- I/O helpers ------------------------------- #


def open_any(path, mode):
    if path == "-" or path is None:
        return sys.stdin if "r" in mode else sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, mode, compresslevel=6)
    return open(path, mode, buffering=1024*1024)

# ---------------- CIGAR-adjusted 5' coordinate (Juicer logic) --------------- #


def adj_pos(pos, cigar, rev):
    # forward: subtract leading S/H
    if not rev:
        i = 0
        while i < len(cigar) and cigar[i].isdigit():
            i += 1
        if i < len(cigar) and i > 0 and cigar[i] in ("S", "H"):
            pos = max(1, pos - int(cigar[:i]))
        return pos
    # reverse: pos + (sum of M/D/N/X/=) - 1 + trailing S/H
    i = 0
    span = 0
    while i < len(cigar):
        j = i
        while j < len(cigar) and cigar[j].isdigit():
            j += 1
        if j == i or j >= len(cigar):
            break
        n = int(cigar[i:j])
        op = cigar[j]
        if op in "MDNX=":
            span += n
        i = j + 1
    pos = pos + span - 1
    # trailing S/H
    k = len(cigar) - 1
    if k >= 0 and cigar[k] in ("S", "H"):
        j = k - 1
        while j >= 0 and cigar[j].isdigit():
            j -= 1
        if j+1 <= k-1:
            pos += int(cigar[j+1:k])
    return pos


def strand_flag(flag):  # returns 0 or 16
    return 16 if (flag & 0x10) else 0


def qname_root(q):
    return q.split("/")[0]

# ----------------------------------- Main ----------------------------------- #


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--norm", required=True)
    ap.add_argument("--sam",  required=True)
    ap.add_argument("--out",  required=True)
    args = ap.parse_args()

    # 1) Load targets from normal file (preserve order).
    # normal columns:
    #  s1 chr1 pos1 s2 chr2 pos2 mq1 cigar1 seq1 mq2 cigar2 seq2 name1 name2
    targets_by_root = {}  # root -> list of two dicts (one per selected end)
    # list of (root, index0/1) in the order we should output
    order = []
    with open_any(args.norm, "rt") as fn:
        for line in fn:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 14:
                continue
            s1, c1, p1, s2, c2, p2, mq1, cg1, sq1, mq2, cg2, sq2, n1, n2 = f[:14]
            root = qname_root(n1)
            pair = [
                {"root": root, "name": n1, "strand": int(s1), "chr": c1, "pos_adj": int(p1),
                 "cigar": cg1, "seq": sq1, "found": False, "samline": None},
                {"root": root, "name": n2, "strand": int(s2), "chr": c2, "pos_adj": int(p2),
                 "cigar": cg2, "seq": sq2, "found": False, "samline": None},
            ]
            # store
            if root in targets_by_root:
                targets_by_root[root].append(pair)
            else:
                targets_by_root[root] = [pair]
            # output order: (root, pair_idx, end_idx)
            idx = len(targets_by_root[root]) - 1
            order.append((root, idx, 0))
            order.append((root, idx, 1))

    # 2) Stream original SAM and capture the exact two records that match each target end.
    headers = []
    with open_any(args.sam, "rt") as fs:
        for line in fs:
            if line.startswith("@"):
                headers.append(line)
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 11:
                continue
            qn = cols[0]
            flag = int(cols[1])
            rname = cols[2]
            pos = int(cols[3])
            mapq = cols[4]  # not used for matching strictly
            cigar = cols[5]
            seq = cols[9]
            # unmapped or no CIGAR
            if cigar == "*" or rname == "*" or (flag & 0x4):
                continue
            root = qname_root(qn)
            if root not in targets_by_root:
                continue
            s = strand_flag(flag)
            p_adj = adj_pos(pos, cigar, rev=(s == 16))
            # try match any still-unfound target end across all pairs for this root
            for pair in targets_by_root[root]:
                for end in pair:
                    if end["found"]:
                        continue
                    # strict match on chr, strand, adj_pos; also require same CIGAR and SEQ to disambiguate
                    if (end["chr"] == rname and end["strand"] == s and end["pos_adj"] == p_adj
                            and end["cigar"] == cigar and end["seq"] == seq):
                        end["found"] = True
                        end["samline"] = line
                        break  # stop looking at other ends
                # optional tiny speedup: if both ends found, skip remaining pairs for this root only if all pairs done

    # 3) Write SAM: headers + matched records in normal-file order.
    with open_any(args.out, "wt") as fo:
        for h in headers:
            fo.write(h)
        fo.write("@PG\tID:normal-restore\tPN:normal-restore\tVN:1.0\n")
        # emit in the same order the normal file listed them (two lines per pair)
        for root, pair_idx, end_idx in order:
            pair = targets_by_root.get(root, [])
            if pair_idx >= len(pair):
                continue
            end = pair[pair_idx][end_idx]
            if end["samline"] is not None:
                fo.write(end["samline"])
            else:
                # not found: silently skip or write a warning to stderr
                sys.stderr.write(
                    f"[warn] no SAM match for {root} end{end_idx+1} at {end['chr']}:{end['pos_adj']}\n")


if __name__ == "__main__":
    main()
