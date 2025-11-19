#!/usr/bin/env python

# Generate site positions in genome from given restriction enzyme
# Juicer 1.5

from __future__ import print_function

import argparse
import os
import re
import sys


def usage():
    print('Usage: {} <restriction enzyme> <genome> [location]'.format(
        sys.argv[0]), file=sys.stderr)
    sys.exit(1)

# ------------------------------------------------------------------------------


DEFAULT_FASTA = {
    'hg19': 'references/Homo_sapiens_assembly19.fasta',
    'mm9': 'references/Mus_musculus_assembly9.fasta',
    'mm10': 'references/Mus_musculus_assembly10.fasta',
    'hg18': 'references/Homo_sapiens_assembly18.fasta',
}

ENZYME_PATTERNS = {
    'HindIII': 'AAGCTT',
    'DpnII': 'GATC',
    'MboI': 'GATC',
    'Sau3AI': 'GATC',
    'Arima': ['GATC', 'GANTC'],
}


def parse_args():
    parser = argparse.ArgumentParser(
        description='Generate restriction enzyme site positions in genome.')
    parser.add_argument(
        '--enzyme', type=str, help='Restriction enzyme name.', choices=ENZYME_PATTERNS.keys(), required=True)
    parser.add_argument(
        '--genome', type=str, help='Genome assembly name.', required=True)
    parser.add_argument(
        '--location', type=str, required=False, default=None,
        help='Path to genome FASTA file (optional).')
    return parser.parse_args()


def process_input(input_file, output_file, pattern):

    f = open(input_file, 'r')
    g = open(output_file, 'w')

    minsize = min([len(p) for p in pattern])
    maxsize = max([len(p) for p in pattern])
    matches = get_match_func(pattern)

    segment = ''
    counter = 0
    endl = ''

    for line in f:

        line = line.strip()

        if line.startswith('>'):

            # This is the beginning of a new sequence, but before starting it we must
            # finish processing of the remaining segment of the previous sequence.

            while len(segment) > minsize:
                segment = segment[1:]
                if matches(segment):
                    g.write(' ' + str(counter - len(segment) + 1))

            if counter > 0:
                # Close the previous sequence here.
                g.write(' ' + str(counter))

            firststr = re.split(r'\s+', line[1:])
            g.write(endl+firststr[0])

            segment = ''
            counter = 0
            endl = '\n'

            continue

        # Process next line of the sequence.

        line = line.upper()

        for symbol in line:

            counter += 1
            segment += symbol

            while len(segment) > maxsize:
                segment = segment[1:]

            # Do pattern matching only if segment size equals maxsize.

            if len(segment) == maxsize:
                if matches(segment):
                    # maxsize == len(segment)
                    g.write(' ' + str(counter - maxsize + 1))

    # Finish the last sequence.

    while len(segment) > minsize:
        segment = segment[1:]
        if matches(segment):
            g.write(' ' + str(counter - len(segment) + 1))

    if counter > 0:
        g.write(' ' + str(counter))

    g.write('\n')  # End the output file with a newline.

    # Close files.

    g.close()
    f.close()


def main():
    args = parse_args()
    if not args.genome and not args.location:
        raise ValueError("Either genome or location must be specified.")
    if args.genome and args.location:
        print(
            "Warning: Both genome and location specified. Using location.",
            file=sys.stderr)
        genome = args.genome
        reference = args.location
    elif args.genome not in DEFAULT_FASTA and not args.location:
        raise ValueError(
            f"Genome '{args.genome}' not found in default references. Please provide location.")
    elif args.genome not in DEFAULT_FASTA and args.location:
        genome = args.genome
        reference = args.location
    else:
        genome = args.genome
        reference = DEFAULT_FASTA[genome]

    if not os.path.isfile(reference):
        raise FileNotFoundError(
            f"Reference file '{reference}' does not exist.")

    output_file = f"{genome}_{args.enzyme}_sites.txt"
    pattern = ENZYME_PATTERNS[args.enzyme]

    process_input(reference, output_file, pattern)


# ------------------------------------------------------------------------------


def has_wildcard(pattern):

    # Input pattern can be a list or a string.

    wildcards = re.compile(r'[NMRWYSKHBVD]')

    if (isinstance(pattern, list)):
        for p in pattern:
            if re.search(wildcards, p):
                return True
    else:
        if re.search(wildcards, pattern):
            return True

    return False

# ------------------------------------------------------------------------------


def pattern2regexp(pattern):

    # Input pattern can be a list or a string.

    wildcards = {
        'N': '[ACGT]',
        'M': '[AC]',
        'R': '[AG]',
        'W': '[AT]',
        'Y': '[CT]',
        'S': '[CG]',
        'K': '[GT]',
        'H': '[ACT]',
        'B': '[CGT]',
        'V': '[ACG]',
        'D': '[AGT]',
    }

    if isinstance(pattern, list):
        return [pattern2regexp(p) for p in pattern]

    pattern = pattern.upper()
    for p, r in wildcards.items():
        pattern = re.sub(p, r, pattern)

    return re.compile(pattern.upper())

# ------------------------------------------------------------------------------


def get_match_func(pattern):

    # Input pattern can be a list or a string.

    if not isinstance(pattern, list):
        pattern = [pattern]

    if has_wildcard(pattern):

        pattern = pattern2regexp(pattern)

        if len(pattern) == 1:  # There is only a single pattern.

            # Use the only element from the list as a single regexp.
            pattern = pattern[0]

            def match_single_regexp(segment):
                if re.match(pattern, segment):
                    return True
                return False

            return match_single_regexp

        else:  # There are multiple patterns.

            def match_multi_regexp(segment):
                for p in pattern:
                    if re.match(p, segment):
                        return True
                return False

            return match_multi_regexp

    else:  # No wildcard in any of the patterns.

        if len(pattern) == 1:  # There is only a single pattern.

            # Use the only element from the list as a single string.
            pattern = pattern[0]

            def match_single_string(segment):
                if segment.startswith(pattern):
                    return True
                return False

            return match_single_string

        else:  # There are multiple patterns.

            def match_multi_string(segment):
                for p in pattern:
                    if segment.startswith(p):
                        return True
                return False

            return match_multi_string

# ------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
