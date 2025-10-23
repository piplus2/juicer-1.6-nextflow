#!/bin/bash

# ##########################################################################
# IMPORTANT:
# This is the custom script for the MegaMap pipeline that as been tested
# on the IIT HPC cluster.
# ##########################################################################

##########
#The MIT License (MIT)
#
# Copyright (c) 2015 Aiden Lab
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
##########
# MegaMap script.
#
#
# [topDir]    - Should contain the results of all base experiments
#
# From the top-level directory, the following two directories are created:
#
# [topDir]/mega     - Location of result of processing the mega map

# Juicer version 1.5
juicer_version="1.5.7"
# top level directory, can also be set in options
topDir=$(pwd)

# Juicer directory, contains scripts/, references/, and restriction_sites/
# can also be set in options via -D

# Use my source dir
juiceDir="/work/pinglese/tools/juicer-1.6"

usageHelp="Usage: ${0##*/} -g genomeID [-d topDir] [-s site] [-S stage] [-b ligation] [-D Juicer scripts directory] [-f] [-h]"
genomeHelp="   genomeID is either defined in the script, e.g. \"hg19\" or \"mm10\" or the path to the chrom.sizes file"
dirHelp="   [topDir] is the top level directory (default \"$topDir\") and must contain links to all merged_nodups files underneath it"
siteHelp="   [site] must be defined in the script, e.g.  \"HindIII\" or \"MboI\" (default \"$site\"); alternatively, this can be the restriction site file"
stageHelp="* [stage]: must be one of \"final\", \"postproc\", or \"early\".\n    -Use \"final\" when the reads have been combined into merged_nodups but the\n     final stats and hic files have not yet been created.\n    -Use \"postproc\" when the hic files have been created and only\n     postprocessing feature annotation remains to be completed.\n    -Use \"early\" for an early exit, before the final creation of the stats and\n     hic files"
ligationHelp="* [ligation junction]: use this string when counting ligation junctions"
scriptDirHelp="* [Juicer scripts directory]: set the Juicer directory,\n  which should have scripts/ references/ and restriction_sites/ underneath it\n  (default ${juiceDir})"
excludeHelp="   -f: include fragment-delimited maps in Hi-C mega map (will run slower)"
helpHelp="   -h: print this help and exit"

printHelpAndExit() {
    echo "$usageHelp"
    echo "$genomeHelp"
    echo "$dirHelp"
    echo "$siteHelp"
    echo "$stageHelp"
    echo "$ligationHelp"
    echo "$excludeHelp"
    echo "$helpHelp"
    exit "$1"
}

# Set defaults
nofrag="1"
# restriction enzyme, can also be set in options
site="none"
# genome ID, default to human, can also be set in options
genomeID="mm10"

while getopts "d:g:hfs:S:b:D:y:" opt; do
    case $opt in
    g) genomeID=$OPTARG ;;
    h) printHelpAndExit 0 ;;
    d) topDir=$OPTARG ;;
    s) site=$OPTARG ;;
    b) ligation=$OPTARG ;;
    y) site_file=$OPTARG ;;
    D) juiceDir=$OPTARG ;;
    f) nofrag="0" ;;
    S) stage=$OPTARG ;;
    [?]) printHelpAndExit 1 ;;
    esac
done

# Add logging
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>"${juiceDir}/mega_map.log" 2>&1

# resolutionsToBuildString="-r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000,500,200,100,50,20,10"

## Set ligation junction based on restriction enzyme
if [ -z "$ligation" ]; then
    case $site in
    HindIII) ligation="AAGCTAGCTT" ;;
    DpnII) ligation="GATCGATC" ;;
    MboI) ligation="GATCGATC" ;;
    Arima) ligation="'(GAATAATC|GAATACTC|GAATAGTC|GAATATTC|GAATGATC|GACTAATC|GACTACTC|GACTAGTC|GACTATTC|GACTGATC|GAGTAATC|GAGTACTC|GAGTAGTC|GAGTATTC|GAGTGATC|GATCAATC|GATCACTC|GATCAGTC|GATCATTC|GATCGATC|GATTAATC|GATTACTC|GATTAGTC|GATTATTC|GATTGATC)'" ;;
    none) ligation="XXXX" ;;
    *)
        ligation="XXXX"
        site_file=$site
        echo "$site not listed as recognized enzyme, so trying it as site file."
        echo "Ligation junction is undefined"
        ;;
    esac
fi

## If DNAse-type experiment, no fragment maps
if [ "$site" == "none" ]; then
    nofrag=1
fi

# In the case of Arima, you must provide the site file defined as <genome>_<site>.txt
if [ -z "$site_file" ]; then
    site_file="${juiceDir}/restriction_sites/${genomeID}_${site}.txt"
fi

## Check that site file exists, needed for fragment number for merged_nodups
if [ ! -e "$site_file" ] && [ "$site" != "none" ]; then
    echo "***! $site_file does not exist. It must be created before running this script."
    echo "The site file is used for statistics even if fragment delimited maps are excluded"
    exit 1
fi

if [ ! -z "$stage" ]; then
    case $stage in
    final) final=1 ;;
    early) early=1 ;;
    hic) hic=1 ;;
    postproc) postproc=1 ;;
    *)
        echo "$usageHelp"
        echo "$stageHelp"
        exit 1
        ;;
    esac
fi

## Directories to be created and regex strings for listing files
megadir=${topDir}"/mega"
outputdir=${megadir}"/aligned"
tmpdir=${megadir}"/HIC_tmp"
export TMPDIR=${tmpdir}
outfile=${megadir}/lsf.out
#output messages
logdir="$megadir/debug"

## Check for existing merged_nodups files:
merged_count=$(find -L ${topDir} | grep merged_nodups.txt | wc -l)
if [ "$merged_count" -lt "1" ]; then
    echo "***! Failed to find at least one merged_nodups files under ${topDir}"
    exit 1
fi

merged_names=$(find -L ${topDir} | grep merged_nodups.txt.gz | awk '{print "<(gunzip -c",$1")"}' | tr '\n' ' ')
if [ ${#merged_names} -eq 0 ]; then
    merged_names=$(find -L ${topDir} | grep merged_nodups.txt | tr '\n' ' ')
fi
inter_names=$(find -L ${topDir} | grep inter.txt | tr '\n' ' ')
if [[ $(echo ${inter_names} | wc -l) -ne "${merged_count}" ]]; then
    echo "***! Number of inter.txt files different from merged_nodups.txt"
fi

## Create output directory, exit if already exists
if [[ -d "${outputdir}" ]] && [ -z $final ] && [ -z $postproc ] && [ -z $hic ]; then
    echo "***! Move or remove directory \"${outputdir}\" before proceeding."
    exit 1
else
    mkdir -p ${outputdir}
fi

## Create temporary directory
if [ ! -d "$tmpdir" ]; then
    mkdir $tmpdir
    chmod 777 $tmpdir
fi

## Create output directory, used for reporting commands output
if [ ! -d "$logdir" ]; then
    mkdir "$logdir"
    chmod 777 "$logdir"
fi

## Arguments have been checked and directories created. Now begins
## the real work of the pipeline

export _JAVA_OPTIONS=-Xmx120g
export LC_ALL=en_US.UTF-8
NTHREADS=20

# Not in final or postproc
if [ -z $final ] && [ -z $postproc ]; then

    # Create top statistics file from all inter.txt files found under current dir
    echo "Creating top stats files from from inter.txt"
    awk -f ${juiceDir}/scripts/common/makemega_addstats.awk ${inter_names} >${outputdir}/inter.txt
    echo "(-: Finished creating top stats files."

    cp ${outputdir}/inter.txt ${outputdir}/inter_30.txt

    echo "Merging all merged_nodups files into a single file"
    sort --parallel=${NTHREADS} -T ${tmpdir} -m -k2,2d -k6,6d ${merged_names} >${outputdir}/merged_nodups.txt
    rm -r ${tmpdir}
    echo "(-: Finished sorting all merged_nodups files into a single merge."

    echo "Running statistics on merged_nodups file:"
    echo "- inter.txt"
    ${juiceDir}/scripts/common/statistics.pl \
        -q 1 \
        -o ${outputdir}/inter.txt \
        -s $site_file \
        -l $ligation \
        ${outputdir}/merged_nodups.txt
    echo "done"

    echo "- inter_30.txt"
    ${juiceDir}/scripts/common/statistics.pl \
        -q 30 \
        -o ${outputdir}/inter_30.txt \
        -s $site_file \
        -l $ligation \
        ${outputdir}/merged_nodups.txt
    echo "done"

    if [ "$nofrag" -eq "1" ]; then
        echo "Juicer Preprocessing without fragments:"

        echo "- inter.hic"
        ${juiceDir}/scripts/common/juicer_tools pre \
            -s ${outputdir}/inter.txt \
            -g ${outputdir}/inter_hists.m \
            -q 1 \
            ${outputdir}/merged_nodups.txt \
            ${outputdir}/inter.hic \
            ${genomeID}
        echo "done"

        echo "- inter_30.hic"
        ${juiceDir}/scripts/common/juicer_tools pre \
            -s ${outputdir}/inter_30.txt \
            -g ${outputdir}/inter_30_hists.m \
            -q 30 \
            ${outputdir}/merged_nodups.txt \
            ${outputdir}/inter_30.hic \
            ${genomeID}
        echo "done"
    else
        echo "Juicer Preprocessing with fragments:"

        echo "- inter.hic"
        ${juiceDir}/scripts/common/juicer_tools pre \
            -f ${site_file} \
            -s ${outputdir}/inter.txt \
            -g ${outputdir}/inter_hists.m \
            -q 1 \
            ${outputdir}/merged_nodups.txt \
            ${outputdir}/inter.hic \
            ${genomeID}
        echo "done"

        echo "- inter_30.hic"
        ${juiceDir}/scripts/common/juicer_tools pre \
            -f ${site_file} \
            -s ${outputdir}/inter_30.txt \
            -g ${outputdir}/inter_30_hists.m \
            -q 30 \
            ${outputdir}/merged_nodups.txt \
            ${outputdir}/inter_30.hic \
            ${genomeID}
        echo "done"
    fi
fi

if [ -z $early ]; then

    echo "Running postprocessing with MAPQ > 30"

    ${juiceDir}/scripts/common/juicer_postprocessing.sh \
        -j ${juiceDir}/scripts/common/juicer_tools \
        -i ${outputdir}/inter_30.hic \
        -m ${juiceDir}/references/motif \
        -t ${NTHREADS} \
        -g ${genomeID}

fi

echo "(-: Successfully completed making mega map. Done. :-)"
