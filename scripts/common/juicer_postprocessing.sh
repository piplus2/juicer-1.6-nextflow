#!/bin/bash
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
# Juicer postprocessing script.
# This will run the major post-processing on the HiC file, including finding
# loops with HiCCUPS, finding motifs of these loops with MotifFinder, and
# finding contact domains with Arrowhead.
# Juicer version 1.5

## Read arguments
usageHelp="Usage: ${0} [-h] -j <juicer_tools_file_path> -i <hic_file_path> -m <bed_file_dir> -g <genome ID> -t <num_threads>"

printHelpAndExit() {
  echo "$usageHelp"
  exit $1
}

#set defaults
n_threads=1
genomeID="mm10"
hic_file_path="$(pwd)/aligned/inter_30.hic"
juicer_tools_path="/opt/juicer/scripts/juicer_tools"
bed_file_dir="/opt/juicer/references/motif"

while getopts "h:g:j:i:m:t:" opt; do
  case $opt in
  h) printHelpAndExit 0 ;;
  j) juicer_tools_path=$OPTARG ;;
  i) hic_file_path=$OPTARG ;;
  m) bed_file_dir=$OPTARG ;;
  g) genomeID=$OPTARG ;;
  t) n_threads=$OPTARG ;;
  [?]) printHelpAndExit 1 ;;
  esac
done

aligned_dir="$(dirname ${hic_file_path})"
log_dir="${aligned_dir}/logs"
mkdir -p "${log_dir}"
# Set the logging
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>"${log_dir}/juicer_postprocessing.log" 2>&1

## Check that juicer_tools exists
if [ ! -e "${juicer_tools_path}" ]; then
  echo "***! Can't find juicer tools in ${juicer_tools_path}"
  exit 1
fi

## Check that hic file exists
if [ ! -e "${hic_file_path}" ]; then
  echo "***! Can't find inter.hic in ${hic_file_path}"
  exit 1
fi

echo "${juicer_tools_path} is post-processing Hi-C for ${genomeID}"
echo "Data read from ${hic_file_path}."
echo "Motifs read from ${bed_file_dir}"

# Arrowhead

echo "##############"
echo "# ARROWHEAD: #"
echo "##############"

# Set ignore_sparsity to true so that Arrowhead will run even if the map is sparse
${juicer_tools_path} arrowhead --ignore_sparsity --threads ${n_threads} ${hic_file_path} ${hic_file_path%.*}"_contact_domains"

if [ $? -ne 0 ]; then
  echo "***! Problem while running Arrowhead"
  exit 1
fi

# HiCCUPS

echo "############"
echo "# HiCCUPS: #"
echo "############"

if hash nvcc 2>/dev/null; then
  # Use 1 thread otherwise it will give CUDA_OUT_OF_MEMORY error
  # Set ignore_sparsity to true so that HiCCUPS will run even if the map is sparse
  ${juicer_tools_path} hiccups --ignore_sparsity ${hic_file_path} ${hic_file_path%.*}"_loops"
  if [ $? -ne 0 ]; then
    echo "***! Problem while running HiCCUPS"
    exit 1
  fi
else
  echo "GPUs are not installed so HiCCUPs cannot be run"
fi

# Annotations

if [ -f ${hic_file_path%.*}"_loops/merged_loops.bedpe" ]; then

  echo "########"
  echo "# APA: #"
  echo "########"

  apa_dir="$(dirname ${hic_file_path})/apa_results"
  ${juicer_tools_path} apa ${hic_file_path} ${hic_file_path%.*}"_loops/merged_loops.bedpe" ${apa_dir}

  ## Check that bed folder exists
  if [ ! -e "${bed_file_dir}" ]; then
    echo "***! Can't find folder ${bed_file_dir}"
    echo "***! Not running motif finder"
  else

    echo "#################"
    echo "# MOTIF FINDER: #"
    echo "#################"

    cp ${hic_file_path%.*}"_loops/merged_loops.bedpe" ${hic_file_path%.*}"_loops.txt"
    ${juicer_tools_path} motifs ${genomeID} ${bed_file_dir} ${hic_file_path%.*}"_loops.txt"
  fi

  echo -e "\n(-: Feature annotation successfully completed (-:"
else
  # if loop lists do not exist but Juicer Tools didn't return an error, likely
  # too sparse
  echo -e "\n(-: Postprocessing successfully completed, maps too sparse to annotate or GPUs unavailable (-:"
fi
