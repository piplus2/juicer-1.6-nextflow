#!/bin/bash

RUN_ID=241023_A00558_0323_AHC2G3DRX5
RAW_DIR=/fastqstore_gefa/BaseCalls/${RUN_ID}
OUTPUT_DIR=/results_gefa/${RUN_ID}
SCRATCH_DIR=/work/pinglese/scratch/${RUN_ID}

for dd in $(find ${RAW_DIR} -name *_R1_001.fastq.gz -type f); do
	fq=$(basename ${dd})
	sample=$(echo ${fq} | sed -e 's/_R1_001.fastq.gz//g')
        
	if [[ -d ${OUTPUT_DIR}/${sample} ]]; then
		echo "Skipping ${sample}"
		continue
	fi
	
	echo "Processing sample ${sample}"
	if [[ -f rync_files.txt ]]; then
		rm rsync_files.txt
	fi
	touch rsync_files.txt
	echo ${fq} > rsync_files.txt
	echo ${fq//_R1_001/_R2_001} >> rsync_files.txt

	mkdir -p ${SCRATCH_DIR}/${sample}/fastq

	rsync -ahP --files-from=rsync_files.txt --no-relative ${RAW_DIR}/ ${SCRATCH_DIR}/${sample}/fastq

	job_id=$(qsub -v RUN_ID=${RUN_ID},SAMPLE_ID=${sample} submit_juicer_nxf.sh)

	echo "Launched job ${job_id}"

	while true; do
		job_state=$(qstat -fx ${job_id} | grep -i "job_state =" | awk '{print $3}')
		if [[ ${job_state} == 'F' ]]; then
			echo "Job ${job_id} finished"
			break
		fi
		sleep 30
	done

	job_exit_code=$(qstat -fx ${job_id} | grep -i "Exit_Status" | awk '{print $3}')
	echo "Exit code: ${job_exit_code})"
	if [[ ${job_exit_code} -ne 0 ]]; then
		echo "Job failed. Exiting"
		exit 1
	fi

	# If everything is ok, sync with results_gefa
	echo "Deleting temp dirs"
	rm -rf ${SCRATCH_DIR}/${sample}/fastq ${SCRATCH_DIR}/${sample}/work 

        cd ${SCRATCH_DIR}
	rsync -ahP ${sample} ${OUTPUT_DIR}

	if [[ $? -ne 0 ]]; then
		echo "Error with rsync"
		exit 1
	fi

	rm -rf ${SCRATCH_DIR}/${sample}

	echo "Processing complete."
done
