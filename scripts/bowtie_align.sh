s3Path=$1
ref=$2

get_adapters() {
    local sample_num=$1
    local config_file="/home/ssm-user/scripts/adapter_config.txt"
    
    # Check if config file exists
    if [[ ! -f ${config_file} ]]; then
        echo "Error: Configuration file $config_file not found!"
        exit 1
    fi
    
    # Search for the folder number in the config file and extract adapters
    local line=$(grep "^${sample_num}," "${config_file}")
    
    if [[ -n "$line" ]]; then
        # Parse the line to extract adapter sequences
        IFS=',' read -r folder_number adapter1 adapter2 <<< "$line"
        echo "$adapter1,$adapter2"
    else
        echo "Error: No adapter configuration found for folder ${sample_num}"
        exit 1
    fi
}

for folderNumber in {640..663}; do	

	#define sample number
	value=$((folderNumber - 559))

	#define local path
	localData=/home/ssm-user/data/${folderNumber}

	mkdir -p ${localData}

	# get adapter sequences specific to this folder
	echo "Pulling adapter sequences for sample ${value}"
	adapter_info=$(get_adapters ${value})

	if [[ $? -ne 0 ]]; then
		echo "Failed to get adapters for folder ${folderNumber}. Skipping..."
		continue
	fi

	# Extract individual adapter sequences
	IFS=',' read -r adapter1 adapter2 <<< "$adapter_info"
	echo "Using adapters for folder ${folderNumber}: Adapter1=${adapter1}, Adapter2=${adapter2}"

	#download data from s3 bucket
	s3Folder=$(aws s3 ls "$s3Path" | awk '{print $2}' | grep "${folderNumber}")
	s3NestFolder=$(aws s3 ls "${s3Path}${s3Folder}" | awk '{print $2}' | egrep "_730?${value}_L001")

	echo "Downloading from "${s3Path}${s3Folder}${s3NestFolder}" to ${localData}"
	aws s3 cp "${s3Path}${s3Folder}${s3NestFolder}" "${localData}" --recursive

	echo "Trimming adapters of sample ${value} using CutAdapt."
	cutadapt -a ${adapter1} -A ${adapter2} -o ${localData}/${value}.R1.trimmed.fastq.gz -p ${localData}/${value}.R2.trimmed.fastq.gz ${localData}/*_R1_*.fastq.gz ${localData}/*_R2_*.fastq.gz --minimum-length 20 --cores 4

	echo "Alignment of ${value} to HG38"
	bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 4 -x ${ref} -1 ${localData}/${value}.R1.trimmed.fastq.gz -2 ${localData}/${value}.R2.trimmed.fastq.gz -S /home/ssm-user/alignment/sam/${value}_bowtie2.sam &> /home/ssm-user/alignment/sam/bowtie2_summary/${value}_bowtie2.txt

	#upload trimmed fastq files to s3 bucket
	echo "Uploading trimmed files to ${s3Path}${s3Folder}${s3NestFolder}"
	aws s3 cp "${localData}/*trimmed*.fastq.gz" "${s3Path}${s3Folder}${s3NestFolder}"  --recursive

	#upload alignment results to s3 bucket
	echo "Uploading SAM files to /alignment/sam in s3 bucket"
	aws s3 cp "/home/ssm-user/alignment/sam/${value}_bowtie2.sam" "${s3Path}alignment/sam/"
	echo "Uploading TXT files to /alignment/sam/bowtie2_summary in s3 bucket"
	aws s3 cp "/home/ssm-user/alignment/sam/bowtie2_summary/${value}_bowtie2.txt" "${s3Path}alignment/sam/bowtie2_summary/"

	#cleaning folders
	echo "Removing raw data"
	rm -rf ${localData} --recursive
	echo "Removing SAM and TXT files"
	rm -rf "/home/ssm-user/alignment/sam/${value}_bowtie2.sam"
	rm -rf "/home/ssm-user/alignment/sam/bowtie2_summary/${value}_bowtie2.txt"

	echo "Bowtie alignment complete for sample ${value}"
	echo "--------------------------------------"
done

echo "All samples processed successfully."
