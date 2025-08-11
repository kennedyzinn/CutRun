s3Path=$1
ref=$2

for folderNumber in {640..663}; do	

	#define sample number
	value=$((folderNumber - 559))

	#define local path
	localData=/home/ssm-user/data/${folderNumber}

	mkdir -p ${localData}

	#download data from s3 bucket
	s3Folder=$(aws s3 ls "$s3Path" | awk '{print $2}' | grep "${folderNumber}")
	s3NestFolder=$(aws s3 ls "${s3Path}${s3Folder}" | awk '{print $2}' | egrep "_730?${value}_L001")

	echo "Downloading from "${s3Path}${s3Folder}${s3NestFolder}" to ${localData}"
	aws s3 cp "${s3Path}${s3Folder}${s3NestFolder}" "${localData}" --recursive

	echo "Alignment of ${value} to HG38"
	bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 16 -x ${ref} -1 ${localData}/*_R1_*.fastq.gz -2 ${localData}/*_R2_*.fastq.gz -S /home/ssm-user/alignment/sam/${value}_bowtie2.sam &> /home/ssm-user/alignment/sam/bowtie2_summary/${value}_bowtie2.txt

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
