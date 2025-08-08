# define base path to s3 bucket

s3Path=$1

for folderNumber in {640..663}; do
	
	# define file identifier
	value=$((folderNumber - 559))

	# define local paths
	localData=./data/${folderNumber}
	qcOutput=./fastqFileQC/${folderNumber}

	# make destination directory for fastqc files
	mkdir -p "$qcOutput"
	mkdir -p "$localData"

	#download data from s3 bucket
	s3Folder=$(aws s3 ls "$s3Path" | awk '{print $2}' | grep "${folderNumber}")
	s3NestFolder=$(aws s3 ls "${s3Path}${s3Folder}" | awk '{print $2}' | egrep "_730?${value}_L001")

	echo "Downloading from ${s3Path}${s3Folder}${s3NestFolder} to ${localData}"
	aws s3 cp "${s3Path}${s3Folder}${s3NestFolder}" "${localData}" --recursive

	echo "Running FastQC on ${value}..."

	fastqc -o "$qcOutput" -f fastq ${localData}/*_R1_*.fastq.gz
	fastqc -o "$qcOutput" -f fastq ${localData}/*_R2_*.fastq.gz

	# Delete local data
	echo "Removing folder ${folderNumber}"

	# Upload fastQC results to S3 bucket
	echo "Uploading FastQC results to ${s3Path}fastqc/${folderNumber}/"
	aws s3 cp "$qcOutput" "${s3Path}fastqc/${folderNumber}/" --recursive

	# Clean up folder ${folderNumber}
	rm -rf "$localData"
	rm -rf "$qcOutput"

	echo "FastQC processing for folder ${folderNumber} complete."
	echo "------------------------------------"
done

echo "All folders processed successfully."
