# define base path to s3 bucket

s3Path=$1

for folderNumber in {81..104}; do

	# define local paths
	localData=~/data/${folderNumber}
	qcOutput=~/fastqFileQC/${folderNumber}

	# make destination directory for fastqc files
	mkdir -p "$qcOutput"
	mkdir -p "$localData"

	#upload data from s3 bucket
	s3Folder=$(aws s3 ls "${s3Path}" | awk '{print $2}' | grep "730${folderNumber}")
	s3Folder=${s3Folder%/}
	fileList=$(aws s3 ls "${s3Path}${s3Folder}/" | awk '{print $4}' | grep -v '^$')

	for file in $fileList; do
		echo "Downloading from ${s3Path}${folderNumber}/ to ${localData}"
		aws s3 cp "${s3Path}${s3Folder}/${file}" "$localData/${file}"

		echo "Running FastQC on ${file}..."

		fastqc -o "$qcOutput" -f fastq "$localData"/*_R1_*.fastq.gz
		fastqc -o "$qcOutput" -f fastq "$localData"/*_R2_*.fastq.gz

		# Delete local data
		echo "Removing folder ${folderNumber}"
	done

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
