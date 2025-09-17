
for folderNumber in {640..663}; do
	s3Path=s3://kbzcutrun/

	value=$((folderNumber - 559))

	s3Folder=$(aws s3 ls "$s3Path" | awk '{print $2}' | grep "${folderNumber}")
	s3NestFolder=$(aws s3 ls "${s3Path}${s3Folder}" | awk '{print $2}' | egrep "_730?${value}_L001")	

	aws s3 sync "${s3Path}${s3Folder}${s3NestFolder}" "${s3Path}"
done
