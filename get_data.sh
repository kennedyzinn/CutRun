# First, see what samples/datasets you have
PROJECT_ID=446978745

# Download one sample at a time
for sample_id in $(bs list biosamples --project-id=${PROJECT_ID} --format=json | jq -r '.[].Id'); do
    echo "Processing sample: $sample_id"
    
    # Create temp directory
    mkdir -p $HOME/basespace_data
    
    # Download this sample only
    bs download biosample --id=$sample_id --output=$HOME/basespace_data/
    
    # Check size
    du -sh $HOME/basespace_data/
    
    # Upload to S3 with checksum validation
    aws s3 sync $HOME/basespace_data/ s3://kbzcutrun/sample_$sample_id/ --checksum-algorithm SHA256
    
    # Clean up immediately
    rm -rf $HOME/basespace_data
    
    echo "Completed sample: $sample_id"
done
