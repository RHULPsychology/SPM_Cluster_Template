#!/bin/sh


# Script to take list of subject names and performs 
# first level processing on Target Task files


### Select the proper version of matlab

source /usr/local/apps/psycapps/config/matlab_bash

# Make sure to change this according to your account and that this folder exists
OUTPUT_LOG_DIR=$(pwd)/log
mkdir -p $OUTPUT_LOG_DIR

script_folder=$(pwd)

for subject_number in {2..3}
do
    echo "Processing Subject $subject_number"
    qsub    -l h_rss=4G \
            -o ${OUTPUT_LOG_DIR}/matlab_${subject_number}.out \
            -e ${OUTPUT_LOG_DIR}/matlab_${subject_number}.err \
            $script_folder/call_matlab_function.sh $subject_number;       
done
