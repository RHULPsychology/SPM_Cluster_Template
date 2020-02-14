#!/bin/sh


# Script to take list of subject names and performs 
# first level processing on Target Task files


### Select the proper version of matlab

source /usr/local/apps/psycapps/config/matlab_bash

# Make sure to change this according to your account and that this folder exists
OUTDIR=/MRIWork/MRIWork03/nf/forida_begum/logs
mkdir -p $OUTDIR

for subject_number in {2..3}
do
    echo "Processing Subject $subject_number"
    qsub    -l h_cpu=00:02:00,h_rss=8G \
            -o ${OUTDIR}/matlab_${subject_number}.out \
            -e ${OUTDIR}/matlab_${subject_number}.err \
            /MRIWork/MRI-Scratch/nick/call_matlab_function.sh $subject_number;       
done
