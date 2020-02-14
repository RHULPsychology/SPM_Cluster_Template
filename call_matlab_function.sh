#!/bin/sh

subject_number=$1
echo $subject_number
script_folder="/MRIWork/MRI-Scratch/nick"
cd $script_folder
matlab -nodisplay -nosplash -nodesktop -r "Cluster_Template_cars_rsa_v8_2019_rstab('$subject_number');exit;"
