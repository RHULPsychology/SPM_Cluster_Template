#!/bin/sh

subject_number=$1
echo $subject_number

### in this case the script folder is same as pwd but in other occasions change the path to your script directory
script_folder=$(pwd)
cd $script_folder
matlab -nodisplay -nosplash -nodesktop -r "Cluster_Template('$subject_number');exit;"
