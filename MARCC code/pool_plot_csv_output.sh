#!/bin/bash
#This script is for the purpose of directing relevant outputs to important locations
#directory_from= ~/data2/gtex_v8/diff_samps
#directory_to= ~/data2/gtex_v8/analysis_output
#mkdir $directory_to
while read -r general
do
        cp      ~/work2/tnieuwe1/data/gtex_v8/diff_samps/${general}/*.pdf ~/work2/tnieuwe1/data/gtex_v8/analysis_output/high_var_output/.
        cp      ~/work2/tnieuwe1/data/gtex_v8/diff_samps/${general}/*.csv ~/work2/tnieuwe1/data/gtex_v8/analysis_output/high_var_output/.
done < ~/work2/tnieuwe1/data/gtex_v8/gen_tiss_lists/general_list_test_v8.txt
