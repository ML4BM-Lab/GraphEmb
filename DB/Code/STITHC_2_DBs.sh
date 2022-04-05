# !/bin/bash

for db_var in 'ATC' 'BindingDB' 'ChEBI' 'ChEMBL' 'KEGG' 'PS' 'PC'
do 
	echo "Processing $db_var"
	grep "$db_var" ./../Data/cross_side_information_DB/STITCH/chemical.sources.v5.0.tsv > ./../Data/cross_side_information_DB/STITCH/dict_STITCH_$db_var.tsv &
done