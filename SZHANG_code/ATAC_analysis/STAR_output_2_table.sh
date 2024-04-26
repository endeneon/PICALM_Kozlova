#!/bin/bash
# Siwei 23 Mar 2021
# convert STAR ReadsPerGene.out to table format

# use inter_df.txt to store the intermediate df
# use current_col.txt to store the current column
# note paste use \t as the default delimiter
# use file_names.txt to store all file names, remove the ReadsPerGene.out.tab part

# init
rm output/inter_df.txt
rm output/file_names.txt
rm output/ReadsPerGene_STAR.txt

# collect .tab files only
# need to add a \n to the end of file_names.txt later or upon assembling the final output

for eachfile in *.tab
do
	echo $eachfile
	if [ ! -f output/inter_df.txt ]; then
		cat $eachfile | awk -F '\t' '{print $4}' > output/inter_df.txt
	else
#		cat $eachfile | awk -F '\t' '{print $4}' > output/current_col.txt
#		paste output/current_col.txt output/inter_df.txt > output/inter_df_1.txt
		paste output/inter_df.txt <(cat $eachfile | awk -F '\t' '{print $4}') > output/inter_df_1.txt

#		paste <(cat $eachfile | awk -F '\t' '{print $4}') output/inter_df.txt > output/inter_df_1.txt
		mv output/inter_df_1.txt output/inter_df.txt
#		echo -n $eachfile | sed 's/ReadsPerGene\.out\.tab/\t/g' >> output/file_names.txt
	fi
	echo -n -e "\t$eachfile" | sed 's/ReadsPerGene\.out\.tab//g' >> output/file_names.txt

done

# conc. file header and intermediate df as file body
# skip the first 4 lines of inter_df.txt since it does not contain gene counts
# tail -n +x starts printing at line number x

cat <(echo -n -e "Geneid") \
	output/file_names.txt \
	<(echo) \
	<(paste <(cat $eachfile | awk -F '\t' '{print $1}' | tail -n +5) \
		<(cat output/inter_df.txt | tail -n +5)) \
	> output/ReadsPerGene_STAR.txt

# cat <(echo -n "Geneid\t") \
#        output/file_names.txt \
#        <(echo -n "\n") \
#        <(paste <(cat $eachfile | awk -F '\t' '{print $1}') <(cat output/inter_df.txt | tail -n +5)) \
#        > output/ReadsPerGene_STAR.txt



# cleanup on-the-fly files

rm output/file_names.txt
rm output/inter_df.txt


