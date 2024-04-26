#! Siwei 09 Dec 2021

for eachfile in *.bed
do
	cat $eachfile \
		| grep -E '^chr[0-9]{1,2}' \
		| grep -v 'random' \
		| sed 's/chr//g' \
		> ${eachfile/%.bed/_ENSEMBL_chr.bed}
done
