#!/bin/sh


#Siwei 14 Feb 2018



date > 23Jun2020_genotyping.txt

echo "All genomic coordinates are based on GRCh38p7." >> 23Jun2020_genotyping.txt

for EACHFILE in *WASPed.bam
do
	#samtools index -@ 24 $EACHFILE
	#samtools index -@ 24 $EACHFILE
	echo $EACHFILE
	#samtools index -@ 24 $EACHFILE
	echo -n $EACHFILE >> 23Jun2020_genotyping.txt
	#echo $EACHFILE
	#samtools index -@ 24 $EACHFILE
	printf "\trs72986630\t" >> 23Jun2020_genotyping.txt
	samtools mpileup -r chr19:11849736-11849736 $EACHFILE >> 23Jun2020_genotyping.txt
	printf "\n" >> 23Jun2020_genotyping.txt

	echo -n $EACHFILE >> 23Jun2020_genotyping.txt
	printf "\trs10933\t" >> 23Jun2020_genotyping.txt
	samtools mpileup -r chr3:52719816-52719816 $EACHFILE >> 23Jun2020_genotyping.txt
	printf "\n" >> 23Jun2020_genotyping.txt

	echo -n $EACHFILE >> 23Jun2020_genotyping.txt
        printf "\trs7148456\t" >> 23Jun2020_genotyping.txt
        samtools mpileup -r chr14:104028270-104028270 $EACHFILE >> 23Jun2020_genotyping.txt
	printf "\n" >> 23Jun2020_genotyping.txt

        echo -n $EACHFILE >> 23Jun2020_genotyping.txt
        printf "\trs56205728\t" >> 23Jun2020_genotyping.txt
        samtools mpileup -r chr15:40275036-40275036 $EACHFILE >> 23Jun2020_genotyping.txt
        printf "\n" >> 23Jun2020_genotyping.txt

        #printf "rs1532276\t" >> 23Mar2018_genotyping.txt
        #samtools mpileup -r chr8:27608640-27608640 $EACHFILE >> 23Mar2018_genotyping.txt
        #printf "rs1532277\t" >> 23Mar2018_genotyping.txt
        #samtools mpileup -r chr8:27608664-27608664 $EACHFILE >> 23Mar2018_genotyping.txt
        #printf "rs1532278\t" >> 23Mar2018_genotyping.txt
        #samtools mpileup -r chr8:27608798-27608798 $EACHFILE >> 23Mar2018_genotyping.txt
        #printf "rs2027349\t" >> genotyping.txt
        #samtools mpileup -r chr1:150067621-150067621 $EACHFILE >> genotyping.txt

done


