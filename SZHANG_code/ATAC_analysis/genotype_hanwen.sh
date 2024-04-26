#!/bin/sh


#Siwei 14 Feb 2018

#Genotype rs1748456 and rs10933 for Hanwen
#Genotype rs78710909, etc.


date > genotyping.txt

echo "All genomic coordinates are based on GRCh38p7." >> genotyping.txt

for EACHFILE in *.bam
do
	echo $EACHFILE >> genotyping.txt
	echo $EACHFILE
#	samtools index -@ 20 $EACHFILE
	#printf "rs78710909\t" >> 23Mar2018_genotyping.txt
	#samtools mpileup -r chr2:127107346-127107346 $EACHFILE >> 23Mar2018_genotyping.txt
	#printf "rs76516995\t" >> 23Mar2018_genotyping.txt
	#samtools mpileup -r chr2:127107345-127107345 $EACHFILE >> 23Mar2018_genotyping.txt
        #printf "rs1532276\t" >> 23Mar2018_genotyping.txt
        #samtools mpileup -r chr8:27608640-27608640 $EACHFILE >> 23Mar2018_genotyping.txt
        #printf "rs1532277\t" >> 23Mar2018_genotyping.txt
        #samtools mpileup -r chr8:27608664-27608664 $EACHFILE >> 23Mar2018_genotyping.txt
        #printf "rs1532278\t" >> 23Mar2018_genotyping.txt
        #samtools mpileup -r chr8:27608798-27608798 $EACHFILE >> 23Mar2018_genotyping.txt
#       printf "rs2027349\t" >> genotyping.txt
#        samtools mpileup -r chr1:150067621-150067621 $EACHFILE >> genotyping.txt

        printf "rs7148456\t" >> genotyping.txt
        samtools mpileup -r chr14:103561933-103561933 $EACHFILE >> genotyping.txt

        printf "rs10933\t" >> genotyping.txt
        samtools mpileup -r chr3:52685800-52685800 $EACHFILE >> genotyping.txt

        printf "rs56205728\t" >> genotyping.txt
        samtools mpileup -r chr15:40275036-40275036 $EACHFILE >> genotyping.txt



done


