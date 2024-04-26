for EACHFILE in *.bed
do
	echo $EACHFILE
	wc -l $EACHFILE
done

