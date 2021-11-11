#!/bin/bash
 
#parallel_blat.sh target query out.psl 8 "parameters"
 
cpu=$4
parameters=$( echo "-noHead "$5 )
query=$(realpath $2)
target=$(realpath $1)

blat=~/Lab/Assembly_tools/Tools/blat/blat/blat

mkdir tmp_p_b_${3}
cd tmp_p_b_${3}
cat $query | awk '{if ($1~"^>") {print $0"\t"} else {print $0} }' | tr -d "\n"|sed 's/>/\n/g'| sed '/^$/d' |split -a 10 -d -e -l 1 --additional-suffix ".fa" - seq
for seq in seq*fa; do sed -i 's:^:>:;s:\t:\n:' $seq; done
 
parallel -j $cpu $blat $parameters $target {} {}.psl ::: seq*fa

# if blat is run with header
#sed -i '1,5d' *.psl

# Header
#echo "psLayout version 3" >> ../$3
#echo "" >> ../$3
#echo "match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts" >> ../$3
#echo "     	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count" >> ../$3
#echo "---------------------------------------------------------------------------------------------------------------------------------------------------------------" >> ../$3

cat *.psl >> ../$3
 
#pigz -9p $cpu *

cd -
#rm -rf tmp_p_b_${3}
