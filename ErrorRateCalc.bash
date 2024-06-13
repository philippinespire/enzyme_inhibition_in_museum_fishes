#!/bin/bash
#SBATCH --job-name=ErrorRateCalc
#SBATCH --time=96:00:00
#SBATCH -p normal

# script to query error rates in newRAD data French Masters
# dir is expected to be mkBAM
# put fastq files into dir called mkBAM/fastq

# give SubsetData value at command line, which is the number of reads to process

# to run

# ErrorRateCalc_10_11_ceb.bash <subset value> <number threads>

# if subset value is omitted or zero, then no subset will be applied

# User defined variables and packages
AlbatrossPATTERN="*-A-*"
ContempPATTERN="*-C-*"
FqExt=.fq
FqPrefix=2019

module load parallel
module load bio-misc
module load R/gcc7/3.6.1

# Read in values from command line
SubsetReads=100000
THREADS=20

if [ -z $SubsetReads ]; then
	echo; echo `date` Processing all reads in FASTQ files
        SubsetReads=0  #num reads to subset, if 0 then no subset
else
	echo; echo `date` Processing $SubsetReads reads per FASTQ file
        #SubsetReads=$1  #num reads to subset, if 0 then no subset
	echo; echo `date` Deleting Sub*fq files from previous run
	rm Sub*fq
fi

if [ -z $THREADS ]; then
	echo; echo `date` Processing data with 1 thread
	THREADS=1
else
	echo; echo `date` Processing data with $THREADS threads
        #THREADS=$2
fi

# Gunzip everything in mkBAM
ls *$FqExt.gz 2> /dev/null | parallel --no-notice -j $THREADS "gunzip {}"

ErrorRateCalc () {
	fqBase=$1
	for barcode in CGATGCTCTGCA AAGCCGGTTGCA; do
		for read in R1 R2; do		
        		agrep -1 -D100 -I100 "^${barcode}" ${fqBase}.${read}.fq | head -n -1 | sed 's/^/2del,--/g' > ${fqBase}.${read}.${barcode}.NoDels.seq
        		agrep -1 -D100 -I100 "^.${barcode}" ${fqBase}.${read}.fq | head -n -1 | sed 's/^/1del,-/g' >> ${fqBase}.${read}.${barcode}.NoDels.seq
	        	agrep -1 -D100 -I100 "^..${barcode}" ${fqBase}.${read}.fq | head -n -1 | sed 's/^/0Ind,/g' >> ${fqBase}.${read}.${barcode}.NoDels.seq
        		agrep -1 -D100 -I100 "^...${barcode}" ${fqBase}.${read}.fq | head -n -1 | cut -c 2- | sed 's/^/1Ins,/g' >> ${fqBase}.${read}.${barcode}.NoDels.seq
        		agrep -1 -D100 -I100 "^....${barcode}" ${fqBase}.${read}.fq | head -n -1 | cut -c 3- | sed 's/^/2Ins,/g' >> ${fqBase}.${read}.${barcode}.NoDels.seq
	        	agrep -1 -D100 -I100 "^.....${barcode}" ${fqBase}.${read}.fq | head -n -1 | cut -c 4- | sed 's/^/3Ins,/g' >> ${fqBase}.${read}.${barcode}.NoDels.seq
        		agrep -1 -D100 -I100 "^......${barcode}" ${fqBase}.${read}.fq | head -n -1 | cut -c 5- | sed 's/^/4Ins,/g' >> ${fqBase}.${read}.${barcode}.NoDels.seq
        		agrep -1 -D100 -I100 "^.......${barcode}" ${fqBase}.${read}.fq | head -n -1 | cut -c 6- | sed 's/^/5Ins,/g' >> ${fqBase}.${read}.${barcode}.NoDels.seq
	        	agrep -1 -D100 -I100 "^........${barcode}" ${fqBase}.${read}.fq | head -n -1 | cut -c 7- | sed 's/^/6Ins,/g' >> ${fqBase}.${read}.${barcode}.NoDels.seq
        		agrep -1 -D100 -I100 "^.........${barcode}" ${fqBase}.${read}.fq | head -n -1 | cut -c 8- | sed 's/^/7Ins,/g' >> ${fqBase}.${read}.${barcode}.NoDels.seq
        		agrep -1 -D100 -I100 "^..........${barcode}" ${fqBase}.${read}.fq | head -n -1 | cut -c 9- | sed 's/^/8Ins,/g' >> ${fqBase}.${read}.${barcode}.NoDels.seq

		        cut -c 1-21 ${fqBase}.${read}.${barcode}.NoDels.seq | sed -e 's/\([GATCN-]\)/\1,/g' | sed -e 's/,$//g' -e "s/^/$fqBase,$barcode,$read,/g" > ${fqBase}_${barcode}_${read}.csv #save to file
        	done
	done
}
export -f ErrorRateCalc

#run function in parallel
if [ $SubsetReads > 0 ]; then
        SubsetData=$(($SubsetReads * 4))
	echo; echo `date` Retrieving $SubsetData lines from each FASTQ file
	ls ${fqPrefix}*$FqExt | parallel --no-notice -j $THREADS "head -n $SubsetData {} > Sub-$SubsetReads-{}"
	echo; echo `date` Begin processing data
	ls Sub-$SubsetReads*R1$FqExt | sed "s/\.R1$FqExt//g" | parallel -j $THREADS --no-notice ErrorRateCalc {}
else
	echo; echo `date` Begin processing data
	ls ${fqPrefix}*R1$FqExt | grep "^2019" | sed "s/\.R1$FqExt//g" | parallel -j $THREADS --no-notice ErrorRateCalc {}
fi

echo; echo `date` Gathering the following $(ls ${AlbatrossPATTERN}[1-2].csv | wc -l) Albatross files for concatenation:
ls ${AlbatrossPATTERN}[1-2].csv

Albatross_Barcodes=$(ls ${AlbatrossPATTERN}[1-2].csv)

echo; echo `date` Gathering the following $(ls ${ContempPATTERN}[1-2].csv | wc -l) Contemporary files for concatenation:
ls ${ContempPATTERN}[1-2].csv

Contemp_Barcodes=$(ls ${ContempPATTERN}[1-2].csv)

Header=fqBase,Barcode,Read,Indels,BP1,BP2,BP3,BP4,BP5,BP6,BP7,BP8,BP9,BP10,BP11,BP12,BP13,BP14,BP15,BP16

cat <(echo $Header) $Albatross_Barcodes > AlbatrossBarcodes.csv
cat <(echo $Header) $Contemp_Barcodes > ContempBarcodes.csv

#echo; echo `date` Calculate statistics in R
R CMD BATCH PIRE_Stats_forR_HPC_9_30.R

echo; echo `date` Script completed!
