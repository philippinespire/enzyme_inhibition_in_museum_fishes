#!/bin/bash
#SBATCH --job-name=SeqCnt
#SBATCH --time=00:00:00
#SBATCH -p cbirdq

# User Defined Variables
THREADS=10
FW=R1
RV=R2
#FWRV=FR

# Print names of files processed to output file
ls *[$FW$RV].fq 

module load parallel
module load bio-misc

rm *.Temp *.NotInIndelCategory SequenceCounts.tsv

# Define Function to Track RAM usage
reportUsage(){
RAM=$(free -m | awk 'NR==2{printf "%.0f%%\t", $3*100/$2 }')
SWAP=$(free -m | awk 'NR==4{printf "%.0f%%\t", $3*100/$2 }')
DISK=$(df -h | awk '$NF=="/"{printf "%s\t", $5}')
CPU=$(top -bn1 | grep load | awk '{printf "%.0f%%\t\n", $(NF-2)}')
DATETIME=`date`
echo "   $RAM$SWAP$DISK$CPU$DATETIME"
}
export -f reportUsage


# Define Main Function
IndelDiag () {
fqBase=$1

# rm -rf $fqBase
# mkdir $fqBase

echo; echo `date` $fqBase
# Report memory usage header
printf "   RAM\tSwapRAM\tDisk\tCPU\tDateTime\n"
reportUsage

# Count the number of sequences in each file and put into "SequenceCounts.tsv" file.
paste <(echo Total All ${fqBase}.fq | tr " " "\t") <(grep "^\+$" ${fqBase}.fq | wc -l) >> SequenceCounts.tsv
# Grab Sequence and put into "${fqBase}.Seqs.Temp" file.
awk 'NR % 4 == 2' ${fqBase}.fq > ${fqBase}.Seqs.Temp
# Grab Sequence IDs and put into "${fqBase}.IDs.Temp" file.
grep "\@" ${fqBase}.fq | tr " " "\t" | cut -f1 > ${fqBase}.IDs.Temp
# Stick Sequence ID and Sequence together then change tabs to commas
paste <(cat ${fqBase}.Seqs.Temp) <(cat ${fqBase}.IDs.Temp) > ${fqBase}.SeqIDs.Temp

# Search each "${fqBAse}.SeqIDs.Temp" file for each barcode, remove the last line of text "Grand Total..." and sort into indel categories 2Del-8Ins.
# Echo the indel category, the barcode, and the "${fqBase}.fq" file name. Change spaces to tabs. Count the number of sequences in each indel category and put into "SequenceCounts.tsv" file.

for PATTERN in $(cat barcodes.txt); do

agrep -1 -D100 -I100 "^$PATTERN" ${fqBase}.SeqIDs.Temp | head -n -1 > ${fqBase}.$PATTERN.2Del.Temp
paste <(echo 2Del $PATTERN ${fqBase}.fq | tr " " "\t") <(cat ${fqBase}.$PATTERN.2Del.Temp | wc -l | awk '{print $1;}') >> SequenceCounts.tsv
agrep -1 -D100 -I100 "^.$PATTERN" ${fqBase}.SeqIDs.Temp | head -n -1 > ${fqBase}.$PATTERN.1Del.Temp
paste <(echo 1Del $PATTERN ${fqBase}.fq | tr " " "\t") <(cat ${fqBase}.$PATTERN.1Del.Temp | wc -l | awk '{print $1;}') >> SequenceCounts.tsv
agrep -1 -D100 -I100 "^..$PATTERN" ${fqBase}.SeqIDs.Temp | head -n -1 > ${fqBase}.$PATTERN.NoInd.Temp
paste <(echo NoInd $PATTERN ${fqBase}.fq | tr " " "\t") <(cat ${fqBase}.$PATTERN.NoInd.Temp | wc -l | awk '{print $1;}') >> SequenceCounts.tsv
agrep -1 -D100 -I100 "^...$PATTERN" ${fqBase}.SeqIDs.Temp | head -n -1 > ${fqBase}.$PATTERN.1Ins.Temp
paste <(echo 1Ins $PATTERN ${fqBase}.fq | tr " " "\t") <(cat ${fqBase}.$PATTERN.1Ins.Temp | wc -l | awk '{print $1;}') >> SequenceCounts.tsv
agrep -1 -D100 -I100 "^....$PATTERN" ${fqBase}.SeqIDs.Temp | head -n -1 > ${fqBase}.$PATTERN.2Ins.Temp
paste <(echo 2Ins $PATTERN ${fqBase}.fq | tr " " "\t") <(cat ${fqBase}.$PATTERN.2Ins.Temp | wc -l | awk '{print $1;}') >> SequenceCounts.tsv
agrep -1 -D100 -I100 "^.....$PATTERN" ${fqBase}.SeqIDs.Temp | head -n -1 > ${fqBase}.$PATTERN.3Ins.Temp
paste <(echo 3Ins $PATTERN ${fqBase}.fq | tr " " "\t") <(cat ${fqBase}.$PATTERN.3Ins.Temp | wc -l | awk '{print $1;}') >> SequenceCounts.tsv
agrep -1 -D100 -I100 "^......$PATTERN" ${fqBase}.SeqIDs.Temp | head -n -1 > ${fqBase}.$PATTERN.4Ins.Temp
paste <(echo 4Ins $PATTERN ${fqBase}.fq | tr " " "\t") <(cat ${fqBase}.$PATTERN.4Ins.Temp | wc -l | awk '{print $1;}') >> SequenceCounts.tsv
agrep -1 -D100 -I100 "^.......$PATTERN" ${fqBase}.SeqIDs.Temp | head -n -1 > ${fqBase}.$PATTERN.5Ins.Temp
paste <(echo 5Ins $PATTERN ${fqBase}.fq | tr " " "\t") <(cat ${fqBase}.$PATTERN.5Ins.Temp | wc -l | awk '{print $1;}') >> SequenceCounts.tsv
agrep -1 -D100 -I100 "^........$PATTERN" ${fqBase}.SeqIDs.Temp | head -n -1 > ${fqBase}.$PATTERN.6Ins.Temp
paste <(echo 6Ins $PATTERN ${fqBase}.fq | tr " " "\t") <(cat ${fqBase}.$PATTERN.6Ins.Temp | wc -l | awk '{print $1;}') >> SequenceCounts.tsv
agrep -1 -D100 -I100 "^.........$PATTERN" ${fqBase}.SeqIDs.Temp | head -n -1 > ${fqBase}.$PATTERN.7Ins.Temp
paste <(echo 7Ins $PATTERN ${fqBase}.fq | tr " " "\t") <(cat ${fqBase}.$PATTERN.7Ins.Temp | wc -l | awk '{print $1;}') >> SequenceCounts.tsv
agrep -1 -D100 -I100 "^..........$PATTERN" ${fqBase}.SeqIDs.Temp | head -n -1 > ${fqBase}.$PATTERN.8Ins.Temp
paste <(echo 8Ins $PATTERN ${fqBase}.fq | tr " " "\t") <(cat ${fqBase}.$PATTERN.8Ins.Temp | wc -l | awk '{print $1;}') >> SequenceCounts.tsv

done

# Report memory usage
reportUsage

# Stick all of the separate indel category files together keeping just the indel IDs and put them into "${fqBase}.IndelCatIDs.Temp" files.
cat <(cat ${fqBase}.CGATGCTCTGCA.*.Temp) \
	<(cat ${fqBase}.AAGCCGGTTGCA.*.Temp) \
	<(cat ${fqBase}.AAGCGACCTGCA.*.Temp) \
	<(cat ${fqBase}.AGTGGTGCTGCA.*.Temp) \
	<(cat ${fqBase}.CAGGCTTGTGCA.*.Temp) \
	<(cat ${fqBase}.CCAGCGCGTGCA.*.Temp) \
	<(cat ${fqBase}.CCGTCAGATGCA.*.Temp) \
	<(cat ${fqBase}.CGGCACCATGCA.*.Temp) \
	<(cat ${fqBase}.CGGTCGCCTGCA.*.Temp) \
	<(cat ${fqBase}.CGTCTCGTTGCA.*.Temp) \
	<(cat ${fqBase}.CTCGCCGATGCA.*.Temp) \
	<(cat ${fqBase}.CTGAGGCATGCA.*.Temp) \
	<(cat ${fqBase}.CTGCGTAGTGCA.*.Temp) \
	<(cat ${fqBase}.GACTGCCATGCA.*.Temp) \
	<(cat ${fqBase}.GCACCATGTGCA.*.Temp) \
	<(cat ${fqBase}.GCGTTCTCTGCA.*.Temp) \
	<(cat ${fqBase}.GCTAGTCGTGCA.*.Temp) \
	<(cat ${fqBase}.GCTCGCTATGCA.*.Temp) \
	<(cat ${fqBase}.GCTGAGACTGCA.*.Temp) \
	<(cat ${fqBase}.GGCGCTCATGCA.*.Temp) \
	<(cat ${fqBase}.GTCGGCATTGCA.*.Temp) \
	<(cat ${fqBase}.TACACCGGTGCA.*.Temp) \
	<(cat ${fqBase}.TCCGTGAGTGCA.*.Temp) \
	<(cat ${fqBase}.TGACTGCGTGCA.*.Temp) | cut -f2 > ${fqBase}.IndelCatIDs.Temp


# Report memory usage
reportUsage

paste <(echo IndelCat All ${fqBase}.fq | tr " " "\t") <(wc -l ${fqBase}.IndelCatIDs.Temp | awk '{print $1;}') >> SequenceCounts.tsv

# sed -i 's/\t/,/g' SequenceCounts.tsv

# Report memory usage
reportUsage

cat ${fqBase}.fq | paste - - - - | grep -vFf ${fqBase}.IndelCatIDs.Temp | tr "\t" "\n" > ${fqBase}.NotInIndelCategory.Temp

# Report memory usage
reportUsage
}

export -f IndelDiag


echo; echo `date` "---------------Running IndelDiag()---------------"
printf "   RAM\tSwapRAM\tDisk\tCPU\tDateTime\n"
reportUsage
ls *[$FW$RV].fq | sed "s/.fq//g" | parallel --no-notice -k -j ${THREADS} IndelDiag {}
reportUsage

# echo; echo `date` "---------------Creating All.IndelCatIDs.Temp---------------"
# printf "   RAM\tSwapRAM\tDisk\tCPU\tDateTime\n"
# reportUsage
# cat *R[12].IndelCatIDs.Temp > All.IndelCatIDs.Temp 
# reportUsage

echo; echo `date` "---------------Creating All.NotInIndelCategory---------------"
# will return double counts
printf "   RAM\tSwapRAM\tDisk\tCPU\tDateTime\n"
reportUsage
# ls *.R[12].fq | parallel --no-notice -kj $THREADS "cat {} | paste - - - - | grep -vFf All.IndelCatIDs.Temp | tr '\t' '\n'" > All.NotInIndelCategory.fq
ls *[$FW$RV].fq | sed 's/\.fq//' | parallel --no-notice -j $THREADS "cat {}.fq | paste - - - - | grep -vFf {}.IndelCatIDs.Temp | tr '\t' '\n' > {}.NonIndelCat.fq"

reportUsage
