for R1 in /staton/projects/peach/raw_data/*_R1_001.fastq.gz
do
    R2=$(echo $R1 | sed 's/R1/R2/g')
    BASE=$( basename $R1 | sed 's/_R1_001.fastq.gz//g')
    echo "R1 $R1"
    echo "R2 $R2"
    echo "BASE $BASE"
    /staton/software/skewer/skewer \
	-x /staton/software/Trimmomatic-0.35/adapters/illuminaClipping.fa \
	-l 50 \
	$R1 \
	$R2 \
	-t 15 \
	-o $BASE \
	>& $BASE.trim_output
done