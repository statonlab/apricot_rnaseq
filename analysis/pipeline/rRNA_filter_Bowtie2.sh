for f1 in 1_skewer/*-pair1.fastq
do
    f2=$(echo $f1 | sed 's/pair1/pair2/g')
    prefix=$( basename $f1 | sed 's/-pair1.fastq//g')
    echo "f1 $f1"
    echo "f2 $f2"
    echo "base $prefix"

    /staton/software/bowtie2-2.2.7/bowtie2 \
        -p 28 \
        -x ./rRNA/rDNA \
        -1 $f1 \
        -2 $f2 \
        -S 6_Aligned_rRNA/$prefix-bowtie_out.sam \
	--un-conc 6_Aligned_rRNA/$prefix.not_aligned.fq \
        >& 6_Aligned_rRNA/$prefix.bowtie2.output
done