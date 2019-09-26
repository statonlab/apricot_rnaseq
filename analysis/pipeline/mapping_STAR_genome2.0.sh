for f1 in filtered_data/*.1.fq
do
    f2=$(echo $f1 | sed 's/.1.fq/.2.fq/g')
    prefix=$( basename $f1 | sed 's/-trimmed.not_aligned.1.fq//g')
    echo "f1 $f1"
    echo "f2 $f2"
    echo "base $prefix"

    /staton/software/STAR-2.5.3a/STAR \
	--genomeDir /staton/projects/apricot/peach_genome/Peach2.0_STAR/ \
	--readFilesIn $f1 $f2 \
	--runThreadN 25 \
	--outFileNamePrefix /staton/projects/peach/7_mapping_filtered/$prefix.
done