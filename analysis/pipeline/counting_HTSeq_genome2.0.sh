BAS_DIR='/staton/projects/peach'
for batch in {0..2}
do
for id in {1..23}
do n=$((23*$batch + $id))
sample=$(ls $BAS_DIR/7_mapping_filtered/genome2.0/*.sam | sed -n "$n"p)
base=$( basename $sample | sed 's/.out.sam//g')
echo "base $base"
echo "sample $sample"

/staton/software/htseq-count \
--format=sam \
--order=name \
--stranded=no \
--type=gene \
--idattr=ID \
$sample \
/staton/projects/apricot/peach_genome/Ppersica_298_v2.1.gene_exons.gff3 \
>& $BAS_DIR/9_HTseq/genome2.0/$base.counts.txt \
2> $BAS_DIR/9_HTseq/genome2.0/$base.out &

done

wait
done