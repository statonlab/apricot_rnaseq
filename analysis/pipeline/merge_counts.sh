files=$(ls *.counts.txt | tr '\n' '\t')
echo -e "gene_id\t" "$files" | sed 's/.counts.txt//g' > gene_counts.txt

awk '{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}' <( for f in $files; do grep -v processed $f | grep -v Warning | head -n -5; done ) | sort \
>> gene_counts.txt