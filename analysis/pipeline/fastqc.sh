for f in ../raw_data/*.fastq
do
    filename=$(basename "$f")
    base="${filename%%.fastq*}"
    echo "filename $filename base $base"
    mkdir $base.fastQC

    /staton/software/FastQC-v0.11.5/FastQC/fastqc -t 28 -o $base.fastQC $f
done