## Scripts for RNASeq analysis
`Pipeline`      
The following scripts are used:
- `fastqc.sh`: quality control    
- `Trim_skewer.sh`: trim adaptors and low quality reads       
- `mapping_STAR_genome1.0.sh`: align reads to peach genome v1.0      
- `mapping_STAR_genome2.0.sh`: align reads to peach genome v2.0       
- `counting_HTSeq_genome1.0.sh`: count reads aligned to genome v1.0      
- `counting_HTSeq_genome2.0.sh`: count reads aligned to genome v2.0      
- `rRNA_filter_Bowtie2.sh`: filter out the rRNA reads      
- `merge_counts.sh`: create a table of read-counts for all samples 

`DEanalysis`     
R scripts for differential expression analyses. The following scripts are:
- `Apricot_DE.R`: DESeq2 scripts on apricot data. 
- `Peach_DEseq.R`: DESeq2 scritps on peach data.
- `compare_apri_peachDEG.R`: Compare apricot DEGs and peach DEGs, and significant GO terms.
- `drawplots.R`: plot gene expression profiles by ggplot2 using the normalized counts.           
- `circlize_circosplot.R`: generate ciros plot.
- `extractQTL_gene.R`: extract genes on the QTL region.

`network_analysis`    
R scripts for co-expression analyses. The following scripts are:
- `Apricot_WGCNA.R`: WGCNA scripts on apricot data.
- `Peach_WGCNA.R`: WGCNA scripts on peach data.
- 'combiWGCNA.R`: WGCNA scritps integrating apricot and peach data.
- 'coNetworkGO.R': scripts for plotting GO enrichment results.

     
