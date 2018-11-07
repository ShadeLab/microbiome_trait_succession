#Mouse project script: making a phylogenetic tree with samples

# we will use living tree project -- very well annotated and small (which is a feature)
wget "https://www.arb-silva.de/fileadmin/silva_databases/living_tree/LTP_release_132/LTPs132_datasets.fasta.tar.gz"
tar -xvf LTPs132_datasets.fasta.tar.gz

#standardize the format with output from the usearch pipeline
#remove blank spaces
sed '/^[^>]/ 's/[[:space:]]//g'' LTPs132_SSU_compressed.fasta > LTP.fasta

#remove line breaks...
sed -i ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' LTP.fasta

#now convert rna to dna (even though it doesn't matter for cluster -- i checked)
sed -i '/^[^>]/ y/uU/tT/' LTP.fasta

#now pick out just id from the header
sed -i 's/\t.*//' LTP.fasta

#now use cutadapt to select just the v4 region
#i got these sequences from the GLBRC pipeline: https://github.com/GLBRC-TeamMicrobiome/TagAnalysis/blob/master/16S_USEARCH_Pipeline.md
#note that we are using the first one in the 3' direction, and the second in the 5' direction
module load cutadapt/1.8.1
cutadapt -a ATTAGAWACCCBDGTAGTCC -o tmp.fasta LTP132.fasta > cut_Fadpt_results_LTP.txt
cutadapt -g GTGCCAGCMGCCGCGGTAA -o LTP132.fasta tmp.fasta > cut_Radpt_results_LTP.txt

#remove sequences that are less than 250, then trim sequences to 250
awk 'BEGIN{RS=">";ORS=""}length($0)>249{print ">"$0}' LTP132.fasta > tmp.fasta
cut -c 1-250 tmp.fasta > LTP132.fasta

#combine LTP and succession data
cat LTP.fasta merged_combined_cut_filtered_uniques_otus.fa > LTP_succession.fa

#calculate distance matrix (~30min)
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -calc_distmx LTP_succession.fa -distmxout LTP_succession_mx.txt

#turn distance matrix into a tree (~1hr)
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -cluster_aggd LTP_succession_mx.txt -treeout LTP_succession.tree -clusterout LTP_succession_clusters.txt -id 0.80 -linkage min

