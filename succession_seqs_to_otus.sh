#Infant gut succession project script: going from raw seqs to an OTU table


#############################################
#############################################
#############################################
####helpful lines to remember
#converting text from windows to unix format
awk '{ sub("\r$", ""); print }' tmp.qsub > unix.txt

#converting text from unix to windows format
awk 'sub("$", "\r")' tmp.qsub > windows.txt
#############################################
#############################################
#############################################

#make new directory and go there
mkdir succession
cd succession/

#download all the sequences
#first Abx FASTA files
wget -r -np -nd https://pubs.broadinstitute.org/diabimmune/data/15

#remove the first period in filenames because they throw off usearch
#rename the '.' to a 'dot'
rename .fna fna *fna
rename \. dot *fna
rename fna \.fna *fna

#rename Abx fasta headers to match filenames (normally usearch does this, but we got unpaired ends)
#note that awk doesn't allow to directly save files in place -- must use the temporary dummy variable tmp.fastq

for i in *.fna; do
	awk '/^>/{print ">" substr(FILENAME,1,length(FILENAME)-4) "." (1+i++); next} 1' $i > tmp.fna && mv tmp.fna $i
done

#merge Abx fna into a single fna using cat
cat *.fna >Abx.fna

#convert fna to fastq using a custom perl script
#(`wget https://code.google.com/archive/p/fasta-to-fastq/`)
perl ../fasta_to_fastq.pl Abx.fna > Abx.fastq

#clean up fna header names so they only include sample name, time, and sequence number
#explanation: for each line that starts with a '@', substitute the first '_' with 'time', and then substitute the first string that starts with an '_', is followed by any character for any number of times, and ends with a(n escaped) '.' with a single '.'
# must again use the dummy variable with awk
awk '/^@/{sub("_","time");sub("_.*\\.",".")}1' Abx.fastq > tmp.fastq && tmp.fastq > Abx.fastq

### now downlaod and prep T1D FASTQ files
wget -r -np -nd https://pubs.broadinstitute.org/diabimmune/data/9
gunzip *gz

#note I don't have to rename headers.... that's what the -relabel is for
#get consensus sequences by merging pairs of sequences from the T1D fastq files
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_mergepairs *R1*.fastq -relabel @ -fastqout T1D.fastq -tabbedout T1D_mergepairs_report.txt -alnout T1D_aln.txt

#merge the two merged fastq files from each study
mkdir merged
cat Abx.fastq T1D.fastq >seqs.fastq

#Now go to GLBTR USEARCH pipeline and proceed from step 2b (but keep in mind that the quality scores from the Abx sequences are FAKE (set to 40, i think))
# the following code pipeline is taken from 
# https://github.com/GLBRC-TeamMicrobiome/TagAnalysis/blob/master/16S_USEARCH_Pipeline.md

#check quality scores
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_eestats2 seqs.fastq -output fastq_info/qualitystats.txt

#filter and truncate [~1hr]
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_filter seqs.fastq -fastq_maxee 1 -fastq_trunclen 250 -fastaout seqs_filtered.fa

#filter for unique [~1 hr]
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_uniques seqs_filtered.fa -fastaout seqs_filtered_uniques.fa -sizeout

#cluster into 97% OTUs and remove singletons [~30min]
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -cluster_otus seqs_filtered_uniques.fa -otus seqs_filtered_uniques_otus.fa -uparseout seqs_filtered_uniques_otus_uparse.txt -relabel OTU

#Identify ZOTUs and remove 8x singletons [~30min]
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -unoise3 seqs_filtered_uniques.fa -zotus seqs_filtered_uniques_zotus.fa -tabbedout seqs_filtered_uniques_zotus_report.txt

#rename representative sequence names from “Zotu” to “ZOTU” for the mapping back to the raw reads step to work correctly
sed -i 's/Zotu/ZOTU/g' seqs_filtered_uniques_zotus.fa

#NOT done
#map reads back to OTUs [~2+ hr]
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -otutab seqs.fastq -otus seqs_filtered_uniques_otus.fa -uc seqs_filtered_uniques_OTU_map.uc -otutabout seqs_filtered_uniques_OTU_table.txt -biomout seqs_filtered_uniques_OTU_jsn.biom -notmatchedfq seqs_filtered_uniques_otu_unmapped.fq

#NOT done
#map reads to zotus
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -otutab seqs.fastq -zotus seqs_filtered_uniques_zotus.fa -uc seqs_filtered_uniques_ZOTU_map.uc -otutabout seqs_filtered_uniques_ZOTU_table.txt -biomout seqs_filtered_uniques_ZOTU_jsn.biom -notmatchedfq seqs_filtered_uniques_otu_unmapped.fq

#map OTU to SILVA
wget http://drive5.com/sintax/silva_16s_v123.fa.gz
gunzip silva_16s_v123.fa.gz
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -sintax seqs_filtered_uniques_otus.fa -db silva_16s_v123.fa -tabbedout seqs_filtered_uniques_otus_taxonomy_SILVA.sintax -strand both

#map ZOTU to SILVA (haven't run this yet)
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -sintax seqs_filtered_uniques_zotus.fa -db silva_16s_v123.fa -tabbedout seqs_filtered_uniques_zotus_taxonomy_SILVA.sintax -strand both

#map OTU to LTP (not sure if necessary, but it's fast)
wget http://drive5.com/sintax/ltp_16s_v123.fa.gz
gunzip ltp_16s_v123.fa.gz
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -sintax seqs_filtered_uniques_otus.fa -db ltp_16s_v123.fa -tabbedout seqs_filtered_uniques_otus_taxonomy_LTP.sintax -strand both

#map OTU to gg (only necessary if i want to reverse map back to gg for bugbase predictions)
wget http://drive5.com/sintax/gg_16s_13.5.fa.gz
gunzip gg_16s_13.5.fa.gz
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -sintax seqs_filtered_uniques_otus.fa -db gg_16s_13.5.fa -tabbedout seqs_filtered_uniques_otus_taxonomy_GG.sintax -strand both

#map ZOTU to LTP (not sure if necessary yet)
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -sintax seqs_filtered_uniques_zotus.fa -db ltp_16s_v123.fa -tabbedout seqs_filtered_uniques_zotus_taxonomy_LTP.sintax -strand both