# In this script, I use two different methods to estimate missing trait data.
# In the first, I use genus-level means if there are any species data available,
# In the second, I use PICRUSt-based inference

library(data.table)
library(tidyverse)

##############################################################
##################  METHOD 1: GENUS-EVEL MEANS   #############
##############################################################

wd <- 'C:\\Users\\John\\Documents\\msu\\microbiome_trait_succession\\'
setwd(wd)
x <- readRDS('data\\traits_sparse.RDS')
tax <- readRDS('data\\ggtax.RDS') 

#calculate genus-level means
x_genus_means <- x %>%
  gather(trait, val, -otu) %>%
  left_join(tax[, c('Genus','otu')], by = c('otu')) %>%
  group_by(trait, Genus) %>%
  filter(sum(!is.na(val)) > 0) %>%
  mutate(val = ifelse(is.na(val), mean(val, na.rm = T), val)) %>%
  spread(trait, val) %>%
  ungroup() %>%
  select(-Genus)

saveRDS(x_genus_means, file = 'data\\traits_genus_means.RDS')





##############################################################
##################  METHOD 2: PICRUSt   ######################
##############################################################


# Finally figured out how to use the picrust algorithm to estimate unknown trait values
# I manually ran dataframes with known trait values through picrust on the hpcc. First, create the trait input files
#
if (FALSE) {
  x_tmp <- x %>%
    gather(trait, val, -otu) %>%
    filter(!is.na(val))
  
  for (i in unique(x_tmp$trait)) {
    tmp <- x_tmp %>%
      filter(trait == i) %>%
      mutate(
        otu = as.numeric(gsub("otu", "", otu)), 
        trait = NULL)
    colnames(tmp) <- c('GenomeID', i)
    write.table(tmp, file = paste0('data\\picrust_input_traits\\', i, '.tab'), row.names = FALSE, sep = '\t', quote = FALSE)
  }
}

##then, on the HPCC:
## open hpcc terminal
#ssh guittarj@hpcc.msu.edu

## to move files between local and remote locations
# the r is if i want to move the whole directory. or I can just omit it. 
#scp -r input.file guittarj@hpcc.msu.edu:/mnt/home/guittarj/input/

######################email from HPCC
##As recommended by the software's guidance, the installation would be best done via anaconda. On the HPCC, anaconda is installed to a user's own home directory, not system wide.

##If you haven't used anaconda before, please go to https://urldefense.proofpoint.com/v2/url?u=https-3A__docs.anaconda.com_anaconda_install_linux&d=DwIDaQ&c=nE__W8dFE-shTxStwXtp0A&r=sp5b5l5mBaeg053XJnK_DQ&m=OJvbiW1gzIux-cnqcJZ9rLJtGS6BNZcxb0lcwByXjTw&s=1RusCwkH7IKuteDxeEY6bMyGGwrjDvQmHsi7FbAhofs&e= to download, install and learn about it. Note that when you choose Anaconda for Linux installer from https://urldefense.proofpoint.com/v2/url?u=https-3A__www.anaconda.com_download_-23linux&d=DwIDaQ&c=nE__W8dFE-shTxStwXtp0A&r=sp5b5l5mBaeg053XJnK_DQ&m=OJvbiW1gzIux-cnqcJZ9rLJtGS6BNZcxb0lcwByXjTw&s=3-UGwQgvetkerzfr0dnECLwICd1WYX_dMm4VumuEV6k&e=, it's better to go with the Python 3.6 version.  Also, we have a short conda tutorial at https://wiki.hpcc.msu.edu/display/ITH/Using+conda in case you find the full documentation (https://urldefense.proofpoint.com/v2/url?u=https-3A__conda.io_docs_user-2Dguide_getting-2Dstarted.html&d=DwIDaQ&c=nE__W8dFE-shTxStwXtp0A&r=sp5b5l5mBaeg053XJnK_DQ&m=OJvbiW1gzIux-cnqcJZ9rLJtGS6BNZcxb0lcwByXjTw&s=pTI8yNeEjgngZ8H_1GQ0PnGZ2r4a4zkNBTF7y4GKoi8&e=) is too much.

##Once anaconda has been installed successfully in your home directory, you can run the following commands to complete the installation of picrust, according to https://urldefense.proofpoint.com/v2/url?u=http-3A__picrust.github.io_picrust_install.html-23install&d=DwIDaQ&c=nE__W8dFE-shTxStwXtp0A&r=sp5b5l5mBaeg053XJnK_DQ&m=OJvbiW1gzIux-cnqcJZ9rLJtGS6BNZcxb0lcwByXjTw&s=lWJWjXPLEMBI2aF2IMAwIWqfNqlkJEWkDkRtknK_yBA&e=:

#conda create --name guittarj_picrust python=2.7.6
#source activate guittarj_picrust
#conda install -c bioconda picrust
#download_picrust_files.py # this creates a dir $HOME/anaconda3/envs/guittarj_picrust/lib/python2.7/site-packages/picrust/data/ to store downloaded data
#################################


##Then, run picrust's predict genotype for each trait*taxa. For example:
#mkdir output/Aggregation_score
#format_tree_and_trait_table.py -t 99_otus.tree -i input/traits_Aggregation_score.tab -o output/Aggregation_score/
#ancestral_state_reconstruction.py -i output/Aggregation_score/trait_table.tab -t output/Aggregation_score/pruned_tree.newick -o output/Aggregation_score/asr_counts.tab
#predict_traits.py --no_round -i output/Aggregation_score/trait_table.tab -t output/Aggregation_score/reference_tree.newick -r output/Aggregation_score/asr_counts.tab -o output/Aggregation_score/predictions_Aggregation_score.tab
#source deactivate guittarj_picrust

# After copying output to a local directory, parse and combine
if (TRUE) {
  wd_tmp <- "C:/Users/John/Documents/msu/analysis/picrust-1.1.2/traits/output/"
  dirs <- list.files(wd_tmp)
  
  pidat <- data.frame(trait = character(0), otu = numeric(0), val = numeric(0))
  for (i in dirs) {
    tmp <- data.frame(i, read.table(paste0(wd_tmp,dirs[dirs == i],"/predictions.tab"), sep = '\t'), 
                      stringsAsFactors = FALSE)
    colnames(tmp) <- colnames(pidat)
    pidat <- rbind(pidat, tmp)
  }
  
  #cleanup
  pidat$otu <- paste0('otu', pidat$otu)
  pidat$trait[pidat$trait == 'iga'] <- 'IgA'
  pidat$trait <- factor(pidat$trait)
  pidat$otu <- factor(pidat$otu)
  saveRDS(pidat, file = 'data\\traits_picrust.RDS')
}




##############################################################
##################  METHOD 3: CASTOR   ######################
##############################################################

require(ape)
require(castor)
mytree <- read.tree("C:/Users/John/Documents/msu/analysis/picrust-1.1.2/99_otus.tree")
traits <- as.data.frame(readRDS('C:/Users/John/Documents/msu/microbiome_trait_succession/data/traits_sparse.RDS'))
traits$otu <- gsub('otu','',traits$otu)
traits <- filter(traits, traits$otu %in% mytree$tip.label)

dat <- data.frame(
  otu = character(), 
  trait = character(), 
  val = numeric(), stringsAsFactors = FALSE)

for(i in names(traits)[names(traits) != 'otu']) {
  
  mytips <- traits[[i]][match(mytree$tip.label, traits$otu)]
  #out <- hsp_independent_contrasts(mytree, mytips)
  out <- hsp_subtree_averaging(mytree, mytips)
  tmp <- data.frame(
    otu = mytree$tip.label, 
    trait = i, 
    val = out$states[1:length(mytips)], stringsAsFactors = FALSE)
  dat <- bind_rows(dat, tmp)

}

dat <- dat %>% 
  mutate(
    otu = factor(paste0('otu',otu)),
    trait = factor(trait))

#saveRDS(dat, file = 'data\\traits_castor_independent_contrasts.RDS')
saveRDS(dat, file = 'data\\traits_castor_subtree_averaging.RDS')



##############################################################
##################  METHOD 4: Bug-Base   #####################
##############################################################


if (FALSE) {
  otus_wide <- readRDS('data\\otus_wide.RDS')
  tmp <- otus_wide[, c(4:ncol(otus_wide))]
  x <- data.frame(otu = gsub("otu", "", colnames(tmp)), t(tmp))
  colnames(x) <- c('#OTU ID',otus_wide$sampleID)
  write.table(x, 'data\\otu_table.txt', sep = '\t', row.names = FALSE, quote = FALSE)
  system('biom convert -i C:/Users/John/Documents/msu/microbiome_trait_succession/data/otu_table.txt -o C:/Users/John/Documents/msu/microbiome_trait_succession/data/otu_table.biom --table-type="OTU table" --to-json')
  # then, upload the resulting biom file to https://bugbase.cs.umn.edu/upload.html and download and unzip results
  
}

#processing Bugdat output
j <- read.table('data\\BugBase_Analyses\\predicted_phenotypes\\predictions.txt')
j <- data.frame(sampleID = row.names(j), j, row.names = NULL)
saveRDS(j, 'data\\traits_bugbase.RDS')
