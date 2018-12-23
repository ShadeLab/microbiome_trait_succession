# Custom functions for microbiome analyses.
# Created by John Guittar

scaleFUN <- function(x) sprintf("%.1f", x)

loadpax <- function(pkg){
  # (1) checks package installation, (2) installs them if not, then (3) loads them
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}

swan <- function (x) {
	# suppress warnings when converting to numeric
	suppressWarnings(as.numeric(x))
}

bac_search <- function (x) {
  # needs list of genera to search for
  
  bacdat <- list()
  #just use the first page of results
  page <- '1'
  taxon <- 'https://bacdive.dsmz.de/api/bacdive/taxon/'
  
  for (i in 1:length(x)) {
    
    url_genus <- URLencode(paste0(taxon, x[i], '/?page=', page, '&format=json'))
    response <- getURL(url_genus, userpwd="guittarj@msu.edu:twkYkcJbxCQzEvwkUi", httpauth = 1L)
    jsondata <- fromJSON(response)
    urls <- unlist(jsondata$results)
    
    if (length(urls) > 0) {
      genus <- x[i] 
      bacdat[[genus]] <- list()
      for (j in 1:length(urls)) {
        urlx <- URLencode(paste0(urls[j],'?&format=json'))
        response <- getURL(urlx, userpwd="guittarj@msu.edu:twkYkcJbxCQzEvwkUi", httpauth = 1L)
        jsondata <- fromJSON(response)
        bacdat[[genus]][[j]] <- jsondata
      }
    }
    print(paste("Downloaded", length(urls), "entries for", x[i]))
  }
  return(bacdat)
}

# the ggplot theme that I like for easy loading...
library(ggplot2)
th <- theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank())

meltdist <- function(x, y, tree = NULL, method = c('bray','euclidean','unifrac')) {
  #meltdist takes a list of samples (x) and an otu table (y) and calculates intersample dissimilarity using (meth)
  
  #prune the original otu table
  otumat <- y[y$sampleID %in% x, !names(y) %in% c('sampleID','subject','t','delivery')] %>% as.data.frame()
  row.names(otumat) <- y$sampleID[y$sampleID %in% x]
  
  if (method == 'unifrac') {
  	
  	#create phyloseq object
  	j <- phyloseq(otu_table(otumat, taxa_are_rows = FALSE),
                    drop.tip(tree, tree$tip.label[!tree$tip.label %in% names(otumat)]))

  	#calculate UniFrac
    j <- UniFrac(j, weighted = TRUE)

  } else {

  	#calculate bray-curtis or euclidean dissimilarity
  	j <- vegdist(otumat, method = method)

  }

  #melt into a dataframe, append times
  j <- data.frame(melt(as.matrix(j)), stringsAsFactors = FALSE) %>%
    transmute(
    	sample1 = Var1,
    	sample2 = Var2,
    	t1 = y$t[match(sample1, y$sampleID)],
    	t2 = y$t[match(sample2, y$sampleID)],
    	dist = value) %>%
    mutate_if(is.factor, as.character)

  return(j)

}

