library(dplyr)
library(ggplot2)
library(vegan)
library(poilog)

set.seed(7)

#in this script we explore how euclidean trait dissimilarity changes as bray-curtis community distance increases

#create an artificial community with a randomly generated Poisson lognormal species abundance distribution
#replicate it (this is the community that will diverge in composition)
nspp <- 1000
c0 <- c(rpoilog(nspp, mu=1.0, sig=1.0, keep0 = TRUE), rep(0, times = nspp))
c1 <- c0
n_individuals <- sum(c0)

#load real traits from the succession datasets. Use these to randomly asign trait values to artificial species
traits <- readRDS('C:\\Users\\John\\Documents\\msu\\microbiome_trait_succession\\data\\traits_sparse.RDS')
traits <- select(traits, -Genus, -Species)

#the strategy here is to have two sets of 1000 species (2000 total potential species across both community). At first, the first 1000 species of both communities have equal abundances (derived by randomly sampling using rpoilog()) and 1000 species with zero individuals.
traits <- sapply(traits, function(x) sample(x[!is.na(x)], 2*nspp, replace = TRUE))
traits <- apply(traits, 2, scale)

#Now, during each time step in the following loop, one individual will leave the first 1000 species, and an new individual will appear in the second set of 1000 species.

#create empty dataframe for results
dists <- data.frame(
    shared = c(sum(c0), rep(NA, sum(c0))),
    EU = c(0, rep(NA, sum(c0))),
    BC = c(0, rep(NA, sum(c0))),
    EU_traits = c(0, rep(NA, sum(c0))))

for (i in 2:(n_individuals + 1)) {
  
  #randomly select the species that will lose an individual
  losers <- seq_along(c1)[seq_along(c1) <= nspp & c1 > 0]
  loser <- if (length(losers) > 1) {
    sample(losers, 1, prob = c1[seq_along(c1) <= nspp & c1 > 0])
  } else {
    losers
  }
  
  #randomly select the species that will gain an individual
  winner <- nspp + sample(seq_along(c0), 1, prob = c0)
  c1[loser] <- c1[loser] - 1
  c1[winner] <- c1[winner] + 1
  
  #calculate distances and populate dataframe with EU and BC results
  dists$shared[i] = sum(c1[c(1:nspp)])
  dists$EU[i] <- vegdist(rbind(c0, c1), method = 'euclidean')
  dists$BC[i] <- vegdist(rbind(c0, c1), method = 'bray')
  
  #calculate trait-based euclidean distance and populate dataframe with results
  tmp <- rbind(
    apply(traits, 2, function(x) weighted.mean(x, c0)),
    apply(traits, 2, function(x) weighted.mean(x, c1)))
  dists$EU_traits[i] <- vegdist(tmp, method = 'euclidean')
  
}

dists %>%
  ggplot(aes(x = BC, y = EU)) +
  geom_point() +
  stat_smooth(method = 'lm', color = 'red') +
  labs(x = "Bray-Curtis distance", y = 'OTU-based Euclidean distance') +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank())

dists %>%
  ggplot(aes(x = BC, y = EU_traits)) +
    geom_point() +
    stat_smooth(method = 'lm', color = 'red') +
    labs(x = "Bray-Curtis distance", y = 'Trait-based Euclidean distance') +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank())





