# Custom functions for microbiome analyses. Created Feb 28, 2017.

loadpax <- function(pkg){
  # (1) checks package installation, (2) installs them if not, then (3) loads them
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}


grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  # A function that lets you put multiple plots together and have a single shared legend (from the first plot)
  # from hadley on the internet...
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)

}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

pick_best <- function(x) {
  
  #to ensure the same random sampling.
  set.seed(7)
  
  # is it numeric?
  if (is.numeric(x) & sum(!is.na(x)) > 0) {
    
    val <- mean(x, na.rm = T)
    
    # if not...
  } else if (!is.numeric(x) & sum(!is.na(x)) > 0) {
    
    x <- table(as.character(x))
    
    # if values are of equal abundance, random sample
    if (length(unique(x)) == 1) val <- sample(names(x), 1)
    # otherwise select maximum value
    if (length(unique(x)) > 1) val <- names(x)[which.max(x)]
    
  } else {
    
    val <- NA
    
  }
  
  return(val)
  
}

bac_extract <- function(x, trait) {
  
  # for each species query result, filter blanks and pick best value
  for(i in 1:length(tax_bacdat)) {
    names(x[[i]]) <- tax_bacdat[[i]]
    x[[i]] <- x[[i]][sapply(x[[i]], length) > 0]
    x[[i]] <- sapply(x[[i]], pick_best)
  }
  
  # unlist, clean, return df
  x <- unlist(x)
  x[x == 'NULL'] <- NA
  x_names <- strsplit(names(x), '.', fixed = TRUE)
  x <- data.frame(
    Genus = sapply(x_names, function(y) y[[1]]),
    Species = sapply(x_names, function(y) y[[2]]),
    trait = trait,
    val = x
  )
  x <- x[!is.na(x$val), ]
  x$val <- as.character(x$val)
  return(x)
}

#suppress warnings as numeric
swan <- function (x) suppressWarnings(as.numeric(x))

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
}

# the ggplot theme that I like for easy loading...
library(ggplot2)
th <- theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank())


library(ggplot2)

#define a helper function (borrowed from the "ez" package)
ezLev=function(x,new_order){
  for(i in rev(new_order)){
    x=relevel(x,ref=i)
  }
  return(x)
}

ggcorplot = function(data,var_text_size,cor_text_limits){
  # normalize data
  for(i in 1:length(data)){
    data[,i]=(data[,i]-mean(data[,i]))/sd(data[,i])
  }
  # obtain new data frame
  z=data.frame()
  i = 1
  j = i
  while(i<=length(data)){
    if(j>length(data)){
      i=i+1
      j=i
    }else{
      x = data[,i]
      y = data[,j]
      temp=as.data.frame(cbind(x,y))
      temp=cbind(temp,names(data)[i],names(data)[j])
      z=rbind(z,temp)
      j=j+1
    }
  }
  names(z)=c('x','y','x_lab','y_lab')
  z$x_lab = ezLev(factor(z$x_lab),names(data))
  z$y_lab = ezLev(factor(z$y_lab),names(data))
  z=z[z$x_lab!=z$y_lab,]
  #obtain correlation values
  z_cor = data.frame()
  i = 1
  j = i
  while(i<=length(data)){
    if(j>length(data)){
      i=i+1
      j=i
    }else{
      x = data[,i]
      y = data[,j]
      x_mid = min(x)+diff(range(x))/2
      y_mid = min(y)+diff(range(y))/2
      this_cor = cor(x,y)
      this_cor.test = cor.test(x,y)
      this_col = ifelse(this_cor.test$p.value<.05,'<.05','>.05')
      this_size = (this_cor)^2
      cor_text = ifelse(
        this_cor>0
        ,substr(format(c(this_cor,.123456789),digits=2)[1],2,4)
        ,paste('-',substr(format(c(this_cor,.123456789),digits=2)[1],3,5),sep='')
      )
      b=as.data.frame(cor_text)
      b=cbind(b,x_mid,y_mid,this_col,this_size,names(data)[j],names(data)[i])
      z_cor=rbind(z_cor,b)
      j=j+1
    }
  }
  names(z_cor)=c('cor','x_mid','y_mid','p','rsq','x_lab','y_lab')
  z_cor$x_lab = ezLev(factor(z_cor$x_lab),names(data))
  z_cor$y_lab = ezLev(factor(z_cor$y_lab),names(data))
  diag = z_cor[z_cor$x_lab==z_cor$y_lab,]
  z_cor=z_cor[z_cor$x_lab!=z_cor$y_lab,]
  #start creating layers
  points_layer = layer(
    geom = 'point'
    , data = z
    , mapping = aes(
      x = x
      , y = y
    )
  )
  lm_line_layer = layer(
    geom = 'line'
    , geom_params = list(colour = 'red')
    , stat = 'smooth'
    , stat_params = list(method = 'lm')
    , data = z
    , mapping = aes(
      x = x
      , y = y
    )
  )
  lm_ribbon_layer = layer(
    geom = 'ribbon'
    , geom_params = list(fill = 'green', alpha = .5)
    , stat = 'smooth'
    , stat_params = list(method = 'lm')
    , data = z
    , mapping = aes(
      x = x
      , y = y
    )
  )
  cor_text = layer(
    geom = 'text'
    , data = z_cor
    , mapping = aes(
      x=y_mid
      , y=x_mid
      , label=cor
      , size = rsq
      , colour = p
    )
  )
  var_text = layer(
    geom = 'text'
    , geom_params = list(size=var_text_size)
    , data = diag
    , mapping = aes(
      x=y_mid
      , y=x_mid
      , label=x_lab
    )
  )
  f = facet_grid(y_lab~x_lab,scales='free')
  o = opts(
    panel.grid.minor = theme_blank()
    ,panel.grid.major = theme_blank()
    ,axis.ticks = theme_blank()
    ,axis.text.y = theme_blank()
    ,axis.text.x = theme_blank()
    ,axis.title.y = theme_blank()
    ,axis.title.x = theme_blank()
    ,legend.position='none'
  )
  size_scale = scale_size(limits = c(0,1),to=cor_text_limits)
  return(
    ggplot()+
    points_layer+
    lm_ribbon_layer+
    lm_line_layer+
    var_text+
    cor_text+
    f+
    o+
    size_scale
  )
}


# print list of loaded functions
print(data.frame(Custom_Functions = 
  c('loadpax: install+load multiple packages',
    'grid_arrange_shared_legend: Multiple plots, one legend',
    'multiplot: multiple plots',
    'pick_best: if numeric, return mean; if character, return most common',	
    'swan: Suppress Warnings, convert to As Numeric',
    'the ggplot theme that I like',
    'ggcorplot: plot pairwise correlations')))

