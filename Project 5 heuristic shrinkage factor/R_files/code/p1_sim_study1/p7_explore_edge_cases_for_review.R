### Set working directory
rm(list=ls())

### Set working directory
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_8.2/")

### Source functions
R.func.sources <- list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Load functions
source("R/sim_functions.R")

### Call libraries
library(ggplot2)

### Read in data
metadata <- readRDS("data/sim_study1_metadata.rds")
metadata <- subset(metadata, C.pop != 0)
metadata <- subset(metadata, cat == FALSE)

### Add the magnitude of bias/absolute bias
metadata$S.VH.AB <- abs(metadata$S.VH.diff.mean)
metadata$S.boot.AB <- abs(metadata$S.boot.diff.mean)

# str(metadata)
# temp <- subset(metadata, S.boot.mean > 1.5)
# temp <- subset(metadata, S.VH.mean < 0)
# temp <- subset(metadata, corr == TRUE)
# temp <- temp[temp$S.boot.mean > 1 & temp$S.boot.mean < 2, ]
# 
# max(metadata$S.VH.mean)
# max(metadata$S.boot.mean)
# temp <- subset(metadata, corr == TRUE)
# sd(temp$S.boot.diff.mean)
# 
# temp <- temp[temp$S.boot.mean > 0, ]
# temp4 <- temp[temp$S.boot.mean > 3, ]
# 
# temp4 <- temp[temp$S.boot.mean < 0, ]

df <- subset(metadata, corr == TRUE & S.boot.mean < 2 & S.boot.mean > 0)


### Create colour range for ggplot
fun_color_range <- colorRampPalette(c("green", "blue"))
my_colors <- fun_color_range(20)
colnames(df)

### Create plot S.boot
  gg.out <- ggplot(aes(x = S.boot.mean, y = S.pop.mean, color = C.pop), data = df) +
    scale_colour_gradientn(colors = my_colors) + xlim(c(0,1.5)) + ylim(c(0,1.5))


  gg.out <- gg.out + geom_point(size = 0.5) + geom_abline(intercept = 0, slope = 1)

  Cairo::CairoPNG("figures/TEMP_bigger_range_axes_S_boot.png", width = 7, height = 7, unit = "in", dpi = 300)
  plot(gg.out)
  dev.off()
  gg.out

  ### Create plot S.VH
  gg.out <- ggplot(aes(x = S.VH.mean, y = S.pop.mean, color = C.pop), data = df) +
    scale_colour_gradientn(colors = my_colors) + xlim(c(0,1.5)) + ylim(c(0,1.5))
  
  
  gg.out <- gg.out + geom_point(size = 0.5) + geom_abline(intercept = 0, slope = 1)
  
  Cairo::CairoPNG("figures/TEMP_bigger_range_axes_S_VH.png", width = 7, height = 7, unit = "in", dpi = 300)
  plot(gg.out)
  dev.off()
  gg.out

  ### which are the flower bits
  df_temp <- subset(df, S.pop.mean > 1.1 & (S.boot.mean > 1.15 | S.boot.mean < 0.9))
  df_temp
  gg.out <- ggplot(aes(x = S.boot.mean, y = S.pop.mean, color = C.pop), data = df_temp) +
    scale_colour_gradientn(colors = my_colors) 
  
  
  gg.out <- gg.out + geom_point(size = 0.5) + geom_abline(intercept = 0, slope = 1)
  
  gg.out
  
  
  
  ### Add extra bits
gg.out <- gg.out + 
  tune::coord_obs_pred() +
  theme(legend.position = "bottom", 
        text = element_text(size = 8), 
        legend.text = element_text(size = 8), 
        legend.key.height = unit(0.3, "cm"),
        legend.margin = unit(c(-0.2,0,0,0), "cm"),
        legend.box.margin = unit(c(-0.2,0,0,0), "cm"),
        plot.margin = grid::unit(c(0,0,0,0), "cm")) +
  labs(color = by.label, x = x.label, y = y.label)


library(ggplot2)
ggplot(data = temp) + geom_point(aes(x = S.boot.mean, y = S.pop.mean)) + geom_abline(slope = 1, intercept = 0)
plot(temp$S.boot.mean, temp$S.pop.mean)


temp2 <- subset(metadata, S.pop.mean > 1.2)

### Write a function to get the quantiles of the bias for a given range of C-statistics
get_quantiles_by_C <- function(C1,C2,correlated = TRUE, magnitude = FALSE){
  
  # Subset dataset to be correlated/uncorrelated
  df <- subset(metadata, C.pop > C1 & C.pop <= C2 & corr == correlated)
  
  if (magnitude == FALSE){
    # Get quantiles for bias
    quantiles.S.VH <- quantile(as.numeric(df$S.VH.diff.mean), p = c(0.025,0.25,0.5,0.75,0.975))
    quantiles.S.boot <- quantile(as.numeric(df$S.boot.diff.mean), p = c(0.025,0.25,0.5,0.75,0.975))
    
    # Get mean for bias
    mean.S.VH <- mean(df$S.VH.diff.mean)
    mean.S.boot <- mean(df$S.boot.diff.mean)
    
    # Get sd for bias
    sd.S.VH <- sd(df$S.VH.diff.mean)
    sd.S.boot <- sd(df$S.boot.diff.mean)
    
  } else if (magnitude == TRUE){
    # Get quantiles for absolute bias
    quantiles.S.VH <- quantile(as.numeric(df$S.VH.AB), p = c(0.025,0.25,0.5,0.75,0.975))
    quantiles.S.boot <- quantile(as.numeric(df$S.boot.AB), p = c(0.025,0.25,0.5,0.75,0.975))
    
    # Get mean for absolute bias
    mean.S.VH <- mean(df$S.VH.AB)
    mean.S.boot <- mean(df$S.boot.AB)
    
    # Get sd for absolute bias
    sd.S.VH <- sd(df$S.VH.AB)
    sd.S.boot <- sd(df$S.boot.AB)
  }
  
  
  # Concatenate into string
  mean.S.VH <- paste(sprintf("%.3f", mean.S.VH), " (", sprintf("%.3f", sd.S.VH), ")", sep = "")
  mean.S.boot <- paste(sprintf("%.3f", mean.S.boot), " (", sprintf("%.3f", sd.S.boot), ")", sep = "")
  names(mean.S.VH) <- "mean"
  names(mean.S.boot) <- "mean"
  
  # Concatenate
  quantiles.S.VH <- c(sprintf("%.3f", quantiles.S.VH), mean.S.VH)
  quantiles.S.boot <- c(sprintf("%.3f", quantiles.S.boot), mean.S.boot)
  
  return(list("S.VH" = quantiles.S.VH,
              "S.boot" = quantiles.S.boot))
}

### Write a function to get the number of observations with a given range of C-statistics
get_N <- function(C1,C2,correlated = TRUE){
  
  df <- subset(metadata, C.pop > C1 & C.pop <= C2 & corr == correlated)
  n <- nrow(df)
  paste("N=", n)
  
}

### Create a list of C-statistic ranges
quantiles.list <- list(c(0.6, 0.65),
                       c(0.65, 0.7),
                       c(0.7, 0.75),
                       c(0.75, 0.8),
                       c(0.8, 0.85),
                       c(0.85, 0.9))


# Bias

## Get results for correlated data

### Get quantiles
S.by.quantiles <- lapply(quantiles.list, function(x) {get_quantiles_by_C(x[1],x[2])})
S.all <- get_quantiles_by_C(0,1)

### Put into a table
S.by.quantiles <- lapply(S.by.quantiles, function(x) {
  out <- rbind(x[[1]], x[[2]])
  rownames(out) <- c("S_VH", "S_boot")
  return(out)
})
S.all <- rbind(S.all[[1]], S.all[[2]])


### Print S_pop by quantiles
knitr::kable(S.all)
knitr::kable(S.by.quantiles[[1]])
knitr::kable(S.by.quantiles[[2]])
knitr::kable(S.by.quantiles[[3]])
knitr::kable(S.by.quantiles[[4]])
knitr::kable(S.by.quantiles[[5]])
knitr::kable(S.by.quantiles[[6]])

### Get N values
get_N(0,1)
lapply(quantiles.list, function(x) {get_N(x[1],x[2])})

```

# Get results for uncorrelated data

```{r}
### Get quantiles
S.by.quantiles <- lapply(quantiles.list, function(x) {get_quantiles_by_C(x[1],x[2], correlated = FALSE)})
S.all <- get_quantiles_by_C(0,1, correlated = FALSE)

### Put into a table
S.by.quantiles <- lapply(S.by.quantiles, function(x) {
  out <- rbind(x[[1]], x[[2]])
  rownames(out) <- c("S_VH", "S_boot")
  return(out)
})
S.all <- rbind(S.all[[1]], S.all[[2]])


### Print S_pop by quantiles
knitr::kable(S.all)
knitr::kable(S.by.quantiles[[1]])
knitr::kable(S.by.quantiles[[2]])
knitr::kable(S.by.quantiles[[3]])
knitr::kable(S.by.quantiles[[4]])
knitr::kable(S.by.quantiles[[5]])
knitr::kable(S.by.quantiles[[6]])

### Get N values
get_N(0,1, correlated = FALSE)
lapply(quantiles.list, function(x) {get_N(x[1],x[2], correlated = FALSE)})
```

# Magnatude of Bias

## Get results for correlated data

```{r}
### Get quantiles
S.by.quantiles <- lapply(quantiles.list, function(x) {get_quantiles_by_C(x[1],x[2], magnitude = TRUE)})
S.all <- get_quantiles_by_C(0,1, magnitude = TRUE)

### Put into a table
S.by.quantiles <- lapply(S.by.quantiles, function(x) {
  out <- rbind(x[[1]], x[[2]])
  rownames(out) <- c("S_VH", "S_boot")
  return(out)
})
S.all <- rbind(S.all[[1]], S.all[[2]])


### Print S_pop by quantiles
knitr::kable(S.all)
knitr::kable(S.by.quantiles[[1]])
knitr::kable(S.by.quantiles[[2]])
knitr::kable(S.by.quantiles[[3]])
knitr::kable(S.by.quantiles[[4]])
knitr::kable(S.by.quantiles[[5]])
knitr::kable(S.by.quantiles[[6]])

### Get N values
get_N(0,1)
lapply(quantiles.list, function(x) {get_N(x[1],x[2])})

```

## Get results for uncorrelated data

```{r}
### Get quantiles
S.by.quantiles <- lapply(quantiles.list, function(x) {get_quantiles_by_C(x[1],x[2], correlated = FALSE, magnitude = TRUE)})
S.all <- get_quantiles_by_C(0,1, correlated = FALSE, magnitude = TRUE)

### Put into a table
S.by.quantiles <- lapply(S.by.quantiles, function(x) {
  out <- rbind(x[[1]], x[[2]])
  rownames(out) <- c("S_VH", "S_boot")
  return(out)
})
S.all <- rbind(S.all[[1]], S.all[[2]])


### Print S_pop by quantiles
knitr::kable(S.all)
knitr::kable(S.by.quantiles[[1]])
knitr::kable(S.by.quantiles[[2]])
knitr::kable(S.by.quantiles[[3]])
knitr::kable(S.by.quantiles[[4]])
knitr::kable(S.by.quantiles[[5]])
knitr::kable(S.by.quantiles[[6]])

### Get N values
get_N(0,1, correlated = FALSE)
lapply(quantiles.list, function(x) {get_N(x[1],x[2], correlated = FALSE)})
```