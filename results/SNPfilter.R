library(ggplot2)
library(reshape2)
library(gridExtra)
library(dplyr)

df <- read.table("heterozygous.txt", header = FALSE)

colnames(df) <- c('sample', 'xvalue', 'allowed', 'repeated', 'heterozygous', 'lowCoverage')

dfm <- melt(df, id.vars = c('xvalue', 'sample'))

colnames(dfm) <- c('value', 'sample', 'filter', 'numSNPs')

ggplot(dfm, aes(x = value, y = numSNPs)) + 
  
  geom_line(aes(color = filter), size = 1) +
  
  coord_trans(y = "log") +
  
  scale_y_continuous(breaks = c(100, 200, 500, 1000, 2000, 4000, 8000)) +
  
  facet_wrap(sample ~ ., ncol = 4 ) +
  
  ggtitle ("Effect of Different Alt:Ref Ratio Values on Filtering")


df <- read.table("low.txt", header = FALSE)

colnames(df) <- c('sample', 'xvalue', 'allowed', 'repeated', 'heterozygous', 'lowCoverage')

dfm <- melt(df, id.vars = c('xvalue', 'sample'))

colnames(dfm) <- c('value', 'sample', 'filter', 'numSNPs')

ggplot(dfm, aes(x = value, y = numSNPs)) + 
  
  geom_line(aes(color = filter), size = 1) +
  
  coord_trans(y = "log") +
  
  scale_y_continuous(breaks = c(100, 200, 500, 1000, 2000, 4000, 8000)) +
  
  facet_wrap(sample ~ ., ncol = 4 ) +
  
  ggtitle ("Effect of Different LowCoverage Values on Filtering")
