library(gatom) 
library(sybilSBML)
library(R.utils)
library(sybil)
library(gtools)
library(slam)
library(Matrix)
library(foreach)
library(data.table)
library(dplyr)
library(plyr)
library(org.Mm.eg.db)
library(DESeq2)
library(igraph)
library(pheatmap)
library(devtools)
library(limma)
library(cluster)
library(ggplot2)
library(ggrepel)
library(knitr)

#read FVA result
df <- read.csv('out.csv', header = T)

#read annotations
reac <- read_delim("bigg_models_reactions.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
reac$model_list <- NULL

metab <- read_delim("bigg_models_metabolites.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
metab$model_list <- NULL

annot <- read.table('annot.txt', "\t", header = T)
annot$Source <- as.character(annot$Source)
annot$Sid <- as.character(annot$Sid)
annot$Equa <- as.character(annot$Equa)

#read analysis from MetaNetX
ana <- read.table('analysis.txt', "\t", header = T)
ana$Scope <- as.character(ana$Scope)
anl$ID <- as.character(ana$ID)
ana$Key <- as.character(ana$Key)
str(ana)

df$range <- df$maximum - df$minimum #range
df$X <- as.character(df$X)
str(df)
df1 <- df[order(df$range),] #sorting
xdf <- merge(df1, reac, by.x = "X", by.y = "bigg_id", all.x = TRUE) #annotation
xdf <- xdf[order(xdf$range),] #sorting
xdf$database_links <- NULL
xdf$old_bigg_ids <- NULL

annot1 <- annot[,c(1,2)]
analy <- merge(ana, annot1, by.x = "ID", by.y = "Sid", all.x = TRUE) #annotation

analy1 <- analy
analy1$Scope <- NULL

#subset results of reaction knockout_____________________
analyRKO_c <- analy1[ which(analy1$Key=='RKO-c'),]
analyRKO_c$Key <- NULL
analyRKO_c$RKO_c <- analyRKO_c$Value
analyRKO_c$Value  <- NULL

analyRKO_r <- analy1[ which(analy1$Key=='RKO-r'),]
analyRKO_r$Key <- NULL
analyRKO_r$RKO_r <- analyRKO_r$Value
analyRKO_r$Value  <- NULL

analyRKO <- merge(analyRKO_c, analyRKO_r, by = "Source")

head(xdf)
head(analyRKO)


#subset results of gene knockout_____________________

analyPKO_c <- ana[ which(ana$Key=='PKO-c'),]
analyPKO_c$Scope <- NULL
analyPKO_c$Key <- NULL
analyPKO_c$PKO_c <- analyPKO_c$Value
analyPKO_c$Value  <- NULL

analyPKO_r <- ana[ which(ana$Key=='PKO-r'),]
analyPKO_r$Scope <- NULL
analyPKO_r$Key <- NULL
analyPKO_r$PKO_r <- analyPKO_r$Value
analyPKO_r$Value  <- NULL

analyPKO <- merge(analyPKO_c, analyPKO_r, by = "ID")

head(xdf)
head(analyPKO)

gen <-read.table('gene_reaction.txt', "\t", header = T) # read extra annotation
genn <- merge(gen, annot1, by.x = "Reaction", by.y = "Sid", all.x = TRUE)
genn$Reaction <- NULL

genne <- merge(analyPKO, genn, by.x = "ID", by.y = "Complex", all.x = TRUE) #annotation
#genne1 <- genne[which(genne$PKO_c == 'impa' | genne$PKO_c == 'leth'),]
genne2 <- genne[ which(genne$Source != "NA"),]

#subset FBA results_____________________
analyFBA <- ana[ which(ana$Key=='FBA-f'),]
analyFBA$Scope <- NULL
analyFBA$Key <- NULL
analyFBA$FBA <- analyFBA$Value
analyFBA$Value  <- NULL

# merging FVA results with results of reaction knockout, gene knockout and FBA 
xdfx <- merge(xdf, analyRKO, by.x = "X", by.y = "Source", all.x = TRUE)
xdfx <- merge(xdfx, genne2, by.x = "X", by.y = "Source", all.x = TRUE)
xdfx$ID.y <- NULL
xdfx <- merge(xdfx, analyFBA, by.x = "ID.x", by.y = "ID", all.x = TRUE)
xdfx1 <- xdfx
xdfx1$FBA <- as.numeric(xdfx$FBA)

#subset impared and lethal reactions RKO-c<0.95
xdfx1 <- xdfx[which(xdfx$RKO_c == 'impa' | xdfx$RKO_c == 'leth' ),]
#subset only lethal reactions RKO-c<0.05
xdfx1 <- xdfx[which(xdfx$RKO_c == 'leth' ),]

#some graphs
ggplot(xdfx1) +
  geom_point(aes(x=FBA, y=X, color=RKO_c), size=5, shape=20) +
  labs(title="Reaction knockout", x ="Flux", y = "Reactions") +
  theme_bw() 


ggplot(xdfx1) +
  geom_point(aes(x=range, y=FBA, color=RKO_c), size=2) +
  geom_text_repel(mapping=aes(x=range, y=FBA, color=RKO_c,label=X), size=3)+
  labs(title="Reaction knockout", x ="Max Flux - Min Flux", y = "Flux") +
  theme_bw() 

#write only lethal reactions
write.table(xdfx1$reaction_string,"met_lethal.txt",quote=F,sep="\t", row.names = F)


#read table with metabolites, which participate in lethal reactions
che <- read.table('met_lethal_sum.txt', "\t", header = T)
che <- che[order(che$coun),]

#read table with metabolites, which participate in lethal reactions (compartment annotation)
chems <- read.table('met_lethal_sum1.txt', "\t", header = T)
#chems1 <- ddply(chems,"chem",numcolwise(sum))

chemsx <- chems[chems$coun >2,]
ggplot(chemsx, aes(chem, coun)) +
  geom_col(aes(fill = compart))+
  theme_bw() 


#view all data
ggplot(xdfx) +
  geom_point(aes(x=range, y=X, color=RKO_c), size=5, shape=20) +
  theme_bw()

#first group
df_0 <- xdfx[xdfx$range<2999, ]
write.table(df_0$ID,"gene_group_0.txt",quote=F,sep="\t", row.names = F)
write.table(df_0$X,"react_group_0.txt",quote=F,sep="\t", row.names = F)
write.table(df_0[,c(2,14)],"reactF_group_0.txt",quote=F,sep="\t", row.names = F)
ggplot(df_0) +
  geom_point(aes(x=range, y=X, color=RKO_c), size=5, shape=20) +
  theme_bw()

ggplot(df_0) +
  geom_point(aes(x=range, y=FBA, color=RKO_c), size=2) +
  geom_text_repel(mapping=aes(x=range, y=FBA, color=RKO_c,label=X), size=3)+
  labs(title="Reaction knockout", x ="Max Flux - Min Flux", y = "Flux") +
  theme_bw() 

#2 group
df_5 <- xdfx[xdfx$range<7999, ]
df_5 <- df_5[df_5$range>3000, ]
write.table(df_5$ID,"gene_group_5.txt",quote=F,sep="\t", row.names = F)
write.table(df_5$X,"react_group_5.txt",quote=F,sep="\t", row.names = F)
write.table(df_5[,c(2,14)],"reactF_group_5.txt",quote=F,sep="\t", row.names = F)
ggplot(df_5) +
  geom_point(aes(x=range, y=X, color=RKO_c), size=5, shape=20) +
  theme_bw() 

ggplot(df_5) +
  geom_point(aes(x=range, y=FBA, color=RKO_c), size=2) +
  geom_text_repel(mapping=aes(x=range, y=FBA, color=RKO_c,label=X), size=3)+
  labs(title="Reaction knockout", x ="Max Flux - Min Flux", y = "Flux") +
  theme_bw() 

#3 group
df_10 <- xdfx[xdfx$range<13999, ]
df_10 <- df_10[df_10$range>8000, ]
write.table(df_10$ID,"gene_group_10.txt",quote=F,sep="\t", row.names = F)
write.table(df_10$X,"react_group_10.txt",quote=F,sep="\t", row.names = F)
write.table(df_10[,c(2,14)],"reactF_group_10.txt",quote=F,sep="\t", row.names = F)
ggplot(df_10) +
  geom_point(aes(x=range, y=X, color=RKO_c), size=5, shape=20) +
  theme_bw() 

ggplot(df_10) +
  geom_point(aes(x=range, y=FBA, color=RKO_c), size=2) +
  geom_text_repel(mapping=aes(x=range, y=FBA, color=RKO_c,label=X), size=3)+
  labs(title="Reaction knockout", x ="Max Flux - Min Flux", y = "Flux") +
  theme_bw() 

#4 group
df_20 <- xdfx[xdfx$range>14000, ]
write.table(df_20$ID,"gene_group_20.txt",quote=F,sep="\t", row.names = F)
write.table(df_20$X,"react_group_20.txt",quote=F,sep="\t", row.names = F)
write.table(df_20[,c(2,14)],"reactF_group_20.txt",quote=F,sep="\t", row.names = F)
ggplot(df_20) +
  geom_point(aes(x=range, y=X, color=RKO_c), size=5, shape=20) +
  theme_bw() 

ggplot(df_20) +
  geom_point(aes(x=range, y=FBA, color=RKO_c), size=2) +
  geom_text_repel(mapping=aes(x=range, y=FBA, color=RKO_c,label=X), size=3)+
  labs(title="Reaction knockout", x ="Max Flux - Min Flux", y = "Flux") +
  theme_bw() 

#reactions with zero range
df_x <- xdfx[xdfx$minimum == 0, ]
df_x <- df_x[df_x$maximum == 0, ]
write.table(df_x$ID,"gene_group_zero.txt",quote=F,sep="\t", row.names = F)
write.table(df_x$X,"react_group_zero.txt",quote=F,sep="\t", row.names = F)
write.table(df_x[,c(2,14)],"reactF_group_zero.txt",quote=F,sep="\t", row.names = F)
ggplot(df_x) +
  geom_point(aes(x=range, y=X, color=PKO_c), size=5, shape=20) +
  theme_bw() 
ggplot(df_x) +
  geom_point(aes(x=range, y=FBA, color=RKO_c), size=2) +
  geom_text_repel(mapping=aes(x=range, y=FBA, color=RKO_c,label=X), size=3)+
  labs(title="Reaction knockout", x ="Max Flux - Min Flux", y = "Flux") +
  theme_bw() 



#1 group without zero range reactions
df0 <- df_0[!df_0$X %in% df_x$X,]
df0 <- df0[ which(df0$RKO_c != "NA"),]

#write.table(df0$ID,"gene_group_0-zero.txt",quote=F,sep="\t", row.names = F)
#write.table(df0$X,"react_group_0-zero.txt",quote=F,sep="\t", row.names = F)
#write.table(df0[,c(2,14)],"reactF_group_0-zero.txt",quote=F,sep="\t", row.names = F)
ggplot(df0) +
  geom_point(aes(x=range, y=X, color=RKO_c), size=5, shape=20) +
  theme_bw()

df0 <- df0[which(df0$RKO_c == 'leth' &  df0$range <0.1),]

ggplot(df0) +
  geom_point(aes(x=range, y=FBA, color=RKO_c), size=2) +
  geom_text_repel(mapping=aes(x=range, y=FBA, color=RKO_c,label=X), size=3)+
  labs(title="Reaction knockout", x ="Max Flux - Min Flux", y = "Flux") +
  theme_bw() 

df0$delta <- abs(df0$FBA) -abs(df0$range)

#__________________________________________________________________________________________
