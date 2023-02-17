#Oral wash sample workflow of the analysis

#Title: Azithromycin Exposure Induces Transient Microbial Composition Shifts and Decreasing Ability of The Airway Microbiota Resilient from PM2.5 Stress in Healthy Adults: A Randomized, Double-blind, Placebo-controlled Trial

#Authors: Sisi Du

#Date:2022/9/9

##Figure S11 The comparison between the oral wash samples and sputum samples bacterial DNA burden  
####DNA in oral wash samples
#D0-D4
wilcox.test(DNAOWD0D14$D0,DNAOWD0D14$D4,paired = T)
#D0-D14
wilcox.test(DNAOWD0D14$D0,DNAOWD0D14$D14,paired = T)
#D4-D14
wilcox.test(DNAOWD0D14$D4,DNAOWD0D14$D14,paired = T)
#D4-D30
wilcox.test(DNAOWD0D60$D4,DNAOWD0D60$D30,paired = T)
#D4-D60
wilcox.test(DNAOWD0D60$D4,DNAOWD0D60$D60,paired = T)
#D14-D30
wilcox.test(DNAOWD0D60$D14,DNAOWD0D60$D30,paired = T)
#D14-D60
wilcox.test(DNAOWD0D60$D14,DNAOWD0D60$D60,paired = T)
#D0-D30
wilcox.test(DNAOWD0D60$D0,DNAOWD0D60$D30,paired = T)
#D0-D60
wilcox.test(DNAOWD0D60$D0,DNAOWD0D60$D60,paired = T)
#D30-D60
wilcox.test(DNAOWD0D60$D30,DNAOWD0D60$D60,paired = T)
p.adjust(c(0.5271,0.02913,0.4389,0.2455,0.33,0.6742,0.6215,0.00169,0.009436,0.3488),n=10,method = 'bonferroni')

####DNA OW and Sputum
#D0
wilcox.test(DNAOSD0D14$D0S,DNAOSD0D14$D0O,paired = T)
#D4
wilcox.test(DNAOSD0D14$D4S,DNAOSD0D14$D4O,paired = T)
#D14
wilcox.test(DNAOSD0D14$D14S,DNAOSD0D14$D14O,paired = T)
#D30
wilcox.test(DNAOSD0D60$D30S,DNAOSD0D60$D30O,paired = T)
#D60
wilcox.test(DNAOSD0D60$D60S,DNAOSD0D60$D60O,paired = T)

p.adjust(c(0.003901,2.98e-06,0.0002052,9.537e-06,0.0003948),n=5,method = 'bonferroni')

library(ggplot2)
ggplot(data =DNApointc,aes(x=OR,y=DNAme,fill=Nith) )+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=SE1, ymax=SE2),width=0.5)+
  geom_jitter(data=DNApoint, mapping=aes(x=OR,y=DNA),shape=1)+
  scale_y_continuous(limits=c(0,2e+8),expand=c(0,0))+
  scale_x_continuous(limits = c(0,14),breaks = seq(0,14,by=1))+
  scale_fill_manual(values=c('#CC3300','#FF9900'))+
  theme_classic()

###Minimal influence of procedural contaminations on oral wash samples sequencing
#Figure S12a PCoA plot based on unweighted UniFrac distance between sputum samples and negative control samples.
#PERMANOVA between sputum samples and negative control samples.
library(vegan)
library(ggplot2)
library(permute)
library(lattice)
path="D:/myfunction"
setwd(path)  
source("D:/myfunction/myfuction2.R")
library(lmPerm)
library(ape)
source("D:/myfunction/myfunction.R")
otu_tree <- read.tree("D:/sputum/otus.tree")
####
adonisdata <-read.csv("D:/OW/adonisam.csv",header = T,row.names = 1)
adonisgroup =read.table("D:/OW/adonisgroup.txt", header=T, row.names= 1,sep="\t")
adonis5(adonisdata~SampleControl,data=adonisgroup,permutations = 999, method="du",tree=otu_tree) -> au
au
#PCoA plot
unifracpcoasam <- vegdist2(adonisdata,method = 'du',tree =otu_tree )
unifracpcoasam1 <- as.matrix(unifracpcoasam)
idx =rownames(adonisgroup) %in% colnames(unifracpcoasam1)
sub_design =adonisgroup[idx,]
unifracpcoasam1 =unifracpcoasam1[rownames(sub_design), rownames(sub_design)]
pcoasam =cmdscale(unifracpcoasam1, k=2, eig=T)
pointsams = as.data.frame(pcoasam$points) 
eigsam = pcoasam$eig
levels(sub_design$SampleControl)=c("neg","sam")
pointsams = cbind(pointsams, sub_design$SampleControl)
colnames(pointsams) = c("PC1", "PC2","SampleControl") 
p = ggplot(pointsams, aes(x=PC1, y=PC2, color=SampleControl)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eigsam[1] / sum(eigsam), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigsam[2] / sum(eigsam), digits=4), "%)", sep=""),
       title="UniFrac distance PCoA",color='Samplecontrol') + theme_classic()+
  scale_color_manual(values = c('#003399','#FF3333'))+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15),legend.key.size = unit(1,'cm'))

p + stat_ellipse(level=0.68,size=1)

#Figure S12b PCoA plot based on weighted UniFrac distance between sputum samples and negative control samples.
#PERMANOVA between sputum samples and negative control samples.
adonis5(adonisdata~SampleControl,data=adonisgroup,permutations = 999, method="dw",tree=otu_tree) -> aw
aw
#PCoA plot
wunifracpcoasam <- vegdist2(adonisdata,method = 'dw',tree =otu_tree )
wunifracpcoasam1 <- as.matrix(wunifracpcoasam)
idx =rownames(adonisgroup) %in% colnames(wunifracpcoasam1)
sub_design =adonisgroup[idx,]
wunifracpcoasam1 =wunifracpcoasam1[rownames(sub_design), rownames(sub_design)]
pcoasam1 =cmdscale(wunifracpcoasam1, k=2, eig=T)
pointsam1s = as.data.frame(pcoasam1$points) 
eigsam1 = pcoasam1$eig
levels(sub_design$SampleControl)=c("neg","sam")
pointsam1s = cbind(pointsam1s, sub_design$SampleControl)
colnames(pointsam1s) = c("PC1", "PC2","SampleControl") 
p = ggplot(pointsam1s, aes(x=PC1, y=PC2, color=SampleControl)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eigsam1[1] / sum(eigsam1), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigsam1[2] / sum(eigsam1), digits=4), "%)", sep=""),
       title="Weight UniFrac distance PCoA",color='Samplecontrol') + theme_classic()+
  scale_color_manual(values = c('#003399','#FF3333'))+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15),legend.key.size = unit(1,'cm'))

p + stat_ellipse(level=0.68,size=1) 

#Figure S12c Top 15 ZOTUs are ranked in descending order of mean relative abundance of oral wash samples.
otua <-read.table("D:/OW/otua.txt", header=T,row.names= 1, sep="\t")
otua_relative <- otua/rowSums(otua)
write.table(otua_relative,file='D:/OW/otua_relative.txt',quote=FALSE,sep='\t',row.names=T,col.names = T)

tapply(X = sdsamow$OTU_1,INDEX = sdsamow$Group,FUN = sd)
tapply(X = sdsamow$OTU_2,INDEX = sdsamow$Group,FUN = sd )
tapply(X = sdsamow$OTU_3,INDEX = sdsamow$Group,FUN = sd )
tapply(X = sdsamow$OTU_5,INDEX = sdsamow$Group,FUN = sd )
tapply(X = sdsamow$OTU_6,INDEX = sdsamow$Group,FUN = sd )
tapply(X = sdsamow$OTU_10,INDEX = sdsamow$Group,FUN = sd )
tapply(X = sdsamow$OTU_8,INDEX = sdsamow$Group,FUN = sd )
tapply(X = sdsamow$OTU_7,INDEX = sdsamow$Group,FUN = sd )
tapply(X = sdsamow$OTU_11,INDEX = sdsamow$Group,FUN = sd )
tapply(X = sdsamow$OTU_9,INDEX = sdsamow$Group,FUN = sd )
tapply(X = sdsamow$OTU_12,INDEX = sdsamow$Group,FUN = sd )
tapply(X = sdsamow$OTU_13,INDEX = sdsamow$Group,FUN = sd )
tapply(X = sdsamow$OTU_14,INDEX = sdsamow$Group,FUN = sd )
tapply(X = sdsamow$OTU_17,INDEX = sdsamow$Group,FUN = sd )
tapply(X = sdsamow$OTU_15,INDEX = sdsamow$Group,FUN = sd )

samconow$Sam <- factor(samconow$Sam)
samconow$Sam <- relevel(samconow$Sam,ref = "Sputum samples")
ggplot(samconow,aes(OR,R,fill=Phylum,ymin=R,ymax=R+SD))+
  geom_bar(stat = 'identity',width = 0.45,position = 'identity')+
  geom_errorbar(width = 0.45,color='black')+
  scale_x_continuous(limits=c(0,8),breaks=seq(0.5,7.5,by=0.5),expand=c(0,0))+
  facet_grid(Sam~.)+
  theme_classic()+
  scale_y_continuous(labels = scales::percent,limits=c(0,0.20),expand=c(0,0))+
  theme(axis.text.x=element_text(size=6,angle=45,color="black",hjust = 1))+
  ylab("Relative abundance (%of all sequences)")+
  xlab('')+
  theme(axis.text.y=element_text(size=6), axis.title.y=element_text(size=6))+
  theme(strip.text = element_text(colour = 'white'), strip.background = element_rect(fill = 'white', colour = 'white'))+
  theme(legend.text = element_text(size=6),legend.title = element_text(size = 6))+
  theme(panel.spacing = unit(2,'lines'))+
  scale_fill_manual(values = c('#FF3333','#330000','#660000','#990000'))

#Figure S12d  Top 15 ZOTUs are ranked in descending order of mean relative abundance of control samples.
tapply(X = sdnegow$OTU_4,INDEX = sdnegow$Group,FUN = sd)
tapply(X = sdnegow$OTU_32,INDEX = sdnegow$Group,FUN = sd )
tapply(X = sdnegow$OTU_41,INDEX = sdnegow$Group,FUN = sd )
tapply(X = sdnegow$OTU_54,INDEX = sdnegow$Group,FUN = sd )
tapply(X = sdnegow$OTU_107,INDEX = sdnegow$Group,FUN = sd )
tapply(X = sdnegow$OTU_64,INDEX = sdnegow$Group,FUN = sd )
tapply(X = sdnegow$OTU_131,INDEX = sdnegow$Group,FUN = sd )
tapply(X = sdnegow$OTU_139,INDEX = sdnegow$Group,FUN = sd )
tapply(X = sdnegow$OTU_1,INDEX = sdnegow$Group,FUN = sd )
tapply(X = sdnegow$OTU_414,INDEX = sdnegow$Group,FUN = sd )
tapply(X = sdnegow$OTU_592,INDEX = sdnegow$Group,FUN = sd )
tapply(X = sdnegow$OTU_3,INDEX = sdnegow$Group,FUN = sd )
tapply(X = sdnegow$OTU_926,INDEX = sdnegow$Group,FUN = sd )
tapply(X = sdnegow$OTU_12,INDEX = sdnegow$Group,FUN = sd )
tapply(X = sdnegow$OTU_851,INDEX = sdnegow$Group,FUN = sd )

ggplot(consamow,aes(OR,R,fill=Phylum,ymin=R,ymax=R+SD))+
  geom_bar(stat = 'identity',width = 0.45,position = 'identity')+
  geom_errorbar(width = 0.45,color='black')+
  scale_x_continuous(limits=c(0,8),breaks=seq(0.5,7.5,by=0.5),expand=c(0,0))+
  facet_grid(Sam~.)+
  theme_classic()+
  scale_y_continuous(labels = scales::percent,limits=c(0,1.2),expand=c(0,0))+
  theme(axis.text.x=element_text(size=6,angle=45,color="black",hjust = 1))+
  ylab("Relative abundance (%of all sequences)")+
  xlab('')+
  theme(axis.text.y=element_text(size=6), axis.title.y=element_text(size=6))+
  theme(strip.text = element_text(colour = 'white'), strip.background = element_rect(fill = 'white', colour = 'white'))+
  theme(legend.text = element_text(size=6),legend.title = element_text(size = 6))+
  theme(panel.spacing = unit(2,'lines'))+
  scale_fill_manual(values = c('#FF9966','#330000','#990000'))

#Figure S12e-f Decontam package confirms contaminants ZOTUs in our data using both the frequency and prevalence methods with threshold 0.1 and 0.5 respectively.
install.packages("phyloseqGraphTest")
BiocManager::install("phyloseq",ask = F,update = F)
BiocManager::install("decontam",ask = F,update = F)
library(phyloseq)
library(decontam)
otudec<-read.csv("D:/OW/otu_table1.csv",header = T,row.names = 1)
meta<-read.csv("D:/OW/sampledata1.csv",header = T,row.names = 1)
OTUDEC <-otu_table(otudec,taxa_are_rows = F)
META <-sample_data(meta)
phy <-phyloseq(OTUDEC,META)
sample_data(phy)$is.neg <- sample_data(phy)$SampleControl== "neg"
contamdf.both <- isContaminant(phy, conc="Quan", neg="is.neg", method="both", threshold=c(0.1,0.5))
table(contamdf.both$contaminant)
which(contamdf.both$contaminant)
head(contamdf.both)
write.table(contamdf.both,file='D:/OW/contamdf.both.txt',quote=FALSE,sep='\t',row.names=F,col.names = T)
#Figure S12e
plot_frequency(phy, taxa_names(phy)[c(78,96,423,483,1186,1241,1658,1664,1)], conc="Quan") + 
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")+
  theme_classic()+
  theme(axis.text.y=element_text(size=15), axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15), axis.title.x=element_text(size=15))+
  theme(strip.text = element_text(colour = 'black',size = 15))

#Figure S12f
phy.pa <- transform_sample_counts(phy, function(abund) 1*(abund>0))
phy.pa.neg <- prune_samples(sample_data(phy.pa)$SampleControl== "neg", phy.pa)
phy.pa.pos <- prune_samples(sample_data(phy.pa)$SampleControl== "sam", phy.pa)
df.pb <- data.frame(pa.pos=taxa_sums(phy.pa.pos), pa.neg=taxa_sums(phy.pa.neg),contaminant=contamdf.both$contaminant)
write.table(df.pb,file='D:/OW/df.pb.txt',quote=FALSE,sep='\t',row.names=F,col.names = T)
df.pb =read.table("D:/OW/df.pb.txt", header=T, row.names= 1,sep="\t")

View(df.pb)
ggplot(data=df_pb, aes(x=pa.neg, y=pa.pos)) + 
  geom_jitter(aes(size=contamdf.both$freq,color=contaminant),alpha=1/3)+
  geom_text(aes(y = pa.pos+5, label = Name))+
  scale_size_continuous(name='Relative abundance',breaks=c(0.000,0.025,0.050,0.075),labels=c('0%','2.5%','5%','7.5%'))+
  theme_classic()+
  xlab("Prevalence (Control samples)") + ylab("Prevalence (Sputum samples)")+
  theme(axis.text.y=element_text(size=15), axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15), axis.title.x=element_text(size=15))+
  scale_color_manual(values = c('#FF3333','#330000'))+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15))

#The shifts in oral wash phylogenetic diversity after antibiotic exposure
#Alpha diversity
#Figure S2a-b Richness and Shannon index in oral cavity microbiota
## Alpha diversity calculation
library(picante) 
alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')    
  Pielou <- Shannon / log(Richness, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
  }
  result
}
tree <- read.tree('D:/sputum/otus.tree')
otu <- read.delim('D:/OW/otutab_norm.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- t(otu)
alpha_all <- alpha(otu, base = 2)
alpha_allt <- alpha(otu, tree, base = 2)
write.csv(alpha_allt, 'D:/OW/otu_alphat.csv', quote = FALSE)
#Difference test within oral wash samples between the 5 timepoints
###Richness
sd(pairzeroD0D4D14$D0Richness)
#D0 D4
shapiro.test(pairzeroD0D4D14$D0D4richness)
t.test(pairzeroD0D4D14$D0Richness,pairzeroD0D4D14$D4Richness,paired=T)
#D0 D14
shapiro.test(pairzeroD0D4D14$D0D14richness)
t.test(pairzeroD0D4D14$D0Richness,pairzeroD0D4D14$D14Richness,paired=T)
#D4 D14
shapiro.test(pairzeroD0D4D14$D4D14richness)
t.test(pairzeroD0D4D14$D4Richness,pairzeroD0D4D14$D14Richness,paired=T)
#D0 D30
shapiro.test(pairzeroD0D4D14D30D60$D0D30richness)
t.test(pairzeroD0D4D14D30D60$D0Richness,pairzeroD0D4D14D30D60$D30Richness,paired = T)
#D0 D60
shapiro.test(pairzeroD0D4D14D30D60$D0D60richness)
t.test(pairzeroD0D4D14D30D60$D0Richness,pairzeroD0D4D14D30D60$D60Richness,paired = T)
#D4 D30
shapiro.test(pairzeroD0D4D14D30D60$D4D30richness)
t.test(pairzeroD0D4D14D30D60$D4Richness,pairzeroD0D4D14D30D60$D30Richness,paired = T)
#D4 D60
shapiro.test(pairzeroD0D4D14D30D60$D4D60richness)
t.test(pairzeroD0D4D14D30D60$D4Richness,pairzeroD0D4D14D30D60$D60Richness,paired = T)
#D14 D30
shapiro.test(pairzeroD0D4D14D30D60$D14D30richness)
t.test(pairzeroD0D4D14D30D60$D14Richness,pairzeroD0D4D14D30D60$D30Richness,paired = T)
#D14 D60
shapiro.test(pairzeroD0D4D14D30D60$D14D60richness)
t.test(pairzeroD0D4D14D30D60$D14Richness,pairzeroD0D4D14D30D60$D60Richness,paired = T)
#D30 D60
shapiro.test(pairzeroD0D4D14D30D60$D30D60richness)
wilcox.test(pairzeroD0D4D14D30D60$D30Richness,pairzeroD0D4D14D30D60$D60Richness,paired = T)

p.adjust(c(7.077e-07,1.036e-08,0.1924,4.767e-06,0.001468,0.06951,0.0001461,0.0003804,3.046e-08,0.0001114),method = 'bonferroni',n=10)

###Shannon
#D0 D4
shapiro.test(pairzeroD0D4D14$D0D4shannon)
t.test(pairzeroD0D4D14$D0Shannon,pairzeroD0D4D14$D4Shannon,paired=T)
#D0 D14
shapiro.test(pairzeroD0D4D14$D0D14shannon)
t.test(pairzeroD0D4D14$D0Shannon,pairzeroD0D4D14$D14Shannon,paired=T)
#D4 D14
shapiro.test(pairzeroD0D4D14$D4D14shannon)
t.test(pairzeroD0D4D14$D4Shannon,pairzeroD0D4D14$D14Shannon,paired=T)
#D0 D30
shapiro.test(pairzeroD0D4D14D30D60$D0D30shannon)
t.test(pairzeroD0D4D14D30D60$D0Shannon,pairzeroD0D4D14D30D60$D30Shannon,paired = T)
#D0 D60
shapiro.test(pairzeroD0D4D14D30D60$D0D60shannon)
t.test(pairzeroD0D4D14D30D60$D0Shannon,pairzeroD0D4D14D30D60$D60Shannon,paired = T)
#D4 D30
shapiro.test(pairzeroD0D4D14D30D60$D4D30shannon)
t.test(pairzeroD0D4D14D30D60$D4Shannon,pairzeroD0D4D14D30D60$D30Shannon,paired = T)
#D4 D60
shapiro.test(pairzeroD0D4D14D30D60$D4D60shannon)
t.test(pairzeroD0D4D14D30D60$D4Shannon,pairzeroD0D4D14D30D60$D60Shannon,paired = T)
#D14 D30
shapiro.test(pairzeroD0D4D14D30D60$D14D30shannon)
t.test(pairzeroD0D4D14D30D60$D14Shannon,pairzeroD0D4D14D30D60$D30Shannon,paired = T)
#D14 D60
shapiro.test(pairzeroD0D4D14D30D60$D14D60shannon)
t.test(pairzeroD0D4D14D30D60$D14Shannon,pairzeroD0D4D14D30D60$D60Shannon,paired = T)
#D30 D60
shapiro.test(pairzeroD0D4D14D30D60$D30D60shannon)
t.test(pairzeroD0D4D14D30D60$D30Shannon,pairzeroD0D4D14D30D60$D60Shannon,paired = T)
p.adjust(c(1.679e-07,0.0003327,0.03564,0.1411,0.9179,0.000122,9.05e-06,0.004822,0.0001311,0.04247),method = 'bonferroni',n=10)

#Difference test
#one
###RICHNESS
#D0 D4
shapiro.test(paironeD0D4D14$D0D4richness)
t.test(paironeD0D4D14$D0Richness,paironeD0D4D14$D4Richness,paired=T)
#D0 D14
shapiro.test(paironeD0D4D14$D0D14richness)
t.test(paironeD0D4D14$D0Richness,paironeD0D4D14$D14Richness,paired=T)
#D4 D14
shapiro.test(paironeD0D4D14$D4D14richness)
t.test(paironeD0D4D14$D4Richness,paironeD0D4D14$D14Richness,paired=T)
#D0 D30
shapiro.test(paironeD0D4D14D30$D0D30richness)
t.test(paironeD0D4D14D30$D0Richness,paironeD0D4D14D30$D30Richness,paired = T)
#D4 D30
shapiro.test(paironeD0D4D14D30$D4D30richness)
t.test(paironeD0D4D14D30$D4Richness,paironeD0D4D14D30$D30Richness,paired = T)
#D14 D30
shapiro.test(paironeD0D4D14D30$D14D30richness)
wilcox.test(paironeD0D4D14D30$D14Richness,paironeD0D4D14D30$D30Richness,paired = T)
#D0 D60
shapiro.test(paironeD0D4D14D30D60$D0D60richness)
t.test(paironeD0D4D14D30D60$D0Richness,paironeD0D4D14D30D60$D60Richness,paired = T)
#D4 D60
shapiro.test(paironeD0D4D14D30D60$D4D60richness)
t.test(paironeD0D4D14D30D60$D4Richness,paironeD0D4D14D30D60$D60Richness,paired = T)
#D14 D60
shapiro.test(paironeD0D4D14D30D60$D14D60richness)
t.test(paironeD0D4D14D30D60$D14Richness,paironeD0D4D14D30D60$D60Richness,paired = T)
#D30 D60
shapiro.test(paironeD0D4D14D30D60$D30D60richness)
t.test(paironeD0D4D14D30D60$D30Richness,paironeD0D4D14D30D60$D60Richness,paired = T)

p.adjust(c(0.3679,0.2409,0.6921,0.6083,0.9102,0.6705,0.3025,0.2946,0.09726,0.2449),method = 'bonferroni',n=10)

###shannon
#D0 D4
shapiro.test(paironeD0D4D14$D0D4shannon)
t.test(paironeD0D4D14$D0Shannon,paironeD0D4D14$D4Shannon,paired=T)
#D0 D14
shapiro.test(paironeD0D4D14$D0D14shannon)
wilcox.test(paironeD0D4D14$D0Shannon,paironeD0D4D14$D14Shannon,paired=T)
#D4 D14
shapiro.test(paironeD0D4D14$D4D14shannon)
t.test(paironeD0D4D14$D4Shannon,paironeD0D4D14$D14Shannon,paired=T)
#D0 D30
shapiro.test(paironeD0D4D14D30$D0D30shannon)
t.test(paironeD0D4D14D30$D0Shannon,paironeD0D4D14D30$D30Shannon,paired = T)
#D4 D30
shapiro.test(paironeD0D4D14D30$D4D30shannon)
t.test(paironeD0D4D14D30$D4Shannon,paironeD0D4D14D30$D30Shannon,paired = T)
#D14 D30
shapiro.test(paironeD0D4D14D30$D14D30shannon)
t.test(paironeD0D4D14D30$D14Shannon,paironeD0D4D14D30$D30Shannon,paired = T)
#D0 D60
shapiro.test(paironeD0D4D14D30D60$D0D60shannon)
t.test(paironeD0D4D14D30D60$D0Shannon,paironeD0D4D14D30D60$D60Shannon,paired = T)
#D4 D60
shapiro.test(paironeD0D4D14D30D60$D4D60shannon)
t.test(paironeD0D4D14D30D60$D4Shannon,paironeD0D4D14D30D60$D60Shannon,paired = T)
#D14 D60
shapiro.test(paironeD0D4D14D30D60$D14D60shannon)
t.test(paironeD0D4D14D30D60$D14Shannon,paironeD0D4D14D30D60$D60Shannon,paired = T)
#D30 D60
shapiro.test(paironeD0D4D14D30D60$D30D60shannon)
t.test(paironeD0D4D14D30D60$D30Shannon,paironeD0D4D14D30D60$D60Shannon,paired = T)

p.adjust(c(0.3126,0.06342, 0.2804,0.6669,0.838,0.3518,0.921,0.5311,0.2917,0.7984),method = 'bonferroni',n=10)

#Difference test between the oral wash samples and sputum samples
###ZERO
###Richness
t.test(OWtestD0D14$D0OWRichness,OWtestD0D14$D0SRichness,paired=T)
t.test(OWtestD0D14$D4OWRichness,OWtestD0D14$D4SRichness,paired=T)
t.test(OWtestD0D14$D14OWRichness,OWtestD0D14$D14SRichness,paired=T)
t.test(OWtestD0D60$D30OWRichness,OWtestD0D60$D30SRichness,paired=T)
t.test(OWtestD0D60$D60OWRichness,OWtestD0D60$D60SRichness,paired=T)
p.adjust(c(0.000538,0.01133,0.03322,0.04879,0.009991),method = 'bonferroni',n=5)

###Shannon
t.test(OWtestD0D14$D0OWShannon,OWtestD0D14$D0SShannon,paired=T)
t.test(OWtestD0D14$D4OWShannon,OWtestD0D14$D4SShannon,paired=T)
t.test(OWtestD0D14$D14OWShannon,OWtestD0D14$D14SShannon,paired=T)
t.test(OWtestD0D60$D30OWShannon,OWtestD0D60$D30SShannon,paired=T)
t.test(OWtestD0D60$D60OWShannon,OWtestD0D60$D60SShannon,paired=T)
p.adjust(c(0.4409,0.7454,0.7173,0.08832,0.2243),method = 'bonferroni',n=5)

###Figure S14
###Richness
t.test(OWtestD0D14one$D0OWRichness,OWtestD0D14one$D0SRichness,paired=T)
t.test(OWtestD0D14one$D4OWRichness,OWtestD0D14one$D4SRichness,paired=T)
t.test(OWtestD0D14one$D14OWRichness,OWtestD0D14one$D14SRichness,paired=T)
t.test(OWtestD0D30one$D30OWRichness,OWtestD0D30one$D30SRichness,paired=T)
t.test(OWtestD0D60one$D60OWRichness,OWtestD0D60one$D60SRichness,paired=T)
p.adjust(c(6.007e-07,9.623e-06,2.311e-05,2.241e-05,0.0002262),method = 'bonferroni',n=5)

###Shannon
t.test(OWtestD0D14one$D0OWShannon,OWtestD0D14one$D0SShannon,paired=T)
t.test(OWtestD0D14one$D4OWShannon,OWtestD0D14one$D4SShannon,paired=T)
t.test(OWtestD0D14one$D14OWShannon,OWtestD0D14one$D14SShannon,paired=T)
t.test(OWtestD0D30one$D30OWShannon,OWtestD0D30one$D30SShannon,paired=T)
t.test(OWtestD0D60one$D60OWShannon,OWtestD0D60one$D60SShannon,paired=T)
p.adjust(c(0.008403,0.01022,0.001705,0.004433,0.0006031),method = 'bonferroni',n=5)

tapply(X=OWap$Richness,INDEX = OWap$Date1,FUN = mean)
tapply(X=OWap$Richness,INDEX = OWap$Date1,FUN = sd)
tapply(X=OWap$Shannon,INDEX = OWap$Date1,FUN = mean)
tapply(X=OWap$Shannon,INDEX = OWap$Date1,FUN = sd)

###ZERO
###Richness plot
ggplot()+
  geom_jitter(data=OWap, mapping=aes(x=OR,y=Richness,shape=Nith,color=Date1),size=5)+
  scale_y_continuous(limits=c(0,800),expand=c(0,0))+
  scale_color_manual(values=c('#003399','#003399','#FF3333','#FF3333','#FF9966','#FF9966','#CC6600','#CC6600','#330000','#330000'))+
  geom_crossbar(data=OWapt,aes(x=OR,y=richmea,ymin=richmea-Sdrich, ymax=richmea+Sdrich),width=0.3)+
  scale_x_continuous(limits = c(0,12),breaks = seq(0,12,by=1),expand = c(0,0))+
  theme_classic()

###Shannon plot
ggplot()+
  geom_jitter(data=OWap, mapping=aes(x=OR,y=Shannon,shape=Nith,color=Date1),size=5)+
  scale_y_continuous(limits=c(0,8),expand=c(0,0))+
  scale_color_manual(values=c('#003399','#003399','#FF3333','#FF3333','#FF9966','#FF9966','#CC6600','#CC6600','#330000','#330000'))+
  geom_crossbar(data=OWapt,aes(x=OR,y=shamea,ymin=shamea-Sdsha, ymax=shamea+Sdsha),width=0.3)+
  scale_x_continuous(limits = c(0,12),breaks = seq(0,12,by=1),expand = c(0,0))+
  theme_classic()

tapply(X=OWapone$Richness,INDEX = OWapone$Date1,FUN = mean)
tapply(X=OWapone$Richness,INDEX = OWapone$Date1,FUN = sd)
tapply(X=OWapone$Shannon,INDEX = OWapone$Date1,FUN = mean)
tapply(X=OWapone$Shannon,INDEX = OWapone$Date1,FUN = sd)

###ONE
###Richness plot
ggplot()+
  geom_jitter(data=OWapone, mapping=aes(x=OR,y=Richness,shape=Nith,color=Date1),size=5)+
  scale_y_continuous(limits=c(0,800),expand=c(0,0))+
  scale_color_manual(values=c('#003399','#003399','#FF3333','#FF3333','#FF9966','#FF9966','#CC6600','#CC6600','#330000','#330000'))+
  geom_crossbar(data=OWaptone,aes(x=OR,y=richmea,ymin=richmea-Sdrich, ymax=richmea+Sdrich),width=0.3)+
  scale_x_continuous(limits = c(0,12),breaks = seq(0,12,by=1),expand = c(0,0))+
  theme_classic()

###Shannon plot
ggplot()+
  geom_jitter(data=OWapone, mapping=aes(x=OR,y=Shannon,shape=Nith,color=Date1),size=5)+
  scale_y_continuous(limits=c(0,8),expand=c(0,0))+
  scale_color_manual(values=c('#003399','#003399','#FF3333','#FF3333','#FF9966','#FF9966','#CC6600','#CC6600','#330000','#330000'))+
  geom_crossbar(data=OWaptone,aes(x=OR,y=shamea,ymin=shamea-Sdsha, ymax=shamea+Sdsha),width=0.3)+
  scale_x_continuous(limits = c(0,12),breaks = seq(0,12,by=1),expand = c(0,0))+
  theme_classic()

adoniszero <- read.csv("D:/OW/adoniszero.csv",header = T,row.names = 1)
samplezero =read.table("D:/OW/samplezero.txt", header=T, row.names= 1,sep="\t")
adonis5(adoniszero~Date1,data=samplezero,permutations = 999, method="du",tree=otu_tree) -> adt
adt
a <- pairwise.adonis5(adoniszero,samplezero$Date1, sim.method ="dw",tree=otu_tree, p.adjust.m= "bonferroni")

#Beta diversity
##Figure S2c PCoA plot based on unweighted UniFrac distance in oral wash samples in azithromycin group
otu_tree <- read.tree("D:/sputum/otus.tree")
adoniszero <- read.csv("D:/OW/adoniszero.csv",header = T,row.names = 1)
samplezero =read.table("D:/OW/samplezero.txt", header=T, row.names= 1,sep="\t")
adonis5(adoniszero~Date1,data=samplezero,permutations = 999, method="du",tree=otu_tree) -> adt
adt
a <- pairwise.adonis5(adoniszero,samplezero$Date1, sim.method ="du",tree=otu_tree, p.adjust.m= "bonferroni")

unifracpcoa <- vegdist2(adoniszero,method = 'du',tree =otu_tree )
unifracpcoa1 <- as.matrix(unifracpcoa)
idx =rownames(samplezero) %in% colnames(unifracpcoa1)
sub_design =samplezero[idx,]
unifracpcoa1 =unifracpcoa1[rownames(sub_design), rownames(sub_design)]
pcoa =cmdscale(unifracpcoa1, k=2, eig=T)
points = as.data.frame(pcoa$points) 
eig = pcoa$eig
levels(sub_design$Date1)=c("D0","D14","D30","D4","D60")
points = cbind(points, sub_design$Date1)
colnames(points) = c("PC1", "PC2","Date1") 
p = ggplot(points, aes(x=PC1, y=PC2, color=Date1)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title="UniFrac distance PCoA",color='Timepoints') + theme_classic()+
  scale_color_manual(values = c('#003399','#FF3333','#FF9966','#CC6600','#330000'))+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15),legend.key.size = unit(1,'cm'))

p + stat_ellipse(level=0.68,size=1)

#Figure S14c PCoA plot based on unweighted UniFrac distance in placebo group
adonisone <- read.csv("D:/sputum/adonisone.csv",header = T,row.names = 1)
sampleone =read.table("D:/sputum/sampleone.txt", header=T, row.names= 1,sep="\t")
adonis5(adonisone~Date1,data=sampleone,permutations = 999, method="du",tree=otu_tree) -> adl
adl
pairwise.adonis5(adonisone,sampleone$Date1, sim.method ="du",tree=otu_tree, p.adjust.m= "bonferroni")

unifracpcoat <- vegdist2(adonisone,method = 'du',tree =otu_tree )
unifracpcoa1t <- as.matrix(unifracpcoat)
idx =rownames(sampleone) %in% colnames(unifracpcoa1t)
sub_designt =sampleone[idx,]
unifracpcoa1t =unifracpcoa1t[rownames(sub_designt), rownames(sub_designt)]
pcoat =cmdscale(unifracpcoa1t, k=2, eig=T)
pointst = as.data.frame(pcoat$points) 
eigt = pcoat$eig
levels(sub_designt$Date1)=c("D0","D14","D30","D4","D60")
pointst = cbind(pointst, sub_designt$Date1)
colnames(pointst) = c("PC1", "PC2","Date1") 
p = ggplot(pointst, aes(x=PC1, y=PC2, color=Date1)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eigt[1] / sum(eigt), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigt[2] / sum(eigt), digits=4), "%)", sep=""),
       title="UniFrac distance PCoA",color='Timepoints') + theme_classic()+
  scale_color_manual(values = c('#003399','#FF3333','#FF9966','#CC6600','#330000'))+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15),legend.key.size = unit(1,'cm'))

p + stat_ellipse(level=0.68,size=1)

#Figure S2d PCoA plot based on weighted UniFrac distance in azithromycin group
adonis5(adoniszero~Date1,data=samplezero,permutations = 999, method="dw",tree=otu_tree) -> adt
adt
pairwise.adonis5(adoniszero,samplezero$Date1, sim.method ="dw",tree=otu_tree, p.adjust.m= "bonferroni")

wunifracpcoa <- vegdist2(adoniszero,method = 'dw',tree =otu_tree )
wunifracpcoa1 <- as.matrix(wunifracpcoa)
idx =rownames(samplezero) %in% colnames(wunifracpcoa1)
sub_design =samplezero[idx,]
wunifracpcoa1 =wunifracpcoa1[rownames(sub_design), rownames(sub_design)]
pcoa =cmdscale(wunifracpcoa1, k=2, eig=T)
points = as.data.frame(pcoa$points) 
eig = pcoa$eig
levels(sub_design$Date1)=c("D0","D14","D30","D4","D60")
points = cbind(points, sub_design$Date1)
colnames(points) = c("PC1", "PC2","Date1") 
p = ggplot(points, aes(x=PC1, y=PC2, color=Date1)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title="Weight UniFrac distance PCoA",color='Timepoints') + theme_classic()+
  scale_color_manual(values = c('#003399','#FF3333','#FF9966','#CC6600','#330000'))+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15),legend.key.size = unit(1,'cm'))

p + stat_ellipse(level=0.68,size=1)

#Figure S14d PCoA plot based on weighted UniFrac distance in placebo group
adonis5(adonisone~Date1,data=sampleone,permutations = 999, method="dw",tree=otu_tree) -> adt
adt
pairwise.adonis5(adonisone,sampleone$Date1, sim.method ="dw",tree=otu_tree, p.adjust.m= "bonferroni")
wunifracpcoat <- vegdist2(adonisone,method = 'dw',tree =otu_tree )
wunifracpcoa1t <- as.matrix(wunifracpcoat)
idx =rownames(sampleone) %in% colnames(wunifracpcoa1t)
sub_designt =sampleone[idx,]
wunifracpcoa1t =wunifracpcoa1t[rownames(sub_designt), rownames(sub_designt)]
pcoat =cmdscale(wunifracpcoa1t, k=2, eig=T)
pointst = as.data.frame(pcoat$points) 
eigt = pcoat$eig
levels(sub_designt$Date1)=c("D0","D14","D30","D4","D60")
pointst = cbind(pointst, sub_designt$Date1)
colnames(pointst) = c("PC1", "PC2","Date1") 
p = ggplot(pointst, aes(x=PC1, y=PC2, color=Date1)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eigt[1] / sum(eigt), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigt[2] / sum(eigt), digits=4), "%)", sep=""),
       title="Weight UniFrac distance PCoA",color='Timepoints') + theme_classic()+
  scale_color_manual(values = c('#003399','#FF3333','#FF9966','#CC6600','#330000'))+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15),legend.key.size = unit(1,'cm'))

p + stat_ellipse(level=0.68,size=1)




#Adonis test
##Adonis test within the oral wash samples
adonis5(adoniszero~Date1,data=samplezero,permutations = 999, method="dw",tree=otu_tree) -> adt
adt
pairwise.adonis5(adoniszero,samplezero$Date1, sim.method ="dw",tree=otu_tree, p.adjust.m= "bonferroni")

##Adonis test between the oral wash samples and sputum samples
##D0
ADOOWD0 <- read.csv("D:/OW/ADOOWD0.csv",header = T,row.names = 1)
OWgroupD0 =read.table("D:/OW/OWgroupD0.txt", header=T, row.names= 1,sep="\t")
adonis5(ADOOWD0~Nith,data=OWgroupD0,permutations = 999, method="du",tree=otu_tree) -> adt
adt

##D4
ADOOWD4 <- read.csv("D:/OW/ADOOWD4.csv",header = T,row.names = 1)
OWgroupD4 =read.table("D:/OW/OWgroupD4.txt", header=T, row.names= 1,sep="\t")
adonis5(ADOOWD4~Nith,data=OWgroupD4,permutations = 999, method="du",tree=otu_tree) -> adt
adt

##D14
ADOOWD14 <- read.csv("D:/OW/ADOOWD14.csv",header = T,row.names = 1)
OWgroupD14 =read.table("D:/OW/OWgroupD14.txt", header=T, row.names= 1,sep="\t")
adonis5(ADOOWD14~Nith,data=OWgroupD14,permutations = 999, method="du",tree=otu_tree) -> adt
adt

##D30
ADOOWD30 <- read.csv("D:/OW/ADOOWD30.csv",header = T,row.names = 1)
OWgroupD30 =read.table("D:/OW/OWgroupD30.txt", header=T, row.names= 1,sep="\t")
adonis5(ADOOWD30~Nith,data=OWgroupD30,permutations = 999, method="du",tree=otu_tree) -> adt
adt

##D60
ADOOWD60 <- read.csv("D:/OW/ADOOWD60.csv",header = T,row.names = 1)
OWgroupD60 =read.table("D:/OW/OWgroupD60.txt", header=T, row.names= 1,sep="\t")
adonis5(ADOOWD60~Nith,data=OWgroupD60,permutations = 999, method="du",tree=otu_tree) -> adt
adt


####one
##D0
ADOOWD0 <- read.csv("D:/OW/ADOOWD0one.csv",header = T,row.names = 1)
OWgroupD0 =read.table("D:/OW/OWgroupD0one.txt", header=T, row.names= 1,sep="\t")
adonis5(ADOOWD0~Nith,data=OWgroupD0,permutations = 999, method="dw",tree=otu_tree) -> adt
adt

##D4
ADOOWD4 <- read.csv("D:/OW/ADOOWD4one.csv",header = T,row.names = 1)
OWgroupD4 =read.table("D:/OW/OWgroupD4one.txt", header=T, row.names= 1,sep="\t")
adonis5(ADOOWD4~Nith,data=OWgroupD4,permutations = 999, method="dw",tree=otu_tree) -> adt
adt

##D14
ADOOWD14 <- read.csv("D:/OW/ADOOWD14one.csv",header = T,row.names = 1)
OWgroupD14 =read.table("D:/OW/OWgroupD14one.txt", header=T, row.names= 1,sep="\t")
adonis5(ADOOWD14~Nith,data=OWgroupD14,permutations = 999, method="dw",tree=otu_tree) -> adt
adt

##D30
ADOOWD30 <- read.csv("D:/OW/ADOOWD30one.csv",header = T,row.names = 1)
OWgroupD30 =read.table("D:/OW/OWgroupD30one.txt", header=T, row.names= 1,sep="\t")
adonis5(ADOOWD30~Nith,data=OWgroupD30,permutations = 999, method="dw",tree=otu_tree) -> adt
adt

##D60
ADOOWD60 <- read.csv("D:/OW/ADOOWD60one.csv",header = T,row.names = 1)
OWgroupD60 =read.table("D:/OW/OWgroupD60one.txt", header=T, row.names= 1,sep="\t")
adonis5(ADOOWD60~Nith,data=OWgroupD60,permutations = 999, method="dw",tree=otu_tree) -> adt
adt
p.adjust(c(0.028,0.126,0.044,0.002,0.001),method = 'bonferroni',n=5)

##Figure S4 Similar but not identical microbial profile between oral cavity and airway microbiota  
###Figure S12a Clustering of the oral cavity microbiota at ZOTUs according to the detection rate
otur <-read.table("D:/OW/otur.txt", header=T,row.names= 1, sep="\t")
otur_relative <- otur/rowSums(otur)
write.table(otur_relative,file='D:/OW/otur_relative.txt',quote=FALSE,sep='\t',row.names=T,col.names = T)
View(otur_relative)
mapdd = read.table("D:/OW/mapdd.txt",header=T,sep='\t',row.names=1)
mapdd <- as.matrix(mapdd)
Heatmap(mapdd,col = colorRamp2(c(0,1), c("white","#0000CC")),show_column_names=T, show_row_names = F,cluster_rows = T, cluster_columns = F,width =  unit(150, "mm"),height = unit(200,'mm'),heatmap_legend_param=list(legend_direction='horizontal'))
maprr = read.table("D:/OW/maprr.txt",header=T,sep='\t',row.names=1)
library(ggplot2)
###Figure S12b Clustering of the oral cavity microbiota at ZOTUs according to the mean relative abundance
maprr$logD0O = log2(maprr$D0O+0.00001)
maprr$logD4O = log2(maprr$D4O+0.00001)
maprr$logD14O = log2(maprr$D14O+0.00001)
maprr$logD30O = log2(maprr$D30O+0.00001)
maprr$logD60O = log2(maprr$D60O+0.00001)

maprr$logD0S = log2(maprr$D0S+0.00001)
maprr$logD4S = log2(maprr$D4S+0.00001)
maprr$logD14S = log2(maprr$D14S+0.00001)
maprr$logD30S = log2(maprr$D30S+0.00001)
maprr$logD60S = log2(maprr$D60S+0.00001)

write.table(maprr,file='D:/OW/maprr1.txt',quote=FALSE,sep='\t',row.names=F,col.names = T)
maprr1 = read.table("D:/OW/maprr1.txt",header=T,sep='\t',row.names=1)

maprr1 <- as.matrix(maprr1)
Heatmap(maprr1, col = colorRamp2(c(-12,-3), c("white","#FF0000")),show_column_names=T, show_row_names = F,cluster_rows = T, cluster_columns = F,width =  unit(150, "mm"),height = unit(200,'mm'),heatmap_legend_param=list(legend_direction='horizontal'))

##Figure S13 Similar but not identical microbial profile between oral cavity and airway microbiota  
###Figure S13a Clustering of the oral cavity microbiota at ZOTUs according to the detection rate

mapddone = read.table("D:/OW/mapddone.txt",header=T,sep='\t',row.names=1)
mapddone <- as.matrix(mapddone)
Heatmap(mapddone,col = colorRamp2(c(0,1), c("white","#0000CC")),show_column_names=T, show_row_names = F,cluster_rows = T, cluster_columns = F,width =  unit(150, "mm"),height = unit(200,'mm'),heatmap_legend_param=list(legend_direction='horizontal'))
maprrone = read.table("D:/OW/maprrone.txt",header=T,sep='\t',row.names=1)
library(ggplot2)
###Figure S13b Clustering of the oral cavity microbiota at ZOTUs according to the mean relative abundance
maprrone$logD0O = log2(maprrone$D0O+0.00001)
maprrone$logD4O = log2(maprrone$D4O+0.00001)
maprrone$logD14O = log2(maprrone$D14O+0.00001)
maprrone$logD30O = log2(maprrone$D30O+0.00001)
maprrone$logD60O = log2(maprrone$D60O+0.00001)

maprrone$logD0S = log2(maprrone$D0S+0.00001)
maprrone$logD4S = log2(maprrone$D4S+0.00001)
maprrone$logD14S = log2(maprrone$D14S+0.00001)
maprrone$logD30S = log2(maprrone$D30S+0.00001)
maprrone$logD60S = log2(maprrone$D60S+0.00001)

write.table(maprrone,file='D:/OW/maprrone1.txt',quote=FALSE,sep='\t',row.names=F,col.names = T)
maprrone1 = read.table("D:/OW/maprrone1.txt",header=T,sep='\t',row.names=1)

maprrone1 <- as.matrix(maprrone1)
Heatmap(maprrone1, col = colorRamp2(c(-14,-3), c("white","#FF0000")),show_column_names=T, show_row_names = F,cluster_rows = T, cluster_columns = F,width =  unit(150, "mm"),height = unit(200,'mm'),heatmap_legend_param=list(legend_direction='horizontal'))

###Figure S12c The distribution of the taxa of the oral cavity microbiota at family level
ggplot(topownew,aes(reorder(Date1,OR),RR,fill=Color))+
  geom_bar(stat="identity",position="stack")+
  facet_grid(.~Nich)+
  theme_classic()+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values = c('grey','#990033','#FF3399','#FF0033','#CC6699','#FF9999',"#FF00CC",'#FF6699','#FF99CC','#CC3399','#666699','#3333CC','#3399FF','#9999FF','#66CC00','#66CC99'))

###Figure S12d The distribution of the taxa of the oral cavity microbiota at family level
ggplot(topownew,aes(reorder(Date1,OR),RR,fill=Color1))+
  geom_bar(stat="identity",position="stack")+
  facet_grid(.~Nich)+
  theme_classic()+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values = c('grey','#990033','#FF3399','#FF0033','#CC6699','#FF9999',"#FF00CC",'#FF6699','#FF99CC','#CC3399','#666699','#3333CC','#3399FF','#9999FF','#66CC00','#66CC99','#339966'))

###Figure S13c The distribution of the taxa of the oral cavity microbiota at family level
ggplot(topowonenew,aes(reorder(Date1,OR),RR,fill=Color))+
  geom_bar(stat="identity",position="stack")+
  facet_grid(.~Nich)+
  theme_classic()+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values = c('grey','#990033','#FF3399','#FF0033','#CC6699','#FF9999',"#FF00CC",'#FF6699','#FF99CC','#CC3399','#666699','#3333CC','#3399FF','#9999FF','#66CC00','#66CC99'))

###Figure S13d The distribution of the taxa of the oral cavity microbiota at family level
ggplot(topowonenew,aes(reorder(Date1,OR),RR,fill=Color1))+
  geom_bar(stat="identity",position="stack")+
  facet_grid(.~Nich)+
  theme_classic()+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values = c('grey','#990033','#FF3399','#FF0033','#CC6699','#FF9999',"#FF00CC",'#FF6699','#FF99CC','#CC3399','#666699','#3333CC','#3399FF','#9999FF','#66CC00','#66CC99','#339966'))

###Difference test
median(testowsone$OWPrevotellaceae)
median(testowsone$Prevotellaceae)
IQR(testows$OWPrevotellaceae)
IQR(testows$Prevotellaceae)
median(testows$OWGemella)
median(testows$Gemella)
IQR(testows$OWVeillonellaceae)
IQR(testows$Veillonellaceae)

wilcox.test(testows$OWGranulicatella,testows$Granulicatella,paired =T)
p.adjust(c(3.04E-05,7.17E-05,1.67E-01,2.87E-01,1.47E-05,1.53E-01
,7.92E-01,7.88E-02,9.66E-15,8.94E-05,9.98E-08,1.62E-02,2.38E-03
,4.23E-08,2.68E-06),method = "bonferroni",n=15)

#Microbial taxonomic variation during 60 days’ follow-up 
#Difference test in detection rate
#McNemar-Bowker test 
#D0 VS D4
D01 <- read_excel("D:/OW/D01.xlsx")
D4 <- read_excel("D:/OW/D4.xlsx")
D14 <- read_excel("D:/OW/D14.xlsx")
D30 <- read_excel("D:/OW/D30.xlsx")
D60 <- read_excel("D:/OW/D60.xlsx")
P <- c(rep(0,ncol(D01)))
for(i in 5:ncol(D01))
{ matrix1 <- matrix(data=0,nrow=2,ncol=2,byrow=FALSE,dimnames=NULL)
matrix1[1,1]=length(which(D01[,i]==1&D4[,i]==1))
matrix1[1,2]=length(which(D01[,i]==1&D4[,i]==0))
matrix1[2,1]=length(which(D01[,i]==0&D4[,i]==1))
matrix1[2,2]=length(which(D01[,i]==0&D4[,i]==0))
y= mcnemar.test(matrix1)
P[i] <- y$p.value}
write.table(P,file='D:/OW/P1.txt',quote=FALSE,sep='\t',row.names=F,col.names = T)

#D0 VS D14
View(D01)

P <- c(rep(0,ncol(D01)))
for(i in 5:ncol(D01))
{ matrix1 <- matrix(data=0,nrow=2,ncol=2,byrow=FALSE,dimnames=NULL)
matrix1[1,1]=length(which(D01[,i]==1&D14[,i]==1))
matrix1[1,2]=length(which(D01[,i]==1&D14[,i]==0))
matrix1[2,1]=length(which(D01[,i]==0&D14[,i]==1))
matrix1[2,2]=length(which(D01[,i]==0&D14[,i]==0))
y= mcnemar.test(matrix1)
P[i] <- y$p.value}
write.table(P,file='D:/OW/P2.txt',quote=FALSE,sep='\t',row.names=F,col.names = T)


#D0 VS D30
P <- c(rep(0,ncol(D03)))
for(i in 5:ncol(D03))
{ matrix1 <- matrix(data=0,nrow=2,ncol=2,byrow=FALSE,dimnames=NULL)
matrix1[1,1]=length(which(D03[,i]==1&D30[,i]==1))
matrix1[1,2]=length(which(D03[,i]==1&D30[,i]==0))
matrix1[2,1]=length(which(D03[,i]==0&D30[,i]==1))
matrix1[2,2]=length(which(D03[,i]==0&D30[,i]==0))
y= mcnemar.test(matrix1)
P[i] <- y$p.value}
write.table(P,file='D:/OW/P3.txt',quote=FALSE,sep='\t',row.names=F,col.names = T)

#D0 VS D60
P <- c(rep(0,ncol(D03)))
for(i in 5:ncol(D03))
{ matrix1 <- matrix(data=0,nrow=2,ncol=2,byrow=FALSE,dimnames=NULL)
matrix1[1,1]=length(which(D03[,i]==1&D60[,i]==1))
matrix1[1,2]=length(which(D03[,i]==1&D60[,i]==0))
matrix1[2,1]=length(which(D03[,i]==0&D60[,i]==1))
matrix1[2,2]=length(which(D03[,i]==0&D60[,i]==0))
y= mcnemar.test(matrix1)
P[i] <- y$p.value}
write.table(P,file='D:/OW/P4.txt',quote=FALSE,sep='\t',row.names=F,col.names = T)


library(lmerTest)
glmotuzero$Date1 <- factor(glmotuzero$Date1)
glmotuzero$Date1 <- relevel(glmotuzero$Date1,ref = 'D30')
fit1=lmer(OTU_99~Date1+(1|Name),data = glmotuzero)
summary(fit1)
glmgenzero$Date1 <- relevel(glmgenzero$Date1,ref = 'D30')

#Figure S15 Heat map
library(ComplexHeatmap)
library(GetoptLong)
library(dendextend)
library(circlize)
library(grid)
egow = read.table("D:/OW/owpvalue.txt",header=T,sep='\t',row.names=1)
owSSP= read.table("D:/OW/owSSP.txt",header=T,sep='\t',row.names=1)
owdet = read.table("D:/OW/owdet.txt",header=T,sep='\t',row.names=1)
owdet <- as.matrix(owdet)
owdetz = apply(owdet, 1, scale)
rownames(owdetz) = colnames(owdet)
write.table(owdetz,file='D:/OW/owdetz.txt',quote=FALSE,sep='\t',row.names=F,col.names = T)
owdetz = read.table("D:/OW/owdetz.txt",header=T,sep='\t',row.names=1)
owdetzt <- t(owdetz)
H1=Heatmap(owdetzt,name="Detection rate\n     Z-score", col = colorRamp2(c(-1,1), c("white","#0000CC")),show_column_names=T, row_names_side='left',show_row_names = T,cluster_rows = FALSE, cluster_columns = FALSE,width =  unit(40, "mm"),height = unit(250,'mm'),heatmap_legend_param=list(legend_direction='horizontal'),split=owSSP$Family)
owrel = read.table("D:/OW/owrel.txt",header=T,sep='\t',row.names=1)
owrel <- as.matrix(owrel)
owrelz = apply(owrel, 1, scale)
rownames(owrelz) = colnames(owrel)
owrelzt <- t(owrelz)
pvalue = row_anno_points(egow$D0_D14)
ha_mix_right = HeatmapAnnotation(p = pvalue,which = 'row',width  = unit( 2, "cm"))
H2=Heatmap(owrelzt,name="Relative abundance\n     Z-score", col = colorRamp2(c(-1,1), c("white","#FF0000")),right_annotation = ha_mix_right,show_column_names=T, show_row_names = F,cluster_rows = FALSE, cluster_columns = FALSE,width =  unit(40, "mm"),height = unit(250,'mm'),heatmap_legend_param=list(legend_direction='horizontal'))

H3=add_heatmap(H1,H2,direction =  "horizontal")

###Figure S16a Sankey plots to describe the number of ZOTUs whose relative abundance shift or return to the baseline level across the five timepoints. 
library(riverplot)
library(RColorBrewer) 
edges = data.frame(ID =c('F','G','H','I','J','K','L','M','N','O','P','Q','A1','A2','A3','A4'),N1 = c('A','A','Cno','A','Dno','A','A','B','B','C','C','Cno','Cno','D','D','Dno'),  
                   N2 = c('Eno','E','D','D','E','C','B','Cno','C','Dno','D','Eno','E','Eno','E','Eno'),  
                   Value = c(21.26,17.86,18.71,31.47,5.1,40.82,218.58,30.62,187.96,9.36,219.43,4.25,7.65,9.36,260.26,4.25),  
                   stringsAsFactors = F)  
nodes = data.frame(ID = unique(c(edges$N1, edges$N2)), stringsAsFactors = FALSE)  
nodes$x = c(1,5,7,3,5,7,9,9)  
nodes$y = c(3,-2,-2,9,7,5,0,4)  
rownames(nodes) = nodes$ID  
head(nodes) 
palette = paste0(brewer.pal(9,"Blues"), "70") 
palette1 = paste0(brewer.pal(9,"BuPu"), "70")
j <- makeRiver(nodes, edges,node_labels= c( A= "D0", Cno='No return',Dno='No return',B= "D4", C= "D14",D= "D30",Eno='No return',E= "D60" ),
               node_styles= list( A= list( col= '#0066CC'),Cno=list(col='#990066'),Dno=list(col='#990066'),B=list(col= "#0066CC"),C=list(col= "#0066CC"),D=list(col= "#0066CC"),Eno=list(col='#990066'),E=list(col= "#0066CC")), 
               edge_styles = list(F=list(col=palette1[6]),G=list(col=palette[5]),H=list(col=palette[5]),I=list(col=palette[5]),J=list(col=palette[5]),
                                  K=list(col=palette[5]), L=list(col=palette[5]),M=list(col=palette1[6]),N=list(col=palette[5]),O=list(col=palette1[6]),
                                  P=list(col=palette[5]),Q=list(col=palette1[6]),A1=list(col=palette[5]),A2=list(col=palette1[6]),A3=list(col=palette[5]),A4=list(col=palette1[6])))

ds <- default.style()
ds[["edgecol"]] <- "col"
ds[['col']]<-"white"

plot(j,yscale=0.02,srt='0',default_style=ds)


#Figure S16b Sankey plots to describe the number of ZOTUs whose relative abundance shift or return to the baseline level across the five timepoints. 

edges1 = data.frame(ID =c('F','G','H','I','J','K','L','M','N','O','P','Q','A1','A2','A3','A4'),N1 = c('A','A','Cno','A','Dno','A','A','B','B','C','C','Cno','Cno','D','D','Dno'),  
                    N2 = c('Eno','E','D','D','E','C','B','Cno','C','Dno','D','Eno','E','Eno','E','Eno'),  
                    Value = c(14.46,11.91,22.11,14.46,14.46,32.32,256.86,38.27,218.58,22.11,228.79,9.36,6.8,44.23,221.13,7.65),  
                    stringsAsFactors = F)  
nodes1 = data.frame(ID = unique(c(edges1$N1, edges1$N2)), stringsAsFactors = FALSE) 
nodes1$x = c(1,5,7,3,5,7,9,9)  
nodes1$y = c(3,-2,-2,9,8,7,0,4)  
rownames(nodes1) = nodes1$ID  
head(nodes1)  
palette2 = paste0(brewer.pal(9,"Reds"), "70") #####
palette1 = paste0(brewer.pal(9,"BuPu"), "70")####
j1 <- makeRiver(nodes1, edges1,node_labels= c( A= "D0", Cno='No return',Dno='No return',B= "D4", C= "D14",D= "D30",Eno='No return',E= "D60" ),
                node_styles= list( A= list( col= '#FF3333'),Cno=list(col='#990066'),Dno=list(col='#990066'),B=list(col= "#FF3333"),C=list(col= "#FF3333"),D=list(col= "#FF3333"),Eno=list(col='#990066'),E=list(col= "#FF3333")), 
                edge_styles = list(F=list(col=palette1[6]),G=list(col=palette2[5]),H=list(col=palette2[5]),I=list(col=palette2[5]),J=list(col=palette2[5]),
                                   K=list(col=palette2[5]), L=list(col=palette2[5]),M=list(col=palette1[6]),N=list(col=palette2[5]),O=list(col=palette1[6]),
                                   P=list(col=palette2[5]),Q=list(col=palette1[6]),A1=list(col=palette2[5]),A2=list(col=palette1[6]),A3=list(col=palette2[5]),A4=list(col=palette1[6])))

ds <- default.style()
ds[["edgecol"]] <- "col"
ds[['col']]<-'white'
plot(j1,yscale=0.02,srt='0',default_style=ds)

####Figure S16c The overlap between the oral cavity and sputum microbiota at ZOTU level
owdd = read.table("D:/OW/owdd.txt",header=T,sep='\t',row.names=1)
owdd <- as.matrix(owdd)
Heatmap(owdd,col = colorRamp2(c(-0.05,0.05,0.1), c("yellow","green","#0000CC")),show_column_names=T, show_row_names = F,cluster_rows = F, cluster_columns = F,width =  unit(150, "mm"),height = unit(200,'mm'),heatmap_legend_param=list(legend_direction='horizontal'))

####Figure S16d The overlap between the oral cavity and sputum microbiota at ZOTU level
owrr = read.table("D:/OW/owrr.txt",header=T,sep='\t',row.names=1)
owrr <- as.matrix(owrr)
Heatmap(owrr,col = colorRamp2(c(-0.05,0.05,0.1), c("yellow","green","#FF0000")),show_column_names=T, show_row_names = F,cluster_rows = F, cluster_columns = F,width =  unit(150, "mm"),height = unit(200,'mm'),heatmap_legend_param=list(legend_direction='horizontal'))

adonisdata <-read.csv("D:/OW/UadonisDOD60.csv",header = T,row.names = 1)
adonisgroup =read.table("D:/OW/UadonisgroupDOD60.txt", header=T, row.names= 1,sep="\t")
adonis5(otu1~Date1,data=RDAfactor1,permutations = 999, method="dw",tree=otu_tree) -> au
au

###Multivariable PERMANOVA：The relation of the environmental factors, antibiotic event and the sputum microbiota variation at D4，D30 in azithromycin in oral wash samples
#D4
otu1 <- read.csv("D:/OW/RDAD0D4.csv", row.names=1)
RDAfactor1 <- read.csv("D:/OW/factorD0D4.csv", row.names=1)
adonis2N(otu1~Date1+Humidity+PM2.5+PM10,data=RDAfactor1,permutations = 999,tree=otu_tree, method="dw",by = "margin") -> au
au
#D30
otu <- read.csv("D:/OW/RDAD0D30.csv", row.names=1)
RDAfactor <- read.csv("D:/OW/factorD0D30.csv", row.names=1)
adonis2N(otu~Date1+Humidity+PM2.5+PM10,data=RDAfactor,permutations = 999,tree=otu_tree, method="dw",by = "margin") -> au
au
##in placebo in oral wash samples
#D4
otu1 <- read.csv("D:/OW/RDAD0D4one.csv", row.names=1)
RDAfactor1 <- read.csv("D:/OW/factorD0D4one.csv", row.names=1)
adonis2N(otu1~Temperature+Humidity+PM2.5+PM10,data=RDAfactor1,permutations = 999,tree=otu_tree, method="dw",by = "margin") -> au
au
#D4
otu <- read.csv("D:/OW/RDAD0D30one.csv", row.names=1)
RDAfactor <- read.csv("D:/OW/factorD0D30one.csv", row.names=1)
adonis2N(otu~Temperature+Humidity+PM2.5+PM10,data=RDAfactor,permutations = 999,tree=otu_tree, method="dw",by = "margin") -> au
au

#The effects of azithromycin on airway microbial interactions within the microbial community network
  #Figure S9a The network analysis in azithromycin group. 
  library(igraph)
  library(reshape2)
  library(bindrcpp)
  library(VennDiagram)
  library(grid)
  library(futile.logger)
  library(ggthemes)
  library(shiny)
  library(agricolae)
  library(vegan)
  library(LearnBayes)
  library(dplyr)
  library(psych)
  library(sqldf)
  library(gsubfn)
  library(proto)
  library(RSQLite)
  library(digest)
  library(AnnotationDbi)
  library(stats4)
  library(BiocGenerics)
  library(parallel)
  library(Biobase)
  library(IRanges)
  library(S4Vectors)
  library(impute)
  library(GO.db)
  library(preprocessCore)
  library(WGCNA)
  library(dynamicTreeCut)
  library(fastcluster)
  library(multtest)
  source("D:/myfunction/matrix2igraph.R")
  source("D:/myfunction/net_pro.R")
  source("D:/myfunction/node_pro.R")
  #D0 network
  owotu_sample_file <- "D:/OW/ownetpre0.txt"
  owotu_tax_file<-"D:/OW/ownetpretax0.txt"
  r.threshold=0.6
  p.threshold=0.05
  size=3
  gcol=2
  glab=3
  owotunetpre0 <- read.table(owotu_sample_file, header=T, row.names= 1, sep="\t", comment.char="") 
  owotunetpre0<-t(owotunetpre0)
  owotunetpretax0 <- read.table(owotu_tax_file, header=T, row.names= 1, sep="\t", comment.char="")
  owotunetpretax0<-as.data.frame(owotunetpretax0[colnames(owotunetpre0),])
  owotu_abundancepre0 <- colSums(owotunetpre0)
  owotu_pronetpre0 <- cbind(owotunetpretax0,owotu_abundancepre0)
  owigraphpre0<-matrix2igraph(owotunetpre0,r.threshold,p.threshold)
  owigraphpre0.weight <- E(owigraphpre0)$weight
  E(owigraphpre0)$weight <- NA
  sum(owigraphpre0.weight>0)# number of postive correlation
  sum(owigraphpre0.weight<0)# number of negative correlation
  
  dir.create("D:/OW/network_results0")
  #The topological properties of microbial networks
  ownetpro_result<-net_pro(owigraphpre0)
  write.csv(ownetpro_result,"D:/OW/network_results0/igraphpre.network.pro.csv")
  #The topological properties of microbial notes
  ownodepro_result<-node_pro(owigraphpre0)
  write.csv(ownodepro_result,"D:/OW/network_results0/igraphpre.node.pro.csv")
  
  ####D0ZERO
  occor<-corAndPvalue(owotunetpre0,method = c( "spearman"))
  # multiple test the p values
  mtadj<-mt.rawp2adjp(unlist(occor$p),proc="bonferroni")
  adpcor<-mtadj$adjp[order(mtadj$index),2]
  occor.p<-matrix(adpcor,dim(owotunetpre0)[2])
  ##R值
  occor.r<-occor$cor
  # 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
  occor.r[occor.p>p.threshold|abs(occor.r)<r.threshold] = 0 
  write.csv(occor.r,"D:/OW/cor.csv")
  
  #D4 network
  owotu_sample_file <- "D:/OW/ownet0D4.txt"
  owotu_tax_file<-"D:/OW/ownet0D4tax.txt"
  r.threshold=0.6
  p.threshold=0.05
  size=3
  gcol=2
  glab=3
  owotunetpre <- read.table(owotu_sample_file, header=T, row.names= 1, sep="\t", comment.char="") 
  owotunetpre<-t(owotunetpre)
  owotunetpretax <- read.table(owotu_tax_file, header=T, row.names= 1, sep="\t", comment.char="")
  owotunetpretax<-as.data.frame(owotunetpretax[colnames(owotunetpre),])
  owotu_abundancepre <- colSums(owotunetpre)
  owotu_pronetpre <- cbind(owotunetpretax,owotu_abundancepre)
  owigraph0D4<-matrix2igraph(owotunetpre,r.threshold,p.threshold)
  owigraph0D4.weight <- E(owigraph0D4)$weight
  E(owigraph0D4)$weight <- NA
  sum(owigraph0D4.weight>0)# number of postive correlation
  sum(owigraph0D4.weight<0)# number of negative correlation
  
  dir.create("D:/OW/network_results0D4")
  #The topological properties of microbial networks
  ownetpro_result<-net_pro(owigraph0D4)
  write.csv(ownetpro_result,"D:/OW/network_results0D4/igraphpre.network.pro.csv")
  #The topological properties of microbial notes
  ownodepro_result<-node_pro(owigraph0D4)
  write.csv(ownodepro_result,"D:/OW/network_results0D4/igraphpre.node.pro.csv")
  
  ####D4ZERO
  occor<-corAndPvalue(owotunetpre,method = c( "spearman"))
  # multiple test the p values
  mtadj<-mt.rawp2adjp(unlist(occor$p),proc="bonferroni")
  adpcor<-mtadj$adjp[order(mtadj$index),2]
  occor.p<-matrix(adpcor,dim(owotunetpre)[2])
  ##R值
  occor.r<-occor$cor
  # 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
  occor.r[occor.p>p.threshold|abs(occor.r)<r.threshold] = 0 
  write.csv(occor.r,"D:/OW/corD4zero.csv")
  
  
  #D14 network
  owotu_sample_file <- "D:/OW/ownet0D14.txt"
  owotu_tax_file<-"D:/OW/ownet0D14tax.txt"
  r.threshold=0.6
  p.threshold=0.05
  size=3
  gcol=2
  glab=3
  owotunetpre <- read.table(owotu_sample_file, header=T, row.names= 1, sep="\t", comment.char="") 
  owotunetpre<-t(owotunetpre)
  owotunetpretax <- read.table(owotu_tax_file, header=T, row.names= 1, sep="\t", comment.char="")
  owotunetpretax<-as.data.frame(owotunetpretax[colnames(owotunetpre),])
  owotu_abundancepre <- colSums(owotunetpre)
  owotu_pronetpre <- cbind(owotunetpretax,owotu_abundancepre)
  owigraph0D14<-matrix2igraph(owotunetpre,r.threshold,p.threshold)
  owigraph0D14.weight <- E(owigraph0D14)$weight
  E(owigraph0D14)$weight <- NA
  sum(owigraph0D14.weight>0)# number of postive correlation
  sum(owigraph0D14.weight<0)# number of negative correlation
  
  dir.create("D:/OW/network_results0D14")
  #The topological properties of microbial networks
  ownetpro_result<-net_pro(owigraph0D14)
  write.csv(ownetpro_result,"D:/OW/network_results0D14/igraphpre.network.pro.csv")
  
  #The topological properties of microbial notes
  ownodepro_result<-node_pro(owigraph0D14)
  write.csv(ownodepro_result,"D:/OW/network_results0D14/igraphpre.node.pro.csv")
  
  ####D14ZERO
  occor<-corAndPvalue(owotunetpre,method = c( "spearman"))
  # multiple test the p values
  mtadj<-mt.rawp2adjp(unlist(occor$p),proc="bonferroni")
  adpcor<-mtadj$adjp[order(mtadj$index),2]
  occor.p<-matrix(adpcor,dim(owotunetpre)[2])
  ##R值
  occor.r<-occor$cor
  # 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
  occor.r[occor.p>p.threshold|abs(occor.r)<r.threshold] = 0 
  write.csv(occor.r,"D:/OW/corD14zero.csv")
  
  
  #D30 network
  owotu_sample_file <- "D:/OW/ownet0D30.txt"
  owotu_tax_file<-"D:/OW/ownet0D30tax.txt"
  r.threshold=0.6
  p.threshold=0.05
  size=3
  gcol=2
  glab=3
  owotunetpre <- read.table(owotu_sample_file, header=T, row.names= 1, sep="\t", comment.char="") 
  owotunetpre<-t(owotunetpre)
  owotunetpretax <- read.table(owotu_tax_file, header=T, row.names= 1, sep="\t", comment.char="")
  owotunetpretax<-as.data.frame(owotunetpretax[colnames(owotunetpre),])
  owotu_abundancepre <- colSums(owotunetpre)
  owotu_pronetpre <- cbind(owotunetpretax,owotu_abundancepre)
  owigraph0D30<-matrix2igraph(owotunetpre,r.threshold,p.threshold)
  owigraph0D30.weight <- E(owigraph0D30)$weight
  E(owigraph0D30)$weight <- NA
  sum(owigraph0D30.weight>0)# number of postive correlation
  sum(owigraph0D30.weight<0)# number of negative correlation
  
  dir.create("D:/OW/network_results0D30")
  #The topological properties of microbial networks
  ownetpro_result<-net_pro(owigraph0D30)
  write.csv(ownetpro_result,"D:/OW/network_results0D30/igraphpre.network.pro.csv")
  #The topological properties of microbial notes
  ownodepro_result<-node_pro(owigraph0D30)
  write.csv(ownodepro_result,"D:/OW/network_results0D30/igraphpre.node.pro.csv")
  
  ####D30ZERO
  occor<-corAndPvalue(owotunetpre,method = c( "spearman"))
  # multiple test the p values
  mtadj<-mt.rawp2adjp(unlist(occor$p),proc="bonferroni")
  adpcor<-mtadj$adjp[order(mtadj$index),2]
  occor.p<-matrix(adpcor,dim(owotunetpre)[2])
  ##R值
  occor.r<-occor$cor
  # 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
  occor.r[occor.p>p.threshold|abs(occor.r)<r.threshold] = 0 
  write.csv(occor.r,"D:/OW/corD30zero.csv")
  
  
  #D60 network
  owotu_sample_file <- "D:/OW/ownet0D60.txt"
  owotu_tax_file<-"D:/OW/ownet0D60tax.txt"
  r.threshold=0.6
  p.threshold=0.05
  size=3
  gcol=2
  glab=3
  owotunetpre <- read.table(owotu_sample_file, header=T, row.names= 1, sep="\t", comment.char="") 
  owotunetpre<-t(owotunetpre)
  owotunetpretax <- read.table(owotu_tax_file, header=T, row.names= 1, sep="\t", comment.char="")
  owotunetpretax<-as.data.frame(owotunetpretax[colnames(owotunetpre),])
  owotu_abundancepre <- colSums(owotunetpre)
  owotu_pronetpre <- cbind(owotunetpretax,owotu_abundancepre)
  owigraph0D60<-matrix2igraph(owotunetpre,r.threshold,p.threshold)
  owigraph0D60.weight <- E(owigraph0D60)$weight
  E(owigraph0D60)$weight <- NA
  sum(owigraph0D60.weight>0)# number of postive correlation
  sum(owigraph0D60.weight<0)# number of negative correlation
  
  
  dir.create("D:/OW/network_results0D60")
  #The topological properties of microbial networks
  ownetpro_result<-net_pro(owigraph0D60)
  write.csv(ownetpro_result,"D:/OW/network_results0D60/igraphpre.network.pro.csv")
  #The topological properties of microbial notes
  ownodepro_result<-node_pro(owigraph0D60)
  write.csv(ownodepro_result,"D:/OW/network_results0D60/igraphpre.node.pro.csv")
  
  ####D60ZERO
  occor<-corAndPvalue(owotunetpre,method = c( "spearman"))
  # multiple test the p values
  mtadj<-mt.rawp2adjp(unlist(occor$p),proc="bonferroni")
  adpcor<-mtadj$adjp[order(mtadj$index),2]
  occor.p<-matrix(adpcor,dim(owotunetpre)[2])
  ##R值
  occor.r<-occor$cor
  # 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
  occor.r[occor.p>p.threshold|abs(occor.r)<r.threshold] = 0 
  write.csv(occor.r,"D:/OW/corD60zero.csv")
  
  #Figure 9b The number of shared edges between oral cavity and sputum network at timepoints D0,D4, D14, D30, D60 in azithromycin
  owD00 <- length(E(intersection(owigraphpre0,igraphpre0)))
  owD44 <- length(E(intersection(owigraph0D4,igraph0D4)))
  owD1414 <- length(E(intersection(owigraph0D14,igraph0D14)))
  owD3030 <- length(E(intersection(owigraph0D30,igraph0D30)))
  owD6060 <- length(E(intersection(owigraph0D60,igraph0D60)))
  
  owedgenumpre0 <-length(E(owigraphpre0))
  owedgenum0D4<-length(E(owigraph0D4))
  owedgenum0D14<-length(E(owigraph0D14))
  owedgenum0D30<-length(E(owigraph0D30))
  owedgenum0D60<-length(E(owigraph0D60))
  
  edgenumpre0 <- length(E(igraphpre0))
  edgenum0D4<-length(E(igraph0D4))
  edgenum0D14<-length(E(igraph0D14))
  edgenum0D30<-length(E(igraph0D30))
  edgenum0D60<-length(E(igraph0D60))
  
  ###D0 
  grid.newpage()
  draw.pairwise.venn(area1=owedgenumpre0,area2=edgenumpre0,cross.area=owD00
                     ,lwd=1,lty=1
                     ,col=c('#FF3399','#009999'),fill=c('#FF3399','#33FF33')
                     ,cat.col=c('#FF3399','#009999')
                     ,rotation.degree=360, scaled = FALSE,cex = 2,cat.cex =1, cat.pos = c(400, 200))
  
  ###D4 
  grid.newpage()
  draw.pairwise.venn(area1=owedgenum0D4,area2=edgenum0D4,cross.area=owD44
                     ,lwd=1,lty=1
                     ,col=c('#FF3399','#009999'),fill=c('#FF3399','#33FF33')
                     ,cat.col=c('#FF3399','#009999')
                     ,rotation.degree=360, scaled = FALSE,cex = 2,cat.cex =1, cat.pos = c(400, 200))
  
  ####D14
  grid.newpage()
  draw.pairwise.venn(area1=owedgenum0D14,area2=edgenum0D14,cross.area=owD1414
                     ,lwd=1,lty=1
                     ,col=c('#FF3399','#009999'),fill=c('#FF3399','#33FF33')
                     ,cat.col=c('#FF3399','#009999')
                     ,rotation.degree=360, scaled = FALSE,cex = 2,cat.cex =1, cat.pos = c(400, 200))
  
  ####D30
  grid.newpage()
  draw.pairwise.venn(area1=owedgenum0D30,area2=edgenum0D30,cross.area=owD3030
                     ,lwd=1,lty=1
                     ,col=c('#FF3399','#009999'),fill=c('#FF3399','#33FF33')
                     ,cat.col=c('#FF3399','#009999')
                     ,rotation.degree=360, scaled = FALSE,cex = 2,cat.cex =1, cat.pos = c(400, 200))
  
  ####D60
  grid.newpage()
  draw.pairwise.venn(area1=owedgenum0D60,area2=edgenum0D60,cross.area=owD6060
                     ,lwd=1,lty=1
                     ,col=c('#FF3399','#009999'),fill=c('#FF3399','#33FF33')
                     ,cat.col=c('#FF3399','#009999')
                     ,rotation.degree=360, scaled = FALSE,cex = 2,cat.cex =1, cat.pos = c(400, 200))
  
  #Figure 9c The closeness centralization of shared nodes between oral cavity and sputum network at timepoints D0,D4, D14, D30, D60 in azithromycin
  library(ComplexHeatmap)
  library(circlize)
  OWD0D0 = read.table("D:/OW/OWD0D0.txt",header=T,sep='\t',row.names=1)
  OWD0D0 <- as.matrix(OWD0D0)
  Heatmap(OWD0D0,name='Rank of the closeness',col = colorRamp2(c(0,0.05), c("#FFFF99",'#3366FF')), show_column_names = F,show_row_names = T,cluster_rows = FALSE, cluster_columns = FALSE,width =  unit(150, "mm"),height = unit(40,'mm'),heatmap_legend_param=list(legend_direction='horizontal'))
  
  OWD4D4 = read.table("D:/OW/OWD4D4.txt",header=T,sep='\t',row.names=1)
  OWD4D4 <- as.matrix(OWD4D4)
  Heatmap(OWD4D4,name='Rank of the closeness',col = colorRamp2(c(0,0.05), c("#FFFF99",'#3366FF')), show_column_names = F,show_row_names = T,cluster_rows = FALSE, cluster_columns = FALSE,width =  unit(150, "mm"),height = unit(40,'mm'),heatmap_legend_param=list(legend_direction='horizontal'))
  
  OWD14D14 = read.table("D:/OW/OWD14D14.txt",header=T,sep='\t',row.names=1)
  OWD14D14 <- as.matrix(OWD14D14)
  Heatmap(OWD14D14,name='Rank of the closeness',col = colorRamp2(c(0,0.05), c("#FFFF99",'#3366FF')), show_column_names = F,show_row_names = T,cluster_rows = FALSE, cluster_columns = FALSE,width =  unit(150, "mm"),height = unit(40,'mm'),heatmap_legend_param=list(legend_direction='horizontal'))
  
  OWD30D30 = read.table("D:/OW/OWD30D30.txt",header=T,sep='\t',row.names=1)
  OWD30D30 <- as.matrix(OWD30D30)
  Heatmap(OWD30D30,name='Rank of the closeness',col = colorRamp2(c(0,0.05), c("#FFFF99",'#3366FF')), show_column_names = F,show_row_names = T,cluster_rows = FALSE, cluster_columns = FALSE,width =  unit(150, "mm"),height = unit(40,'mm'),heatmap_legend_param=list(legend_direction='horizontal'))
  
  OWD60D60 = read.table("D:/OW/OWD60D60.txt",header=T,sep='\t',row.names=1)
  OWD60D60 <- as.matrix(OWD60D60)
  Heatmap(OWD60D60,name='Rank of the closeness',col = colorRamp2(c(0,0.05), c("#FFFF99",'#3366FF')), show_column_names = F,show_row_names = T,cluster_rows = FALSE, cluster_columns = FALSE,width =  unit(150, "mm"),height = unit(40,'mm'),heatmap_legend_param=list(legend_direction='horizontal'))
  
  #Figure S8 The topological properties of microbial networks in azithromycin.
  library(ggplot2)
  ggplot()+
    geom_point(owtra1,mapping = aes(reorder(Date1,OR),Parameter,colour=Group),size=5)+
    geom_line(owtra1,mapping=aes(reorder(Date1,OR),Parameter,group=Group,colour=Group))+
    scale_color_manual(values = c('#CC3300','#FF9900'))+
    theme_classic()+
    xlab("Timepoints") + ylab("Vertices (n)")
  
  ggplot()+
    geom_point(owtra2,mapping = aes(reorder(Date1,OR),Parameter,colour=Group),size=5)+
    geom_line(owtra2,mapping=aes(reorder(Date1,OR),Parameter,group=Group,colour=Group))+
    scale_color_manual(values = c('#CC3300','#FF9900'))+
    theme_classic()+
    xlab("Timepoints") + ylab("Edges (n)")
  
  ggplot()+
    geom_point(owtra3,mapping = aes(reorder(Date1,OR),Parameter,colour=Group),size=5)+
    geom_line(owtra3,mapping=aes(reorder(Date1,OR),Parameter,group=Group,colour=Group))+
    scale_color_manual(values = c('#CC3300','#FF9900'))+
    theme_classic()+
    xlab("Timepoints") + ylab("Positive edges (n)")
  
  ggplot()+
    geom_point(owtra4,mapping = aes(reorder(Date1,OR),Parameter,colour=Group),size=5)+
    geom_line(owtra4,mapping=aes(reorder(Date1,OR),Parameter,group=Group,colour=Group))+
    scale_color_manual(values = c('#CC3300','#FF9900'))+
    theme_classic()+
    xlab("Timepoints") + ylab("Negative edges (n)")
  
  ggplot()+
    geom_point(owtra5,mapping = aes(reorder(Date1,OR),Parameter,colour=Group),size=5)+
    geom_line(owtra5,mapping=aes(reorder(Date1,OR),Parameter,group=Group,colour=Group))+
    scale_color_manual(values = c('#CC3300','#FF9900'))+
    theme_classic()+
    xlab("Timepoints") + ylab("Average degree")
  
  ggplot()+
    geom_point(owtra6,mapping = aes(reorder(Date1,OR),Parameter,colour=Group),size=5)+
    geom_line(owtra6,mapping=aes(reorder(Date1,OR),Parameter,group=Group,colour=Group))+
    scale_color_manual(values = c('#CC3300','#FF9900'))+
    theme_classic()+
    xlab("Timepoints") + ylab("Connectance")
  
  ggplot()+
    geom_point(owtra7,mapping = aes(reorder(Date1,OR),Parameter,colour=Group),size=5)+
    geom_line(owtra7,mapping=aes(reorder(Date1,OR),Parameter,group=Group,colour=Group))+
    scale_color_manual(values = c('#CC3300','#FF9900'))+
    theme_classic()+
    xlab("Timepoints") + ylab("Average clustering coefficient")
  
  ggplot()+
    geom_point(owtra8,mapping = aes(reorder(Date1,OR),Parameter,colour=Group),size=5)+
    geom_line(owtra8,mapping=aes(reorder(Date1,OR),Parameter,group=Group,colour=Group))+
    scale_color_manual(values = c('#CC3300','#FF9900'))+
    theme_classic()+
    xlab("Timepoints") + ylab("Average centralization closeness")
  
  
  
 