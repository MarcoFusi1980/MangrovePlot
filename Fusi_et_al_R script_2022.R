# Packages needed to run this code
library(vegan)
library(tidyverse)
library('BiodiversityR')
library(devtools)
library(pairwiseAdonis)
library(mvabund)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(rmarkdown)
library(RVAideMemoire)
library(OTUtable)
library(compositions)
library(randomForest)
library(dplyr)

setwd("~/Dropbox/Articoli/PLOT/0_PLOT_MS_NEW/0_Mature/SpectrumMicrobiology/00_Markdown/")
# LOAD DATA AND PREPARE IT FOR ANALISIS
otu_mat<-read.table("OTU_mat.txt",header=T,sep="\t",row.names = 1)
map_mat<-read.table("map_mat.txt",header=T,sep="\t")
otu_tax<-read.table("taxononomy.txt",header=T,sep="\t")

otu_mat_abund<-mvabund(otu_mat_filtered)
otu_mglm<-manyglm(otu_mat_abund~map_mat$MONTH*map_mat$DEPTH*map_mat$CRAB, offset = otu_mat[rowSums()])
otu_mglm<-anova(otu_mglm, block = 'plotty')


matrix_clr<-compositions::clr(matrice)
## export in Primer for PcoA analysis
write.table(matrie_clr,"~/clrotutable.txt",sep='\t')


otu_rf_test_d<-otu_rf_test %>% dplyr::filter(otu_rf_test$depth=="D")
otu_rf_test_s<-otu_rf_test %>% dplyr::filter(otu_rf_test$depth=="S")
otu_rf_test_ss<-otu_rf_test %>% dplyr::filter(otu_rf_test$depth=="SS")


otu_rf_test_s<-otu_rf_test_s[,colSums(otu_rf_test_s[1:1241])>0]
otu_rf_test_ss<-otu_rf_test_ss[,colSums(otu_rf_test_ss[1:1241])>0]
otu_rf_test_d<-otu_rf_test_d[,colSums(otu_rf_test_d[1:1241])>0]


rf_d<-randomForest(factor(otu_rf_test_d$full_factor)~.,data=otu_rf_test_d[,c(1:1241)],proximity=TRUE,ntree=5000,importance=TRUE)
rf_s<-randomForest(factor(otu_rf_test_s$full_factor)~.,data=otu_rf_test_s[,c(1:1208)],proximity=TRUE,ntree=5000,importance=TRUE)
rf_ss<-randomForest(factor(otu_rf_test_ss$full_factor)~.,data=otu_rf_test_ss[,c(1:1235)],proximity=TRUE,ntree=5000,importance=TRUE)

plot(rf)
varImpPlot(rf,
           sort = T,
           n.var = 50,
           main = "Top 50 - Variable Importance")

otu_imp<-as.data.frame(importance(rf))
otu_imp_s<-as.data.frame(importance(rf_s))
otu_imp_ss<-as.data.frame(importance(rf_ss))
otu_imp_d<-as.data.frame(importance(rf_d))


otu_imp<-otu_imp[order(-otu_imp$MeanDecreaseAccuracy),c(35,36)]
otu_imp_s<-otu_imp_s[order(-otu_imp_s$MeanDecreaseAccuracy),c(11,12)]
otu_imp_ss<-otu_imp_ss[order(-otu_imp_ss$MeanDecreaseAccuracy),c(13,14)]
otu_imp_d<-otu_imp_d[order(-otu_imp_d$MeanDecreaseAccuracy),c(13,14)]



otu_imp <- tibble::rownames_to_column(otu_imp, "OTUID")
otu_imp_s <- tibble::rownames_to_column(otu_imp_s, "OTUID")
otu_imp_ss <- tibble::rownames_to_column(otu_imp_ss, "OTUID")
otu_imp_d <- tibble::rownames_to_column(otu_imp_d, "OTUID")

names(otu_imp)
otus<-as.data.frame(t(otu_mat))
otus <- tibble::rownames_to_column(otus, "OTUID")

otu_selected<-otus[ otus$OTUID %in% otu_imp$OTUID[c(1:100)],]
otu_selected_s<-otus[ otus$OTUID %in% otu_imp_s$OTUID[c(1:100)],]
otu_selected_ss<-otus[ otus$OTUID %in% otu_imp_ss$OTUID[c(1:100)],]
otu_selected_d<-otus[ otus$OTUID %in% otu_imp_d$OTUID[c(1:100)],]

taxo<-otu_tax[otu_tax$OTUID %in% otu_imp$OTUID[c(1:100)],]
taxo_s<-otu_tax[otu_tax$OTUID %in% otu_imp_s$OTUID[c(1:100)],]
taxo_ss<-otu_tax[otu_tax$OTUID %in% otu_imp_ss$OTUID[c(1:100)],]
taxo_d<-otu_tax[otu_tax$OTUID %in% otu_imp_d$OTUID[c(1:100)],]

otu_selected<-otu_selected[order(otu_selected$OTUID),]
otu_selected_s<-otu_selected[order(otu_selected_s$OTUID),]
otu_selected_ss<-otu_selected[order(otu_selected_ss$OTUID),]
otu_selected_d<-otu_selected[order(otu_selected_d$OTUID),]

taxo<-taxo[order(taxo$OTUID),]
taxo_s<-taxo[order(taxo_s$OTUID),]
taxo_ss<-taxo[order(taxo_ss$OTUID),]
taxo_d<-taxo[order(taxo_d$OTUID),]

ot<-otu_selected[otu_selected$OTUID%in%taxo$OTUID,]
ot_s<-otu_selected_s[otu_selected_s$OTUID%in%taxo_s$OTUID,]
ot_ss<-otu_selected_ss[otu_selected_ss$OTUID%in%taxo_ss$OTUID,]
ot_d<-otu_selected_d[otu_selected_d$OTUID%in%taxo_d$OTUID,]

ot<-ot[order(ot$OTUID),]
ot_s<-ot_s[order(ot_s$OTUID),]
ot_ss<-ot_ss[order(ot_ss$OTUID),]
ot_d<-ot_d[order(ot_d$OTUID),]

ot2<-cbind(ot,taxo)

ot2_s<-cbind(ot_s,taxo_s[1:97,])
ot2_ss<-cbind(ot_ss,taxo_ss[1:97,])
ot2_d<-cbind(ot_d,taxo_d)

otu_imp2<-otu_imp[otu_imp$OTUID%in%ot2$OTUID,]

otu_imp2_s<-otu_imp_s[otu_imp_s$OTUID%in%ot2_s$OTUID,]
otu_imp2_ss<-otu_imp_ss[otu_imp_ss$OTUID%in%ot2_ss$OTUID,]
otu_imp2_d<-otu_imp_d[otu_imp_d$OTUID%in%ot2_d$OTUID,]

table_selection<-read.table("../01_Taxonomy/importance_taxonomomy_all_values.csv",header=T,sep=",")
table_selected_s <- table_selection[table_selection$OUT %in% otu_imp_s$OTUID[c(1:100)],]
table_selected_ss <- table_selection[table_selection$OUT %in% otu_imp_ss$OTUID[c(1:100)],]
table_selected_d <- table_selection[table_selection$OUT %in% otu_imp_d$OTUID[c(1:100)],]


write.table(table_selected,"../01_Taxonomy/taxonomy_averaged_importance.txt",sep='\t')

write.table(table_selected_s,"../01_Taxonomy/taxonomy_averaged_importance_s.txt",sep='\t')
write.table(table_selected_ss,"../01_Taxonomy/taxonomy_averaged_importance_ss.txt",sep='\t')
write.table(table_selected_d,"../01_Taxonomy/taxonomy_averaged_importance_d.txt",sep='\t')


model_otu<-manyglm(otu_mat_abund~MONTH*DEPTH*CRAB,data = map_mat)
anova(model_otu, show.time = "all",pairwise.com=MONTH:DEPTH:CRAB,block = PLOT)

## analysis of pH

pH<-read.table("~/Dropbox/Articoli/PLOT/NEW_analysis/biochemistry/pH.txt",header=T,sep="\t")
pH$DEPTH<-factor(pH$DEPTH,level=c("S","SS","D"), labels=c("S","SS","D"))
pH$SITE<-factor(pH$SITE,level=c("MATURE","JUVENILE"), labels=c("MATURE","JUVENILE"))

pH<-na.omit(pH)
pH_m<-subset(pH,SITE=='MATURE')
pH_m<-subset(pH,BIOTURBATION!='NONE')
pH_m$BIOTURBATION<-factor(pH_m$BIOTURBATION,level=c("HIGH","LOW"), labels=c("HIGH","NORMAL"))


p <- ggplot(pH_m, aes(MONTH, pH,linetype=BIOTURBATION))+ 
  geom_jitter(aes(shape=BIOTURBATION,colour = DEPTH), size = 3,width = 0.15) +
  geom_smooth(aes(linetype=BIOTURBATION, color=DEPTH,group = paste(BIOTURBATION,DEPTH)), method = "loess",se=FALSE)+
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 18), 
        axis.title.x = element_text(size = 18), axis.text.y = element_text(size = 18), 
        axis.title.y = element_text(size = 18)) +  
  scale_x_discrete(limits=c("Aug_16","Nov_16","Feb_17","May_17","Aug_17"))+
  xlab('Month')+
  ylab('pH')+
  ylim(8,10.5)+
  scale_color_manual(values=c("snow4","dodgerblue1"))+
  scale_shape_manual(values=c(19, 1))

p

summary(pH_m)
model2<-manylm(pH~DEPTH*factor(MONTH)*BIOTURBATION, pH_m)
plotty<-factor(pH_m$PLOT)
modello<-anova(model2, block = 'plotty', resamp='case', test="F")

## analysis of Salinity
sal<-read.table("~/Dropbox/Articoli/PLOT/NEW_analysis/biochemistry/SALINITY.txt",header=T,sep="\t")
sal$DEPTH<-factor(sal$DEPTH,level=c("S","SS","D"), labels=c("S","SS","D"))
sal$SITE<-factor(sal$SITE,level=c("MATURE","JUVENILE"), labels=c("MATURE","JUVENILE"))
sal<-na.omit(sal)
sal_m<-subset(sal,SITE=='MATURE')
sal_m<-subset(sal_m,CRAB!='NONE')
sal_m$CRAB<-factor(sal_m$CRAB,level=c("HIGH","MEDIUM"), labels=c("HIGH","NORMAL"))

unique(sal_m$MONTH)

p <- ggplot(sal_m, aes(MONTH, SALINITY,linetype=CRAB))+ 
  geom_jitter(aes(shape=CRAB,colour = DEPTH), size = 3,width = 0.15) +
  geom_smooth(aes(linetype=CRAB, color=DEPTH,group = paste(CRAB,DEPTH)), method = "loess",se=FALSE)+
  theme_bw()+
  scale_x_discrete(limits=c("Aug_16","Nov_16","Feb_17","May_17","Aug_17"))+
  theme(legend.position = "none",
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 18), 
        axis.title.x = element_text(size = 18), axis.text.y = element_text(size = 18), 
        axis.title.y = element_text(size = 18)) +
  xlab('Month')+
  ylab('Salinity %')+
  scale_x_discrete(limits=c("Aug_16","Nov_16","Feb_17","May_17","Aug_17"))+
  scale_color_manual(values=c("snow4","dodgerblue1"))+
  scale_shape_manual(values=c(19, 1))

p


summary(sal_m)
model2<-manylm(SALINITY~DEPTH*factor(MONTH)*CRAB, sal_m)
plotty<-factor(sal_m$PLOT)
modello<-anova(model2, block = 'plotty', resamp='case', test="F")


## Multivariate analysis of the biochemical parameters
biochem<-read.table("~/Dropbox/Articoli/PLOT/NEW_analysis/biochemistry/BIOCHEMISTRY.txt",header=T,sep="\t")

biochem1<-subset(biochem, biochem$SITE=="MATURE")
biochem2<-subset(biochem1, biochem$BIOTURBATION!="NONE")

biochem.mod <- manylm(cbind(pH,Nitrite,
                            Phosphate,Silicate,Nitrate,POC,PIC,PIN,PON) ~ 
                            DEPTH * MONTH2 * BIOTURBATION,
                  data=biochem2)
plotty<-factor(biochem2$PLOT)
summary(biochem.mod, block = 'plotty', resamp='case', test="LR")
anova(biochem.mod)

biochem2<-na.omit(biochem2)
biochem2$DEPTH<-factor(biochem2$DEPTH,level=c("S","SS","D"), labels=c("Surface","SubSurface","Deep"))
biochem2$MONTH2<-factor(biochem2$MONTH2,levels = c("May-16","Aug-16","Nov-16","Feb-17","May-17","Aug-17"))


Sulphate<-ggplot(biochem2, aes(x=BIOTURBATION, y=Sulphate, fill=DEPTH)) + 
  geom_boxplot(notch=TRUE)+ 
  facet_grid(~MONTH2)+
  scale_fill_manual(values=c("dodgerblue1","orange1","snow4"))+
  xlab("") + 
  ylab(expression(paste("Sulphate (", mu, "g/mg dry sediment)")))+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), axis.text.y = element_text(size = 6), 
        axis.title.y = element_text(size = 6), strip.text = element_text(size = 6)) 

Nitrite<-ggplot(biochem2, aes(x=BIOTURBATION, y=Nitrite, fill=DEPTH)) + 
  geom_boxplot(notch=TRUE)+ 
  facet_grid(~MONTH2)+
  scale_fill_manual(values=c("dodgerblue1","orange1","snow4"))+
  xlab("") + 
  ylab(expression(paste("Nitrite (", mu, "g/mg dry sediment)")))+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 6), strip.text = element_text(size = 6)) 

Phosphate<-ggplot(biochem2, aes(x=BIOTURBATION, y=Phosphate, fill=DEPTH)) + 
  geom_boxplot(notch=TRUE)+ 
  facet_grid(~MONTH2)+
  scale_fill_manual(values=c("dodgerblue1","orange1","snow4"))+
  xlab("") + 
  ylab(expression(paste("Phosphate (", mu, "g/mg dry sediment)")))+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), axis.text.y = element_text(size = 6), 
        axis.title.y = element_text(size = 6), strip.text = element_text(size = 6)) 


Silicate<-ggplot(biochem2, aes(x=BIOTURBATION, y=Silicate, fill=DEPTH)) + 
  geom_boxplot(notch=TRUE)+ 
  facet_grid(~MONTH2)+
  scale_fill_manual(values=c("dodgerblue1","orange1","snow4"))+
  xlab("Bioturbation intensity") + 
  ylab(expression(paste("Silicate (", mu, "g/mg dry sediment)")))+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6), 
        axis.title.x = element_text(size = 6), axis.text.y = element_text(size = 6), 
        axis.title.y = element_text(size = 6), strip.text = element_text(size = 6))

Nitrate<-ggplot(biochem2, aes(x=BIOTURBATION, y=Nitrate, fill=DEPTH)) + 
  geom_boxplot(notch=TRUE)+ 
  facet_grid(~MONTH2)+
  scale_fill_manual(values=c("dodgerblue1","orange1","snow4"))+
  xlab("") + 
  ylab(expression(paste("Nitrate (", mu, "g/mg dry sediment)")))+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), axis.text.y = element_text(size = 6), 
        axis.title.y = element_text(size = 6), strip.text = element_text(size = 6)) 

POC<-ggplot(biochem2, aes(x=BIOTURBATION, y=POC, fill=DEPTH)) + 
  geom_boxplot(notch=TRUE)+ 
  facet_grid(~MONTH2)+
  scale_fill_manual(values=c("dodgerblue1","orange1","snow4"))+
  xlab("") + 
  ylab(expression(paste("POC (", mu, "g/mg dry sediment)")))+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), axis.text.y = element_text(size = 6), 
        axis.title.y = element_text(size = 6), strip.text = element_text(size = 6)) 

PIC<-ggplot(biochem2, aes(x=BIOTURBATION, y=PIC, fill=DEPTH)) + 
  geom_boxplot(notch=TRUE)+ 
  facet_grid(~MONTH2)+
  scale_fill_manual(values=c("dodgerblue1","orange1","snow4"))+
  xlab("") + 
  ylab(expression(paste("PIC (", mu, "g/mg dry sediment)")))+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(), 
        axis.title.x =element_blank(), axis.text.y = element_text(size = 6), 
        axis.title.y = element_text(size = 6), strip.text = element_text(size = 6)) 


PIN<-ggplot(biochem2, aes(x=BIOTURBATION, y=PIN, fill=DEPTH)) + 
  geom_boxplot(notch=TRUE)+ 
  facet_grid(~MONTH2)+
  scale_fill_manual(values=c("dodgerblue1","orange1","snow4"))+
  xlab("") + 
  ylab(expression(paste("PIN (", mu, "g/mg dry sediment)")))+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), axis.text.y = element_text(size = 6), 
        axis.title.y = element_text(size = 6), strip.text = element_text(size = 6)) 

PON<-ggplot(biochem2, aes(x=BIOTURBATION, y=PON, fill=DEPTH)) + 
  geom_boxplot(notch=TRUE)+ 
  facet_grid(~MONTH2)+
  scale_fill_manual(values=c("dodgerblue1","orange1","snow4"))+
  ylab(expression(paste("PON (", mu, "g/mg dry sediment)")))+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), axis.text.y = element_text(size = 6), 
        axis.title.y = element_text(size = 6), strip.text = element_text(size = 6)) 

PH<-ggplot(biochem2, aes(x=BIOTURBATION, y=pH, fill=DEPTH)) + 
  geom_boxplot(notch=TRUE)+ 
  facet_grid(~MONTH2)+
  scale_fill_manual(values=c("dodgerblue1","orange1","snow4"))+
  xlab("Bioturbation intensity") + 
  ylab("pH") +
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6), 
        axis.title.x = element_text(size = 6), axis.text.y = element_text(size = 6), 
        axis.title.y = element_text(size = 6), strip.text = element_text(size = 6))

plot_grid(POC,PON,PIC,PIN,Phosphate,Nitrite, Nitrate,Sulphate,Silicate,PH, ncol=2, nrow=5, align="v")


##  analysis of alpha-diversity indexes
#Preparation of the data
divers<-read.table("~/Dropbox/Articoli/PLOT/NEW_analysis/divrich/divrich.txt",header=T,sep="\t")
divers$DEPTH<-factor(divers$DEPTH,level=c("S","SS","D"), labels=c("S","SS","D"))
divers$SITE<-factor(divers$SITE,level=c("MATURE","JUVENILE"), labels=c("MATURE","JUVENILE"))
divers<-na.omit(divers)
divers_m<-subset(divers,SITE=='MATURE')
divers_m<-subset(divers_m,BIOTURBATION!='NONE')
unique(divers_m$MONTH2)
divers_m$MONTH2<-factor(divers_m$MONTH2,levels = c("May-16","Aug-16","Nov-16","Feb-17","May-17","Aug-17"))
divers_m$BIOTURBATION<-factor(divers_m$BIOTURBATION,levels = c("HIGH","MEDIUM"), labels = c("HIGH","NORMAL"))

#tests
summary(divers_m)
model2<-manylm(log(H)~BIOTURBATION*DEPTH*factor(MONTH), divers_m,offset=divers_m[rowSums()])
plotty<-factor(divers_m$PLOT)
anova(model2, block = 'plotty', resamp='case', test="F")

model2<-manylm(log(N)~BIOTURBATION*DEPTH*factor(MONTH), divers_m,offset=divers_m[rowSums()])
plotty<-factor(divers_m$PLOT)
anova(model2, block = 'plotty', resamp='case', test="F")

#graphs
Shannon<-ggplot(divers_m, aes(x=BIOTURBATION, y=H, fill=DEPTH)) + 
  geom_boxplot(notch=TRUE)+ 
  facet_grid(~MONTH2)+
  scale_fill_manual(values=c("dodgerblue1","orange1","snow4"))+
  xlab("Bioturbation intensity") + 
  ylab("Diversity (H)") +
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6), 
        axis.title.x = element_text(size = 6), axis.text.y = element_text(size = 6), 
        axis.title.y = element_text(size = 6), strip.text = element_text(size = 6))


Richness<-ggplot(divers_m, aes(x=BIOTURBATION, y=N, fill=DEPTH)) + 
  geom_boxplot(notch=TRUE)+ 
  facet_grid(~MONTH2)+
  scale_fill_manual(values=c("dodgerblue1","orange1","snow4"))+
  xlab("Bioturbation intensity") + 
  ylab("Richness (J)") +
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6), 
        axis.title.x = element_text(size = 6), axis.text.y = element_text(size = 6), 
        axis.title.y = element_text(size = 6), strip.text = element_text(size = 6))


##  analysis of qPCR
#import and preparation of the data
qpcr<-read.table("~/Dropbox/Articoli/PLOT/NEW_analysis/qpcr/qpcr_mature.txt",header=T,sep="\t")
qpcr$Depth<-factor(qpcr$Depth,level=c("S","SS","D"), labels=c("Surface","SubSurface","Deep"))
qpcr$Stage<-factor(qpcr$Stage,level=c("MATURE","JUVENILE"), labels=c("MATURE","JUVENILE"))
qpcr<-na.omit(qpcr)
qpcr_m<-subset(qpcr,Stage=='MATURE')
qpcr_m<-subset(qpcr_m,BioT!='NONE')
qpcr_m$MONTH<-factor(qpcr_m$MONTH,levels = c("May-16","Aug-16","Nov-16","Feb-17","May-17","Aug-17"))
qpcr_m$BioT<-factor(qpcr_m$BioT,levels = c("HIGH","LOW"), labels = c("HIGH","NORMAL"))

#tests
summary(qpcr_m)
model2<-manylm(log(Copies)~BioT*Depth*MONTH, qpcr_m)
plotty<-factor(qpcr_m$PLOT)
anova(model2, block = 'plotty', resamp='case', test="LR")

#graph
bacterio_a<-ggplot(qpcr_m,aes(BioT,log(Copies),fill=Depth))+ 
  geom_boxplot(notch=TRUE)+
  scale_fill_manual(values=c("dodgerblue1","orange1","snow4"))+
  ylab("Log Bacterial 16S Copies / gr sediment")+ 
  facet_grid(~MONTH) +
  xlab("Bioturbation intensity") +
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6), 
        axis.title.x = element_text(size = 6), axis.text.y = element_text(size = 6), 
        axis.title.y = element_text(size = 6), strip.text = element_text(size = 6))

##  analysis of Fluorescine diacetate
#importing and formatting data
fda<-read.table("~/Dropbox/Articoli/PLOT/NEW_analysis/fda/FDA_mature.txt",header=T,sep="\t")
fda$Depth<-factor(fda$Depth,level=c("S","SS","D"), labels=c("Surface","SubSurface","Deep"))
fda$MONTH<-factor(fda$MONTH,levels = c("May-16","Aug-16","Nov-16","Feb-17","May-17","Aug-17"),labels=c("May-16","Aug-16","Nov-16","Feb-17","May-17","Aug-17"))
fda$SITE<-factor(fda$SITE,levels = c("MATURE","JUVENILE"), labels=c("MATURE","JUVENILE"))

fda$Crab<-factor(fda$Crab,levels = c("HIGH","MEDIUM"), labels = c("HIGH","NORMAL"))

fda<-na.omit(fda)

fda_m<-subset(fda,SITE=='MATURE')
fda_m<-subset(fda_m,Crab!="NONE")
#testing
summary(fda_m)
model2<-manylm(log(Fluorescein)~Crab*Depth*MONTH, fda_m,offset = fda_m[rowSums()])
plotty<-factor(fda_m$Plot)
anova(model2, block = 'plotty',resamp='case', test="LR")

#graph
fda_plot<-ggplot(fda_m, aes(x=Crab, y=Fluorescein, fill=Depth)) + 
  geom_boxplot(notch=TRUE)+
  facet_grid(~MONTH)+
  scale_fill_manual(values = c("dodgerblue1","orange1","snow4")) +
  xlab("Bioturbation intensity") +
  ylab(bquote("Fluorescein (mg "*g^-1*"sediment 3"*h^-1*")")) +
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6), 
        axis.title.x = element_text(size = 6), axis.text.y = element_text(size = 6), 
        axis.title.y = element_text(size = 6), strip.text = element_text(size = 6))


##  Assemblage Figure S5

plot_grid(Shannon,Richness,bacterio_a,fda_plot, ncol=1, nrow=4, align="v",axis = 'l')




##  analysis of bacterial function

funcmat<-read.table("~/Dropbox/Articoli/PLOT/NEW_analysis/function/func_table_reduced_MATURE.txt",header=T,sep="\t")
funcmat<-na.omit(funcmat)
funcmat<-subset(funcmat, funcmat$BIOTURBATION!="NONE")
funcmat$DEPTH<-factor(funcmat$DEPTH,level=c("S","SS","D"), labels=c("Surface","SubSurface","Deep"))
funcmat$BIOTURBATION<-factor(funcmat$BIOTURBATION,levels = c("HIGH","LOW"), labels = c("HIGH","NORMAL"))

attach(funcmat)

str(funcmat)

func.matrix.distance.bray.mat<- vegdist(funcmat[, c(8:41)], method="bray")
func_test<-adonis2(log10(funcmat[, c(8:41)]+1) ~ DEPTH* MONTH *BIOTURBATION,data=funcmat, block=PLOT,na.rm=TRUE)

funcmat$MONTH<-factor(funcmat$MONTH,levels = c("May-16","Aug-16","Nov-16","Feb-17","May-17","Aug-17"))

photoauto<-ggplot(funcmat, aes(x=BIOTURBATION, y=photoautotrophy, fill=DEPTH)) + 
  geom_boxplot(notch=TRUE)+
  facet_grid(~MONTH)+
  scale_fill_manual(values=c("dodgerblue1","orange1","snow4"))+
  xlab("Bioturbation Intensity") +
  ylab(bquote('Photoautotrophy')) + 
  theme_classic()+
   theme(legend.position = "none",
        axis.text.x = element_text(size = 6), 
        axis.title.x = element_text(size = 6), axis.text.y = element_text(size = 6), 
        axis.title.y = element_text(size = 6), strip.text = element_text(size = 6))

sulfate_resp<-ggplot(funcmat, aes(x=BIOTURBATION, y=sulfate_respiration, fill=DEPTH)) + 
  geom_boxplot(notch=TRUE)+
  facet_grid(~MONTH)+
  scale_fill_manual(values=c("dodgerblue1","orange1","snow4"))+
  xlab("Bioturbation Intensity") +
  ylab(bquote('Sulfate respiration')) +
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6), 
        axis.title.x = element_text(size = 6), axis.text.y = element_text(size = 6), 
        axis.title.y = element_text(size = 6), strip.text = element_text(size = 6))

Cell<-ggplot(funcmat, aes(x=BIOTURBATION, y=cellulolysis, fill=DEPTH)) + 
  geom_boxplot(notch=TRUE)+
  facet_grid(~MONTH)+
  scale_fill_manual(values=c("dodgerblue1","orange1","snow4"))+
  xlab("Bioturbation Intensity") +
  ylab(bquote('Cellulolysis')) + 
  theme_classic()+
 theme(legend.position = "none",
        axis.text.x = element_text(size = 6), 
        axis.title.x = element_text(size = 6), axis.text.y = element_text(size = 6), 
        axis.title.y = element_text(size = 6), strip.text = element_text(size = 6))


Nitri<-ggplot(funcmat, aes(x=BIOTURBATION, y=nitrification, fill=DEPTH)) + 
  geom_boxplot(notch=TRUE)+
  facet_grid(~MONTH)+
  scale_fill_manual(values=c("dodgerblue1","orange1","snow4"))+
  xlab("Bioturbation Intensity") +
  ylab(bquote('Nitrification')) + 
  theme_classic()+
   theme(legend.position = "none",
        axis.text.x = element_text(size = 6), 
        axis.title.x = element_text(size = 6), axis.text.y = element_text(size = 6), 
        axis.title.y = element_text(size = 6), strip.text = element_text(size = 6))

nitrate_red<-ggplot(funcmat, aes(x=BIOTURBATION, y=nitrate_reduction, fill=DEPTH)) + 
  geom_boxplot(notch=TRUE)+
  facet_grid(~MONTH)+
  scale_fill_manual(values=c("dodgerblue1","orange1","snow4"))+
  xlab("Bioturbation Intensity") +
  ylab(bquote('Nitrate reduction')) + 
  theme_classic()+
   theme(legend.position = "none",
        axis.text.x = element_text(size = 6), 
        axis.title.x = element_text(size = 6), axis.text.y = element_text(size = 6), 
        axis.title.y = element_text(size = 6), strip.text = element_text(size = 6))

photoheterotrophy<-ggplot(funcmat, aes(x=BIOTURBATION, y=photoheterotrophy, fill=DEPTH)) + 
  geom_boxplot(notch=TRUE)+
  facet_grid(~MONTH)+
  scale_fill_manual(values=c("dodgerblue1","orange1","snow4"))+
  xlab("Bioturbation Intensity") +
  ylab(bquote('Photoheterotrophy')) + 
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6), 
        axis.title.x = element_text(size = 6), axis.text.y = element_text(size = 6), 
        axis.title.y = element_text(size = 6), strip.text = element_text(size = 6))

Nfixation<-ggplot(funcmat, aes(x=BIOTURBATION, y=nitrogen_fixation, fill=DEPTH)) + 
  geom_boxplot(notch=TRUE)+
  facet_grid(~MONTH)+
  scale_fill_manual(values=c("dodgerblue1","orange1","snow4"))+
  xlab("Bioturbation Intensity") +
  ylab(bquote('Nitrogen fixation')) + 
  theme_classic()+
   theme(legend.position = "none",
        axis.text.x = element_text(size = 6), 
        axis.title.x = element_text(size = 6), axis.text.y = element_text(size = 6), 
        axis.title.y = element_text(size = 6), strip.text = element_text(size = 6))

Sulphi_respiration<-ggplot(funcmat, aes(x=BIOTURBATION, y=sulfite_respiration, fill=DEPTH)) + 
  geom_boxplot(notch=TRUE)+
  facet_grid(~MONTH)+
  scale_fill_manual(values=c("dodgerblue1","orange1","snow4"))+
  xlab("Bioturbation Intensity") +
  ylab(bquote('Sulphite respiration')) + 
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6), 
        axis.title.x = element_text(size = 6), axis.text.y = element_text(size = 6), 
        axis.title.y = element_text(size = 6), strip.text = element_text(size = 6))


plot_grid(photoauto,sulfate_resp,Cell,Nitri,nitrate_red,photoheterotrophy ,Nfixation,Sulphi_respiration, ncol=2, nrow=4, align="v",axis = 'l')






##  Plant performance

setwd("~/Dropbox/Articoli/PLOT/NEW_analysis/plantgrowth")


matplant1<-read.table("~/Dropbox/Articoli/PLOT/NEW_analysis/plantgrowth/plant_mat1.txt",header=T,sep="\t")
matplant2<-read.table("~/Dropbox/Articoli/PLOT/NEW_analysis/plantgrowth/plant_mat2.txt",header=T,sep="\t")


matplant1$MONTH2<-factor(matplant1$MONTH2,level=c("Feb-16","May-16","Aug-16","Nov-16","Feb-17","May-17","Aug-17"),labels=c("Feb-16","May-16","Aug-16","Nov-16","Feb-17","May-17","Aug-17"))



#ancova with biot as continuous

model.pneu<-aov(log(PNEUM+1)~BIOT*MONTH, matplant1)
anova(model.pneu)

model.pneu2<-lm(log(PNEUM+1)~BIOT*tempo, matplant1)
anova(model.pneu2)
summary(model.pneu2)
plot(model.pneu2)
model.height<-aov(log(HEIGHT+1)~BIOT*MONTH, matplant1)
anova(model.height)

model.height2<-lm(log(HEIGHT+1)~BIOT*tempo, matplant1)
summary(model.height2)

matplant2$MONTH2<-factor(matplant2$MONTH2,level=c("Feb-16","May-16","Aug-16","Nov-16","Feb-17","May-17","Aug-17"),labels=c("Feb-16","May-16","Aug-16","Nov-16","Feb-17","May-17","Aug-17"))

matplant2<-na.omit(matplant2)

#ancova with biot as continuous
summary(matplant2)
model.branch<-aov(log(BRANCH+1)~BIOTURBATION*MONTH2, matplant2)
anova(model.branch)

model.branch2<-lm(log(BRANCH+1)~BIOTURBATION*tempo, matplant2)
anova(model.branch2)
summary(model.branch2)

##plot

pneu<-ggplot(matplant1, aes(x=MONTH2, y=PNEUM, fill=BIOT)) + 
  geom_boxplot(notch=TRUE)+
  xlab("")+
  ylab(bquote('Number of pneumatophores\n'))+
  theme_bw()+
  scale_fill_manual(values = c("red4", "red"), name = "Bioturbation intensity")+
  theme(panel.border = element_blank(), legend.position = "none",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 10), 
        axis.title.x = element_text(size = 10), axis.text.y = element_text(size = 10), 
        axis.title.y = element_text(size = 10))

pneu2<-ggplot(matplant1, aes(x=tempo, y=PNEUM,color=BIOT)) + 
  scale_color_manual(values = c("red4", "red"), name = "Bioturbation intensity")+
  geom_point()+
  xlab("Days")+
  ylab(bquote('Number of pneumatophores\n'))+
  theme_bw()+
  geom_smooth(method = "glm")+
  theme(panel.border = element_blank(), legend.position = "none",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 10), 
        axis.title.x = element_text(size = 10), axis.text.y = element_text(size = 10), 
        axis.title.y = element_text(size = 10))


height<-ggplot(matplant1, aes(x=MONTH2, y=HEIGHT, fill=BIOT)) + 
  geom_boxplot(notch=TRUE)+
  xlab("")+
  ylab(bquote('Height (cm)'))+
  theme_bw()+
  scale_fill_manual(values = c("red4", "red"), name = "Bioturbation intensity")+
  theme(panel.border = element_blank(), legend.position = "none",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 10), 
        axis.title.x = element_text(size = 10), axis.text.y = element_text(size = 10), 
        axis.title.y = element_text(size = 10))

height2<-ggplot(matplant1, aes(x=tempo, y=HEIGHT, color=BIOT)) + 
  geom_point()+
  xlab("Days")+
  ylab(bquote('Height (cm)'))+
  theme_bw()+
  geom_smooth(method = "glm")+
  scale_color_manual(values = c("red4", "red"), name = "Bioturbation intensity")+
  theme(panel.border = element_blank(), legend.position = "none",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 10), 
        axis.title.x = element_text(size = 10), axis.text.y = element_text(size = 10), 
        axis.title.y = element_text(size = 10))

branch<-ggplot(matplant2, aes(x=MONTH2, y=BRANCH, fill=BIOTURBATION)) + 
  geom_boxplot(notch=TRUE)+
  xlab("Time")+
  ylab(bquote('Branch diameter (mm)'))+
  theme_bw()+
  scale_fill_manual(values = c("red4", "red"), name = "Bioturbation intensity")+
  theme(panel.border = element_blank(), legend.position = "none",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 10), 
        axis.title.x = element_text(size = 10), axis.text.y = element_text(size = 10), 
        axis.title.y = element_text(size = 10))
branch2<-ggplot(matplant2, aes(x=tempo, y=BRANCH, color=BIOTURBATION)) + 
  geom_point()+
  xlab("Days")+
  ylab(bquote('Branch diameter (mm)'))+
  theme_bw()+
  geom_smooth(method = "glm")+
  scale_color_manual(values = c("red4", "red"), name = "Bioturbation intensity")+
  theme(panel.border = element_blank(), legend.position = "none",
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 10), 
        axis.title.x = element_text(size = 10), axis.text.y = element_text(size = 10), 
        axis.title.y = element_text(size = 10))

plot_grid(pneu,pneu2,height,height2,branch,branch2, align="v",axis = 'l',ncol = 2, labels = c("A","B","C","D","E","F"))


##  Correlation test

library(vegan)
library(ecodist)
library(readr)
library(dplyr)
library(phyloseq)
library(vegan)
library(ade4)
library(ggplot2)
library(lavaan)


# Mature
maturo_map<-read.table("../06_SEM/map_mat3.txt",header=T,sep='\t')
maturo<-read.table("../06_SEM/OTU_mat2.txt",header=T,sep='\t',row.names = 1)
nrow(maturo)

maturo.na<-na.omit(maturo)
map_veg<-maturo_map[,10:12]
map_crab<-maturo_map[,8]
map_chemical<-maturo_map[,14:27]
map_pcr<-maturo_map[,13]

maturo.bray <- vegdist(maturo, method = "bray") # dissimilarity matrix using bray-curtis distance indices on the varespec dataset native to vegan

pcoaVS <- pco(maturo.bray, negvals = "zero", dround = 0) # if negvals = 0 sets all negative eigenvalues to zero; if = "rm" corrects for negative eigenvalues using method 1 of Legendre and Anderson 1999
componentPCO<-as.data.frame(pcoaVS$vectors)

maturo2<-cbind(maturo_map,componentPCO$V1,componentPCO$V2,componentPCO$V3)

names(maturo2)
colnames(maturo2)[30]<-"Component1"
colnames(maturo2)[31]<-"Component2"
colnames(maturo2)[32]<-"Component3"

attach(maturo2)
names(maturo2)
summary(maturo2$Nitrite_mat)

x<-prcomp(maturo2[c(14,27,31,33:43)])
pca_maturo<-as.data.frame(x$x)
maturo2$Nitrite2<-scale(maturo2$Nitrite)
maturo2$Phosphate2<-scale(maturo2$Phosphate)
maturo2$Silicate2<-scale(maturo2$Silicate)
maturo2$Nitrate2<-scale(maturo2$Nitrate)
maturo2$Chloride2<-scale(maturo2$Chloride)
maturo2$Sulphate2<-scale(maturo2$Sulphate)
maturo2$TPC2<-scale(maturo2$TPC)
maturo2$POC2<-scale(maturo2$POC)
maturo2$PIC2<-scale(maturo2$PIC)
maturo2$PIN2<-scale(maturo2$PIN)
maturo2$PON2<-scale(maturo2$PON)
maturo2$TPN2<-scale(maturo2$TPN)
maturo2$Component12<-scale(maturo2$Component1)
maturo2$Component22<-scale(maturo2$Component2)
maturo2$Component32<-scale(maturo2$Component3)
maturo2$Salinity2<-scale(maturo2$Salinity)
maturo2$pH2<-scale(maturo2$pH)


maturo2$PCAChem<-scale(pca_maturo$PC1)



maturo2$qpcr2<-scale(maturo2$qpcr)

model <-'
# LATENT
Plant=~MATHEIGHT+PNEUMATOPHORE+BRANCHDIAM
# STRUCTURAL
Plant~TPC2+PON2+PIC2
Plant~CRABNO+Component1+H+qpcr2
Component1~CRABNO+pH2+TPC2+PON2+PIC2+POC2+TPC2+Silicate2+Phosphate2
H~CRABNO+Component1+TPC2+PIN2+PIC2+POC2+TPC2+Sulphate2+Nitrite2
qpcr2~H+CRABNO+Component1+Salinity2+TPC2+PON2+PIC2+POC2+TPC2+Silicate2+Phosphate2
TPC2~CRABNO
PON2~CRABNO
PIN2~CRABNO
PIC2~CRABNO
POC2~CRABNO
Nitrate2~CRABNO
Phosphate2~CRABNO
Nitrite2~CRABNO'

fit <- sem(model, data = maturo2)
summary(fit, standardized=TRUE, fit.measures=TRUE, rsquare=TRUE)


maturo.dist <- distance(maturo,method = "bray-curtis") # Bray-Curtis
map_veg.dist <- distance(map_veg, "euclidean")
c.dist <- vegdist(map_crab, "euclidean")
pcr.dist<-vegdist(map_pcr, "euclidean")
map_chem.dist <- distance(map_veg, "euclidean")

vegan::mantel(c.dist, maturo.dist, permutations = 999)    
vegan::mantel(maturo.dist, map_chem.dist, permutations = 999)
vegan::mantel(maturo.dist, map_veg.dist, permutations = 999)
vegan::mantel(maturo.dist, pcr.dist, permutations = 999)
vegan::mantel(map_chem.dist, pcr.dist, permutations = 999)
vegan::mantel(map_veg.dist, pcr.dist, permutations = 999)
vegan::mantel(c.dist, pcr.dist, permutations = 999)
vegan::mantel(c.dist, map_veg.dist, permutations = 999)
vegan::mantel(map_chem.dist, map_veg.dist, permutations = 999)
vegan::mantel(c.dist, map_chem.dist, permutations = 999)






