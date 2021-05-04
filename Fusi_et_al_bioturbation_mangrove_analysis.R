# Packages needed to run this code
library(vegan)
library('BiodiversityR')
library(devtools)
library(pairwiseAdonis)
library(devtools)
library(ggord)

# LOAD DATA AND PREPARE IT FOR ANALISIS
otu_mat<-read.table("path/OTU_mat.txt",header=T,sep="\t",row.names = 1)
map_mat<-read.table("path/map_mat.txt",header=T,sep="\t")

# Transformation of the data for homogeneity
otu_mat<-log10(otu_mat+1)

#This is for testing the full set of explanatory variable
your.matrix.distance.BC.mat<- vegdist(otu_mat, method="bray")
mat<-adonis2(your.matrix.distance.BC.mat ~ DEPTH*MONTH*CRABNO, data = map_mat, block=PLOT)

# For strata, extract factors into a new data.frame
fac <- data.frame(Depth=map_mat$DEPTH, MONTH = map_mat$MONTH)
pairwise.adonis2(your.matrix.distance.BC.mat~MONTH,data=fac,strata='MONTH')
pairwise.adonis2(your.matrix.distance.BC.mat~Depth,data=fac,strata='Depth')






your.matrix.distance.BC <- vegdist(otu_mat, method="bray")
adonis(your.matrix.distance.BC ~ SITE, data = map_mat)

source("~/Downloads/parwise.adonis.r")
pairwise.adonis(x, factors, sim.method = "bray", p.adjust.m = "none")



# Homogeneity of dispersion test
beta <- betadisper(your.matrix.distance.BC.juv, map_juv$CRABNO)
permutest(beta)
plot(beta)
boxplot(beta)


#analysis of pH
pH<-read.table("~/Dropbox/Articoli/PLOT/NEW_analysis/biochemistry/pH.txt",header=T,sep="\t")
pH$DEPTH<-factor(pH$DEPTH,level=c("S","SS","D"), labels=c("S","SS","D"))
pH$SITE<-factor(pH$SITE,level=c("MATURE","JUVENILE"), labels=c("MATURE","JUVENILE"))
pH<-na.omit(pH)
pH_m<-subset(pH,SITE=='MATURE')


summary(pH_m)
model2<-manylm(pH~DEPTH*factor(MONTH), pH_m)
plot(model2,1)
plotty<-factor(pH_m$PLOT)
anova(model2, block = 'plotty', resamp='case', test="F")
summary(model2, block = 'plotty', resamp='case', test="F")

#analysis of Salinity
sal<-read.table("~/Dropbox/Articoli/PLOT/NEW_analysis/biochemistry/SALINITY.txt",header=T,sep="\t")
sal$DEPTH<-factor(sal$DEPTH,level=c("S","SS","D"), labels=c("S","SS","D"))
sal$SITE<-factor(sal$SITE,level=c("MATURE","JUVENILE"), labels=c("MATURE","JUVENILE"))
sal<-na.omit(sal)
sal_m<-subset(sal,SITE=='MATURE')


summary(sal_m)
model2<-manylm(SALINITY~DEPTH*factor(MONTH), sal_m)
plot(model2,1)
plotty<-factor(sal_m$PLOT)
anova(model2, block = 'plotty', resamp='case', test="F")
summary(model2, block = 'plotty', resamp='case', test="F")


# Multivariate analysis of the biochemical parameters
biochem<-read.table("~/Dropbox/Articoli/PLOT/NEW_analysis/biochemistry/BIOCHEMISTRY.txt",header=T,sep="\t")

biochem1<-subset(biochem, biochem$SITE=="MATURE")
attach(biochem1)

library(mvabund)
biochem.mod <- manylm(cbind(pH,Nitrite,Phosphate,Silicate,Nitrate,POC,PIC,PIN,PON) ~ DEPTH * MONTH2 * CRABNUMBER,
                  data=biochem1)
plotty<-factor(PLOT)
summary(biochem.mod, block = 'plotty', resamp='case', test="LR")
plot.manylm(biochem.mod, which = c(1,2))



#  analysis of alpha-diversity indexes
divers<-read.table("~/Dropbox/Articoli/PLOT/NEW_analysis/divrich/divrich.txt",header=T,sep="\t")
divers$DEPTH<-factor(divers$DEPTH,level=c("S","SS","D"), labels=c("S","SS","D"))
divers$SITE<-factor(divers$SITE,level=c("MATURE","JUVENILE"), labels=c("MATURE","JUVENILE"))
divers<-na.omit(divers)
divers_m<-subset(divers,SITE=='MATURE')


summary(divers_m)
model2<-manylm(log(H)~CRABNO*DEPTH*factor(MONTH), divers_m)
plot(model2,1)
plotty<-factor(divers_m$PLOT)
anova(model2, block = 'plotty', resamp='case', test="F")
summary(model2, block = 'plotty', resamp='case', test="F")

model2<-manylm(log(J)~CRABNO*DEPTH*factor(MONTH), divers_m)
plot(model2,1)
plotty<-factor(divers_m$PLOT)
anova(model2, block = 'plotty', resamp='case', test="F")
summary(model2, block = 'plotty', resamp='case', test="F")


#  analysis of qPCR
qpcr<-read.table("~/Dropbox/Articoli/PLOT/NEW_analysis/qpcr/qpcr_mature.txt",header=T,sep="\t")
qpcr$Depth<-factor(qpcr$Depth,level=c("S","SS","D"), labels=c("S","SS","D"))
qpcr$Stage<-factor(qpcr$Stage,level=c("MATURE","JUVENILE"), labels=c("MATURE","JUVENILE"))
qpcr<-na.omit(qpcr)
qpcr_m<-subset(qpcr,Stage=='MATURE')


summary(qpcr_m)
model2<-manylm(log(Copies)~CRABNO*Depth*MONTH, qpcr_m)
plot(model2,1)
plotty<-factor(divers_m$PLOT)
anova(model2, block = 'plotty', resamp='case', test="LR")
summary(model2, block = 'plotty', resamp='case', test="LR")

model2<-manylm(log(J)~CRABNO*DEPTH*factor(MONTH), divers_m)
plot(model2,1)
plotty<-factor(divers_m$PLOT)
anova(model2, block = 'plotty', resamp='case', test="F")
summary(model2, block = 'plotty', resamp='case', test="F")



#  analysis of bacterial function
funcmat<-read.table("~/Dropbox/Articoli/PLOT/NEW_analysis/function/func_table_reduced_MATURE.txt",header=T,sep="\t")
funcmat<-na.omit(funcmat)
attach(funcmat)

str(funcmat)
aa<-cbind(methylotrophy,methanotrophy,aerobic_ammonia_oxidation,nitrification,sulfate_respiration,
          sulfur_respiration,sulfite_respiration,denitrification,chitinolysis,dark_hydrogen_oxidation,
          nitrogen_fixation,cellulolysis,xylanolysis,dark_sulfide_oxidation,dark_thiosulfate_oxidation,
          fermentation,animal_parasites_or_symbionts,aromatic_hydrocarbon_degradation,hydrocarbon_degradation,
          nitrate_reduction,predatory_or_exoparasitic,cyanobacteria,anoxygenic_photoautotrophy,photoheterotrophy,
          ureolysis,chemoheterotrophy)


is.numeric(photoautotrophy)
b<-as.numeric(CRABNO)
is.numeric(b)
library(mvabund)
func.mat.mod <- manylm( aa~ DEPTH* MONTH *b,data=funcmat)

anova(func.mat.mod, block = 'plotty', resamp='case', test="F")
summary(func.mat.mod, block = 'plotty', resamp='case', test="F")
plot(func.mat.mod, which=1)


#  analysis of CO2 fluxes
setwd("~/Dropbox/Articoli/PLOT/NEW_analysis/gasfluxes/")
CO2<-read.table("INSITUCO2.txt",header = T,sep="\t")
CO2_m<-subset(CO2,Site =='Mature')
CO2_m$MONTH2<-factor(CO2_m$Month2,level=c("May, 16","Aug, 16","Nov, 16","Feb, 17","May, 17","Aug, 17"),labels=c("May, 16","Aug, 16","Nov, 16","Feb, 17","May, 17","Aug, 17"))


summary(CO2_m)
model2<-manylm(CO2~Bioturbation+Month2+Bioturbation:Month2, CO2_m)
plot(model2)
summary(model2, block = 'plotty', resamp='case', test="F")
anova(model2, block = 'plotty', resamp='case', test="F")

CO2_m$Bioturbation<-factor(CO2_m$Bioturbation,level=c("B","NB"),labels=c("High Bioturbated","Not Bioturbated"))


p<-ggplot(CO2_m,aes(x=Month2,y=CO2,color=Bioturbation)) +
  facet_wrap(~Site) + theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 14), 
        axis.title.x = element_text(size = 14), axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14))+
  scale_x_discrete(limits=c("May, 16","Aug, 16","Nov, 16","Feb, 17","May, 17","Aug, 17"))+
  ylab(bquote('Flux of mmol of' ~CO[2]~ m^-2~d^-1))+
  xlab('Month')+
  ylim(-100,300)+
  geom_boxplot(notch = T)
p







#  analysis of Fluorescine diacetate
fda<-read.table("~/Dropbox/Articoli/PLOT/NEW_analysis/fda/FDA_mature.txt",header=T,sep="\t")
fda$Depth<-factor(fda$Depth,level=c("S","SS","D"), labels=c("S","SS","D"))
fda$MONTH<-factor(fda$MONTH,levels = c("May-16","Aug-16","Nov-16","Feb-17","May-17","Aug-17"),labels=c("May-16","Aug-16","Nov-16","Feb-17","May-17","Aug-17"))
fda$SITE<-factor(fda$SITE,levels = c("MATURE","JUVENILE"), labels=c("MATURE","JUVENILE"))
fda<-na.omit(fda)

fda_m<-subset(fda,SITE=='MATURE')

summary(fda_m)
model2<-manylm(log(Fluorescein)~CRABNO*Depth*MONTH, fda_m)
plot(model2,2)
summary(model2)
anova(model2)




