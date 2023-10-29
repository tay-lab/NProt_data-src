#Load libraries

library(dplyr)
library(Matrix)
library(ggplot2)
library(patchwork) 
library(bioanalyzeR)


my.theme <- theme(panel.background = element_rect(fill='transparent'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title = element_text(size = 7), axis.text.x = element_text(angle = 0, hjust=0.5),
                  axis.text = element_text(size=7, color='black'), plot.title = element_text(size=7, face="plain", hjust = 0.5), plot.background = element_rect(fill='transparent'),
                 legend.position = "none")


#Load Nanodrop data

nanodrop_DBCO <- read.csv(data.dir = "DBCO.csv")
nanodrop_CD147 <- read.csv(data.dir = "CD147-DBCO.csv")


#Load Tapestation data for Smart-seq PLA

Plate <- read.electrophoresis(data.dir = "Tapestation_Smart-seq_PLA.xml")
Plate$assay.info
Plate$samples


#Load Tapestation data for 10X 

cDNApreAmp <- read.electrophoresis(data.dir = "Tapestation_10X_cDNApreAmp.xml")
cDNApreAmp$assay.info
cDNApreAmp$samples

Final10X <- read.electrophoresis(data.dir = "Tapestation_10X_FInal.xml")
Final10X$assay.info
Final10X$samples

PLApreAmp <- read.electrophoresis(data.dir = "Tapestation_10X_PLApreAmp.xml")
PLApreAmp$assay.info
PLApreAmp$samples


#Identify the required samples

cDNA <- subset(cDNApreAmp, sample.name == "Nicky 2.3A cDNA 10F")
cDNA_PLC <- subset(Final10X, sample.name == "Nicky_cDNA_PLC_10F")
PLA <- subset(PLApreAmp, sample.name == "Nicky_PLA2.3b_10F")
PLA_PLC <- subset(Final10X, sample.name == "Nicky_PLA_PLC_10F")
Plate_PLA <- subset(Plate, sample.name == "P8 10F")



#Plots

margin0 <- theme(plot.margin = unit(c(1,1,1,1), "mm"))

p1 <- ggplot(nanodrop_DBCO, aes(x= Wavelength, y = Absorbance),geom = c("line"))+ ggtitle("")  + theme_bw() + my.theme +geom_line(aes(x= Wavelength, y = Absorbance), colour='#00BF7D',size=1)+ xlim (NA, 800) + ylim (NA, 8) 

p2 <- ggplot(nanodrop_CD147, aes(x= Wavelength, y = Absorbance),geom = c("line"))+ ggtitle("") + theme_bw() + my.theme +geom_line(aes(x= Wavelength, y = Absorbance), colour='#00BF7D',size=1)+ xlim (NA, 800) + ylim (-0.2, 6)


p3 <- qplot.electrophoresis(cDNA, x= "length", y = "fluorescence", xlim = c(200, 5000),ylim = c(NA, 400),geom = c("line"),facets = NULL, title = "10X cDNA after pre-amplification")+ theme_bw() + my.theme +geom_line(aes(x= length, y = fluorescence), colour='#00BF7D',size=1)

p4 <- qplot.electrophoresis(cDNA_PLC, x= "length", y = "fluorescence", xlim = c(NA, 3000),ylim = c(NA, 800),geom = c("line"),facets = NULL, title = "10X cDNA Final")+ theme_bw() + my.theme +geom_line(aes(x= length, y = fluorescence), colour='#00BF7D',size=1)

p5 <- qplot.electrophoresis(PLA, x= "length", y = "fluorescence", xlim = c(NA, 1000),ylim = c(NA, 400),geom = c("line"),facets = NULL, title = "10X PLA after pre-amplification")+ theme_bw() + my.theme +geom_line(aes(x= length, y = fluorescence), colour='#00BF7D',size=1)

p6 <- qplot.electrophoresis(PLA_PLC, x= "length", y = "fluorescence", xlim = c(NA, 1000),ylim = c(NA, 450),geom = c("line"),facets = NULL, title = "10X PLA Final")+ theme_bw() + my.theme +geom_line(aes(x= length, y = fluorescence), colour='#00BF7D',size=1)

p7 <- qplot.electrophoresis(Plate_PLA, x= "length", y = "fluorescence", xlim = c(NA, 1000),ylim = c(NA, 450),geom = c("line"),facets = NULL, title = "Smart-seq PLA Final")+ theme_bw() + my.theme +geom_line(aes(x= length, y = fluorescence), colour='#00BF7D',size=1)


fig <- grid.arrange(p1,p2,p3,p4,p5,p6,p7, ncol=2, nrow=4)


