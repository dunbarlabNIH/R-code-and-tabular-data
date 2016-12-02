#Samson Koelle
#Copywrite Cynthia Dunbar Group
#R code for Blood paper
#This script loads data and calls source code to recreate results from the paper.
#You should be able to copy and paste it into your R console to get results (after installing R packages and setting working directories)

#First, we set some parameters
ppi = 300
size = 8
set.seed(1)

#set the following variables
workingdir = "~/results"
sourcedir = '~/source'

#load required packages
library(pheatmap)
library(DiversitySampler)
library(RColorBrewer)
library(nnet)
library(foreach)
library(stringr)
library(biclust)
library(scales)

#Load functions. Function descriptions are available within the functions themselves
setwd(sourcedir)
#thresholding functions
source('basethresholdinsamples.R')
source('thresholdinsamples.R')

#functions for Figure 2
source('diversitytrack.R')
source('clonenumbertrack.R')
source('GFPtrack.R')

#functions for Figure 3
source("cor_analysisv1.3.R")
source("BCLH v1.3.R")

#functions for Figure 4
source("cluster_curveovertime v1.2.R")
source("autocorrelation.R")

#functions for Figure 6
source('biastracking3.R')
source('rollingvarianceofbias.R')

#functions for Supplemental Figure 1
source('readtrack v1.2.R')
source('threshcomparelatesttimepoint.R')

#functions for Supplemental Figure 2
source('cumulativebarcodefrequency.R')
source('contamcheck v1.3.R')
source('diversitycurve.R')
source('cnattpsensitivitybyanimal.R')
source('diversitysensitivitybyanimal.R')

#functions for Supplemental Figure 3
source('diversityplot.R')

#functions for Supplemental Figure 6
source("TCcontrib v1.2.R")
#functions for Supplemental Figure 7
source('BCLHbinarizeall v1.2.R')

#Response to reviewer
source('compare2bpwith0bp.R')

#other
source('howmuchpolyclonal.R')
# source('clonenumbertrack.R')
# source("stabilityplot.R")
# source("threshold.R")
# source('biasheatmapandtracking.R')
# source('GFPcomparison.R')
# source('dist_analysis.R')
# source('clonetracking.R')

#Load Data
setwd("/Users/samsonkoelle/Desktop/Barcode Projects/DataCodeResults113016final/loaders")
source("zh33 loader v1.6.R")
source("zg66 loader v1.6.R")
source("zj31 loader v1.6.R")
source("zh19 loader v1.6.R")
#load library data
#source("Library loader v1.1.R")

#Make Figures
setwd(workingdir)
#Create images for Figure 2
dirname = paste("Figure 2 Images", Sys.time())
dir.create(dirname)
setwd(dirname)
#Compare Shannon Diversities
diversitytrack(data = zh33alldata, indexmatrix = zh33indexmatrix, timepoints = zh33timepoints, folder = "ZH33")
diversitytrack(data = zg66alldata, indexmatrix = zg66indexmatrix, timepoints = zg66timepoints, folder = "ZG66")
diversitytrack(data = zh19alldata, indexmatrix = zh19indexmatrix, timepoints = zh19timepoints, folder = "ZH19")
diversitytrack(data = zj31alldata, indexmatrix = zj31indexmatrix, timepoints = zj31timepoints, folder = "ZJ31")
#GFP tracking
GFPtrack(data = zh33GFP, timepoints = zh33timepoints, folder = "ZH33")
GFPtrack(data = zg66GFP, timepoints = zg66timepoints, folder = "ZG66")
GFPtrack(data = zh19GFP, timepoints = zh19timepoints, folder = "ZH19")
GFPtrack(data = zj31GFP, timepoints = zj31timepoints[c(1,2,3,4,5,6,8,9,10,11)], folder = "ZJ31")
#Compare clone numbers
clonenumbertrack(data = zh33alldata, indexmatrixnoNA = zh33indexmatrixnoNA, timepoints = zh33timepoints, folder = "ZH33")
clonenumbertrack(data = zg66alldata, indexmatrixnoNA = zg66indexmatrixnoNA, timepoints = zg66timepoints, folder = "ZG66")
clonenumbertrack(data = zh19alldata, indexmatrixnoNA = zh19indexmatrixnoNA, timepoints = zh19timepoints, folder = "ZH19")
clonenumbertrack(data = zj31alldata, indexmatrixnoNA = zj31indexmatrixnoNA, timepoints = zj31timepoints, folder = "ZJ31")
setwd(workingdir)

#Create Images for Figure 3
dirname = paste("Figure 3 Images", Sys.time())
dir.create(dirname)
setwd(dirname)
#analyze Pearson correlations
cor_analysis(data = zh33alldata, indexmatrix = zh33indexmatrix, timepoints = zh33timepoints, folder = "ZH33" , celltypenames = c("T", "B","Mono",'Gr'))
cor_analysis(data = zg66alldata, indexmatrix = zg66indexmatrix, timepoints = zg66timepoints, folder = "ZG66" , celltypenames = c("T", "B","Mono",'Gr'))
#make top clone heatmap
BCLH(data = zh33alldata[,as.vector(zh33indexmatrixnoNA)], names = zh33allnames[as.vector(zh33indexmatrixnoNA)], n_clones = 10, folder = paste("ZH33"))
setwd(workingdir)

#Create Images for Figure 4
dirname = paste("Figure 4 Images", Sys.time())
dir.create(dirname)
setwd(dirname)
#make autocorrelation plots
autocorrelation(data = zh33alldata[,as.vector(zh33indexmatrixnoNA[,4])], names = zh33timepointnames, timepoints = zh33timepoints, folder = "ZH33 Gr")
autocorrelation(data = zg66alldata[,as.vector(zg66indexmatrixnoNA[,4])], names = zg66timepointnames, timepoints = zg66timepoints, folder = "ZG66 Gr")
#make top clone heatmaps
BCLH(data = zh33alldata[,as.vector(zh33indexmatrixnoNA[,4])], names = zh33timepoints,n_clones = 100, folder = "ZH33 Gr")
BCLH(data = zg66alldata[,as.vector(zg66indexmatrixnoNA[,4])], names = zg66timepoints, n_clones = 100, folder = "ZG66 Gr")
#make cumulative distribution plots
cluster_curveovertime(data = zh33alldata[,as.vector(zh33indexmatrixnoNA[,4])],names = zh33timepointnames, folder = "ZH33 Gr", thresh = 2000)
cluster_curveovertime(data = zg66alldata[,as.vector(zg66indexmatrixnoNA[,4])], names = zg66timepointnames, folder = "ZG66 Gr", thresh = 2000)
setwd(workingdir)

#Create Images for Figure 5
dirname = paste("Figure 5 Images", Sys.time())
dir.create(dirname)
setwd(dirname)
#make autocorrelation plots
autocorrelation(data = zh33alldata[,as.vector(zh33indexmatrixnoNA[,1])], names = zh33timepoints, timepoints = zh33timepoints, folder = "ZH33 T")
autocorrelation(data = zh33alldata[,as.vector(zh33indexmatrixnoNA[,2])], names = zh33timepoints, timepoints = zh33timepoints, folder = "ZH33 B")
#make top clone heatmaps
BCLH(data = zh33alldata[,as.vector(zh33indexmatrixnoNA[,1])], names = zh33timepoints,n_clones = 100, folder = "ZH33 T")
BCLH(data = zh33alldata[,as.vector(zh33indexmatrixnoNA[,2])], names = zh33timepoints,n_clones = 100, folder = "ZH33 B")
#make cumulative distribution plots
cluster_curveovertime(data = zh33alldata[,as.vector(zh33indexmatrixnoNA[,1])], names = zh33timepoints, folder = "ZH33 T", thresh = 2000)
cluster_curveovertime(data = zh33alldata[,as.vector(zh33indexmatrixnoNA[,2])],names = zh33timepoints, folder = "ZH33 B", thresh = 2000)
setwd(workingdir)

#Create Images for Figure 6
dirname = paste("Figure 6", Sys.time())
dir.create(dirname)
setwd(dirname)
#make bias tracking plots
biastracking(data = zh33alldata, uppersamples = zh33indexmatrix[,2],lowersamples =  zh33indexmatrix[,4], timepointnames = zh33timepointnames,folder = "ZH33 B v Gr",upperbiaslabel = "B", lowerbiaslabel = 'Gr')
biastracking(data = zh33alldata, uppersamples = zh33indexmatrix[,1],lowersamples =  zh33indexmatrix[,2], timepointnames = zh33timepointnames,folder = "ZH33 T v B",upperbiaslabel = "T", lowerbiaslabel = 'B')
biastracking(data = zh33alldata, uppersamples = zh33indexmatrix[,3], lowersamples =  zh33indexmatrix[,4], timepointnames = zh33timepointnames,folder = "ZH33 Mono v Gr",upperbiaslabel = "Mo", lowerbiaslabel = 'Gr')
#make rolling variance of bias plots
rollingvarianceofbias(data = zh33alldata, uppersamples = zh33indexmatrixnoNA[,3],lowersamples = zh33indexmatrixnoNA[,4], folder = "ZH33 Mono v Gr", timepoints = zh33timepoints) 
rollingvarianceofbias(data = zh33alldata, uppersamples = zh33indexmatrixnoNA[,1],lowersamples = zh33indexmatrixnoNA[,2], folder = "ZH33 T v B", timepoints = zh33timepoints)
rollingvarianceofbias(data = zh33alldata, uppersamples = zh33indexmatrixnoNA[,2],lowersamples = zh33indexmatrixnoNA[,4], folder = "ZH33 B v Gr", timepoints = zh33timepoints) 
setwd(workingdir)

#Create Images for Supplemental Figure 1
dirname = paste("Figure S1", Sys.time())
dir.create(dirname)
setwd(dirname)
#make figures showing effect of 100 thresh on long term granulocyte samples
threshcomparelatesttimepoint()
#track read abundances at different processing steps.
readtrack(indexmatrix = zh33indexmatrix, folder = 'ZH33', timepoints = zh33timepoints, rawreads = zh33readmeordered$READS, readswithLibID = zh33readmeordered$MAPPED, readswithLibIDover100 = apply(FUN = sum, MARGIN = 2, X = zh33alldatanothresh),readswithLibIDover100overthresh = apply(FUN = sum, MARGIN = 2, X = zh33alldata),celltypenames = c('T',"B",'Gr',"Mono"))
readtrack(indexmatrix = zg66indexmatrix, folder = 'ZG66',timepoints = zg66timepoints,rawreads = zg66readmeordered$READS,readswithLibID = zg66readmeordered$MAPPED, readswithLibIDover100 = apply(FUN = sum, MARGIN = 2, X = zg66alldatanothresh),readswithLibIDover100overthresh = apply(FUN = sum, MARGIN = 2, X = zg66alldata),celltypenames = c('T',"B",'Gr',"Mono"))
readtrack(indexmatrix = zh19indexmatrix, folder = 'ZH19',timepoints = zh19timepoints,rawreads = zh19readmeordered$READS,readswithLibID = zh19readmeordered$MAPPED,readswithLibIDover100 = apply(FUN = sum, MARGIN = 2, X = zh19alldatanothresh),readswithLibIDover100overthresh = apply(FUN = sum, MARGIN = 2, X = zh19alldata),celltypenames = c('T',"B",'Gr',"Mono"))
readtrack(indexmatrix = zj31indexmatrix, folder = 'ZJ31',timepoints = zj31timepoints,rawreads = zj31readmeordered$READS,readswithLibID = zj31readmeordered$MAPPED,readswithLibIDover100 = apply(FUN = sum, MARGIN = 2, X = zj31alldatanothresh),readswithLibIDover100overthresh = apply(FUN = sum, MARGIN = 2, X = zj31alldata),celltypenames = c('T',"B",'Gr',"Mono"))
setwd(workingdir)

#Create Images for Supplemental Figure S2
dirname = paste("Figure S2", Sys.time())
dir.create(dirname)
setwd(dirname)
#plot overall clonal distributions
diversitycurve(datalist = list(zh33alldata,zg66alldata,zh19alldata,zj31alldata),datalistnothresh = list(zh33alldatanothresh, zg66alldatanothresh,zh19alldatanothresh, zj31alldatanothresh),indexmatrices = list(zh33indexmatrix, zg66indexmatrix,zh19indexmatrix, zj31indexmatrix),folder = "DIVCURVE")
#plot barcodes in all animals
contamcheck(data = list(zh33alldata,zg66alldata,zh19alldata,zj31alldata), folder = "post-threshold",indexmatricesnoNA = list(zh33indexmatrixnoNA, zg66indexmatrixnoNA,zh19indexmatrixnoNA, zj31indexmatrixnoNA))
#plot the sensitivity of clone number to threshold in each animal
cnattpsensitivitybyanimal(data = list(zh33alldatanothresh,zg66alldatanothresh,zh19alldatanothresh,zj31alldatanothresh),indexmatricesnoNA = list(zh33indexmatrixnoNA, zg66indexmatrixnoNA,zh19indexmatrixnoNA, zj31indexmatrixnoNA),folder = "Clone number at timepoint sensitivity", threshseq = seq(from = 0, to = 10000, by = 1000))
#plot the sensitivity of Shannon diversity to threshold in each animal
diversitysensitivitybyanimal(data = list(zh33alldatanothresh,zg66alldatanothresh,zh19alldatanothresh,zj31alldatanothresh),indexmatricesnoNA = list(zh33indexmatrixnoNA, zg66indexmatrixnoNA,zh19indexmatrixnoNA, zj31indexmatrixnoNA),folder = "Diversity Sensitivity", threshseq = seq(from = 0, to = 10000, by = 1000))
#plot the cumulative frequencies of the barcodes
cumulativebarcodefrequency(data = zh33alldatanothresh[,as.vector(zh33indexmatrix[,4])], folder = 'ZH33 Gr')
cumulativebarcodefrequency(data = zh19alldatanothresh[,as.vector(zh19indexmatrix[,4])], folder = 'ZH19 Gr')
cumulativebarcodefrequency(data = zg66alldatanothresh[,as.vector(zg66indexmatrix[,4])], folder = 'ZG66 Gr')
cumulativebarcodefrequency(data = zj31alldatanothresh[,as.vector(zj31indexmatrix[,4])], folder = 'ZJ31 Gr')
setwd(workingdir)

#Create Images for Figure S3
dirname = paste("Figure S3", Sys.time())
dir.create(dirname)
setwd(dirname)
#show and top 10 contribution
diversityplot(data = zh33alldata[,as.vector(zh33indexmatrixnoNA[,4])], folder = "ZH33 Gr", thresh = 2000, names = zh33timepoints)
diversityplot(data = zh33alldata[,as.vector(zh33indexmatrixnoNA[,1])], folder = "ZH33 T", thresh = 2000, names = zh33timepoints)
diversityplot(data = zh33alldata[,as.vector(zh33indexmatrixnoNA[,2])], folder = "ZH33 B", thresh = 2000, names = zh33timepoints)
diversityplot(data = zh33alldata[,as.vector(zh33indexmatrixnoNA[,3])], folder = "ZH33 Mono", thresh = 2000, names = zh33timepoints)
diversityplot(data = zg66alldata[,as.vector(zg66indexmatrixnoNA[,4])], folder = "ZG66 Gr", thresh = 2000, names = zg66timepoints)
diversityplot(data = zg66alldata[,as.vector(zg66indexmatrixnoNA[,1])], folder = "ZG66 T", thresh = 2000, names = zg66timepoints)
diversityplot(data = zg66alldata[,as.vector(zg66indexmatrixnoNA[,2])], folder = "ZG66 B", thresh = 2000, names = zg66timepoints)
diversityplot(data = zg66alldata[,as.vector(zg66indexmatrixnoNA[,3])], folder = "ZG66 Mono", thresh = 2000, names = zg66timepoints)
diversityplot(data = zh19alldata[,as.vector(zh19indexmatrixnoNA[,4])], folder = "ZH19 Gr", thresh = 2000, names = zh19timepoints)
diversityplot(data = zh19alldata[,as.vector(zh19indexmatrixnoNA[,1])], folder = "ZH19 T", thresh = 2000, names = zh19timepoints)
diversityplot(data = zh19alldata[,as.vector(zh19indexmatrixnoNA[,2])], folder = "ZH19 B", thresh = 2000, names = zh19timepoints)
diversityplot(data = zh19alldata[,as.vector(zh19indexmatrixnoNA[,3])], folder = "ZH19 Mono", thresh = 2000, names = zh19timepoints)
diversityplot(data = zj31alldata[,as.vector(zj31indexmatrixnoNA[,4])], folder = "ZJ31 Gr", thresh = 2000, names = zj31timepoints)
diversityplot(data = zj31alldata[,as.vector(zj31indexmatrixnoNA[,1])], folder = "ZJ31 T", thresh = 2000, names = zj31timepoints)
diversityplot(data = zj31alldata[,as.vector(zj31indexmatrixnoNA[,2])], folder = "ZJ31 B", thresh = 2000, names = zj31timepoints)
diversityplot(data = zj31alldata[,as.vector(zj31indexmatrixnoNA[,3])], folder = "ZJ31 Mono", thresh = 2000, names = zj31timepoints)
setwd(workingdir)

#Create Images for Figure S4
dirname = paste("Figure S4", Sys.time())
dir.create(dirname)
setwd(dirname)
#make correlation tracking plots in other animals
cor_analysis(data = zh19alldata, indexmatrix = zh19indexmatrix, timepoints = zh19timepoints,folder = paste("ZH19"))
cor_analysis(data = zj31alldata, indexmatrix = zj31indexmatrix, timepoints = zj31timepoints,folder = paste("ZJ31"))
#make top clone heatmaps in other animals
BCLH(data = zh19alldata[,as.vector(zh19indexmatrix)], names = zh19allnames[as.vector(zh19indexmatrix)],n_clones = 10,folder = paste("ZH19"))
BCLH(data = zg66alldata[,as.vector(zg66indexmatrixnoNA)], names = zg66allnames[as.vector(zg66indexmatrixnoNA)], n_clones = 10,folder = paste("ZG66"))
BCLH(data = zj31alldata[,as.vector(zj31indexmatrixnoNA)], names = zj31allnames[as.vector(zj31indexmatrixnoNA)],n_clones = 10,folder = paste("ZJ31"))
setwd(workingdir)

#Create Images for Figure S5
dirname = paste("Figure S5", Sys.time())
dir.create(dirname)
setwd(dirname)
#make top 100 heatmaps
BCLH(data = zh19alldata[,as.vector(zh19indexmatrix)], names = zh19allnames[as.vector(zh19indexmatrix)],n_clones = 100,folder = paste("ZH19 top 100"), cexset = .3)
BCLH(data = zg66alldata[,as.vector(zg66indexmatrixnoNA)], names = zg66allnames[as.vector(zg66indexmatrixnoNA)], n_clones = 100,folder = paste("ZG66 top 100"), cexset = .3)
BCLH(data = zj31alldata[,as.vector(zj31indexmatrixnoNA)], names = zj31allnames[as.vector(zj31indexmatrixnoNA)],n_clones = 100,folder = paste("ZJ31 top 100"), cexset = .3)
BCLH(data = zh33alldata[,as.vector(zh33indexmatrixnoNA)], names = zh33allnames[as.vector(zh33indexmatrixnoNA)],n_clones = 100,folder = paste("ZH33 top 100"), cexset = .3)
setwd(workingdir)

#Create Images for Figure S6
dirname = paste("Figure S6", Sys.time())
dir.create(dirname)
setwd(dirname)
#make top clone contribution plots
TCcontrib(data = zh19alldata, indexmatrix = zh19indexmatrix, folder  = paste("Union of top 10 clones in ZH19 samples"),n_clones = 10, timepoints =  zh19timepoints) 
TCcontrib(data = zh33alldata, indexmatrix = zh33indexmatrix, folder  = paste("Union of top 10 clones in ZH33"), n_clones = 10, timepoints =  zh33timepoints) 
TCcontrib(data = zg66alldata, indexmatrix = zg66indexmatrixnoNA, folder  = paste("Union of top 10 clones in ZG66"),n_clones = 10, timepoints =  zg66timepoints) 
TCcontrib(data = zj31alldata, indexmatrix = zj31indexmatrixnoNA, folder  = paste("Union of top 10 clones in ZJ31"), n_clones = 10, timepoints =  zj31timepoints) 
TCcontrib(data = zh19alldata, indexmatrix = zh19indexmatrix, folder  = paste("Union of top 100 clones in ZH19"), n_clones = 100, timepoints =  zh19timepoints) 
TCcontrib(data = zh33alldata, indexmatrix = zh33indexmatrixnoNA, folder  = paste("Union of top 100 clones in ZH33"), n_clones = 100, timepoints =  zh33timepoints) 
TCcontrib(data = zj31alldata, indexmatrix = zj31indexmatrixnoNA, folder  = paste("Union of top 100 clones in ZJ31"), n_clones = 100, timepoints =  zj31timepoints) 
TCcontrib(data = zg66alldata, indexmatrix = zg66indexmatrixnoNA, folder  = paste("Union of top 100 clones in ZG66"),n_clones = 100, timepoints =  zg66timepoints) 

setwd(workingdir)

#Create Images for Figure S7
dirname = paste("Figure S7", Sys.time())
dir.create(dirname)
setwd(dirname)

#Binarized tracking
BCLHbinarizeall(data = zh33alldata[,as.vector(zh33indexmatrixnoNA[,4])],   names = zh33timepoints,folder = paste("ZH33 Gr binarized"))
BCLHbinarizeall(data = zh33alldata[,as.vector(zh33indexmatrixnoNA[,1])],   names = zh33timepoints,folder = paste("ZH33 T binarized"))
BCLHbinarizeall(data = zh33alldata[,as.vector(zh33indexmatrixnoNA[,2])],  names = zh33timepoints, folder = paste("ZH33 B binarized"))
BCLHbinarizeall(data = zh33alldata[,as.vector(zh33indexmatrixnoNA[,3])],  names = zh33timepoints, folder = paste("ZH33 Mono binarized"))
BCLHbinarizeall(data = zh19alldata[,as.vector(zh19indexmatrixnoNA[,4])],  names = zh19timepoints,  folder = paste("ZH19 Gr binarized"))
BCLHbinarizeall(data = zh19alldata[,as.vector(zh19indexmatrixnoNA[,1])],  names = zh19timepoints,  folder = paste("ZH19 T binarized"))
BCLHbinarizeall(data = zh19alldata[,as.vector(zh19indexmatrixnoNA[,2])],   names = zh19timepoints,   folder = paste("ZH19 B binarized"))
BCLHbinarizeall(data = zh19alldata[,as.vector(zh19indexmatrixnoNA[,3])],   names = zh19timepoints,   folder = paste("ZH19 Mono binarized"))
BCLHbinarizeall(data = zg66alldata[,as.vector(zg66indexmatrixnoNA[,4])],   names = zg66timepoints, folder = paste("ZG66 Gr binarized"))
BCLHbinarizeall(data = zg66alldata[,as.vector(zg66indexmatrixnoNA[,1])],  names = zg66timepoints,  folder = paste("ZG66 T binarized"))
BCLHbinarizeall(data = zg66alldata[,as.vector(zg66indexmatrixnoNA[,2])],   names = zg66timepoints, folder = paste("ZG66 B binarized"))
BCLHbinarizeall(data = zg66alldata[,as.vector(zg66indexmatrixnoNA[,3])],  names = zg66timepoints,  folder = paste("ZG66 Mono binarized"))
BCLHbinarizeall(data = zj31alldata[,as.vector(zj31indexmatrixnoNA[,4])],  names = zj31timepoints,  folder = paste("ZJ31 Gr binarized"))
BCLHbinarizeall(data = zj31alldata[,as.vector(zj31indexmatrixnoNA[,1])], names = zj31timepoints,folder = paste("ZJ31 T binarized"))
BCLHbinarizeall(data = zj31alldata[,as.vector(zj31indexmatrixnoNA[,2])],  names = zj31timepoints, folder = paste("ZJ31 B binarized"))
BCLHbinarizeall(data = zj31alldata[,as.vector(zj31indexmatrixnoNA[,3])],  names = zj31timepoints, folder = paste("ZJ31 Mono binarized"))
setwd(workingdir)

#Create Images for Response to Reviewers
dirname = paste("Response to Reviewers", Sys.time())
dir.create(dirname)
setwd(dirname)

compare2bpwith0bp(
  twobpalldatalist = list(zh33alldata, zg66alldata, zh19alldata, zj31alldata),
  zerobpalldata = list(zh330bp100threshalldata, zg660bp100threshalldata, zh190bp100threshalldata, zj310bp100threshalldata),
  indexmatrices = list(zh33indexmatrixnoNA, zg66indexmatrixnoNA, zh19indexmatrixnoNA, zj31indexmatrixnoNA),
  folder = 'twobpzerobpcompare', animalnames = c('ZH33','ZG66','ZH19','ZJ31'))

howmuchpolyclonal(data = zh33alldata, indexmatrix = zh33indexmatrix)
howmuchpolyclonal(data = zg66alldata, indexmatrix = zg66indexmatrix)
howmuchpolyclonal(data = zh19alldata, indexmatrix = zh19indexmatrix)
howmuchpolyclonal(data = zj31alldata, indexmatrix = zj31indexmatrix)

