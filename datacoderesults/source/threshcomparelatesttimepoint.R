#Samson Koelle
#Copywrite Cynthia Dunbar Lab
#11-16-16

#read 2bp mismatch 0 thresh
threshcomparelatesttimepoint = function(){
workingdirtemp = getwd()

setwd('~/data/')
zj31_20m_0thresh_2bp = read.delim('zj31_20m_0_20160829.txt')
zh33_49m_0thresh_2bp= read.delim('zh33_49m_0_20160829.txt')
zh19_36m_0thresh_2bp= read.delim('zh19_36m_0_20160829.txt')
zg66_42m_0thresh_2bp= read.delim('zg66_42m_0_20160829.txt')


#####
# setwd("/Users/samsonkoelle/Desktop/Barcode Projects/DataCodeResults113016final/data")
# 
# zg66 = read.delim("zg66_independent_100_20160829.txt")
# zg66redo = read.delim('zg66_independent_100_20161129.txt')
# zg66latesttimepoint = read.delim('zg66_42m_100_20160829.txt')
# zg66keyfile = read.delim("zg66_keyfile_20160725.txt")
# 
# zg66indexmatrix  = matrix(nrow = 13, ncol = 4, data = c(18,23,28,33, 38,48,53,58, 71,63,257,264,271,
#                                                         19,24,29,34,39,49,54,59,67,64,252,251,265,
#                                                         21,26,31,36,41,51,56,61, 70,NA,256,261,267,
#                                                         22,27,32,37,42,132,57,62, 68,66,255,259,266))
# 
# zg66ordered = zg66[,match(zg66keyfile$FILENAME, colnames(zg66), nomatch = length(zg66))]
# zg66redoordered = zg66redo[,match(zg66keyfile$FILENAME, colnames(zg66redo), nomatch = length(zg66redo))]
# 
# dim(basethresholdinsamples(data = zg66latesttimepoint, samples = 2:5, thresh = 100))
# dim(basethresholdinsamples(data = zg66ordered, samples = zg66indexmatrixnoNA[dim(zg66indexmatrixnoNA)[1],], thresh = 100))
# dim(basethresholdinsamples(data = zg66redoordered, samples = zg66indexmatrixnoNA[dim(zg66indexmatrixnoNA)[1],], thresh = 100))
# 
# 
# zh19 = read.delim("zh19_independent_100_20160829.txt")
# zh19redo = read.delim('zh19_independent_100_20161129.txt')
# zh19latesttimepoint = read.delim('zh19_42m_100_20160829.txt')
# zh19keyfile = read.delim("zh19_keyfile_20160725.txt")
# 
# zh19indexmatrix  = matrix(nrow = 13, ncol = 4, data = c(18,23,28,33, 38,48,53,58, 71,63,257,264,271,
#                                                         19,24,29,34,39,49,54,59,67,64,252,251,265,
#                                                         21,26,31,36,41,51,56,61, 70,NA,256,261,267,
#                                                         22,27,32,37,42,132,57,62, 68,66,255,259,266))
# 
# zh19ordered = zh19[,match(zh19keyfile$FILENAME, colnames(zh19), nomatch = length(zh19))]
# zh19redoordered = zh19redo[,match(zh19keyfile$FILENAME, colnames(zh19redo), nomatch = length(zh19redo))]
# 
# dim(basethresholdinsamples(data = zh19latesttimepoint, samples = 1:5, thresh = 100))
# dim(basethresholdinsamples(data = zh19ordered, samples = zh19indexmatrixnoNA[dim(zh19indexmatrixnoNA)[1],], thresh = 100))
# dim(basethresholdinsamples(data = zh19redoordered, samples = zh19indexmatrixnoNA[dim(zh19indexmatrixnoNA)[1],], thresh = 100))
# 

####

setwd(workingdirtemp)
rownames(zj31_20m_0thresh_2bp) = zj31_20m_0thresh_2bp[,1]
rownames(zh33_49m_0thresh_2bp) = zh33_49m_0thresh_2bp[,1]
rownames(zh19_36m_0thresh_2bp) = zh19_36m_0thresh_2bp[,1]
rownames(zg66_42m_0thresh_2bp) = zg66_42m_0thresh_2bp[,1]
zh19_36m_0thresh_2bp = zh19_36m_0thresh_2bp[-which(rownames(zh19_36m_0thresh_2bp) %in% contam),]

#lots of excluded <100 level clones mean threshold will be higher in non-excluded samples
#note that the applied threshold will then be lower for more diverse libraries (since more <100 clones mean a higher threshold if those clones are included)
#this leads to the counterintuitive result that using the 100 threshold actually leads to MORE clones being included in the final analysis
zj31_20m_100thresh_2bp_gr = basethresholdinsamples(data = zj31ordered, samples = zj31indexmatrixnoNA[dim(zj31indexmatrixnoNA)[1],], thresh = 100)
zj31100threshgr = thresholdinsamples(data = zj31_20m_100thresh_2bp_gr, samples = 256, thresh = 2000)
zj310threshgr = thresholdinsamples(data = zj31_20m_0thresh_2bp, samples = 5, thresh = 2000)
zj31merge= merge(zj310threshgr, zj31100threshgr, all = T, by = 'row.names')
zj31merge[is.na(zj31merge)] = 0

png(paste("ZJ31 preprocess.png"),width=1.4*size*ppi,height=size*ppi,res=ppi)
par(new = F, mar = c(7,6,3,7), xpd = TRUE)
plot(zj31merge[[6]]/sum(zj31merge[[6]]), zj31merge[[262]]/sum(zj31merge[[262]]), xlab = '', ylab = '', cex.axis = 2)
mtext("Preprocessing Threshold = 0", side=1, line=3, cex = 2)
mtext("Preprocessing Threshold = 100", side = 2, line = 3, cex = 2)
mtext("ZJ31 20m Gr", side=3, line=1, cex = 3)
dev.off()

zh33_49m_100thresh_2bp_gr = basethresholdinsamples(data = zh33ordered, samples = zh33indexmatrixnoNA[dim(zh33indexmatrixnoNA)[1],], thresh = 100)
zh33100threshgr = thresholdinsamples(data = zh33_49m_100thresh_2bp_gr, samples = 387, thresh = 2000)
zh330threshgr = thresholdinsamples(data = zh33_49m_0thresh_2bp, samples = 4, thresh = 2000)
zh33merge= merge(zh330threshgr, zh33100threshgr, all = T, by = 'row.names')
zh33merge[is.na(zh33merge)] = 0

png(paste("ZH33 preprocess.png"),width=1.4*size*ppi,height=size*ppi,res=ppi)
par(new = F, mar = c(7,6,3,7), xpd = TRUE)
plot(zh33merge[[5]]/sum(zh33merge[[5]]), zh33merge[[393]]/sum(zh33merge[[393]]), xlab = '', ylab = '')
mtext("Preprocessing Threshold = 0", side=1, line=3, cex = 2)
mtext("Preprocessing Threshold = 100", side = 2, line = 3, cex = 2)
mtext("ZH33 49m Gr", side=3, line=1, cex = 3)
dev.off()

zh19_36m_100thresh_2bp_gr = basethresholdinsamples(data = zh19ordered, samples = zh19indexmatrixnoNA[dim(zh19indexmatrixnoNA)[1],], thresh = 100)
zh19100threshgr = thresholdinsamples(data = zh19_36m_100thresh_2bp_gr, samples = 161, thresh = 2000)
zh190threshgr = thresholdinsamples(data = zh19_36m_0thresh_2bp, samples = 5, thresh = 2000)
zh19merge= merge(zh190threshgr, zh19100threshgr, all = T, by = 'row.names')
zh19merge[is.na(zh19merge)] = 0

png(paste("ZH19 preprocess.png"),width=1.4*size*ppi,height=size*ppi,res=ppi)
par(new = F, mar = c(7,6,3,7), xpd = TRUE)
plot(zh19merge[[6]]/sum(zh19merge[[6]]), zh19merge[[167]]/sum(zh19merge[[167]]), xlab = '', ylab = '')
mtext("Preprocessing Threshold = 0", side=1, line=3, cex = 2)
mtext("Preprocessing Threshold = 100", side = 2, line = 3, cex = 2)
mtext("ZH19 36m Gr", side=3, line=1, cex = 3)
dev.off()

zg66_42m_100thresh_2bp_gr = basethresholdinsamples(data = zg66ordered, samples = zg66indexmatrixnoNA[dim(zg66indexmatrixnoNA)[1],], thresh = 100)
zg66100threshgr = thresholdinsamples(data = zg66_42m_100thresh_2bp_gr, samples = 266, thresh = 2000)
zg660threshgr = thresholdinsamples(data = zg66_42m_0thresh_2bp, samples = 2, thresh = 2000)
zg66merge= merge(zg660threshgr, zg66100threshgr, all = T, by = 'row.names')
zg66merge[is.na(zg66merge)] = 0

png(paste("ZG66 preprocess.png"),width=1.4*size*ppi,height=size*ppi,res=ppi)
par(new = F, mar = c(7,6,3,7), xpd = TRUE)
plot(zg66merge[[3]]/sum(zg66merge[[3]]), zg66merge[[272]]/sum(zg66merge[[272]]), xlab = '', ylab = '')
mtext("Preprocessing Threshold = 0", side=1, line=3, cex = 2)
mtext("Preprocessing Threshold = 100", side = 2, line = 3, cex = 2)
mtext("ZG66 42m Gr", side=3, line=1, cex = 3)
dev.off()
}
