#Samson Koelle
#copywrite Cynthia Dunbar Lab
#10-18-16
#This script loads data from animal ZJ31
wkdirsave = getwd()

#set this directory to data containing folder
setwd("~/data")

#read main data
zj31 = read.delim("zj31_independent_100_forpub.txt")
#read readme
zj31readme = read.delim('zj31readmeforpub.txt')
#read 0 bp mismatch 100 thresh data
zj31100thresh0bp = read.delim('zj310bp100threshforpub.txt')
#read key file containing proper names and sampleIDs
zj31keyfile = read.delim("zj31keyfileforpub.txt")

#order data and readme
tomatch = match(zj31keyfile$FILENAME, colnames(zj31))
fillervalue  = min(setdiff(seq(from = 1, to = 100, by = 1), tomatch[!is.na(tomatch)]))
tomatch[is.na(tomatch)] = fillervalue
zj31ordered = zj31[,tomatch]

tomatch = match(zj31keyfile$FILENAME, colnames(zj31100thresh0bp))
fillervalue = min(setdiff(seq(from = 1, to = 100, by = 1), tomatch[!is.na(tomatch)]))
tomatch[is.na(tomatch)] = fillervalue
zj310bp100threshordered = zj31100thresh0bp[,tomatch]

tomatch = match(zj31keyfile$FILENAME, zj31readme$FILENAME)
fillervalue = min(setdiff(seq(from = 1, to = 100, by = 1), tomatch[!is.na(tomatch)]))
tomatch[is.na(tomatch)] = fillervalue
zj31readmeordered = zj31readme[tomatch,]

#add null column to data
zj31ordered = cbind(zj31ordered, list(rep(0, length(zj31[[1]]))))
zj310bp100threshordered = cbind(zj310bp100threshordered, list(rep(0, length(zj310bp100threshordered[[1]]))))

#doesnt work in publication version (data masking)
#set rownames of barcode abundance matrix equal to barcode names
#rownames(zj31ordered) = zj31[,1]
#rownames(zj310bp100threshordered) = zj310bp100thresh[,1]

#prepublication index configuration
zj31indexmatrix = matrix(nrow = 11, ncol = 4, data = c(32,49,93, 144,153,183,208,241,13,252,260,
                                                       28,34,84,136,145,157,204,218,1,249,255,
                                                       #31,43,97,142,149,173,206,227,10,251,261,
                                                       NA,40, 95,141,152,170,205,225, 9,250,257,
                                                       30, 39,131,139, 151,166,198,223,7,14,256))
#create index matrix with 0 column instead of NA
zj31indexmatrixnoNA = zj31indexmatrix
zj31indexmatrixnoNA[is.na(zj31indexmatrix)] = length(zj31)
#save data with no threshold applied
zj31alldatanothresh <- zj31ordered
zj310bp100threshalldatanothresh <- zj310bp100threshordered
#apply base threshold which replicates preselection of samples in fastq pipeline (prepublication only)
zj31alldatanothresh = basethresholdinsamples(data = zj31alldatanothresh, samples = zj31indexmatrixnoNA, thresh = 100)
zj310bp100threshalldatanothresh = basethresholdinsamples(data = zj310bp100threshalldatanothresh,
                                                         samples = zj31indexmatrixnoNA, thresh = 100)
#apply threshold
zj31alldata = thresholdinsamples(data = zj31alldatanothresh, samples = zj31indexmatrixnoNA, thresh = 2000)
zj310bp100threshalldata = thresholdinsamples(data = zj310bp100threshalldatanothresh, samples = zj31indexmatrixnoNA, thresh = 2000)
#get names
zj31allnames = zj31keyfile$GIVENNAME
#create vector of timepoints
zj31timepointnames = c("1","2","3.5", "4",'5',"6",'8.5',"9.5","12",'17.5', '20')
zj31timepoints = c(1,2,3.5,4,5,6,8.5,9.5,12,17.5, 20)
#import GFP levels
zj31GFPraw = read.delim("zj31 GFP levels 050316.txt")
rownames(zj31GFPraw) = zj31GFPraw[[1]]
zj31GFP = (zj31GFPraw[c(1,2,3,4,5,6,7,8,9,10),c(4,5,2,3)])
colnames(zj31GFP) = c("T","B",'Mono',"Gr")

setwd(wkdirsave)

