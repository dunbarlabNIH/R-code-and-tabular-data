#Samson Koelle
#copywrite Cynthia Dunbar Lab
#10-18-16
#This script loads data from animal ZG66

wkdirsave = getwd()
#set this directory to data containing folder
setwd("~/data")
#read main data
zg66 = read.delim("zg66_independent_100_forpub.txt")

#read readme
zg66readme = read.delim('zg66readmeforpub.txt')

#read 0 bp mismatch 100 thresh data
zg66100thresh0bp = read.delim('zg660bp100threshforpub.txt')

#read 2bp mismatch 0 thresh latest timepoint
#zg660thresh2bplatesttimepoints = read.delim('zg66_42m_0_20160829.txt')

#read key file containing proper names and sampleIDs
zg66keyfile = read.delim("zg66keyfileforpub.txt")

#order data and readme
#this method is adhoc for publication
tomatch = match(zg66keyfile$FILENAME, colnames(zg66))
fillervalue  = min(setdiff(seq(from = 1, to = 100, by = 1), tomatch[!is.na(tomatch)]))
tomatch[is.na(tomatch)] = fillervalue
zg66ordered = zg66[,tomatch]

tomatch = match(zg66keyfile$FILENAME, colnames(zg66100thresh0bp))
fillervalue = min(setdiff(seq(from = 1, to = 100, by = 1), tomatch[!is.na(tomatch)]))
tomatch[is.na(tomatch)] = fillervalue
zg660bp100threshordered = zg66100thresh0bp[,tomatch]

tomatch = match(zg66keyfile$FILENAME, zg66readme$FILENAME)
fillervalue = min(setdiff(seq(from = 1, to = 100, by = 1), tomatch[!is.na(tomatch)]))
tomatch[is.na(tomatch)] = fillervalue
zg66readmeordered = zg66readme[tomatch,]

#add null column to data
zg66ordered = cbind(zg66ordered, list(rep(0, length(zg66[[1]]))))
zg660bp100threshordered = cbind(zg660bp100threshordered , list(rep(0, length(zg660bp100threshordered[[1]]))))

#set rownames of barcode abundance matrix equal to barcode names
#doesn't work in publication version
#rownames(zg66ordered) = zg66[,1]
#rownames(zg660bp100threshordered) = zg66100thresh0bp[,1]
#zg66100thresh0bpordered = zg66100thresh0bp[,match(zg66keyfile$FILENAME, colnames(zg66))]

#index configuration
zg66indexmatrix  = matrix(nrow = 13, ncol = 4, data = c(18,23,28,33, 38,48,53,58, 71,63,257,264,271,
                                                        19,24,29,34,39,49,54,59,67,64,252,251,265,
                                                        21,26,31,36,41,51,56,61, 70,NA,256,261,267,
                                                        22,27,32,37,42,132,57,62, 68,66,255,259,266))
#create index matrix with 0 column instead of NA
zg66indexmatrixnoNA = zg66indexmatrix
zg66indexmatrixnoNA[is.na(zg66indexmatrix)] = length(zg66)

#save data with no threshold applied
zg66alldatanothresh = zg66ordered
zg660bp100threshalldatanothresh <- zg660bp100threshordered

#apply base threshold which replicates preselection of samples in fastq pipeline (prepublication only)
zg66alldatanothresh = basethresholdinsamples(data=  zg66alldatanothresh, samples = zg66indexmatrixnoNA, thresh = 100)
zg660bp100threshalldatanothresh = basethresholdinsamples(data = zg660bp100threshalldatanothresh,
                                                         samples = zg66indexmatrixnoNA, thresh = 100)
#apply threshold
zg66alldata = thresholdinsamples(data = zg66alldatanothresh, samples = zg66indexmatrixnoNA, thresh = 2000)
zg660bp100threshalldata = thresholdinsamples(data = zg660bp100threshalldatanothresh, samples = zg66indexmatrixnoNA, thresh = 2000)

#get names
zg66allnames = zg66keyfile$GIVENNAME
zg66allnames = append(as.vector(zg66allnames), "NULL SAMPLE")

#create vector of timepoints
zg66timepointnames = c("1","2","3","4.5","6.5","12","14.5","17","22","24","27", "36", '42')
zg66timepoints = c(1,2,3,4.5,6.5,12,14.5,17,22,24,27, 36, 42)

#import GFP levels
zg66GFPraw = read.delim("zg66 GFP levels 050316.txt")
rownames(zg66GFPraw) = zg66GFPraw[[1]]
zg66GFP = (zg66GFPraw[c(2,3,4,5,6,8,9,10,12,13,14,15,16),c(4,5,2,3)])
colnames(zg66GFP) = c("T","B",'Mono',"Gr")

setwd(wkdirsave)
