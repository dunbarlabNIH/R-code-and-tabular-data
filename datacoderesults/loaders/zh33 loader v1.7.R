#Samson Koelle
#copywrite Cynthia Dunbar Lab
#10-18-16
#This script loads data from animal ZH33

wkdirsave = getwd()
#set this directory to data containing folder
setwd("~/data")
#read main data
zh33 = read.delim("zh33_independent_100_forpub.txt")

#read readme
zh33readme = read.delim('zh33readmeforpub.txt')

#read 0 bp mismatch 100 thresh data
#zh330bp100thresh = read.delim('zh33_0bp_extracted_20160725.txt')
zh33100thresh0bp = read.delim("zh330bp100threshforpub.txt")
#read key file containing proper names and sampleIDs
#zh33keyfile = read.delim("zh33_keyfile_20160725.txt")
zh33keyfile = read.delim("zh33keyfileforpub.txt")

#order data and readme
tomatch = match(zh33keyfile$FILENAME, colnames(zh33))
fillervalue  = min(setdiff(seq(from = 1, to = 100, by = 1), tomatch[!is.na(tomatch)]))
tomatch[is.na(tomatch)] = fillervalue
zh33ordered = zh33[,tomatch]

tomatch = match(zh33keyfile$FILENAME, colnames(zh33100thresh0bp))
fillervalue = min(setdiff(seq(from = 1, to = 100, by = 1), tomatch[!is.na(tomatch)]))
tomatch[is.na(tomatch)] = fillervalue
zh330bp100threshordered = zh33100thresh0bp[,tomatch]

tomatch = match(zh33keyfile$FILENAME, zh33readme$FILENAME)
fillervalue = min(setdiff(seq(from = 1, to = 100, by = 1), tomatch[!is.na(tomatch)]))
tomatch[is.na(tomatch)] = fillervalue
zh33readmeordered = zh33readme[tomatch,]

#add null column to data
zh33ordered = cbind(zh33ordered, list(rep(0, length(zh33[[1]]))))
zh330bp100threshordered = cbind(zh330bp100threshordered , list(rep(0, length(zh330bp100threshordered[[1]]))))

#doesn't work for publication because of data masking
#set rownames of barcode abundance matrix equal to barcode names
#rownames(zh33ordered) = zh33[,1]
#rownames(zh330bp100threshordered) = zh330bp100thresh[,1]

#prepublication index configuration
zh33indexmatrix = matrix(nrow = 15, ncol = 4, data = c(1,6,11,20,41,61, 75,90,129, 205, 236,351, 360,385,391,
                                                       2,7,12,21,42,62, 76, 91,130,  206, 237, 350,357,372,386,
                                                       #3,8,13,22, 43,63,77,92, 131,  207, 238, 343,NA,379,389,
                                                       4,9,14,23,44, 64,78,93, 132, 208, 245,353, 358,378,388,
                                                       5,10,15,24, 45, 65,79,94, 133, 209, 240, 352,362,373,387))

#publication index configuration
#zh33indexmatrix = matrix(nrow = 15, ncol = 4, data = c(1:15,16:30,31:42,44:45,46:60,61:75))

#create index matrix with 0 column instead of NA
zh33indexmatrixnoNA = zh33indexmatrix
zh33indexmatrixnoNA[is.na(zh33indexmatrix)] = length(zh33)

#save data with no threshold applied
zh33alldatanothresh <- zh33ordered
zh330bp100threshalldatanothresh <- zh330bp100threshordered

#apply base threshold which replicates preselection of samples in fastq pipeline (prepublication only)
zh33alldatanothresh = basethresholdinsamples(data = zh33alldatanothresh, samples = zh33indexmatrixnoNA, thresh = 100)
zh330bp100threshalldatanothresh = basethresholdinsamples(data = zh330bp100threshalldatanothresh,
                                                         samples = zh33indexmatrixnoNA, thresh = 100)

#apply threshold
zh33alldata = thresholdinsamples(data = zh33alldatanothresh, samples = zh33indexmatrixnoNA, thresh = 2000)
zh330bp100threshalldata = thresholdinsamples(data = zh330bp100threshalldatanothresh, samples = zh33indexmatrixnoNA, thresh = 2000)

#create vector of sample names
zh33allnames = zh33keyfile$GIVENNAME

#get names
zh33allnames = append(as.vector(zh33allnames), "NULL SAMPLE")
 
#create vector of timepoints
zh33timepointnames = c("1","2","3", "4.5","6.5","9.5","12","14","21","28","30", '38','43','46','49')
zh33timepoints = c(1,2,3,4.5,6.5,9.5,12,14,21,28,30,38,43,46,49)

#import GFP levels
zh33GFPraw = read.delim("zh33 GFP levels 050316.txt")
rownames(zh33GFPraw) = zh33GFPraw[[1]]
zh33GFP = (zh33GFPraw[c(2,3,4,5,6,7,8,9,11,14,16,17,18,19,20),c(4,5,2,3)])
colnames(zh33GFP) = c("T","B",'Mono',"Gr")

setwd(wkdirsave)
