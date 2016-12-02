compare2bpwith0bp = function(#twobpalldata = zj31alldata,
                            #  zerobpalldata = zj310bp100threshalldata, indexmatrix = zj31indexmatrixnoNA, folder = 'test'){
  twobpalldatalist = list(zh33alldata, zg66alldata, zh19alldata, zj31alldata),
  zerobpalldata = list(zh330bp100threshalldata, zg660bp100threshalldata, zh190bp100threshalldata, zj310bp100threshalldata),
  indexmatrices = list(zh33indexmatrixnoNA, zg66indexmatrixnoNA, zh19indexmatrixnoNA, zj31indexmatrixnoNA),
  folder = 'twobpzerobpcompare', animalnames = c('ZH33','ZG66','ZH19','ZJ31')){
 
  twobp = list()
  zerobp = list()
  nsamples = list()
  combineddata = list()
  twobpfull = list()
  zerobpfull = list()
  cors = list()
  
  for(i in 1:length(twobpalldatalist)){
    twobp[[i]] = twobpalldatalist[[i]][,as.vector(indexmatrices[[i]])]
    zerobp[[i]] = zerobpalldata[[i]][,as.vector(indexmatrices[[i]])]
    nsamples[[i]] = dim(twobp[[i]])[2]
    combineddata[[i]] = merge(twobp[[i]],zerobp[[i]],all=TRUE,by='row.names')
    combineddata[[i]] = combineddata[[i]][,2:dim(combineddata[[i]])[2]]
    combineddata[[i]][is.na(combineddata[[i]])] = 0
    
   # plot(combineddata[[i]][[1]], combineddata[[i]][[45]])
    twobpfull[[i]] = combineddata[[i]][,1:nsamples[[i]]]
    zerobpfull[[i]] = combineddata[[i]][,(nsamples[[i]] + 1):(2*nsamples[[i]])]
    cors[[i]] = mapply(FUN = cor, x = split(t(twobpfull[[i]]), 1:ncol(twobpfull[[i]])), y = split(t(zerobpfull[[i]]), 1:ncol(twobpfull[[i]])))
  }
  png(paste("compare two bp with zero bp", folder,".png"),width=size*ppi,height=size*ppi,res=ppi)
  boxplot(cors, xlab = '', xaxt = 'n', cex.axis = 1.5)
  axis(side = 1, at = 1:length(twobp), labels = animalnames, las = 2, cex.axis= 1.25)
  mtext(paste("Pearson correlation"), side=2, line=2.5, cex = 2)
  mtext(paste("2bp mismatch/indel matching versus 0 bp"), side=3, line=1, cex = 2)
  dev.off()
  
  for(i in 1:4){
    png(paste("compare two bp with zero bp scatter", animalnames[i],".png"),width=size*ppi,height=size*ppi,res=ppi)
    par(mar = c(7,7,4,4))
   plot(combineddata[[i]][,nsamples[[i]]], combineddata[[i]][,2*nsamples[[i]]], xlab= '2 bp mismatch/indel matching read count',
        ylab= 'No matching read count', main = animalnames[i], cex.lab = 2, cex.main = 2,cex.axis = 1.5)
   dev.off()
  }
 
}
