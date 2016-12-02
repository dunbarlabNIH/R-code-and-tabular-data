#Samson Koelle
#copywrite Cynthia Dunbar Lab
#10-18-16

#this function plots both the cumulative hematopoietic contribution and cumulative clone number against clone size
#the threshold is also plotted
cumulativebarcodefrequency = function(data = zh33alldatanothresh[,as.vector(zh33indexmatrix[,4])], folder = ''){
  
  #split data into individual columns
  dataforapply = split(t(data), f = 1:ncol(data))
  
  #normalize
  dataforapply = lapply(dataforapply, function(x){x/sum(x)})
   
  datasorted = lapply(FUN = sort, X = dataforapply, decreasing  = F)
  datasortednonzero = list()
  cumtotal = lapply(FUN = cumsum, X = datasorted)
  png(paste("showthreshold", folder,".png"),width=size*ppi,height=size*ppi,res=ppi)
  #par(new = F, mar = c(7,5,3,1), xpd = TRUE)
  par(new = F, mar = c(7,6,3,7))
  plot(1,xlim = log(c(0.00003,.1)), ylim = c(0,1), cex = 0.001, yaxt = 'n', xlab= '', ylab = '')
  par(new = T)
  #for(i in 1:ncol(data)){
  for(i in 1:length(data)){
    nonzero = which(datasorted[[i]] > 0)
    datasortednonzero[[i]] = datasorted[[i]][nonzero]
    lines(log(datasortednonzero[[i]]), seq(from = 0, to = 1 - 1/length(datasortednonzero[[i]]), by = 1/length(datasortednonzero[[i]])))
    lines(log(datasortednonzero[[i]]), cumtotal[[i]][nonzero], xlim = log(c(0.0001,.1)), lwd = .5, col = 'grey')
    par(new = T)
  }
  abline(v = log(0.0005), lwd = 3, col = 1)
  axis(side = 4, at = c(0,.2,.4,.6,.8,1), labels = c('0','0.2','0.4','0.6','0.8','1'), col = 'grey', col.axis = 'grey', col.ticks = 'grey')
  axis(side = 2, at = c(0,.2,.4,.6,.8,1), labels = c('0','0.2','0.4','0.6','0.8','1'), col = 1)
  mtext(side = 4, line = 5, cex = 2, text = 'Fraction of hematopoiesis from clones\nsmaller than size', col = 'grey')
  mtext(side = 2, line = 3, cex = 2, text = 'Fraction of clones smaller than size')
  mtext(side = 1, line = 3, cex = 2, text = 'Log fractional clone size')
  mtext(side = 3, line = 1, cex = 2, text = folder)
  text(x = log(2*0.0005), y = .1, labels = "Threshold")
  dev.off()
}
