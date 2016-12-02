#Samson Koelle
#copywrite Cynthia Dunbar Lab
#10-18-16

#Diversity curve plots the overall size (i.e. sum in all samples) of each clone in all animals.
#Clones are sorted by overall size, both above the threshold and below the threshold
#Clone position is independent for each animal
#Thus contrast to Figures 4 and 5 and also Figure S2, two clones being on the same vertical line does not indicate that they are in fact the same clone
diversitycurve = function(datalist = list(zh33alldata,zg66alldata,zh19alldata,zj31alldata),datalistnothresh = list(zh33alldatanothresh, zg66alldatanothresh,zh19alldatanothresh, zj31alldatanothresh),indexmatrices = list(zh33indexmatrix, zg66indexmatrix,zh19indexmatrix, zj31indexmatrix),folder = "DIVCURVE"){
  
  abovecumsums = list()
  belowcumsums = list()
  for(k in 1:length(datalist)){
    nsamples = as.numeric(table(is.na(indexmatrices[[k]]))[1])
    data = datalistnothresh[[k]]
    indexmatrix = indexmatrices[[k]]
    #sums = apply(FUN = sum, MARGIN = 2, X = data)
    #data[data == 0] = 1
    for(i in 1:length(data)){
      data[[i]]=(data[[i]]/sum(data[[i]]))
    }
    datasum = apply(FUN = sum, MARGIN = 1, X = data[,as.vector(indexmatrix)[!is.na(as.vector(indexmatrix))]])
    datasum = datasum / nsamples
    abovethresh = which(rownames(data) %in% rownames(datalist[[k]]))
    belowthresh = which(!(rownames(data) %in% rownames(datalist[[k]])))
    above = datasum[abovethresh]
    below = datasum[belowthresh]
    abovesorted = sort(datasum[abovethresh], decreasing = TRUE)
    belowsorted = sort(datasum[belowthresh], decreasing = TRUE)
    abovecumsum = cumsum(abovesorted)
    belowcumsum = cumsum(belowsorted)
    belowcumsum = belowcumsum + max(abovecumsum)
    abovecumsums[[k]] = abovecumsum
    belowcumsums[[k]] = belowcumsum
    }
  
    
  png(paste("diversitycurve", folder,".png"),width=1.4*size*ppi,height=size*ppi,res=ppi)
  par(new = F, mar = c(7,6,3,7))
  for(k in 1:length(datalist)){
    
    plot((1:length(abovecumsums[[k]])), abovecumsums[[k]], pch = 16, xlim = c(0,40000), ylim = c(0,1), xaxt = "n", xlab= "", 
         ylab= '', yaxt = 'n')
    par(new = T)
    plot(((1+length(abovecumsums[[k]])):(length(abovecumsums[[k]]) + length(belowcumsums[[k]]))), 
         belowcumsums[[k]], xlim = c(0,40000), ylim = c(0,1), pch = 16, col = 'gray', xaxt = "n", xlab= "", ylab= '', yaxt = 'n')
    par(new = T)
    plot(length(abovecumsums[[k]]), abovecumsums[[k]][length(abovecumsums[[k]])], pch = 14+k,
         col = 1, xaxt = "n", xlab= "", ylab= '', yaxt = 'n',ylim = c(0,1),xlim = c(0,40000),cex = 5)
    par(new = T)
    abline(v = length(abovecumsums[[k]]), col = 'grey')
    par(new = T)
  }
  mtext("Clone Size Rank", side=1, line=3, cex = 2)
  mtext("Cumulative Clone Size", side = 2, line = 3, cex = 2)
  mtext("Clone Sizes", side = 3, line = 1, cex = 3)
  axis(1, at = c(0,10000,20000,30000, 40000), labels = c("0", "10000", "20000", "30000", '40000'), las = 1, cex.axis = 1.5)
  axis(2, at=c(0,0.2,0.4,0.6,0.8,1), labels = c(0,0.2,0.4,0.6,0.8,1), las = 3, cex.axis = 2)
  legend(22000,0.4, c("ZH33",'ZG66','ZH19','ZJ31'), title = "Animal",cex = 1.6, 
         col = c(1), pch = 15:18)
  legend(7000,0.4, c("Threshold"), title = "",cex = 1.6, 
         col = 'grey', pch = '|', text.col = 'grey', bty = "n")
  dev.off()
  
    

}


