#Samson Koelle
#copywrite Cynthia Dunbar Lab
#10-18-16

#This function plots the contribution of each clone in all four animals
contamcheck = function(data = list(zh33alldata,zg66alldata,zh19alldata,zj31alldata), folder = "post-threshold",indexmatricesnoNA = list(zh33indexmatrixnoNA, zg66indexmatrixnoNA,zh19indexmatrixnoNA, zj31indexmatrixnoNA)){
  
  animalnames = c('ZH33','ZG66','ZH19','ZJ31')
  #sum and normalize clone size in all samples
  sumdata = list()
  for(k in 1:length(data)){
    datatemp = data[[k]][,as.vector(indexmatricesnoNA[[k]])]
    sumdata[[k]] = apply(FUN = sum, MARGIN = 1, X = datatemp)
    sumdata[[k]] = sumdata[[k]][sumdata[[k]] > 0]
    sumdata[[k]] = sumdata[[k]] / sum(sumdata[[k]])
  }
  
  #find the set of all barcodes in all animals
  barcodesunified = c()
  for(k in 1:length(data)){
    barcodesunified = union(barcodesunified, names(sumdata[[k]]))
  }
  
  #extend data from each animal to include barcodes from other animals
  extendeddata = list()
  for(k in 1:length(data)){
    barcodesnotfound = setdiff(barcodesunified, names(sumdata[[k]]))
    extendeddata[[k]] = append(sumdata[[k]],rep(0, length(barcodesnotfound)))
    names(extendeddata[[k]]) = append(names(sumdata[[k]]), barcodesnotfound)
    extendeddata[[k]] = extendeddata[[k]][order(names(extendeddata[[k]]))]
  }
  
  #join data
  dataunified = do.call(what = cbind, args = extendeddata)

  #find which animal is each barcode principally found in and order
  maxanimal = apply(FUN = which.max, X = dataunified, MARGIN = 1)
  orderforplotting = c()
  runningindex = c(1)
  for(i in 1:dim(dataunified)[2]){
    whichbarcodeismaxini = which(maxanimal == i)
    orderwithinmaxanimal = whichbarcodeismaxini[order(-dataunified[whichbarcodeismaxini,i])]
    orderforplotting = append(orderforplotting, orderwithinmaxanimal)
    runningindex = append(runningindex, length(whichbarcodeismaxini) + runningindex[i])
  }
  
  #order data
  ordereddata = dataunified[orderforplotting,]
  
  #convert to cumulative distribution
  ordereddatacumulative = apply(FUN = cumsum, MARGIN = 2, X = ordereddata)
  
  #compute sparse indices (for plotting clarity)
  sparse_indices = seq(from = 1, to = dim(dataunified)[1], by = 1)
  
  for(j in 1:length(data)){
  #plot
    print(j)
  png(paste("contamination check", toString(j),folder,".png"),width=size*ppi,height=size*ppi,res=ppi)
  par(new = F, mar = c(7,5,3,1), xpd = TRUE)
  plot(sparse_indices, ordereddatacumulative[sparse_indices,j], col = 1,
       pch = 16, ylim = c(0,1), xlim = c(0,dim(dataunified)[1]), ylab= "", 
       xlab = "", cex.lab = 1.5, cex.axis = 1.7)
  for(i in 1:dim(ordereddata)[2]){
    
    par(new = T)
#     polygon(border = alpha(colour = 1, alpha = 0.1), lty = 1, lwd = 0.1,
#     col = alpha(colour = (rainbow(256, s = 1, v = 1, start = 0, end = 0.75, alpha = 1))[seq(from = (256/dim(dataunified)[2]), to = 256, by = 256/dim(dataunified)[2])[i]], alpha = 0.5), 
#             x = c(runningindex[i], runningindex[i],
#                   runningindex[i+1], runningindex[i+1]), y = c(-0.1,.1,.1,-0.1))
    lines(x = c(rep(runningindex[i],2)), y = c(0,1))
    lines(x = c(rep(runningindex[i+1],2)), y = c(0,1))
    par(new = T)
  }
  
  #add info
  mtext(paste('Barcodes principally found in', animalnames[j]), side=3, line=1, cex = 2)
  mtext("Cumulative Barcode Size (all animals)", side=2, line=3, cex = 2)
  mtext(paste("Barcodes sorted by principal animal"), side=1, line=3, cex = 1.75)
  text(x = runningindex[j] + 1000, y = 0.1, labels = toString(round(sum(ordereddata[runningindex[j]:(runningindex[j+1]-1),j]), digits = 2)), cex =2)
  dev.off()
  }
  
  #make legends
  png(paste("Contamination check legend 1", folder,".png"),width=size*ppi,height=size*ppi,res=ppi)
  plot.new()
  legend(0,1, c('ZH33','ZG66','ZH19','ZJ31'), title = "Animal",cex = 1.4, pch = 15:18)
  dev.off()
  
  png(paste("Contamination check legend 2", folder,".png"),width=size*ppi,height=size*ppi,res=ppi)
  plot.new()
  legend(0,1, c('ZH33','ZG66','ZH19','ZJ31'), title = "Principal animal",cex = 1, pch = 16,
         col = alpha(colour = (rainbow(256, s = 1, v = 1, start = 0, end = 0.75, alpha = 1))[seq(from = (256/dim(dataunified)[2]), to = 256, by = 256/dim(dataunified)[2])]))
  dev.off()
  
}
