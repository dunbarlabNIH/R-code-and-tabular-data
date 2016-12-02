#Samson Koelle
#Copywrite Cynthia Dunbar Lab
#11-16-16

#biastracking tracks the bias of individual clones between two cell types over time
biastracking = function(data = zh33alldata, uppersamples = zh33indexmatrix[,2],lowersamples =  zh33indexmatrix[,4], timepointnames = zh33timepointnames,folder = "ZH33 B Gr",upperbiaslabel = "B", lowerbiaslabel = 'Gr'){
  ntimepoints = length(uppersamples)
  
  #Get rid of NAs
  notNA = c()
  for(i in 1:ntimepoints){
    if(!is.na(uppersamples[i] + lowersamples[i])){
      notNA = append(notNA, i)
    }
  }
  
  #upper samples and lower samples refer to the top and bottom of the plot respectively
  uppersamples = uppersamples[notNA]
  lowersamples = lowersamples[notNA]
  timepointnames = timepointnames[notNA]
  #normalize data
  for(i in 1:length(data)){
    data[[i]]= data[[i]]/sum(data[[i]])
  }
  #get barcodes that aren't all zero
  data = data[(apply(X = data[,union(uppersamples,lowersamples)], FUN = sum, MARGIN = 1) >0),]
  nclones = length(data[,1])
  ratios = array(dim = c(nclones, ntimepoints))
  sizes = array(dim = c(nclones, ntimepoints))
  for(j in 1:ntimepoints){
    print(j)
    ratios[,j] = data[,uppersamples[j]] / (data[,uppersamples[j]] + data[,lowersamples[j]])
    sizes[,j] = (data[,uppersamples[j]] + data[,lowersamples[j]])
  }
  ratios[ratios == "NaN"] = 0
  ratios = data.frame(ratios)
  colnames(ratios) = timepointnames
  sizes = data.frame(sizes)
  colnames(sizes) = timepointnames
  timepointlist = c()
  colourdf = array(0, dim = c(nclones, ntimepoints - 1))
  layerdf = array(0, dim = c(nclones, ntimepoints - 1))
  for(i in 1:(ntimepoints-1)){
    colourdf[,i] = floor(100 * (sizes[,i] + sizes[,i+1]) / (2*max(sizes)))
    colourdf[,i]
    layerdf[,i] = order(sizes[,i+1])
  }
  
  #colour corresponds to the total size of the clone in all plotted samples
  #order of plotting is thee same as order of colour
  #i.e. the darkest lines go last
  #these are the largest clones in the two cell types over the time course shown
  colours = floor(100* apply(X = sizes, FUN = sum, MARGIN = 1) / max(apply(X = sizes, FUN = sum, MARGIN = 1)))
  ratios = ratios[order(colours),]
  colours = colours[order(colours)]
  png(paste("biaslines", folder,".png"),width=size*ppi,height=size*ppi,res=ppi)
  par(new = F, mar = c(8,7,4,2))
  plot(c(1,ntimepoints), c(0,1), xlim = c(1,ntimepoints), ylim= c(0,1), col = 'white',
       xaxt = 'n', yaxt = 'n', xlab= '', ylab ='')
  for(i in 1:nclones){
    par(new = T)
    
    lines(1:ntimepoints, ratios[i,] , ylim = c(0,1), xlim = c(1,ntimepoints), 
    col = colorRampPalette((brewer.pal(n = 7, name = "Greys")))(100)[colours[i]], lwd = 2,
          xaxt = 'n', yaxt = 'n', xlab= '', ylab ='')
  }
  axis(side = 1, at = 1:ntimepoints, labels = timepointnames, las = 2, cex.axis = 1.5)
  axis(side = 2, at = c(0,.33,.5,.67,1), labels = c(paste(lowerbiaslabel, "1:0", upperbiaslabel), "2:1", 
                                                    "1:1", "1:2", paste(lowerbiaslabel,"0:1", upperbiaslabel)), las = 1, cex.axis = 1.5)
  
  mtext("Months post-transplant", side=1, line=4, cex = 2)
  mtext(paste("Bias"), side=2, line=4, cex = 2)
  mtext(paste(folder, "bias"), side=3, line=1, cex = 2.5)
  dev.off()
  
  
  
  
  
  
}