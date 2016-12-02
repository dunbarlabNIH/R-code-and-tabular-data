#Samson Koelle
#Copywrite Cynthia Dunbar Lab
#11-16-16

#rolling variance of bias tracks the variance of the bias in all previous time points
#that is, the variance measurement at 12 months would be the variance of the biases between the two cell types in question
#over the course of all time points up to and including 12 months
rollingvarianceofbias <- function(data = zh33alldata, uppersamples = zh33indexmatrixnoNA[,3],lowersamples = zh33indexmatrixnoNA[,4], folder = "ZH33 Mono v Gr", timepoints = zh33timepoints) 
{
  
    ntimepoints = length(uppersamples)
  #NAs are not okay... neither are zero columns!
  notNA = c()
  for(i in 1:ntimepoints){
    if(!is.na(uppersamples[i] + lowersamples[i])){
      notNA = append(notNA, i)
    }
  }
  uppersamples = uppersamples[notNA]
  lowersamples = lowersamples[notNA]
  timepoints = timepoints[notNA]
  #normalize data
  for(i in 1:length(data)){
    data[[i]]= data[[i]]/sum(data[[i]])
  }
  #get barcodes that aren't all zero
  data = data[(apply(X = data[,union(uppersamples,lowersamples)], FUN = sum, MARGIN = 1) >0),]
  nclones = length(data[,1])
  tans = array(dim = c(nclones, ntimepoints))
  sizes = array(dim = c(nclones, ntimepoints))
  for(j in 1:ntimepoints){
    tans[,j] = data[,uppersamples[j]] / (data[,uppersamples[j]] + data[,lowersamples[j]])
    sizes[,j] = (data[,uppersamples[j]] + data[,lowersamples[j]])
  }
  tans[tans == "NaN"] = 0
  tans = data.frame(tans)
  colnames(tans) = timepoints
  sizes = data.frame(sizes)
  colnames(sizes) = timepoints
  timepointlist = c()
  rollingvar = array(0, dim = c(nclones, ntimepoints))
  for(i in 2:ntimepoints){
    rollingvar[,i] = apply(FUN = var, X = tans[,1:i], MARGIN = 1)
  }
  
  colours = floor(100* apply(X = sizes, FUN = sum, MARGIN = 1) / max(apply(X = sizes, FUN = sum, MARGIN = 1)))
  rollingvar = rollingvar[order(colours),]
  colours = colours[order(colours)]
  png(paste("Bias Variance", folder,".png"),width=size*ppi,height=size*ppi,res=ppi)
  par(new = F, mar = c(8,6,4,2))
  plot(c(1,ntimepoints), c(0,1), xlim = c(1,ntimepoints), ylim= c(0,.5), col = 'white',
       xaxt = 'n', yaxt = 'n', xlab= '', ylab ='')
  for(i in 1:nclones){
    par(new = T)
    
    lines(rollingvar[i,] , ylim = c(0,.5), xlim = c(1,ntimepoints), 
          col = colorRampPalette((brewer.pal(n = 7, name = "Greys")))(100)[colours[i]], lwd = 2,
          xaxt = 'n', yaxt = 'n', xlab= '', ylab ='')
          #col = colorRampPalette(colors = c("gray", "red"))(100)[colours], lwd = 2)
  }
  axis(side = 1, at = 1:ntimepoints, labels = timepoints, las = 2, cex.axis = 1.5)
  axis(side = 2, at = c(0,.25,.5), labels = c(0,.25,.5), las = 1, cex.axis = 1.5)
  
  mtext("Months post-transplant", side=1, line=4, cex = 2)
  mtext(paste("Variance of bias"), side=2, line=4, cex = 2)
  mtext(paste(folder, "bias stability"), side=3, line=1, cex = 2.5)
  dev.off()
  
}