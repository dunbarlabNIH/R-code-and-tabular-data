#Samson Koelle
#copywrite Cynthia Dunbar Lab
#10-18-16

#this function plots the overall detected clone number in each animal as a function of threshold
cnattpsensitivitybyanimal = function(data = list(zh33alldatanothresh,zg66alldatanothresh,zh19alldatanothresh,zj31alldatanothresh),indexmatricesnoNA = list(zh33indexmatrixnoNA, zg66indexmatrixnoNA,zh19indexmatrixnoNA, zj31indexmatrixnoNA),folder = "Clone number at timepoint sensitivity", threshseq = seq(from = 0, to = 10000, by = 1000)){
  
  
  cnattp = list()
  foreach(k=1:length(data)) %do% {
    
    #initialize diversity matrices
    cnattp[[k]] = array(dim = length(threshseq))
    for(h in 1:length(threshseq)){
      print(h)
          cnattp[[k]][h] = dim(thresholdinsamples(data = data[[k]], samples = as.vector(indexmatricesnoNA[[k]]), thresh = threshseq[h]))[1]
      }
  }
  
  png(paste("cnattpsensitivitybyanimal", folder,".png"),width=1.4*size*ppi,height=size*ppi,res=ppi)
  par(new = F, mar = c(7,6,3,7), xpd = TRUE)
  
  for(k in 1:length(data)){
        plot(threshseq, cnattp[[k]],
             cex= 2, pch = 14 + k, xlim = c(min(threshseq),max(threshseq)),
             ylim = c(0,max(unlist(cnattp))), ylab = "", xlab = "", cex.lab = 3,
             cex.axis = 2)
        lines(threshseq, cnattp[[k]], cex= 1,  xlim = c(min(threshseq),max(threshseq)),
              ylim = c(0,max(unlist(cnattp))), ylab = "", xlab = "", cex.lab = 3,
              cex.axis = 2, lwd = 2)
        par(new = T)
  } 
  mtext("Clone Number Sensitivity to Threshold", side=3, line=1, cex = 2)
  mtext("Threshold", side=1, line=3, cex = 2)
  mtext("Clone Number Detected", side = 2, line = 3, cex = 2)
  #mtext(folder, side = 3, line = 1, cex = 3)
  #legend(45,8, c("T",'B','Mono',"Gr"), title = "Cell Type",cex = 1.6, 
  #col = c(1:4), pch = 15)
  dev.off()
  
  png(paste("cnattp sensitivity legend", folder,".png"),width=1.4*size*ppi,height=size*ppi,res=ppi)
  par(new = F, mar = c(7,6,3,7), xpd = TRUE)
  plot(1,1, cex= .02, pch = 15, col = 1, xlim = c(0,10),
       ylim = c(-5,20), ylab = "", xlab = "", cex.lab = 3,
       cex.axis = 2)
  legend(5,8, c("ZH33",'ZG66','ZH19','ZJ31'), title = "Animal",cex = 1.6, 
         col = c(1), pch = 15:18)
  dev.off()
  
}