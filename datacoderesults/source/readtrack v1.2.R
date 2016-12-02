#Samson Koelle
#copywrite Cynthia Dunbar Lab
#10-18-16

#This function tracks the proportion of sequencing reads retained by the various processing steps
readtrack = function(indexmatrix = zh33indexmatrix, folder = 'ZH33', timepoints = zh33timepoints, rawreads = zh33readmeordered$READS, readswithLibID = zh33readmeordered$MAPPED, readswithLibIDover100 = apply(FUN = sum, MARGIN = 2, X = zh33alldatanothresh),readswithLibIDover100overthresh = apply(FUN = sum, MARGIN = 2, X = zh33alldata),celltypenames = c('T',"B",'Gr',"Mono")){
  
  #compare overall effect of various preprocessing steps on number of reads observed
  png(paste("read count boxplot", folder,".png"),width=1.4*size*ppi,height=size*ppi,res=ppi)
  par(new = F, mar = c(10,6,3,7), xpd = TRUE)
  boxplot(rawreads[as.vector(indexmatrix)], 
          readswithLibID[as.vector(indexmatrix)],
          readswithLibIDover100[as.vector(indexmatrix)],
          readswithLibIDover100overthresh[as.vector(indexmatrix)], xaxt = 'n', cex.axis = 2.3)
  axis(side = 1, at = 1:4, labels = c("Raw", "With LibID", "With LibID \n over \n 100 reads", "Post- \n 0.05% \n threshold"), las = 2, cex.axis= 2)
  mtext("Sequencing read count", side = 2, line = 3, cex = 2.5)
  mtext(folder, side = 3, line = 1, cex = 3)
  dev.off()
  
  #analyze one cell type at a time
  for(i in 1:dim(indexmatrix)[2]){
    png(paste('Pipeline analysis', folder, celltypenames[i]), width=1.4*size*ppi,height=size*ppi,res=ppi)
    par(new = F, mar = c(7,6,3,7), xpd = TRUE)
    
    #plot number of raw reads
    plot(rawreads[indexmatrix[,i]], cex= .2, pch = 15, col = 1, xlim = c(0,length(timepoints)),
         ylim = c(0, 4000000), ylab = "", xlab = "", cex.lab = 3,
         cex.axis = 2, xaxt = 'n', yaxt= 'n')
    lines(rawreads[indexmatrix[,i]], col = 1, lwd = 3, lty = 1)
    par(new = T)
    
    #plot number of reads with LibID
    plot(readswithLibID[indexmatrix[,i]], cex= .2, pch = 16, col = 1, xlim = c(0,length(timepoints)),
         ylim = c(0, 4000000), ylab = "", xlab = "", cex.lab = 3,
         cex.axis = 2, xaxt = 'n', yaxt= 'n')
    lines(readswithLibID[indexmatrix[,i]], col = 1, lwd = 3, lty = 2)
    par(new = T)
    
    #plot number of reads with LibID and raw abundance over 100
    plot(readswithLibIDover100[indexmatrix[,i]], cex= .2, pch = 17, col = 1, xlim = c(0,length(timepoints)),
         ylim = c(0, 4000000), ylab = "", xlab = "", cex.lab = 3,
         cex.axis = 2, xaxt = 'n', yaxt= 'n')
    lines(readswithLibIDover100[indexmatrix[,i]], col = 1, lwd = 3, lty = 3)
    par(new = T)
    
    #plot number of reads with LibID and raw abundance over 100 that exceed threshold
    plot(readswithLibIDover100overthresh[indexmatrix[,i]], cex= .2, pch = 18, col = 1, xlim = c(0,length(timepoints)),
         ylim = c(0, 4000000), ylab = "", xlab = "", cex.lab = 3,
         cex.axis = 2, xaxt = 'n', yaxt= 'n')
    lines(readswithLibIDover100overthresh[indexmatrix[,i]], col = 1, lwd = 3, lty = 4)
    axis(side = 1, at = 1:length(timepoints), labels = timepoints, las = 2, cex.axis= 2.25)
    axis(side = 2, at = c(0, 2000000,4000000), labels = c(0, '2E6', '4E6'), cex.axis = 2.25)
    mtext("Months Post-Transplant", side=1, line=5, cex = 2.5)
    mtext("Sequencing read count", side = 2, line = 3, cex = 2.5)
    mtext(paste(folder, celltypenames[i]), side = 3, line = 1, cex = 3)
    
    dev.off()
    
  }
  
  #create legend
  png(paste("read count legend", folder,".png"),width=1.4*size*ppi,height=size*ppi,res=ppi)
  par(new = F, mar = c(7,6,3,7), xpd = TRUE)
  plot(0, 1, xlim = c(0,1), ylim = c(0,1), col = 'white')
  legend(x = 0, y = 1, legend = c("Total reads", "Reads with LibID", "Reads over 100", "Post-Threshold"), lty = c(1:4))
  
  dev.off()
  
}

  
  