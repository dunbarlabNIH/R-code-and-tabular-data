#Samson Koelle
#Copywrite Cynthia Dunbar Lab
#11-16-16

#diversity track tracks the shannon diversity of cell type clonalities over time in one animal
diversitytrack = function(data = zh33alldata, indexmatrix = zh33indexmatrix, timepoints = zh33timepoints, folder = "ZH33"){
  
    diversity = indexmatrix
    sums = apply(FUN = sum, MARGIN = 2, X = data)
    data[data == 0] = 1
    for(i in 1:length(data)){
      data[[i]]=(data[[i]]/sum(data[[i]]))
    }
    for(i in 1:dim(indexmatrix)[1]){
      #print(i)
      for(j in 1:dim(indexmatrix)[2]){
        #print(j)
        #tempdata = zh19alldata[zh19indexmatrix[i,j]]
        if(!is.na(indexmatrix[i,j])){
          if(sums[indexmatrix[i,j]] > 0){
           diversity[i,j] = Hs(data[,indexmatrix[i,j]])
          } else {
            diversity[i,j] = NA
        }
        }
      }
    }
  
    png(paste("diversity", folder,".png"),width=1.4*size*ppi,height=size*ppi,res=ppi)
    par(new = F, mar = c(7,6,3,7), xpd = TRUE)
    for(i in 1:dim(indexmatrix)[2]){
      print(i)
      plot(timepoints, diversity[,i], cex= 2, pch = 15, col = i, xlim = c(0,49),
           ylim = c(3.5,8), ylab = "", xlab = "", cex.lab = 3,
           cex.axis = 2)
      lines(timepoints, diversity[,i], col = i, lwd = 2)
      par(new = T)
    }
    mtext("Months Post-Transplant", side=1, line=3, cex = 2)
    mtext("Shannon Diversity", side = 2, line = 3, cex = 2)
    mtext(folder, side = 3, line = 1, cex = 3)
    dev.off()
    
    #Generate legend (actual plot for convenience only)
    png(paste("diversity legend", folder,".png"),width=1.4*size*ppi,height=size*ppi,res=ppi)
    par(new = F, mar = c(7,6,3,7), xpd = TRUE)
    plot(timepoints[1], diversity[1,1], cex= 2, pch = 15, col = i, xlim = c(0,49),
         ylim = c(3.5,8), ylab = "", xlab = "", cex.lab = 3,
         cex.axis = 2)
    legend(35,8, c("T",'B','Mono',"Gr"), title = "Cell Type",cex = 1.6, 
           col = c(1:4), pch = 15:18)
    dev.off()
}


                              