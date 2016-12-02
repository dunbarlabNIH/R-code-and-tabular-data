#Samson Koelle
#Copywrite Cynthia Dunbar Lab
#11-16-16

#GFPtrack tracks the GFP % positivity of each cell type over time
GFPtrack = function(data = zh33GFP, timepoints = zh33timepoints, folder = "ZH33"){
  
  png(paste("GFP", folder,".png"),width=1.4*size*ppi,height=size*ppi,res=ppi)
  par(new = F, mar = c(7,6,3,7), xpd = TRUE)
    for(i in 1:4){
      print(i)
      plot(timepoints, data[,i], cex= 2, pch = 15, col = i, xlim = c(0,49),
           ylim = c(0,55), ylab = "", xlab = "", cex.lab = 3,
           cex.axis = 2)
      lines(timepoints, data[,i], col = i, lwd = 3)
      par(new = T)
    }
  mtext("Months Post-Transplant", side=1, line=3, cex = 2)
  mtext("GFP+ %", side = 2, line = 3, cex = 2)
  mtext(folder, side = 3, line = 1, cex = 3)
  dev.off()
}


                              