#Samson Koelle
#Copywrite Cynthia Dunbar Lab
#11-16-16

#Clone number track tracks the number of clones appearing in each cell type over time
#Note that this is cumulative, i.e. it plots the # of clones which have appeared in that cell type, either 
#at the time point in question, or previously 
clonenumbertrack = function(data = zh33alldata, indexmatrixnoNA = zh33indexmatrixnoNA, thresh = 2000, timepoints = zh33timepoints, folder = "ZH33"){
  
    ntimepoints = dim(indexmatrixnoNA)[1]
    ncelltypes = dim(indexmatrixnoNA)[2]
    nclones = length(data[[1]])
    
    threshes = array(dim = c(ntimepoints, ncelltypes))
    
    #for timepoint specific
    clonedetectedincelltypetimepoint = array(0, dim = c(nclones,ntimepoints, ncelltypes))

    for(i in 1:ntimepoints){
      for(j in 1:ncelltypes){
        threshes[i,j] = thresh * sum(data[,indexmatrixnoNA[i,j]]) / 4000000
        #clone is detected if it is above the threshold in a given sample
        clonedetectedincelltypetimepoint[(data[,indexmatrixnoNA[i,j]] > threshes[i,j]),i,j] = 1
      }
    }
    
    clonedetectedincelltype = array(0, dim = c(nclones,ntimepoints, ncelltypes))
    clonedetected = array(0, dim = c(nclones, ntimepoints))
    for(i in 1:nclones){
      #print(i)
      for(j in 1:ntimepoints){
        #cumulative clone numbers in celltype
        if(any(clonedetectedincelltypetimepoint[i,1:j,] == 1)){
          clonedetected[i,j] = 1
        }
        #cumulative clone numbers total
        for(h in 1:ncelltypes){
          if(any(clonedetectedincelltypetimepoint[i,1:j,h]) == 1){
            clonedetectedincelltype[i,j,h] = 1
          }
        }
      }
    }
    
    nclonesincelltype = array(0, dim = c(ntimepoints, ncelltypes))
    nclonesincelltypeattimepoint = array(0, dim = c(ntimepoints, ncelltypes))
    for(h in 1:ncelltypes){
      nclonesincelltype[,h] = apply(FUN = sum, MARGIN = 2, X = clonedetectedincelltype[,,h])
      nclonesincelltypeattimepoint[,h] = apply(FUN = sum, MARGIN = 2, X = clonedetectedincelltypetimepoint[,,h])
    }
    
    clonenumberbycelltype = nclonesincelltype
    clonenumberoverall  = apply(FUN = sum, MARGIN = 2, X = clonedetected)
    clonenumberbycelltypeattimepoint = nclonesincelltypeattimepoint
  
  png(paste("clonenumber", folder,".png"),width=1.4*size*ppi,height=size*ppi,res=ppi)
  par(new = F, mar = c(7,6,3,8), xpd = TRUE)
    for(h in 1:ncelltypes){
    plot(timepoints, clonenumberbycelltype[,h], pch = 15, xlim = c(0,49), ylim = c(0,6500), xlab= "", 
         ylab= '',  col = h, cex = 2, cex.axis = 2)
    lines(timepoints, clonenumberbycelltype[,h], pch = 15, xlim = c(0,49), ylim = c(0,6500), xaxt = "n", xlab= "", 
          ylab= '', yaxt = 'n', col = h, lwd = 2)
    par(new = T)
    
    }
    plot(timepoints, 
         clonenumberoverall, xlim = c(0,49), ylim = c(0,6500),
         pch = 15, col = 'gray', xaxt = "n", xlab= "", ylab= '',  yaxt = 'n', cex = 2, cex.axis = 2)
    lines(timepoints, clonenumberoverall,
          pch = 15, xlim = c(0,49), ylim = c(0,6500), xaxt = "n", xlab= "", 
          ylab= '', yaxt = 'n', col = "black", lty = 2, lwd = 2)
    par(new = T)
  
  mtext("Months Post-Transplant", side=1, line=3, cex = 2)
  mtext("Number of clones detected", side = 2, line = 3, cex = 2)
  mtext(folder, side = 3, line = 1, cex = 3)
  #mtext(folder, side = 3, line = 1, cex = 3)
  dev.off()
  
  png(paste("clonenumber legend", folder,".png"),width=1.4*size*ppi,height=size*ppi,res=ppi)
  par(new = F, mar = c(7,6,3,8), xpd = TRUE)
  plot(timepoints, 
       clonenumberoverall, xlim = c(0,49), ylim = c(0,6500),
       pch = 15, col = 'gray', xaxt = "n", xlab= "", ylab= '',  yaxt = 'n', cex = 2, cex.axis = 2)
  #lines(timepoints, clonenumberoverall,
   #     pch = 15, xlim = c(0,43), ylim = c(0,6500), xaxt = "n", xlab= "", 
   #     ylab= '', yaxt = 'n', col = "black", lty = 2, lwd = 2)
  legend(52,6500, c("T",'B','Mono',"Gr","Total"), title = "Cell Type",cex = 1.6, 
         col = append(1:4, "gray"), pch = c(15:18,15))
  dev.off()
  
}
