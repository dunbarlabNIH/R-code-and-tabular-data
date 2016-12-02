#Samson Koelle
#copywrite Cynthia Dunbar Lab
#10-18-16

#diversity plot plots all clones detected above the threshold at a given time point in a certain cell type
diversityplot = function(data = zh33data, thresh = 2000, folder = "ZH33", names = zh33timepoints){

  data <- data[apply(X = data, FUN = sum, MARGIN  = 1) > 0, ]
  data <- thresholdinsamples(data = data, samples = 1:length(data), thresh = thresh)
  for(i in 1:length(data)){
    if(sum(data[[i]]) >0){
    data[,i]=(data[[i]]/sum(data[[i]]))
    }
  }
  data <- data[apply(X = data, FUN = sum, MARGIN  = 1) > 0, ]
  
  counter = c()
  for(i in 1:length(data)){
    ordered_data = data[order(-data[[i]]),i]
    counter = append(counter, length(ordered_data[ordered_data >0]))
  }

  png(paste("diversity curve", folder,".png"),width=size*ppi,height=size*ppi,res=ppi)
  max10 = 0
  par(new = F, mar = c(7,5,3,1), xpd = FALSE)
  for(i in 1:length(data)){
    ordered_data = data[order(-data[[i]]),i]
    vector_tobeplotted = c()
    total_reads = 0
    for(j in 1:length(ordered_data)){
      total_reads = total_reads + ordered_data[j]
      vector_tobeplotted = append(vector_tobeplotted, total_reads)
    }
    if(sum(vector_tobeplotted >0)){
    plot(vector_tobeplotted, 
         col = (rainbow(256, s = 1, v = 1, start = 0, end = 0.75, alpha = 1))[i*256/length(data)],
         pch = 16, ylim = c(0,1), xlim = c(0,max(counter)), xlab = '', ylab = "",
         cex.axis = 1.25, cex.lab = 1.5, xaxt = 'n', yaxt ='n')
    max10 = max(vector_tobeplotted[10], max10)
    #plot(log(ordered_data), col =i, pch = 14, ylim = c(-10,-2), xlim = c(0,max(counter)), xlab = ("Clone rank"), ylab = ("Cumulative fractional read contribution"))
    abline(h = vector_tobeplotted[100], col = (rainbow(256, s = 1, v = 1, start = 0, end = 0.75, alpha = 1))[i*256/length(data)], lwd = 2)
    par(new=T)		
    }
  }
  
  axis(side = 1, at = seq(from = 1, to = max(counter), by = 200), labels = seq(from = 1, to = max(counter), by = 200), las = 1, cex.axis = 1.75)
  axis(side = 2, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 1, by = 0.2), las = 1, cex.axis = 1.75)
  
  mtext("Clone size rank at time point \n (large to small)", side=1, line=4, cex = 1.5)
  mtext("Cumulative Clone Size", side=2, line=3.25, cex = 2)
  mtext(paste(folder, "clone sizes"), side=3, line=1, cex = 2)
  #text(labels = "Total sizes of 100 \n largest clones in each sample", x = (400), y = max10 + 0.05, cex = 2)
  dev.off()
  
  png(paste("diversity curve legend", folder,".png"),width=size*ppi,height=size*ppi,res=ppi)
  plot.new()
  legend(0,1, names, title = "Months \n post-transplant",bty = 'n', cex = 1.4, col = (rainbow(256, s = 1, v = 1, start = 0, end = 0.75, alpha = 1))[seq(from = (256/length(data)), to = 256, by = 256/length(data))], pch = 16)
  dev.off()
  
}
