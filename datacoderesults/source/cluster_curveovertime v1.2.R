#Samson Koelle
#Copywrite Cynthia Dunbar Lab
#11-16-16

#cluster_curveovertime plots the distribution of clone sizes for all time points of a specific cell type
cluster_curveovertime = function(data = zh33alldata, names = zh33names, folder = "", thresh = 2000){

#because the number of clones is implicit in the plot, we reapply the 0.0005 threshold
threshes = array(dim = length(data))
for(i in 1:length(data)){
	threshes[i] = thresh * sum(data[[i]]) / 4000000
}
clone_detected = array(dim = c(length(data[[1]]), length(data)))
clone_detected[is.na(clone_detected)] = 0
#clone has to be above threshold in THAT sample
for(i in 1:length(data)){
	clone_detected[(data[[i]] > threshes[i]),i] = 1
}
overthreshold = c()
for(i in 1:length(data[[1]])){
	if(any(clone_detected[i,]) == 1){
		overthreshold = append(overthreshold, i)
	}
}
#get the clones that are over the threshold in the samples we care about
data = data[overthreshold,]
#normalize
for(i in 1:length(data)){
	data[[i]]=(data[[i]]/sum(data[[i]]))
}

#compute distance matrix
data_distance = dist(data)

e= c()
#which sample is clone maximized in
iclonemax = apply(FUN = which.max, X = data, MARGIN = 1)
for(i in 1:length(data)){
  # print("anotha one")
  datai = data[,i]
  whichiclonemax = which(iclonemax == i)
  #order those clones which are maximized in i
  ordermaxi = order(-datai[whichiclonemax])
  
  e = append(e,whichiclonemax[ordermaxi])
  
  # what are the differences between this sample and the max contribution of each clone
  #maxdiff = data[[i]] - apply(X = data, MARGIN = 1, FUN = max)
  #select clones for which this is the max sample
  #e = append(e,order(-maxdiff)[1:length(which(iclonemax == i))])
}

png(paste("cluster curve", folder,".png"),width=size*ppi,height=size*ppi,res=ppi)
par(new = F, mar = c(7,5,3,1), xpd = TRUE)
for(i in 1:length(data)){
	ordered_data = data[(e),i]
	vector_tobeplotted = cumsum(ordered_data)
  plot(vector_tobeplotted, col = (rainbow(256, s = 1, v = 1, start = 0, end = 0.75, alpha = 1))[i*256/length(data)],  
       pch = 16, ylim = c(0,1), xlim = c(0,length(data[[1]])), ylab= "", 
       xlab = "", cex.lab = 1.5, cex.axis = 1.7)
	par(new=T)
}

#legend(2,1, names, title = "Time (Months)",cex = 1.8, col = (rainbow(256, s = 1, v = 1, start = 0, end = 0.75, alpha = 1))[seq(from = (256/length(data)), to = 256, by = 256/length(data))], pch = 16)
mtext("Clones sorted by time of maximum contribution", side=1, line=3, cex = 2)
mtext("Cumulative Clone Size", side=2, line=3, cex = 2)
mtext(paste(folder, "clone size distribution"), side=3, line=1, cex = 2)
dev.off()

png(paste("cluster curve legend", folder,".png"),width=size*ppi,height=size*ppi,res=ppi)
plot.new()
legend(0,1, names, title = "Months \n post-transplant",bty = 'n', cex = 1.4, col = (rainbow(256, s = 1, v = 1, start = 0, end = 0.75, alpha = 1))[seq(from = (256/length(data)), to = 256, by = 256/length(data))], pch = 16)
dev.off()
}
