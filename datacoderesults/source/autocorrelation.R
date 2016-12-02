#Samson Koelle
#Copywrite Cynthia Dunbar Lab
#11-16-16

#autocorrelation tracks the pearson correlation of the clonality of a cell type at a certain time point all other time points
autocorrelation = function(data = zh33alldata[,as.vector(zh33indexmatrixnoNA[,4])], names = zh33timepointnames, timepoints = zh33timepoints, folder = "ZH33 Gr"){


for(i in 1:length(data)){
	data[[i]]=(data[[i]]/sum(data[[i]]))
}
ntimepoints = length(data)
correlations = cor(data)
corforplot = array(dim = c(ntimepoints, ntimepoints))
for(i in 2:ntimepoints){
  for(j in 1:i){
    corforplot[i,j] = correlations[i,j]
  }
}

corforplot = t(corforplot)
corforplot[1,1] = 1
corforplot[corforplot <= 0] = 0

png(paste("autocorrelation", folder,".png"),width=size*ppi,height=size*ppi,res=ppi)
par(new = F, mar = c(7,5,3,1), xpd = TRUE)
for(i in 1:ntimepoints){
  plot(0:(ntimepoints-1), corforplot[i,], col = (rainbow(256, s = 1, v = 1, start = 0, end = 0.75, alpha = 1))[i*256/length(data)],  
       pch = 16, ylim = c(0,1), xlim = c(0,ntimepoints), ylab= "", 
       xlab = "", yaxt = 'n', xaxt = "n", cex.lab = 1.5, cex.axis = 1.25, cex = 2)
  lines(0:(ntimepoints-1), corforplot[i,], col = (rainbow(256, s = 1, v = 1, start = 0, end = 0.75, alpha = 1))[i*256/length(data)],
        ylim = c(0,1), xlim = c(0,ntimepoints), ylab= "", lwd = 2,
        xlab = "")
	par(new=T)
}
for(i in 1:length(timepoints)){
  axis(side = 1, at = (i-1),  labels = timepoints[i],  las = 2, cex.axis = 1.7
       #, col.axis = (rainbow(256, s = 1, v = 1, start = 0, end = 0.75, alpha = 1))[i*256/length(data)]
  )
}
axis(side = 2,  at= c(0,.2,.4,.6,.8,1), labels = c("<= 0",.2,.4,.6,.8,1), cex.axis = 1.7)
mtext("Months post-transplant", side=1, line=4, cex = 2)
mtext("Pearson Correlation", side=2, line=3, cex = 2)
mtext(paste(folder, "autocorrelation"), side=3, line=1, cex = 2)
dev.off()

#Generate legend
png(paste("autocor legend", folder,".png"),width=1.4*size*ppi,height=size*ppi,res=ppi)
par(new = F, mar = c(7,6,3,7), xpd = TRUE)
plot(1, cex = .001, xlim = c(0,1), ylim = c(0,1))
legend(0,1, timepoints, title = "Months \n post-transplant",cex = .7,
       col = (rainbow(256, s = 1, v = 1, start = 0, end = 0.75, alpha = 1))[seq(from = (256/length(data)), 
                                                                                to = 256, by = 256/length(data))], pch = 16, bty = 'n')
dev.off()

}
