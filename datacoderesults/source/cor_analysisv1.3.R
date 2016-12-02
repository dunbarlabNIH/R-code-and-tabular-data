#Samson Koelle
#Copywrite Cynthia Dunbar Lab
#11-16-16

#cor_analysis tracks the Pearson correlation between pairs of cell types at each timepoint in a single animal
cor_analysis = function(data = zh33alldata, indexmatrix = zh33indexmatrix, timepoints = zh33timepoints, folder = "", celltypenames = c("T", "B","Mono",'Gr')){

n_timepoints = length(timepoints)
for(i in 1:length(data)){
	data[[i]]=(data[[i]]/sum(data[[i]]))
}

cors = array(dim = c(dim(indexmatrix)[2], dim(indexmatrix)[2], dim(indexmatrix)[1]))

for(i in 1:(dim(indexmatrix)[2] - 1)){
	for(j in (i+1):dim(indexmatrix)[2]){
		print(j)
		for(k in 1:dim(indexmatrix)[1]){
	    if(!is.na(indexmatrix[k,i]) && !is.na(indexmatrix[k,j])){
		    cors[i,j,k] = cor(data[,indexmatrix[k,i]], data[,indexmatrix[k,j]])
	    } else {
	      cors[i,j,k] = NA
	    }
		}
	}
}

ymin = min(cors[!is.na(cors)])
ymax = max(cors[!is.na(cors)])

colors = alpha(colour = 1:4, alpha = c(1,.4,.8,1))
png(paste(folder,"cors.png"),width=size*ppi,height=size*ppi,res=ppi)

par(new = F, mar = c(7,5,3,1))
for(i in 1:(dim(indexmatrix)[2] - 1)){
	for(j in i:dim(indexmatrix)[2]){
		plot(cors[i,j,], pch = 15, col = colors[i], ylim = c(ymin, ymax), xaxt = "n", xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 2, cex = 2)
		lines(cors[i,j,], lty = 1, col = colors[i], lwd = 7)
		#lty = 2,1
		par(new = T)
		plot(cors[i,j,], pch = 18, col = 'white', ylim = c(ymin, ymax),  xaxt = "n", xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 2, cex = 2)
		lines(cors[i,j,], lty = 2, col = 'white', lwd = 7)
		par(new = T)
		plot(cors[i,j,], pch = 18, col = colors[j], ylim = c(ymin, ymax),  xaxt = "n", xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 2, cex = 2)
		lines(cors[i,j,], lty = 2, col = colors[j], lwd = 7)
		#lty =4,2
		par(new = T)
	}
}
axis(1, at=1:n_timepoints, labels = timepoints, las = 2, cex.axis = 2)
mtext("Months Post-Transplant", side=1, line=5, cex = 2)
mtext("Pearson Correlation", side = 2, line = 3, cex = 2)
mtext(folder, side = 3, line = 1, cex = 3)
dev.off()

#Generate legend (actual plot for convenience only)
png(paste("cor legend", folder,".png"),width=1.4*size*ppi,height=size*ppi,res=ppi)
par(new = F, mar = c(7,6,3,7), xpd = TRUE)
plot(1, cex = .0001, xlim = c(-1,3), ylim = c(-1,11))
count = 10
for(i in 1:(dim(indexmatrix)[2] - 1)){
  for(j in (i+1):dim(indexmatrix)[2]){
    polygon(x = c(0,0,1), y = c(count-1,count-.1,count-.1), col = colors[i], border = colors[i])
    polygon(x = c(0,1,1), y = c(count-1,count-1,count-.1), col = colors[j], border = colors[j])
    text(x = 1.1, y = count-.7, labels = paste(celltypenames[i], "and", celltypenames[j]), cex = 1.2, adj = c(0,0))
    count = count - 1
  }
}
dev.off()

}


