#Samson Koelle
#copywrite Cynthia Dunbar Lab
#10-18-16

#plots the heatmap of all clones binarized and sorted by time of emergence
BCLHbinarizeall = function(data = zh33alldata, names = zh33names, folder  = "", n_clones = 10) {

names = names[(apply(X = data, FUN = sum, MARGIN = 2) >0)]
data = data[(apply(X = data, FUN = sum, MARGIN = 1) >0),]
#get column sums.  this cuts out null columns
data = data[,(apply(X = data, FUN = sum, MARGIN = 2) >0)]

totalreads = apply(FUN = sum, MARGIN = 2, X = data)
for(i in 1:length(data)){
	#data[[i]]=binarize((data[[i]]/totalreads[i]), threshold = 0)
	data[[i]]=binarize((data[[i]]/totalreads[i]), threshold = 0.0005)
}
data = data[(apply(X = data, FUN = sum, MARGIN = 1) >0),]
data = data[do.call(order, -as.data.frame(data)),]
#e = hclust(dist(data))$order

png(paste("logscale",folder,".png"),width=6*ppi,height=6*ppi,res=ppi)
par(new = F, mar = c(7,5,3,1), xpd = TRUE)
#image(t(as.matrix(data[e,])), col = c("white","gray"), 
image(t(as.matrix(data[,])), col = c("white","gray"), 
       xaxt = 'n', yaxt = 'n')
axis(side = 1, at = seq(from = 0, to = 1, by = (1/(length(names)-1))), labels = names, las = 2)
mtext("Months post-transplant", side=1, line=3, cex = 2)
mtext(paste(dim(data)[1], " detected clones"), side = 2, line = 2, cex = 2)
mtext(folder, side = 3, line = 1, cex = 2)

dev.off()




}
