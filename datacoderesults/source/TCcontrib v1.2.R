#Samson Koelle
#copywrite Cynthia Dunbar Lab
#10-18-16

#plots the summed contribution not of the top 10 clones in a specific sample, but rather the top 10 clones in any sample
TCcontrib = function(data = zh33alldata, indexmatrix = zh33indexmatrix, folder  = "", n_clones = 10, timepoints =  zh33timepoints) {

n_timepoints = length(timepoints)
ntimepoints = n_timepoints
data_save = data

for(i in 1:length(data)){
	data[[i]]=(data[[i]]/sum(data[[i]]))
}
lineages = c(1:length(data))
#findtopclones(data = data, lineages= lineages
# data =  data[apply(data[,-1], 1, function(x) !all(x==0)),]


order.list = c()
for(i in 1:length(as.vector(indexmatrix))){
	order.list = append(order.list, order(-data[[as.vector(indexmatrix)[i]]])[1:n_clones])
}

#contam = order(-data[[1]])[1:6]

top_clones = unique(order.list)
#top_clones = setdiff(top_clones,contam)

data = data[top_clones,]

png(paste("TC contrib", folder,".png"),width=size*ppi,height=size*ppi,res=ppi)
par(new = F, mar = c(7,5,4,1))
for(i in 1:dim(indexmatrix)[2]){
  nonNAtimepoints = !is.na(indexmatrix[,i])
  dataforapply = data[,indexmatrix[!is.na(indexmatrix[,i]),i]]
	datatobeplotted = apply(X = dataforapply, FUN = sum, MARGIN = 2)
	plot((1:ntimepoints)[nonNAtimepoints], datatobeplotted, col = i, ylim = c(0,1), cex.axis = 2, 
	     xlim = c(1,ntimepoints), pch = 16,xlab = "", ylab= "",  cex = 2, xaxt = "n")
	lines((1:ntimepoints)[nonNAtimepoints], datatobeplotted, col = i, ylim = c(0,1), lwd = 7)
	par(new=T)
}
axis(1, at=1:n_timepoints, labels = timepoints, las = 2, cex.axis = 2)
mtext(folder, side=3, line=1.7, cex = 2.3)
mtext("Months Post-Transplant", side=1, line=5, cex = 3)
mtext("Size", side=2, line=3, cex = 3)
dev.off()

}