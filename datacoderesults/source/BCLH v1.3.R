#Samson Koelle
#Copywrite Cynthia Dunbar Lab
#11-16-16

#BCLH stands for "Barcode Celltype Log Heatmap"
#It first selects the clones which are top 10 high contributors to any sample
#It then plots the contribution of those clones to all samples
#clones are organized by heirarchical clustering
BCLH = function(data = zh33alldata[,as.vector(zh33indexmatrixnoNA)], names = zh33allnames[as.vector(zh33indexmatrixnoNA)], n_clones = 10, folder = paste("ZH33"), cexset = .5) {


names = names[(apply(X = data, FUN = sum, MARGIN = 2) >0)]
data = data[(apply(X = data, FUN = sum, MARGIN = 1) >0),]
#cut out null columns
data = data[,(apply(X = data, FUN = sum, MARGIN = 2) >0)]

#normalize
for(i in 1:length(data)){
	data[[i]]=(data[[i]]/sum(data[[i]]))
}

samples = c(1:length(data))

#find the top clones
order.list = c()
for(i in 1:length(samples)){
	order.list = append(order.list, order(-data[[samples[i]]])[1:n_clones])
}
top_clones = unique(order.list)
data = data[top_clones,]

#list_topclones gets the position of the top clones in each sample
list_topclones = data
for(i in 1:length(samples)){
	list_topclones[[i]] = order(-data[[i]])
}

#is_a_topclone is for plotting
is_a_topclone = data.frame(matrix(ncol = length(data), nrow = length(data[[1]])))
for(i in 1:length(data)){
	is_a_topclone[list_topclones[1:n_clones,i],i] = "*"
}
is_a_topclone[is.na(is_a_topclone)] = " "

#the advantage of taking the log this way (as opposed to
#usual +1 before log adjustment is that all non-present clones will appear the same)
data.log = log(data)
data.log[data.log==-Inf] <- (min(data.log[data.log > -Inf]) - 1)

#get the order of the top clones
e = hclust(dist(data.log))$order

png(paste("logscale",folder,"BCLH.png"),width=6*ppi,height=6*ppi,res=ppi)

par(new = F, mar = c(7,5,3,1), xpd = TRUE)
image(t(as.matrix(data.log[e,])), col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), 
       xaxt = 'n', yaxt = 'n')
xpos = seq(from = 0, to = 1, by = 1 / (dim(data.log)[2] - 1))
ypos = seq(from = -.5 / (dim(data.log)[1] - 1), to = 1 - .5 / (dim(data.log)[1] - 1), by = 1 / (dim(data.log)[1] - 1))
for(i in 1:dim(data.log)[1]){
  for(j in 1:dim(data.log)[2]){
    text(x = xpos[j], y = ypos[i], labels = is_a_topclone[e[i],j], adj = c(0.5,0.5), cex = cexset)
  }
}

axis(side = 1, at = seq(from = 0, to = 1, by = (1/(length(names)-1))), labels = names, las = 2, cex.axis = .6, tick = F)
#mtext("Months post-transplant", side=1, line=3, cex = 2)
mtext("High-contributing clones", side = 2, line = 2, cex = 2)
mtext(folder, side = 3, line = 1, cex = 2)
dev.off()

#plot legend markers at 0, 0.0005, .005,.05
range = max(data.log) - min(data.log)
if(max(data.log) > log(0.05)){
markpoints = c(0,1 - ((max(data.log) - log(c(0.0005,0.005,.05))) / range))
markpointlabels = c("0", '0.0005','0.005','0.05')
} else {
  if(max(data.log) > log(0.005)){
    markpoints = c(0,1 - ((max(data.log) - log(c(0.0005,0.005,.05))) / range))
    markpointlabels = c("0", '0.0005','0.005','0.05')
} else {
  print("Max value too low.  Please input color scale manually.")
}
}
png(paste("logscale legend",folder,"BCLH.png"),width=6*ppi,height=6*ppi,res=ppi)
par(new = F, mar = c(7,5,3,1), xpd = TRUE)
image(matrix(c(1:100), nrow = 1), col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), xaxt= 'n', yaxt = 'n')
axis(2, at= markpoints, labels = markpointlabels, cex.axis = 1.75)
mtext(3, text = paste('Log normalized clone size'), line = 1.5, cex = 2)
dev.off()

}
