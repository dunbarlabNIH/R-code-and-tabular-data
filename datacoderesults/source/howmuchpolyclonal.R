howmuchpolyclonal = function(data = zh33alldata, indexmatrix = zh33indexmatrix){
  
  
  for(i in 1:length(data)){
    data[[i]]=(data[[i]]/sum(data[[i]]))
  }
  #not present is in all four lineages
  notpresent = c()
  presinone = c()
  for(i in 1:4){
    notpresent = union(notpresent, which(data[,indexmatrix[dim(indexmatrix)[1],i]] == 0))
    presinone = union(presinone, which(data[,indexmatrix[dim(indexmatrix)[1],i]] != 0))
  }
  
  #ispresent in all four lineages
  present = seq(from = 1, to = dim(data)[1], by = 1)
  ispresent = setdiff(present, notpresent)
  toreturn = array(dim = 3)
  toreturn[1] = sum(data[ispresent,indexmatrix[dim(indexmatrix)[1], 4]])
  toreturn[2] = length(ispresent) / length(presinone)
  
  return(toreturn)
  
  #dim(thresholdinsamples(data = zh33alldata, thresh = 2000, samples = as.vector(zh33indexmatrixnoNA[dim(zh33indexmatrixnoNA)[1],])))
}