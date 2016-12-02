basethresholdinsamples = function(data= zh33alldatanothresh, samples = as.vector(zh33indexmatrixnonNA), thresh = 100){
  
  clone_detected = array(dim = c(length(data[[1]]), length(samples)))
  clone_detected[is.na(clone_detected)] = 0

  for(i in 1:length(samples)){
    clone_detected[(data[[samples[i]]] > 100),i] = 1
  }
  
  overthreshold = c()
  for(i in 1:length(data[[1]])){
    if(any(clone_detected[i,]) == 1){
      overthreshold = append(overthreshold, i)
    } else {
    }
  }
  
  data = data[overthreshold,]
  
  return(data)
  
}

  