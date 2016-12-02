#Samson Koelle
#copywrite Cynthia Dunbar Lab
#10-18-16

#This function plots the effect of the threshold on the diversity of each animal
diversitysensitivitybyanimal = function(data, indexmatricesnoNA, folder = "diversity sensitivity by animal",
                                        threshseq = seq(from = 0, to = 10000, by = 1000)){
  
  datathreshcombo = list()
  for(k in 1:length(data)){
    datathreshcombo[[k]] = list()
    for(h in 1:length(threshseq)){
      #this fun
      print(h)
      datathreshcombo[[k]][[h]] = thresholdinsamples(data = data[[k]], samples = indexmatricesnoNA[[k]], thresh = threshseq[h])
    }
   }
  
  diversity = list()
  for(k in 1:length(data)){
    
    #initialize diversity matrices
    diversity[[k]] = array(dim = length(threshseq))
    
    #cycle through threholds
    for(h in 1:length(threshseq)){
      
      #set 0 values to 1 to permit application of log function in diversity calculation
      datathreshcombo[[k]][[h]][datathreshcombo[[k]][[h]] == 0] = 1
      

            #this is always the case except for the null column in certain animals
            #should be length(data) type thing since 1
              diversity[[k]][h] = Hs(datathreshcombo[[k]][[h]][,as.vector(indexmatricesnoNA[[k]])])

          
        }
      }
  
  
  
  
  png(paste("diversitysensitivitybyanimal", folder,".png"),width=1.4*size*ppi,height=size*ppi,res=ppi)
  par(new = F, mar = c(7,6,3,7), xpd = TRUE)

  for(k in 1:length(data)){
    
        plot(threshseq, diversity[[k]], cex= 2, pch = 14+k, xlim = c(min(threshseq),max(threshseq)),
            ylim = c(0,max(unlist(diversity))), ylab = "", xlab = "", cex.lab = 3,
            cex.axis = 2)
        lines(threshseq, diversity[[k]], cex= 1,   xlim = c(min(threshseq),max(threshseq)),
              ylim = c(0,max(unlist(diversity))), ylab = "", xlab = "", cex.lab = 3,
              cex.axis = 2, lwd = 2)
        par(new = T)
  }
  mtext("Threshold", side=1, line=3, cex = 2)
  mtext("Shannon Diversity", side = 2, line = 3, cex = 2)
  mtext("Diversity Sensitivity to Threshold", side=3, line=1, cex = 2)
  #mtext(folder, side = 3, line = 1, cex = 3)
  #legend(45,8, c("T",'B','Mono',"Gr"), title = "Cell Type",cex = 1.6, 
  #col = c(1:4), pch = 15)
  dev.off()
  
  png(paste("diversity sensitivity legend", folder,".png"),width=1.4*size*ppi,height=size*ppi,res=ppi)
  par(new = F, mar = c(7,6,3,7), xpd = TRUE)
  plot(1,1, cex= .02, pch = 15, col = 1, xlim = c(0,10),
       ylim = c(-5,20), ylab = "", xlab = "", cex.lab = 3,
       cex.axis = 2)
  legend(5,8, c("ZH33",'ZG66','ZH19','ZJ31'), title = "Animal",cex = 1.6, 
         col = c(1), pch = 15:18)
  dev.off()
  
  
}
