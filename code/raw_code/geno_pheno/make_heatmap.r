library("gplots")
library("viridis")

setwd("C:\\Users\\desch\\git\\enhancer-clustering\\data\\mohumaraba\\filt_k7000p06_cosine_norm-strict")
for (j in 1:117) {
  temp = read.csv(paste0("cluster", j, ".csv", collapse=""), sep=",", header=F)
  temp1 <- as.matrix(temp)
  if (nrow(temp1) > 9){
  #Sketchy color scheme: -1 -> gray, 0-1 -> light green-indigo via Viridis
  #Strange scale is to exclude the vellow part of viridis
  colors = c(seq(-1,-0.001,length=18),seq(-0.001,1,length=30))
  my_palette <- colorRampPalette(c(rep("gray", 18), rev(viridis(30))))(n = 47)
  
  png(file = paste0("cluster", j, ".png", collapse=""), width=1500, height=1200)
  #heatmap.2(temp1, col=rev(viridis(50)), Rowv=FALSE, Colv=FALSE, trace='none')
  heatmap.2(temp1, col=my_palette, Rowv=FALSE, Colv=FALSE, trace='none', symm=F,symkey=F,symbreaks=T, scale="none")
  dev.off()
  }
}

library("fitdistrplus")
setwd("C:\\Users\\desch\\git\\enhancer-clustering\\data")
temp = read.csv("k3000_pearson.csv", sep=",", header=F)
temp1 = as.vector(as.matrix(temp))
temp2 = (2 * (temp1 - 1)) + 1
hist(temp1)
fitdist(temp1, "norm", fix.arg = list(mean = 1), method="mle")
plot(res)

setwd("C:\\Users\\desch\\git\\enhancer-clustering\\data")
temp = read.csv("SomeEnhancersUnique.txt", sep="\t", header=F)
temp$V1 = NULL
temp1 = as.vector(as.matrix(temp))
temp2 = temp1[which(temp1 > -0.5)]
png(file = paste0("te.png", collapse=""), width=1500, height=1200)
hist(temp2, col="blue", main = "", )
abline(v = 0.25, col = "black", lty=2, lwd=2)
abline(v = 0.75, col = "black", lty=2, lwd=2)
dev.off()
fitdist(temp1, "norm", fix.arg = list(mean = 1), method="mle")
plot(res)

