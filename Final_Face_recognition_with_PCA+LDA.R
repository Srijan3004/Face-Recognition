#We start with a bunch of facial photographs, 5 photos of 10 persons each. 
#All the images are 200$\times$180 in size. So they are 36000 dimensional 
#vectors. Our aim is to reduce the dimension as much as possible without 
#losing our ability to recognise the faces.

#We start with loading the images.

library(jpeg)
#------
loadImages = function() {
  name = c("male.pacole","male.pspliu","male.sjbeck","male.skumar","male.rsanti",
           "female.anpage","female.asamma","female.klclar","female.ekavaz","female.drbost")
  #------  
  x = matrix(0,nrow=10*5,ncol=200*180)
  #------
  k = 0
  
  for(i in 1:10) {
    for(j in 1:5) {
      k = k + 1
      #The first argument below should be the location of the face folder.
      filename = paste("C:/Users/user/Downloads/facesbwandcol.tar/facesbwandcol.tar/facesbwandcol/grayfaces/",name[i],".",j,".jpg",sep="")
      x[k,] = as.vector(readJPEG(filename))
    }
  }
  #------
  return(x)
}
#---------------
#Next, we perform PCA. The following function does precisely this, but without using R's built-in tools like
#prcomp or princomp. We perform the operations explicitly so that you can appreciate what goes on
#"behind the scene".

process = function(x) {
  meanx = apply(x,2,mean)
  
  y = scale(x,scale=FALSE) #Rows of y are the cases
  A = y %*% t(y)
  eig = eigen(A)
  
  P = t(y) %*% eig$vec[,-50] #Columns of P are e vectors of Y'Y
  
  Q = apply(P,2,function(x) x/sqrt(sum(x*x))) #Columns of Q form onb for rowspace of Y
  
  
  scores = y %*% Q
  #scores[i,] is for i-th image
  return(list(centre=meanx,onb=Q,scores=scores,values=eig$values))
}
#------------------
#Each principal component is again a 200*180 dimensional vector, and hence may
#be considered as an image. It might be instructive to take a look at these.
#The following function helps you to just that.

showFace = function(newCoord, i) {
  plot(1:2,ty='n',main="0")
  y = abs(newCoord$onb[,i])
  extreme = range(y)
  y = (y-extreme[1])/(extreme[2]-extreme[1])
  dim(y) = c(200,180)
  rasterImage(as.raster(y),1,1,2,2)
}
#----------------
#The next function reconstructs the faces from the principal components.
#It starts with the "mean face" and then gradually adds the details.
#You'll need to hit "enter" to step through the process.

showSteps = function(newCoord,i) {
  meanx = newCoord$centre
  Q = newCoord$onb
  scores = newCoord$scores
  values = newCoord$values
  expl = 100*cumsum(newCoord$values)/sum(newCoord$values)
  
  coeff = as.vector(scores[i,])
  
  plot(1:2,ty='n',main="0")
  y = meanx
  dim(y) = c(200,180)
  rasterImage(as.raster(y),1,1,2,2)
  readline()
  for(k in 1:49) {
    if(k==1)
      temp = Q[,1]*coeff[1]
    else
      temp=Q[,1:k] %*% as.vector(coeff[1:k])
    recons = meanx + temp
    recons[recons<0]=0
    recons[recons>1]=1
    dim(recons) = c(200,180)
    plot(1:2,ty='n',main=paste(k,":",values[k],", ",expl[k]))
    rasterImage(as.raster(recons),1,1,2,2)
    readline()
  }
}
#----------------
#Hence we are done with PCA, and we found out that actually the big 36000 dimensional dataset can be mostly approximated by 10 things only.
#So, we will now use our data set 50*10 and use LDA upon it.

Cluster_identify = function(newcoord){
  scores = newcoord$scores
  R=matrix(0, nrow = 50, ncol = 10)
  for(i in 1:10)
  {
    R[,i]=scores[,i]
  }
  X1=R[1:5,]
  X2=R[6:10,]
  X3=R[11:15,]
  X4=R[16:20,]
  X5=R[21:25,]
  X6=R[26:30,]
  X7=R[31:35,]
  X8=R[36:40,]
  X9=R[41:45,]
  X10=R[46:50,]
  Y1=scale(X1,scale=FALSE)
  Y2=scale(X2,scale=FALSE)
  Y3=scale(X3,scale=FALSE)
  Y4=scale(X4,scale=FALSE)
  Y5=scale(X5,scale=FALSE)
  Y6=scale(X6,scale=FALSE)
  Y7=scale(X7,scale=FALSE)
  Y8=scale(X8,scale=FALSE)
  Y9=scale(X9,scale=FALSE)
  Y10=scale(X10,scale=FALSE)
  m1=apply(X1,2,mean)
  m2=apply(X2,2,mean)
  m3=apply(X3,2,mean)
  m4=apply(X4,2,mean)
  m5=apply(X5,2,mean)
  m6=apply(X6,2,mean)
  m7=apply(X7,2,mean)
  m8=apply(X8,2,mean)
  m9=apply(X9,2,mean)
  m10=apply(X10,2,mean)
  Z1=rbind(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10)
  W=4/49*t(Z1)%*%Z1
  B=5*cov(rbind(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10))
  library(matlib)
  H=inv(W)
  G=H%*%B
  eigG=eigen(G)
  #S=t(inv(M))%*%eigG$vectors
  S=eigG$vectors
  library(matlib)
  scores_new=R%*%inv(t(S))
  N=scores_new[,1]
  plot(N)
}

x=loadImages()
newcoord=process(x)
Cluster_identify(newcoord)
