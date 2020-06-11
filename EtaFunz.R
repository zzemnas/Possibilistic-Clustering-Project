
################Function used for the computation of the centroids given by the possibilistic
################clustering algorithm.

possclust <-function (X, k, m, conv, maxit,RS=1,proto){
  # X: data
  # k: number of clusters
  # m: parameter of fuzziness
  # conv: convergence rate
  # maxit: maximum number of iterations
  # U: (possibilistic) membership matrix
  # D: squared distance matrix
  # H: centroid matrix
  # value: optimization function value
  #RS=random starts (defaults to one)
  #proto=initialization of H and U given by the FKM algorithm
  
  n = nrow(X)  # number of units/objects
  p = ncol(X)  # number of variables
  
  value=vector(length(RS),mode="numeric") #initialization of the values of the objective function
  it=vector(length(RS), mode="numeric") #vector of iterations (the length is given by the number of random starts)
  func.opt = 1e+10
  
  for(rs in 1:RS){
    # STEP 0 #
    U=proto$U  # initialization U^(0)
    D=matrix(0,nrow=n,ncol=k)
    E<-rep(0,k)
    # end STEP 0 #
    H=proto$H
    U.old=U+1
    iter=0
    
    while ((sum(abs(U.old-U))>conv) && (iter<maxit))  # covergence criteria
    {
      
      iter=iter+1
      U.old=U
      H.old=H
      # STEP 1 #
      for (c in 1:k) 
        H[c,]=(t(U[,c]^m)%*%X)/sum(U[,c]^m)   # compute matrix H
      #H is the centroid matrix
      
      # end STEP 1 #
      
      # STEP 2: this is the part that separates the FKM algorithm from the possibilistic one #
      for (c in 1:k) 
      {
        for (i in 1:n) 
        {
          D[i,c]=sum((H[c,]-X[i,])^2) # compute distance matrix D
        }
      }
      
      #the cycle that follows computes the values of eta based on the formula (9) of the paper
      #We could have also used the formula with the alfa-cut:in R, E[i]<-mean((D[which(U[,i]>alfa),i]))
      
      for (i in 1:k){
        E[i]<-sum(U[,i]^m*D[,i])/sum(U[,i])
      }
      
      for (c in 1:k)
      {
        for (i in 1:n)
        {
          U[i,c]=1/(1+(D[i,c]/E[c])^(1/(m-1)))   # compute membership matrix U
        }
        
        # end STEP 2 #
      }
    }
    
    
    func=sum((U^m)*D)+sum(E*(1-U)^m)   #possibilistic objective function
    value[rs]=func
    it[rs]=iter
    
    if (func<func.opt)  #check the behaviour of the objective function between random starts
    {
      U.opt=U
      H.opt=H 
      E.opt=E
      D.opt=D
      func.opt=func
    }
  }
  
  
  out=list()
  out$U=U.opt
  out$D=D.opt
  out$H=H.opt
  out$value=value
  out$iter=it
  out$eta=E.opt
  
  return(out)
}
