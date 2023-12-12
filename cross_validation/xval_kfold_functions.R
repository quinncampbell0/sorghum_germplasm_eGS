######################################################################################
######### Cross Validation - K-fold (2,5,10) - Exponential kernel - Minicore #########
######################################################################################
# Define the function for k-fold cross validation (i.e. 10-fold)
k.xval.EXP <- function(g.in = NULL, # A n x m genotype matrix (coded as -1, 0, 1), with the row.names of g.in equal to the line names
                       y.in = NULL, # A n x t matrix of phentypes, with the row.names of y.in equal to the line names
                       k.fold = 10, # The number of "folds" to perfom cross-validation with
                       reps = 25, # The number of iterations to divide the lines into k "folds"
                       model = "Exponential kernel",
                       y.trainset = NULL, # identical columns to y.in, but consisting of only the prediction data
                       k_dist = k_dist
) {
  
  # Error reporting
  if(is.null(g.in)) { stop("The genotype matrix was not specified.") }
  if(is.null(y.in)) { stop("The phenotype data was not specified.") }
  if(nrow(g.in) != nrow(y.in)) { stop("The number of lines in the genotype matrix and in the phenotype data do not match.") }
  
  # Create a function for using rrBLUP
  cv.rrblup <- function(y.train = y.train, k_dist = k_dist) {
    solve.out <- kin.blup(data = y.train,  geno = "gen_id", pheno = trait, GAUSS = T, K = k_dist)
    return(solve.out$g)
  }
  
  # Create a function that detects the model and trait and implement cross-validation
  implement.xval <- function(g.in, y.in, trait, reps, k.fold, model) {
    
    n.lines <- nrow(y.datexist)
    
    # Notify which trait is being cross-validated
    writeLines(paste("\nNow conducting k-fold cross validation on ", trait, " using ", model, ".", sep = ""))
    
    # sapply function of all reps of k.fold cv for a trait
    acc.out <- sapply(X = 1:reps, FUN = function(x) {
      
      # Design a list of folds
      k.it <- split(x = sample(x = 1:n.lines, size = n.lines), f = factor(cut(1:n.lines, breaks = k.fold)))
      # Apply the function of the list of folds
      k.fold.out <- lapply(X = k.it, FUN = function(x) {
        pred <- row.names(y.datexist)[x] # Names of prediction set
        train <- setdiff(row.names(y.datexist), pred) # Names of training set
        #g.train <- g.in[train,] # Set training genos
        #g.pred <- g.in[pred,] # Set prediction genos
        y.train <- y.in[train,]
        # Marker effects
        BLUP <- cv.rrblup(y.train = y.train, k_dist = k_dist)
        val_set <- setdiff(row.names(g.in),train)
        BLUP <- BLUP[val_set]
        BLUP <- as.matrix(BLUP)
        return(BLUP) # correlate GEBVs to actual phenos
      })
      
      # Create a vector of GEBVs
      GEBV <- do.call("rbind", k.fold.out)
      # Sort it by lines
      y.pred <- y.in[row.names(GEBV),]
      # GEBV <- GEBV[order(row.names(GEBV))]
      # Calculate correlation
      acc <- cor(GEBV, y.pred[,trait], use = "complete.obs")
      return(acc)
    })
    
    # Return the mean and sd of the correlations
    r.mean <- mean(acc.out)
    r.sd <- sd(acc.out)
    return(cbind(trait, model, r.mean, r.sd))
    
  } # Close the function
  
  # Create the output matrix
  results.out <- as.data.frame(matrix(nrow = ncol(y.in), ncol = 4))
  colnames(results.out) <- c("trait", "model", "r.mean", "r.sd")
  
  # Iterate through each trait and perform cross-validation
  for (t in 2:length(colnames(y.in))) {
    trait <- colnames(y.in)[t]
    y.datexist <- y.trainset[which(!is.na(y.trainset[,t])),]
    results.out[t,] <- implement.xval(g.in = g.in, y.in = y.in, trait = trait, reps = reps, k.fold = k.fold, model = model)
  }
  
  # Return a list of results
  cross.val.results <- list(xval.result = results.out, folds = k.fold, reps = reps, model = "Gaussian kernel")
  return(cross.val.results)
} # Close the function


Kernel_computation=function(X,name, degree, nL){
  
  p=ncol(X)
  
  x1=X
  
  x2=X
  
  d=degree
  
  ############Polynomial kernel##################
  
  K.Polynomial=function(x1, x2=x1, gamma=1, b=0, d=3)
    
  { (gamma*(as.matrix(x1)%*%t(x2))+b)^d}
  
  ############Sigmoid kernel####################
  
  K.Sigmoid=function(x1,x2=x1, gamma=1, b=0)
    
  { tanh(gamma*(as.matrix(x1)%*%t(x2))+b) }
  
  ############Gaussian kernel##################
  
  l2norm=function(x){sqrt(sum(x^2))}
  
  K.Gaussian=function(x1,x2=x1, gamma=1){
    
    exp(-gamma*outer(1:nrow(x1<- as.matrix(x1)), 1:ncol(x2<- t(x2)),
                     
                     Vectorize(function(i, j) l2norm(x1[i,]-x2[,j])^2)))}
  
  ##########Arc-cosine kernel with 1 hidden layer
  
  K.AK1_Final<-function(x1,x2){
    
    n1<-nrow(x1)
    
    n2<-nrow(x2)
    
    x1tx2<-x1%*%t(x2)
    
    norm1<-sqrt(apply(x1,1,function(x) crossprod(x)))
    
    norm2<-sqrt(apply(x2,1,function(x) crossprod(x)))
    
    costheta = diag(1/norm1)%*%x1tx2%*%diag(1/norm2)
    
    costheta[which(abs(costheta)>1,arr.ind = TRUE)] = 1
    
    theta<-acos(costheta)
    
    normx1x2<-norm1%*%t(norm2)
    
    J = (sin(theta)+(pi-theta)*cos(theta))
    
    AK1 = 1/pi*normx1x2*J
    
    AK1<-AK1/median(AK1)
    
    colnames(AK1)<-rownames(x2)
    
    rownames(AK1)<-rownames(x1)
    
    return(AK1)
    
  }
  
  ####Kernel Arc-Cosine with deep=L#########
  
  AK_L_Final<-function(AK1,nL){
    
    n1<-nrow(AK1)
    
    n2<-ncol(AK1)
    
    AKl1 = AK1
    
    for (l in 1:nL){
      
      AKAK<-tcrossprod(diag(AKl1),diag(AKl1))
      
      costheta<-AKl1*(AKAK^(-1/2))
      
      costheta[which(costheta>1,arr.ind = TRUE)] = 1
      
      theta<-acos(costheta)
      
      AKl<-(1/pi)*(AKAK^(1/2))*(sin(theta)+(pi-theta)*cos(theta))
      
      AKl1 = AKl
      
    }
    
    AKl<-AKl/median(AKl)
    
    rownames(AKl)<-rownames(AK1)
    
    colnames(AKl)<-colnames(AK1)
    
    return(AKl)
    
  }
  
  ########Exponencial Kernel############
  
  K.exponential=function(x1,x2=x1, gamma=1){
    
    exp(-gamma*outer(1:nrow(x1<- as.matrix(x1)), 1:ncol(x2<- t(x2)),
                     
                     Vectorize(function(i, j) l2norm(x1[i,]-x2[,j]))))}
  
  if (name=="Linear") {
    
    K=X%*%t(X)/p
    
  } else if (name=="Polynomial") {
    
    K=K.Polynomial(x1=x1, x2=x1, gamma=1/p, b=0, d=d)
    
  } else if (name=="Sigmoid") {
    
    K=K.Sigmoid(x1=x1, x2=x1, gamma=1/p, b=0)
    
  }else if (name=="Gaussian") {
    
    K=K.Gaussian(x1=x1, x2=x1, gamma=1/p)
    
  } else if (name=="AK1") {
    
    K= K.AK1_Final(x1=x1, x2=x1)
    
  } else if (name=="AKL") {
    
    AK1=K.AK1_Final(x1=x1, x2=x1)
    
    K=AK_L_Final(AK1=AK1,nL=nL)
    
  } else {
    
    K=K.exponential(x1=x1,x2=x1,gamma=1/p)
    
  }
  
}

###########################################################################
######### Cross Validation - K-fold (2,5,10) - Gaussian kernel - Minicore #########
###########################################################################

# Define the function for k-fold cross validation (i.e. 10-fold)
k.xval.GAUSS <- function(g.in = NULL, # A n x m genotype matrix (coded as -1, 0, 1), with the row.names of g.in equal to the line names
                         y.in = NULL, # A n x t matrix of phentypes, with the row.names of y.in equal to the line names
                         k.fold = 10, # The number of "folds" to perfom cross-validation with
                         reps = 25, # The number of iterations to divide the lines into k "folds"
                         model = "Gaussian kernel",
                         y.trainset = NULL, # identical columns to y.in, but consisting of only the prediction data
                         k_dist = k_dist
) {
  
  # Error reporting
  if(is.null(g.in)) { stop("The genotype matrix was not specified.") }
  if(is.null(y.in)) { stop("The phenotype data was not specified.") }
  if(nrow(g.in) != nrow(y.in)) { stop("The number of lines in the genotype matrix and in the phenotype data do not match.") }
  
  # Create a function for using rrBLUP
  cv.rrblup <- function(y.train = y.train, k_dist = k_dist) {
    solve.out <- kin.blup(data = y.train,  geno = "gen_id", pheno = trait, GAUSS = T, K = k_dist)
    return(solve.out$g)
  }
  
  # Create a function that detects the model and trait and implement cross-validation
  implement.xval <- function(g.in, y.in, trait, reps, k.fold, model) {
    
    n.lines <- nrow(y.datexist)
    
    # Notify which trait is being cross-validated
    writeLines(paste("\nNow conducting k-fold cross validation on ", trait, " using ", model, ".", sep = ""))
    
    # sapply function of all reps of k.fold cv for a trait
    acc.out <- sapply(X = 1:reps, FUN = function(x) {
      
      # Design a list of folds
      k.it <- split(x = sample(x = 1:n.lines, size = n.lines), f = factor(cut(1:n.lines, breaks = k.fold)))
      # Apply the function of the list of folds
      k.fold.out <- lapply(X = k.it, FUN = function(x) {
        pred <- row.names(y.datexist)[x] # Names of prediction set
        train <- setdiff(row.names(y.datexist), pred) # Names of training set
        #g.train <- g.in[train,] # Set training genos
        #g.pred <- g.in[pred,] # Set prediction genos
        y.train <- y.in[train,]
        # Marker effects
        BLUP <- cv.rrblup(y.train = y.train, k_dist = k_dist)
        val_set <- setdiff(row.names(g.in), train)
        BLUP <- BLUP[val_set]
        BLUP <- as.matrix(BLUP)
        return(BLUP) # correlate GEBVs to actual phenos
      })
      
      # Create a vector of GEBVs
      GEBV <- do.call("rbind", k.fold.out)
      # Sort it by lines
      y.pred <- y.in[row.names(GEBV),]
      # GEBV <- GEBV[order(row.names(GEBV))]
      # Calculate correlation
      acc <- cor(GEBV, y.pred[,trait], use = "complete.obs")
      return(acc)
    })
    
    # Return the mean and sd of the correlations
    r.mean <- mean(acc.out)
    r.sd <- sd(acc.out)
    return(cbind(trait, model, r.mean, r.sd))
    
  } # Close the function
  
  # Create the output matrix
  results.out <- as.data.frame(matrix(nrow = ncol(y.in), ncol = 4))
  colnames(results.out) <- c("trait", "model", "r.mean", "r.sd")
  
  # Iterate through each trait and perform cross-validation
  for (t in 2:length(colnames(y.in))) {
    trait <- colnames(y.in)[t]
    y.datexist <- y.trainset[which(!is.na(y.trainset[,t])),]
    results.out[t,] <- implement.xval(g.in = g.in, y.in = y.in, trait = trait, reps = reps, k.fold = k.fold, model = model)
  }
  
  # Return a list of results
  cross.val.results <- list(xval.result = results.out, folds = k.fold, reps = reps, model = "Gaussian kernel")
  return(cross.val.results)
} # Close the function


# Define the function for k-fold cross validation (i.e. 10-fold)
k.xval <- function(g.in = NULL, # A n x m genotype matrix (coded as -1, 0, 1), with the row.names of g.in equal to the line names
                   y.in = NULL, # A n x t matrix of phentypes, with the row.names of y.in equal to the line names
                   k.fold = 10, # The number of "folds" to perfom cross-validation with
                   reps = 25, # The number of iterations to divide the lines into k "folds"
                   model = "rrBLUP",
                   y.trainset = NULL # identical columns to y.in, but consisting of only the prediction data
) {
  
  # Error reporting
  if(is.null(g.in)) { stop("The genotype matrix was not specified.") }
  if(is.null(y.in)) { stop("The phenotype data was not specified.") }
  if(nrow(g.in) != nrow(y.in)) { stop("The number of lines in the genotype matrix and in the phenotype data do not match.") }
  
  # Create a function for using rrBLUP
  cv.rrblup <- function(y.train = y.train, g.train = g.train) {
    solve.out <- mixed.solve(y = y.train,  Z = g.train, SE = F, return.Hinv = F)
    return(solve.out$u)
  }
  
  # Create a function that detects the model and trait and implement cross-validation
  implement.xval <- function(g.in, y.in, y.datexist, trait, reps, k.fold, model) {
    
    n.lines <- nrow(y.datexist)
    
    # Notify which trait is being cross-validated
    writeLines(paste("\nNow conducting fractional cross validation on ", trait, " using ", model, ".", sep = ""))
    
    # sapply function of all reps of k.fold cv for a trait
    acc.out <- sapply(X = 1:reps, FUN = function(x) {
      
      # Design a list of folds
      k.it <- split(x = sample(x = 1:n.lines, size = n.lines), f = factor(cut(1:n.lines, breaks = k.fold)))
      # Apply the function of the list of folds
      k.fold.out <-lapply(X = k.it, FUN = function(x) {
        pred <- row.names(y.datexist)[x] # Names of training set
        train <- setdiff(row.names(y.datexist), pred) # Names of prediction set
        g.train <- g.in[train,] # Set training genos
        #g.pred <- g.in[pred,] # Set prediction genos
        y.train <- as.vector(y.in[train,trait])
        
        # Marker effects
        u.hat <- cv.rrblup(y.train = y.train, g.train = g.train)
        
        # GEBVs and correlation
        val_set <- setdiff(row.names(g.in), train)
        g.pred <- g.in[val_set,]
        GEBV <- g.pred %*% u.hat
        #browser()
        return(GEBV) # correlate GEBVs to actual phenos
      })
      # Create a vector of GEBVs
      GEBV <- do.call("rbind", k.fold.out)
      # Sort it by lines
      y.pred <- y.in[row.names(GEBV),]
      #GEBV <- GEBV[order(row.names(GEBV))]
      # Calculate correlation
      acc <- cor(GEBV, y.pred[,trait], use = "complete.obs")
      #browser()
      return(acc)
    })
    
    # Return the mean and sd of the correlations
    r.mean <- mean(acc.out)
    r.sd <- sd(acc.out)
    return(cbind(trait, model, r.mean, r.sd))
    
  } # Close the function
  
  # Create the output matrix
  results.out <- as.data.frame(matrix(nrow = ncol(y.in), ncol = 4))
  colnames(results.out) <- c("trait", "model", "r.mean", "r.sd")
  
  # Iterate through each trait and perform cross-validation
  for (t in 1:length(colnames(y.in))) {
    trait <- colnames(y.in)[t]
    y.datexist <- y.trainset[which(!is.na(y.trainset[,t])),]
    results.out[t,] <- implement.xval(g.in = g.in, y.in = y.in, y.datexist = y.datexist, trait = trait, 
                                      reps = reps, k.fold = k.fold, model = "rrBLUP")
  }
  
  # Return a list of results
  cross.val.results <- list(xval.result = results.out, folds = k.fold, reps = reps, model = model)
  return(cross.val.results)
} # Close the function

# Define the function for k-fold cross validation (i.e. 10-fold)
k.xval.BayesCpi <- function(g.in = NULL, # A n x m genotype matrix (coded as -1, 0, 1), with the row.names of g.in equal to the line names
                            y.in = NULL, # A n x t matrix of phentypes, with the row.names of y.in equal to the line names
                            k.fold = 10, # The number of "folds" to perfom cross-validation with
                            reps = 25, # The number of iterations to divide the lines into k "folds"
                            model = "BayesCpi",
                            y.trainset = NULL, # identical columns to y.in, but consisting of only the prediction data
                            niter = niter,
                            nburn = nburn
) {
  
  # Error reporting
  if(is.null(g.in)) { stop("The genotype matrix was not specified.") }
  if(is.null(y.in)) { stop("The phenotype data was not specified.") }
  if(nrow(g.in) != nrow(y.in)) { stop("The number of lines in the genotype matrix and in the phenotype data do not match.") }
  
  # Create a function for using rrBLUP
  cv.rrblup <- function(y.train = y.train) {
    solve.out <- eval(parse(text = paste("bayes(",trait,'~1, data=y.train, M=g.in, M.id=row.names(g.in), 
                       method="BayesCpi", niter = ', niter,', nburn=',nburn,', thin=5,threads=1, verbose = F)')))
    return(solve.out$g)
  }
  
  # Create a function that detects the model and trait and implement cross-validation
  implement.xval <- function(g.in, y.in, trait, reps, k.fold, model) {
    
    n.lines <- nrow(y.datexist)
    
    # Notify which trait is being cross-validated
    writeLines(paste("\nNow conducting k-fold cross validation on ", trait, " using ", model, ".", sep = ""))
    
    # sapply function of all reps of k.fold cv for a trait
    acc.out <- sapply(X = 1:reps, FUN = function(x) {
      
      # Design a list of folds
      k.it <- split(x = sample(x = 1:n.lines, size = n.lines), f = factor(cut(1:n.lines, breaks = k.fold)))
      # Apply the function of the list of folds
      k.fold.out <- lapply(X = k.it, FUN = function(x) {
        pred <- row.names(y.datexist)[x] # Names of prediction set
        train <- setdiff(row.names(y.datexist), pred) # Names of training set
        #g.train <- g.in[train,] # Set training genos
        #g.pred <- g.in[pred,] # Set prediction genos
        y.train <- y.in
        val_set <- setdiff(row.names(g.in), train)
        y.train[val_set,trait] <- NA
        
        # Marker effects
        BLUP <- cv.rrblup(y.train = y.train)
        row.names(BLUP) <- BLUP$id
        BLUP <- BLUP[val_set,]
        #BLUP <- as.matrix(BLUP)
        return(BLUP) # correlate GEBVs to actual phenos
      })
      
      # Create a vector of GEBVs
      GEBV <- do.call("rbind", k.fold.out)
      
      # Sort it by lines
      #browser()
      y.pred <- y.in[GEBV$id,trait]
      #GEBV <- GEBV[order(row.names(GEBV))]
      # Calculate correlation
      acc <- cor(GEBV$gebv, y.pred, use = "complete.obs")
      return(acc)
    })
    
    # Return the mean and sd of the correlations
    r.mean <- mean(acc.out)
    r.sd <- sd(acc.out)
    return(cbind(trait, model, r.mean, r.sd))
    
  } # Close the function
  
  # Create the output matrix
  results.out <- as.data.frame(matrix(nrow = ncol(y.in), ncol = 4))
  colnames(results.out) <- c("trait", "model", "r.mean", "r.sd")
  
  # Iterate through each trait and perform cross-validation
  for (t in 2:length(colnames(y.in))) {
    trait <- colnames(y.in)[t]
    y.datexist <- y.trainset[which(!is.na(y.trainset[,t])),]
    results.out[t,] <- implement.xval(g.in = g.in, y.in = y.in, trait = trait, reps = reps, k.fold = k.fold, model = model)
  }
  
  # Return a list of results
  cross.val.results <- list(xval.result = results.out, folds = k.fold, reps = reps, model = model)
  return(cross.val.results)
} # Close the function

# Define the function for k-fold cross validation (i.e. 10-fold)
k.xval.BayesLASSO <- function(g.in = NULL, # A n x m genotype matrix (coded as -1, 0, 1), with the row.names of g.in equal to the line names
                              y.in = NULL, # A n x t matrix of phentypes, with the row.names of y.in equal to the line names
                              k.fold = 10, # The number of "folds" to perfom cross-validation with
                              reps = 25, # The number of iterations to divide the lines into k "folds"
                              model = "BayesLASSO",
                              y.trainset = NULL, # identical columns to y.in, but consisting of only the prediction data
                              niter = niter,
                              nburn = nburn
) {
  
  # Error reporting
  if(is.null(g.in)) { stop("The genotype matrix was not specified.") }
  if(is.null(y.in)) { stop("The phenotype data was not specified.") }
  if(nrow(g.in) != nrow(y.in)) { stop("The number of lines in the genotype matrix and in the phenotype data do not match.") }
  
  # Create a function for using rrBLUP
  cv.rrblup <- function(y.train = y.train) {
    solve.out <- eval(parse(text = paste("bayes(",trait,'~1, data=y.train, M=g.in, M.id=row.names(g.in), 
                       method="BayesL", niter = ', niter,', nburn=',nburn,', thin=5,threads=1, verbose = F)')))
    return(solve.out$g)
  }
  
  # Create a function that detects the model and trait and implement cross-validation
  implement.xval <- function(g.in, y.in, trait, reps, k.fold, model) {
    
    n.lines <- nrow(y.datexist)
    
    # Notify which trait is being cross-validated
    writeLines(paste("\nNow conducting k-fold cross validation on ", trait, " using ", model, ".", sep = ""))
    
    # sapply function of all reps of k.fold cv for a trait
    acc.out <- sapply(X = 1:reps, FUN = function(x) {
      
      # Design a list of folds
      k.it <- split(x = sample(x = 1:n.lines, size = n.lines), f = factor(cut(1:n.lines, breaks = k.fold)))
      # Apply the function of the list of folds
      k.fold.out <- lapply(X = k.it, FUN = function(x) {
        pred <- row.names(y.datexist)[x] # Names of prediction set
        train <- setdiff(row.names(y.datexist), pred) # Names of training set
        #g.train <- g.in[train,] # Set training genos
        #g.pred <- g.in[pred,] # Set prediction genos
        y.train <- y.in
        val_set <- setdiff(row.names(g.in), train)
        y.train[val_set,trait] <- NA
        
        # Marker effects
        BLUP <- cv.rrblup(y.train = y.train)
        row.names(BLUP) <- BLUP$id
        BLUP <- BLUP[val_set,]
        #BLUP <- as.matrix(BLUP)
        return(BLUP) # correlate GEBVs to actual phenos
      })
      
      # Create a vector of GEBVs
      GEBV <- do.call("rbind", k.fold.out)
      
      # Sort it by lines
      #browser()
      y.pred <- y.in[GEBV$id,trait]
      #GEBV <- GEBV[order(row.names(GEBV))]
      # Calculate correlation
      acc <- cor(GEBV$gebv, y.pred, use = "complete.obs")
      return(acc)
    })
    
    # Return the mean and sd of the correlations
    r.mean <- mean(acc.out)
    r.sd <- sd(acc.out)
    return(cbind(trait, model, r.mean, r.sd))
    
  } # Close the function
  
  # Create the output matrix
  results.out <- as.data.frame(matrix(nrow = ncol(y.in), ncol = 4))
  colnames(results.out) <- c("trait", "model", "r.mean", "r.sd")
  
  # Iterate through each trait and perform cross-validation
  for (t in 2:length(colnames(y.in))) {
    trait <- colnames(y.in)[t]
    y.datexist <- y.trainset[which(!is.na(y.trainset[,t])),]
    results.out[t,] <- implement.xval(g.in = g.in, y.in = y.in, trait = trait, reps = reps, k.fold = k.fold, model = model)
  }
  
  # Return a list of results
  cross.val.results <- list(xval.result = results.out, folds = k.fold, reps = reps, model = model)
  return(cross.val.results)
} # Close the function


