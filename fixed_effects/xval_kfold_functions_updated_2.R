###########################################################################
######### Cross Validation - K-fold (2,5,10) - rrBLUP with Fixed SNPs #####
###########################################################################
# Define the function for k-fold cross validation (i.e. 10-fold)
# Define the function for k-fold cross validation (i.e. 10-fold)
k.xval <- function(g.in = NULL, # A n x m genotype matrix (coded as -1, 0, 1), with the row.names of g.in equal to the line names
                   y.in = NULL, # A n x t matrix of phenotypes, with the row.names of y.in equal to the line names
                   k.fold = 10, # The number of "folds" to perform cross-validation with
                   reps = 25, # The number of iterations to divide the lines into k "folds"
                   model = "rrBLUP",
                   y.trainset = NULL, # Identical columns to y.in, but consisting of only the prediction data
                   fixed.snps = NULL # Parameter for fixed SNPs
) {
  
  # Error reporting
  if(is.null(g.in)) { stop("The genotype matrix was not specified.") }
  if(is.null(y.in)) { stop("The phenotype data was not specified.") }
  if(nrow(g.in) != nrow(y.in)) { stop("The number of lines in the genotype matrix and in the phenotype data do not match.") }
  
  # Create a function for using rrBLUP
  cv.rrblup <- function(y.train = y.train, g.train = g.train, fixed.train = fixed.train) {
    solve.out <- mixed.solve(y = y.train, Z = g.train, X = fixed.train, SE = F, return.Hinv = F)
    return(solve.out$u)
  }
  
  # Create a function that detects the model and trait and implements cross-validation
  implement.xval <- function(g.in, y.in, y.datexist, trait, reps, k.fold, model, fixed.snps) {
    
    n.lines <- nrow(y.datexist)
    
    # Notify which trait is being cross-validated
    writeLines(paste("\nNow conducting k-fold cross-validation on ", trait, " using ", model, ".", sep = ""))
    
    # sapply function of all reps of k.fold cv for a trait
    acc.out <- sapply(X = 1:reps, FUN = function(x) {
      
      # Design a list of folds
      k.it <- split(x = sample(x = 1:n.lines, size = n.lines), f = factor(cut(1:n.lines, breaks = k.fold)))
      # Apply the function on the list of folds
      k.fold.out <- lapply(X = k.it, FUN = function(x) {
        pred <- row.names(y.datexist)[x] # Names of the prediction set
        train <- setdiff(row.names(y.datexist), pred) # Names of the training set
        g.train <- g.in[train, ] # Set training genotypes
        y.train <- as.vector(y.in[train, trait])
        
        # Extract fixed SNPs for the training set
        if (!is.null(fixed.snps)) {
          fixed.train <- fixed.snps[train, ]
        } else {
          fixed.train <- NULL
        }
        
        # Marker effects
        u.hat <- cv.rrblup(y.train = y.train, g.train = g.train, fixed.train = fixed.train)
        
        # GEBVs and correlation
        val_set <- setdiff(row.names(g.in), train)
        g.pred <- g.in[val_set, ]
        GEBV <- g.pred %*% u.hat
        return(GEBV)
      })
      
      # Create a vector of GEBVs
      GEBV <- do.call("rbind", k.fold.out)
      # Sort it by lines
      y.pred <- y.in[row.names(GEBV), ]
      # Calculate correlation
      acc <- cor(GEBV, y.pred[, trait], use = "complete.obs")
      return(acc)
    })
    
    # Return the mean and standard deviation of the correlations
    r.mean <- mean(acc.out)
    r.sd <- sd(acc.out)
    return(cbind(trait, model, r.mean, r.sd))
    
  } # Close the function
  
  # Create the output matrix
  results.out <- as.data.frame(matrix(nrow = ncol(y.in), ncol = 4))
  colnames(results.out) <- c("trait", "model", "r.mean", "r.sd")
  
  # Ensure the loop only runs for the number of columns in y.trainset
  for (t in 1:min(ncol(y.in), ncol(y.trainset))) {
    trait <- colnames(y.in)[t]
    
    # Filter for non-NA rows in y.trainset for this specific trait
    non_na_indices <- which(!is.na(y.trainset[, t]))
    
    # Only proceed if there are non-NA values for the training set
    if (length(non_na_indices) > 0) {
      y.datexist <- y.trainset[non_na_indices, ]  # Use only non-NA rows for training
      results.out[t, ] <- implement.xval(g.in = g.in, y.in = y.in, y.datexist = y.datexist, 
                                         trait = trait, reps = reps, k.fold = k.fold, 
                                         model = model, fixed.snps = fixed.snps)
    } else {
      print(paste("Skipping trait", trait, "due to all NA values in the training set"))
    }
  }
  
  # Return a list of results
  cross.val.results <- list(xval.result = results.out, folds = k.fold, reps = reps, model = model)
  return(cross.val.results)
} # Close the function

######################################################################################
######### Cross Validation - K-fold (2,5,10) - Exponential kernel
######################################################################################
k.xval.EXP <- function(g.in = NULL, 
                       y.in = NULL, 
                       k.fold = 10, 
                       reps = 25, 
                       snps.fixed = NULL,  # Add snps.fixed argument
                       k_dist = k_dist) {
  
  if(is.null(g.in)) { stop("The genotype matrix was not specified.") }
  if(is.null(y.in)) { stop("The phenotype data was not specified.") }
  if(nrow(g.in) != nrow(y.in)) { stop("The number of lines in the genotype matrix and in the phenotype data do not match.") }
  
  cv.rrblup <- function(y.train, k_dist, snps.fixed = NULL) {
    if (!is.null(snps.fixed)) {
      X_fixed <- g.train[, snps.fixed, drop = FALSE]
      solve.out <- kin.blup(data = y.train, geno = colnames(y.train)[1], pheno = trait, GAUSS = TRUE, K = k_dist, X = X_fixed)
    } else {
      solve.out <- kin.blup(data = y.train, geno = colnames(y.train)[1], pheno = trait, GAUSS = TRUE, K = k_dist)
    }
    return(solve.out$g)
  }
  
  implement.xval <- function(g.in, y.in, trait, reps, k.fold, snps.fixed = NULL) {
    n.lines <- nrow(y.trainset)
    acc.out <- sapply(1:reps, function(rep) {
      folds <- split(sample(1:n.lines), rep(1:k.fold, length.out = n.lines))
      fold.cor <- sapply(folds, function(test_indices) {
        test_names <- row.names(y.trainset)[test_indices]
        train_names <- setdiff(row.names(y.trainset), test_names)
        y.train <- y.in[train_names, ]
        
        # Marker effects with fixed SNPs
        BLUP <- cv.rrblup(y.train = y.train, k_dist = k_dist, snps.fixed = snps.fixed)
        GEBV <- BLUP[setdiff(row.names(g.in), train_names)]
        
        cor(GEBV, y.in[test_names, trait], use = "complete.obs")
      })
      mean(fold.cor)
    })
    return(cbind(trait, model, mean(acc.out), sd(acc.out)))
  }
  
  results.out <- data.frame(trait = character(), model = character(), r.mean = numeric(), r.sd = numeric(), stringsAsFactors = FALSE)
  for (trait in colnames(y.in)) {
    y.train_subset <- y.trainset[!is.na(y.trainset[[trait]]), ]
    cv_result <- implement.xval(g.in, y.in, y.train_subset, trait, reps, k.fold, snps.fixed)
    results.out <- rbind(results.out, cv_result)
  }
  
  return(list(xval.result = results.out, folds = k.fold, reps = reps, model = "Exponential kernel"))
}


#######################################################################################
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
######### Cross Validation - K-fold (2,5,10) - Gaussian kernel
###########################################################################

k.xval.GAUSS <- function(g.in = NULL, 
                         y.in = NULL, 
                         k.fold = 10, 
                         reps = 25, 
                         snps.fixed = NULL,  # Add snps.fixed argument
                         k_dist = k_dist) {
  
  if(is.null(g.in)) { stop("The genotype matrix was not specified.") }
  if(is.null(y.in)) { stop("The phenotype data was not specified.") }
  if(nrow(g.in) != nrow(y.in)) { stop("The number of lines in the genotype matrix and in the phenotype data do not match.") }
  
  cv.rrblup <- function(y.train, k_dist, snps.fixed = NULL) {
    if (!is.null(snps.fixed)) {
      X_fixed <- g.train[, snps.fixed, drop = FALSE]
      solve.out <- kin.blup(data = y.train, geno = colnames(y.train)[1], pheno = trait, GAUSS = TRUE, K = k_dist, X = X_fixed)
    } else {
      solve.out <- kin.blup(data = y.train, geno = colnames(y.train)[1], pheno = trait, GAUSS = TRUE, K = k_dist)
    }
    return(solve.out$g)
  }
  
  implement.xval <- function(g.in, y.in, trait, reps, k.fold, snps.fixed = NULL) {
    n.lines <- nrow(y.trainset)
    acc.out <- sapply(1:reps, function(rep) {
      folds <- split(sample(1:n.lines), rep(1:k.fold, length.out = n.lines))
      fold.cor <- sapply(folds, function(test_indices) {
        test_names <- row.names(y.trainset)[test_indices]
        train_names <- setdiff(row.names(y.trainset), test_names)
        y.train <- y.in[train_names, ]
        
        # Marker effects with fixed SNPs
        BLUP <- cv.rrblup(y.train = y.train, k_dist = k_dist, snps.fixed = snps.fixed)
        GEBV <- BLUP[setdiff(row.names(g.in), train_names)]
        
        cor(GEBV, y.in[test_names, trait], use = "complete.obs")
      })
      mean(fold.cor)
    })
    return(cbind(trait, model, mean(acc.out), sd(acc.out)))
  }
  
  # Return cross-validation results
  results.out <- data.frame(trait = character(), model = character(), r.mean = numeric(), r.sd = numeric(), stringsAsFactors = FALSE)
  for (trait in colnames(y.in)) {
    y.train_subset <- y.trainset[!is.na(y.trainset[[trait]]), ]
    cv_result <- implement.xval(g.in, y.in, y.train_subset, trait, reps, k.fold, snps.fixed)
    results.out <- rbind(results.out, cv_result)
  }
  
  return(list(xval.result = results.out, folds = k.fold, reps = reps, model = "Gaussian kernel"))
}


###########################################################################
######### Cross Validation - K-fold (2,5,10) - Bayes Cpi
###########################################################################
k.xval.BayesCpi <- function(g.in = NULL,  # Genotype matrix (n x m), row.names of g.in should be line names
                            y.in = NULL,  # Phenotype matrix (n x t), row.names should be line names
                            k.fold = 10,  # Number of folds for cross-validation
                            reps = 25,    # Number of repetitions for cross-validation
                            snps.fixed = NULL,  # SNPs to include as fixed effects
                            niter = 3000,  # Number of iterations for BayesCpi
                            nburn = 1200) {  # Burn-in period for BayesCpi
  
  # Error handling
  if(is.null(g.in)) { stop("The genotype matrix was not specified.") }
  if(is.null(y.in)) { stop("The phenotype data was not specified.") }
  if(nrow(g.in) != nrow(y.in)) { stop("The number of lines in the genotype matrix and in the phenotype data do not match.") }
  
  # Define the BayesCpi function with fixed effects
  cv.BayesCpi <- function(y.train, g.train, snps.fixed = NULL) {
    if (!is.null(snps.fixed)) {
      # Extract fixed SNPs
      X_fixed <- g.train[, snps.fixed, drop = FALSE]
      
      # Fit the BayesCpi model with fixed SNPs
      model <- bayes(y = y.train, M = g.train, X = X_fixed, method = "BayesCpi", 
                     niter = niter, nburn = nburn, thin = 5, verbose = FALSE)
    } else {
      # Fit the BayesCpi model without fixed SNPs
      model <- bayes(y = y.train, M = g.train, method = "BayesCpi", 
                     niter = niter, nburn = nburn, thin = 5, verbose = FALSE)
    }
    
    # Return the estimated marker effects
    return(model$u)
  }
  
  # Implement cross-validation for BayesCpi
  implement.xval <- function(g.in, y.in, trait, reps, k.fold, snps.fixed = NULL) {
    n.lines <- nrow(y.in)
    
    # Perform cross-validation repeats
    acc.out <- sapply(1:reps, function(rep) {
      # Create random folds
      folds <- split(sample(1:n.lines), rep(1:k.fold, length.out = n.lines))
      
      fold.cor <- sapply(folds, function(test_indices) {
        test_names <- row.names(y.in)[test_indices]
        train_names <- setdiff(row.names(y.in), test_names)
        
        g.train <- g.in[train_names, , drop = FALSE]
        y.train <- y.in[train_names, trait]
        
        # Fit BayesCpi with or without fixed SNPs
        u.hat <- cv.BayesCpi(y.train = y.train, g.train = g.train, snps.fixed = snps.fixed)
        
        # Predict GEBVs for test set
        g.test <- g.in[test_names, !(colnames(g.in) %in% snps.fixed), drop = FALSE]
        if (!is.null(snps.fixed)) {
          X.test_fixed <- g.in[test_names, snps.fixed, drop = FALSE]
          GEBV <- X.test_fixed %*% u.hat[snps.fixed] + g.test %*% u.hat[!(names(u.hat) %in% snps.fixed)]
        } else {
          GEBV <- g.test %*% u.hat
        }
        
        # Correlate GEBVs with actual phenotypes
        cor(GEBV, y.in[test_names, trait], use = "complete.obs")
      })
      
      # Mean correlation across folds
      mean(fold.cor, na.rm = TRUE)
    })
    
    # Return the mean and SD of correlations
    r.mean <- mean(acc.out)
    r.sd <- sd(acc.out)
    return(cbind(trait, model, r.mean, r.sd))
  }
  
  # Cross-validation for each trait
  results.out <- data.frame(trait = character(), model = character(), r.mean = numeric(), r.sd = numeric(), stringsAsFactors = FALSE)
  
  for (trait in colnames(y.in)) {
    # Subset training data for the trait (remove NAs)
    y.train_subset <- y.in[!is.na(y.in[, trait]), ]
    
    # Perform cross-validation for the trait
    cv_result <- implement.xval(g.in, y.in, trait, reps, k.fold, snps.fixed)
    
    # Append results
    results.out <- rbind(results.out, cv_result)
  }
  
  # Return cross-validation results
  cross.val.results <- list(xval.result = results.out, folds = k.fold, reps = reps, model = "BayesCpi")
  return(cross.val.results)
}




