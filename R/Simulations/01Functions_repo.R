MCmad <- function(vec){
  median(abs(vec), na.rm=TRUE)
}

MCsd <- function(vec){
  sd(vec, na.rm=TRUE)
}

once <- function(h,P){
  s <<- h
  print(paste("Iteration:", s, sep=" "))
  nmed <- 1
  ncols <- 2*P + nmed*(P+2)
  # Oracle, naive, nr methods x (P + 2) 
  MC <- matrix(NA, nrow=m, ncol=ncols)
  
  # (Number of different estimators x Number of regressors P) + 2
  
  alpha <- c(a_all[1:s], rep(0, L-s)) 
  # Setup ----
  
  for(i in 1:m){
    print(i)
    eps <- rnorm(n,0,epsvar)
    yps <- matrix(rnorm(P*n,0,1),nrow=n, ncol=P) + 0.5*eps
    # Error terms
    
    D <- Z%*%gammam + yps 
    # First stage
    
    beta <- rep(0,P)
    beta <- matrix(beta)
    # True beta
    
    Y <- D %*% beta + Z %*% alpha + eps
    # Structural equation
    
    cind <- 1
    sind <- ncols - 2*nmed + 1  # For selection statistics later
    MC[i, cind:(cind+P-1)] <- ivreg(Y ~ D | Z)$coefficients[2:(P+1)] # Naive ----
    cind <- cind + P
    # Naive 2SLS
    
    Z.end <- Z[, 1:s]
    # Endogenous Z matrix
    MC[i, cind:(cind+P-1)] <- ivreg(Y ~ D + Z.end| Z)$coefficients[2:(P+1)]
    cind <- cind + P
    # Oracle
    
    # Vecomed
    bvec <- matrix(NA, nrow=Len, ncol=P)
    
    for(j in 1:Len){
      u <- comb[j,]
      Zt.end <- Z[, -u]
      bvec[j,] <- ivreg(Y ~ D + Zt.end| Z)$coefficients[2:(P+1)]
    }
    if(dim(bvec)[2]==1){
      plot(bvec)
    } else {
      plot(bvec[,1],bvec[,2]) 
    }
    
    # Vmmada ----
    bmed <- apply(bvec, FUN=median, MARGIN=2)
    
    (alpham <-solve(t(Z)%*%Z)%*%t(Z)%*%(Y - D%*%bmed))
    
    QR <- qr(Z)
    Dhat <- qr.fitted(QR, D)
    QR2 <- qr(Dhat)
    Ztilde <- qr.resid(QR2, Z)
    rm(QR, QR2, Dhat)
    
    w3 <- 1/abs(alpham)
    # Weights for adaLasso
    
    # Adaptive Lasso
    out <- glmnet(Ztilde, Y, family='gaussian', alpha=1, standardize=TRUE, intercept=FALSE, penalty.factor = w3)
    alpha.ad <- out$beta
    alpha.ad <- drop(alpha.ad)
    alpha.ad[alpha.ad != 0] <- 1
    alpha.ad <- unique(as.matrix(alpha.ad), MARGIN=2)
    rm(Ztilde)
    
    l <- 1 # < Start at first adaLasso step
    pv <- 0 # < Initial p-value
    broken <- F
    while(pv < 0.1/log(n)){ # While p-value smaller prespecified significance level
      if(l > dim(alpha.ad)[2]){# Tries to exceed last adaLasso step
        print("Tries to exceed last adaLasso step.")
        MC[i, cind:(cind+P-1)] <- NA
        cind <- cind + P
        broken <- T
        break
      }
      
      if(L-sum(as.numeric(as.logical(alpha.ad[,l]))) <= P){# Only P IVs chosen as valid
        print("Reached just-identified model. Return the just-identified model.")
        Ze <- matrix(Z[,which(alpha.ad[,l]!=0)], nrow=n)
        Zp <- Z[,which(drop(alpha.ad[,l])==0)]
        # Valid and invalid instruments
        # Post-adaLasso 2SLS
        MC[i, cind:(cind+P-1)] <- ivreg(Y ~ D | Z)$coeff[2:(P+1)]
        cind <- cind + P
        broken <- T
        break
      }
      
      Ze <- matrix(Z[,which(alpha.ad[,l]!=0)], nrow=n) # IVs chosen as invali
      X <- cbind(1, Ze, D)
      
      # Sargan-test >
      if(sum(alpha.ad[,l])==0){ # all valid
        tempivr <- ivreg(Y ~ D | Z)
      } else {
        tempivr <- ivreg(Y ~ D + Ze | Z)
      }
      
      tempivr<- resid(tempivr)
      Zu <- Zc*tempivr
      lambda <- qr.solve(t(Zu)%*%Zu)
      (bgmm <- qr.solve((t(X)%*%Zc%*%lambda%*%t(Zc)%*%X))%*%t(X)%*%Zc%*%lambda%*%t(Zc)%*%Y)
      rgmm <- Y - X%*%bgmm
      (J <- t(t(Zc)%*%rgmm)%*%lambda%*%(t(Zc)%*%rgmm)) 
      
      (pv <- 1 - pchisq(J, L-ncol(Ze)-P))
      # < pi-value HS
      l <- l + 1
    } # end while
    # < Hansen-Andrews stopping algorithm
    
    if(broken == F){
      # l-1 because index was set + 1 inside while-loop >
      Ze <- matrix(Z[,which(alpha.ad[,l-1]!=0)], nrow=n)
      Zp <- Z[,which(drop(alpha.ad[,l-1])==0)]
      # < Valid and invalid instruments, corrected shift-share IV
      
      # Post-adaLasso 2SLS >
      if(sum(alpha.ad[,l-1])==0){
        MC[i, cind:(cind+P-1)] <- ivreg(Y ~ D | Z)$coeff[2:(P+1)]
      } else {
        MC[i, cind:(cind+P-1)] <- ivreg(Y ~ D + Ze | Z)$coeff[2:(P+1)]
      }
      cind <- cind + P
    }
    # include nInv
    MC[i,sind]  <- length(which(alpha.ad[,l-1]!=0))
    MC[i,sind+1]  <- prod(as.numeric(1:s %in% which(alpha.ad[,l-1] != 0)))
    sind <- sind + 2
    # This ifelse-statement makes sure that when no invalid IVs chosen, no empty Ze is chosen
    
    rm(alpha.ad, Ze, Zp)
    
  } #end for
  
  out <- c(apply(MC[,1:((nmed + 2)*P)], FUN=MCmad, MARGIN=2), apply(MC[,1:((nmed + 2)*P)], FUN=MCsd, MARGIN=2), apply(MC[,-(1:((nmed + 2)*P))], FUN=mean, MARGIN=2, na.rm=TRUE))
  return(out)
} # end function

simulone <- function(M,s, L, alpha, gammaval, title, ylimup, breaks, yticklab, xlab, ylimF, breaksF, yticklabF){
  gamma <- rep(gammaval, L)
  
  standard <- matrix(NA, nrow=15, ncol=4)
  oracle <- matrix(NA, nrow=15, ncol=5)
  ada <- matrix(NA, nrow=15, ncol=4)
  cim <- matrix(NA, nrow=15, ncol=4)
  
  mats <- list(standard, oracle, ada, cim)
  
  nseq <- seq(400, 6000, by=400)
  for(i in 1:15){
    print(paste("MC iteration: ", i, sep=""))
    n <- nseq[i]
    set.seed(42)
    
    registerDoRNG(142, once=FALSE)
    Mc <- foreach(m = 1:M, 
                  .combine = rbind, .packages = c("MASS","dplyr", "AER", "glmnet", "gmm", "ddpcr", "ivpack"))  %dopar% {
                    tryCatch({
                      #Create IVs
                      Z <- matrix(runif(n*L,0, 1/L), ncol=L)
                      Z <- scale(Z,center=TRUE,scale=TRUE)
                      colnames(Z) <- c("Z1", "Z2", "Z3", "Z4", "Z5", "Z6","Z7", "Z8", "Z9", "Z10")
                      # Matrix of instruments
                      tF <- (t(gamma)%*%t(Z)%*%Z%*%gamma)/(10*100)
                      # Concentration parameter
                      # 1. Create data
                      errormat <- mvrnorm(n, mu=c(0,0), Sigma=cbind(c(1,0.5), c(0.5,10)))
                      # Structural error:
                      eps <- errormat[,1]
                      # First-stage error:
                      ny <- errormat[,2]
                      
                      D <- Z%*%gamma + ny
                      # Create first-stage correlation
                      Y <- D*0 + Z%*%alpha + eps
                      # Create dependent variable
                      d <- cbind(Y,D, Z, eps, ny)
                      colnames(d)[1:2] <- c("Y","D")
                      d <- data.frame(d)
                      # All in one matrix
                      
                      Y <- Y-mean(Y)
                      D <- D-mean(D)
                      
                      # -1. naive 2SLS
#                      (nMc <- ivreg(d$Y ~ d$D - 1| Z)$coeff[1])
                      (nMc <- ivreg(d$Y ~ d$D - 1| Z -1))
                      nMc_summary <- summary(nMc)
                      nMc <- nMc_summary$coeff[1,1]
                      nMc_se <- nMc_summary$coeff[1,2]
                      nMc_cover <- as.numeric(dplyr::between(0, nMc - 1.96*nMc_se, nMc + 1.96*nMc_se))
                      
                      # 0. FD- oracle 2SLS
                      Z.end <- Z[, 1:s]
                      Zp <- Z[, (s+1):L]
                      # Endogenous Z matrix
                      orMc_TSLS <- ivreg(d$Y ~ d$D + Z.end -1| Z -1)
                      orMc_summary <- summary(orMc_TSLS)
                      orMc <- orMc_summary$coeff[1,1]
                      orMc_se <- orMc_summary$coeff[1,2]
                      orMc_cover <- as.numeric(dplyr::between(0, orMc - 1.96*orMc_se, orMc + 1.96*orMc_se))

                      rm(Z.end, Zp)
                      
                      # 2. Naive paL
                      bvec <- rep(NA, L)
                      sevec <- rep(NA,L)
                      # Empty vectors for coefficients and standard errors
                      
                      for(ind in 1:L){# this 
                        (ivres <- ivreg(d$Y ~ d$D + Z[,-ind] - 1| Z))
                        bvec[ind] <- ivres$coefficients[1]
                        quiet(sevec[ind] <- robust.se(ivres)[1,2])
                      }
                      # Vector of just-identified estimates
                      
                      (bvec <- unlist(bvec))
                      (bmed <- median(bvec))
                      # Initial median estimate
                      (alpham <-solve(t(Z)%*%Z)%*%t(Z)%*%(Y - D%*%bmed))
                      # Plug-in estimate of alpha
                      
                      QR <- qr(Z)
                      Dhat <- qr.fitted(QR, D)
                      QR2 <- qr(Dhat)
                      Ztilde <- qr.resid(QR2, Z)
                      rm(QR, QR2, Dhat)
                      
                      # Panel
                      w3 <- 1/abs(alpham)
                      # Weights for adaLasso
                      print("ada")
                      # Adaptive Lasso
                      out <- glmnet(Ztilde, d$Y, family='gaussian', alpha=1, standardize=TRUE, intercept=FALSE, penalty.factor = w3)
                      alpha.ad <- out$beta
                      alpha.ad <- drop(alpha.ad)
                      alpha.ad[alpha.ad != 0] <- 1
                      alpha.ad <- unique(as.matrix(alpha.ad), MARGIN=2)
                      rm(Ztilde)
                      
                      l <- 1 # < Start at first adaLasso step
                      pi <- 0 # < Initial p-value
                      while(pi < 0.1/log(n)){ # While p-value smaller prespecified significance level
                        print("Start while.")
                        if(l > dim(alpha.ad)[2]){# Tries to exceed last adaLasso step
                          print("Tries to exceed last adaLasso step.")
                          pMc <- NA
                          pnInv <- NA
                          pfreq <- NA
                          break
                        }
                        
                        if(L-sum(as.numeric(as.logical(alpha.ad[,l]))) == 1){# Only one IV chosen as valid
                          print("Number of overidentifying restrictions")
                          print(L-sum(as.numeric(as.logical(alpha.ad[,l]))) - 1)
                          print("Reached just-identified model. Return the just-identified model.")
                          Ze <- matrix(Z[,which(alpha.ad[,l]!=0)], nrow=n)
                          Zp <- Z[,which(drop(alpha.ad[,l])==0)]
                          # Valid and invalid instruments
                          # Post-adaLasso 2SLS
                          pMc <- ivreg(d$Y ~ d$D + Ze -1 | Z)$coeff[2]
                          pnInv <- length(which(drop(alpha.ad[,l])!=0)) # nr of IVs chosen as invalid
                          pfreq <- prod(as.numeric(1:s %in% which(drop(alpha.ad[,l]) != 0))) #all invalid chosen as invalid
                          break
                        }
                        
                        Ze <- matrix(Z[,which(alpha.ad[,l]!=0)], nrow=n) # IVs chosen as invalid
                        
                        X <- cbind(Ze, d$D)
                        
                        # Sargan-test >
                        if(sum(alpha.ad[,l])==0){
                          iv1 <- ivreg(Y ~ D - 1 | Z)
                        } else {
                          iv1 <- ivreg(Y ~ D + Ze - 1 | Z)
                        }
                        res <- resid(iv1)
                        (J <- (t(res)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%res)/(t(res)%*%res/n))
                        (pi <- 1 - pchisq(J, L-ncol(Ze)-1))
                        # < pi-value HS
                        l <- l + 1
                      } # end while
                      # < Hansen-Andrews stopping algorithm
                      # l-1 because index was set + 1 inside while-loop >
                      Ze <- matrix(Z[,which(alpha.ad[,l-1]!=0)], nrow=n)
                      Zp <- Z[,which(drop(alpha.ad[,l-1])==0)]
                      # < Valid and invalid instruments, corrected shift-share IV
                      
                      # Post-adaLasso 2SLS >
                      if(sum(alpha.ad[,l-1])==0){
                        (pMc <- ivreg(d$Y ~ d$D -1 | Z-1))
                      } else {
                        (pMc <- ivreg(d$Y ~ d$D + Ze -1 | Z-1))
                      }
                      pMc_summary <- summary(pMc)
                      pMc <- pMc_summary$coeff[1,1]
                      pMc_se <- pMc_summary$coeff[1,2]
                      pMc_cover <- as.numeric(dplyr::between(0, pMc - 1.96*pMc_se, pMc + 1.96*pMc_se))
                      # This ifelse-statement makes sure that when no invalid IVs chosen, no empty Ze is chosen
                      pnInv <- length(which(drop(alpha.ad[,l-1])!=0))
                      pfreq <- prod(as.numeric(1:s %in% which(drop(alpha.ad[,l-1]) != 0)))
                      rm(alpha.ad, Ze, Zp)
                      
                      # 3. naive CIM
                      print("naive CIM")
                      Ji <- Inf
                      crit <- 3*sqrt(2.01^2*log(n)) 
                      wvold <- matrix(2, ncol=L)
                      dfZ <- ncol(Z)
                      
                      while(Ji > qchisq(p=1-0.1/log(n), df=dfZ)){
                        votematrix <- matrix(rep(0, L*L), ncol=L)
                        CImat <- matrix(rep(NA,L*3), ncol=3)
                        CImat[,1] <- as.matrix(bvec) - crit*sevec
                        CImat[,2] <- as.matrix(bvec) + crit*sevec
                        CImat[,3] <- 1:L
                        CImat <- data.frame(CImat)
                        CImat <- arrange(CImat, CImat[,1])
                        
                        for(indV in 1:L){
                          for(oth in 1:indV){
                            if(CImat[indV,1] < CImat[oth,2]){
                              votematrix[indV,oth] <- 1
                            } else {
                              votematrix[indV,oth] <- 0 
                            }
                          }
                        }
                        # Fill in votematrix
                        
                        wv <- votematrix[which(rowSums(votematrix) == max(rowSums(votematrix))),]
                        wv <- as.matrix(wv)
                        # wv is row with maximal votes
                        
                        if(length(wv)==L){ # only one plurality group
                          wv <- t(wv)}
                        
                        #if(dim(wvold)[1] ==  dim(wv)[1] && as.logical(prod(wvold != as.vector(wv)))){
                        if(as.logical(prod(wvold==2)) | (dim(wvold)[1] != dim(wv)[1]) || (dim(wvold)[1] == dim(wv)[1] &&  as.logical(prod(wvold!=wv)))){  
                          # wvold is last wv, wv are votematrix row with maximal nr chosen as valid
                          # if wvold is initial one (=2)
                          # OR nr of rows of wvold and wv is unequal
                          # OR (nr of rows is equal AND the entries are not equal)
                          if(sum(wv)==dim(wv)[1]){
                            print("Tried to run HS-test with only one IV. Fall back to last overidentified model.")
                            break 
                          }
                          
                          # HS-test for plurality group, also ties are taken care of
                          Jvals <- rep(NA, dim(wv)[1])
                          # Empty vector for p-values
                          for(indM in 1:dim(wv)[1]){
                            Ze <- Z[,CImat[,3][as.logical(wv[indM,])==FALSE]]
                            Ze <- as.matrix(Ze, nrow=n)
                            X <- cbind(Ze, d$D)
                            
                            # Sargan-test >
                            if(sum(as.logical(wv)==FALSE)==0){
                              iv1 <- ivreg(Y ~ D - 1 | Z-1)
                            } else {
                              iv1 <- ivreg(Y ~ D + Ze - 1 | Z-1)
                            }
                            res <- resid(iv1)
                            (J <- (t(res)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%res)/(t(res)%*%res/n))
                            #  < J-value HS
                            Jvals[indM] <- J
                          } # < end for, wv selection
                          Ji <- min(Jvals)
                          wve <- sort(CImat[,3][as.logical(wv[which(Jvals==Ji),]) == TRUE])
                          wi <- sort(CImat[,3][as.logical(wv[which(Jvals==Ji),]) == FALSE])
                          wvold <- wv
                          dfZ <- L - ncol(Ze) - 1
                          rm(Ze)
                          #print(Jvals)
                          print(length(wi))
                          print(crit)
                        }
                        crit <- 0.99*crit
                      } # < end while
                      
                      Ze <- Z[,wi]
                      Ze <- as.matrix(Ze, nrow=n)
                      Zp <- Z[,wve]
                      Zp <- as.matrix(Zp, nrow=n)
                      if(dim(Ze)[2] == 0){
                        naiveciMc <- ivreg(d$Y ~ d$D -1 | Z -1)
                      } else {
                        naiveciMc <- ivreg(d$Y ~ d$D + Ze -1 | Z-1)
                      }
                      ciMc_summary <- summary(naiveciMc)
                      naiveciMc <- ciMc_summary$coeff[1,1]
                      ciMc_se <- ciMc_summary$coeff[1,2]
                      ciMc_cover <- as.numeric(dplyr::between(0, naiveciMc - 1.96*ciMc_se, naiveciMc + 1.96*ciMc_se))
                      
                      (naivecinInv <- L-length(wve))
                      (naivecifreq <- prod(as.numeric(1:s %in% wi))) # it is important that this is wve! not wv, check also for panel case
                      
                      rm(wi,wve,Zp)
                      list(nMc, nMc_cover, orMc, orMc_cover, tF, pMc, pMc_cover, pnInv, pfreq, naiveciMc, ciMc_cover, naivecinInv, naivecifreq)
                      #####1####2##########3#####4###########5###6####7##########8######9######10#########11##########12###########13#####
                    },  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})   
                  }
    #    oneN(m)
    (Mc <- matrix(unlist(Mc), ncol=13))
    
    # Standard
    (mats[[1]][i,1] <- median(abs(Mc[,1])))
    (mats[[1]][i,2] <- mean(Mc[,2]))
    (mats[[1]][i,3] <- 0)
    (mats[[1]][i,4] <- 0)
    
    # Oracle
    (mats[[2]][i,1] <- median(abs(Mc[,3])))
    (mats[[2]][i,2] <- mean(Mc[,4]))
    (mats[[2]][i,3] <- s)
    (mats[[2]][i,4] <- 1)
    (mats[[2]][i,5] <- mean(Mc[,5]))
    
    # adaLasso
    (mats[[3]][i,1] <- median(abs(Mc[,6]), na.rm=T))
    (mats[[3]][i,2] <- mean(Mc[,7], na.rm=T))
    (mats[[3]][i,3] <- mean(Mc[,8], na.rm=T))
    (mats[[3]][i,4] <- mean(Mc[,9], na.rm=T))
    
    # CIM
    (mats[[4]][i,1] <- median(abs(Mc[,10]), na.rm=T))
    (mats[[4]][i,2] <- mean(Mc[,11], na.rm=T))
    (mats[[4]][i,3] <- mean(Mc[,12], na.rm=T))
    (mats[[4]][i,4] <- mean(Mc[,13], na.rm=T))
  }
  
  standard <- data.frame(mats[[1]], seq(400,6000, by=400), "Standard")
  colnames(standard) <- c("mad", "cover", "nInv", "freq", "ev", "Method") 
  oracle <- data.frame(mats[[2]][,1:4], seq(400,6000, by=400),"Oracle")
  colnames(oracle) <- c("mad", "cover", "nInv", "freq", "ev", "Method") 
  ada <- data.frame(mats[[3]], seq(400,6000, by=400),"AL")
  colnames(ada) <- c("mad", "cover","nInv", "freq", "ev", "Method") 
  cim <- data.frame(mats[[4]],seq(400,6000, by=400),"CIM")
  colnames(cim) <- c("mad", "cover","nInv", "freq", "ev", "Method") 
  
  mnoac <- rbind.data.frame(standard, oracle, ada, cim)
  mnoac_F <- data.frame(mats[[2]][,5], seq(400,6000, by=400),"oracle")
  
  res <- list(mnoac=mnoac,mnoac_F=mnoac_F)
  return(res)
}

# Functions ----
# dataprep
# and.rub
# adamHS
# adamAR

dataprep <- function(wFB){
  # Variables ----
  # IVs
  Z <- analysis_data[,grep(wFB, colnames(analysis_data))]
  Z <- as.matrix(Z)
  n <- dim(Z)[1]
  
  # Endogenous variables
  D <- cbind(d_immi_std, ld_immi_std)
  D1 <- matrix(D[,1])
  D2 <- matrix(D[,2])
  Y1 <- dlweekly
  Y2 <- dlweekly_hskill
  Y3 <- dlweekly_lskill
  
  # Weighted variables
  wY1 <- Y1 * sqrt(lpop)
  wY2 <- Y2*sqrt(lpop)
  wY3 <- Y3*sqrt(lpop)
  wD1 <- D1 * sqrt(lpop) 
  wD2 <- D2 * sqrt(lpop) 
  wD <- sweep(D, MARGIN=1, sqrt(lpop), `*`)
  wZ <- sweep(Z, MARGIN=1, sqrt(lpop), `*`)
  wyear3 <- year3 * sqrt(lpop)
  wyear4 <- year4 * sqrt(lpop) 
  
  #FWL
  wX <- cbind(wyear3, wyear4, sqrt(lpop))
  
  Y1 <- qr.resid(qr(wX), wY1)
  Y2 <- qr.resid(qr(wX), wY2)
  Y3 <- qr.resid(qr(wX), wY3)
  D <- qr.resid(qr(wX), wD)
  D1 <- matrix(D[,1])
  D2 <- matrix(D[,2])
  Z <- qr.resid(qr(wX), wZ)
  
  X <- cbind(1, year3, year4)
  
  NAMES <- c("Z", "n", "D", "D1", "D2", "Y1", "Y2", "Y3", "X",
             "wZ", "wD", "wD1", "wD2", "wY1", "wY2", "wY3", "wX", "wyear3", "wyear4")
  
  lapply(seq_along(NAMES), 
         function(x) {
           assign(NAMES[x], get(NAMES[x]), envir=.GlobalEnv)
         }
  )
}

fnMatSqrtInverse = function(mA) {
  ei = eigen(mA)
  d = ei$values
  d = (d+abs(d))/2
  d2 = 1/sqrt(d)
  d2[d == 0] = 0
  return(ei$vectors %*% diag(d2) %*% t(ei$vectors))
}

# Anderson-Rubin test, allows controls
and.rub <- function(Y, D, Z, X){
  vY <- Y # Dependent variable
  mY2 <- as.matrix(D) # Endogenous regressor
  mZ1 <- X  # Controls
  mZ <- cbind(Z, mZ1[,which(!colnames(mZ1) %in% intersect(colnames(mZ1),colnames(Z)))]) # IVs and controls which are not IVs

  mMZ <- diag(n) - mZ %*% solve(t(mZ) %*% mZ) %*% t(mZ)
  mMZ1 <- diag(n) - mZ1 %*% solve(t(mZ1) %*% mZ1) %*% t(mZ1)
  mYStar <- cbind(vY, mY2)
  
  mYZY <- t(mYStar) %*% mMZ %*% mYStar
  mYZ1Y <- t(mYStar) %*% mMZ1 %*% mYStar
  mLeftRight = fnMatSqrtInverse(mYZY)
  
  mLeftRight%*%mYZY%*%mLeftRight
  
  (kappa.liml = sort(eigen(mLeftRight%*%mYZ1Y%*%mLeftRight, 
                           only.values = TRUE)$values)[1])
  AR <<- n*log(kappa.liml)
  print(AR)
  pAR <<- pchisq(AR, df=ncol(mZ)-ncol(mZ1)-ncol(mY2), lower.tail = F) # df: nr IVs - nr controls - nr regressors
  pAR
}

adamAR <- function(Y,D){
  P <- dim(D)[2]
  # Vecomed
  comb <- combinations(n=length(colnames(Z)), v=colnames(Z), r=P, repeats.allowed=F)
  Len <- choose(length(colnames(Z)), P)
  
  bvec <- matrix(NA, nrow=Len, ncol=P)
  
  for(j in 1:Len){
    u <- comb[j,]
    Zt.end <- Z[, colnames(Z)[which(colnames(Z) %in% u)]]
    if(P==1){
      bvec[j] <- ivreg(Y ~ D + Zt.end -1 | Z -1 )$coefficients[2] 
    } else {
      bvec[j,] <- ivreg(Y ~ D + Zt.end -1 | Z -1 )$coefficients[2:(P+1)]
    }
  }
  plot(bvec)
  
  bmed <- apply(bvec, FUN=median, MARGIN=2)
  bmed
  
  (alpham <-solve(t(Z)%*%Z)%*%t(Z)%*%(Y - D%*%bmed))
  
  QR <- qr(Z)
  Dhat <- qr.fitted(QR, D)
  QR2 <- qr(Dhat)
  Ztilde <- qr.resid(QR2, Z)
  rm(QR, QR2, Dhat)
  
  w3 <- 1/abs(alpham)
  # Weights for adaLasso
  
  # Adaptive Lasso
  out <- glmnet(Ztilde, Y, family='gaussian', nlambda = 1000, alpha=1, standardize=F, intercept=FALSE, penalty.factor = w3)
  alpha.ad <- out$beta
  alpha.ad <- drop(alpha.ad)
  alpha.ad[alpha.ad != 0] <- 1
  alpha.ad <- unique(as.matrix(alpha.ad), MARGIN=2)
  alpha.ad <<- alpha.ad
  print("Number of IVs on path")
  print(colSums(alpha.ad))
  rm(Ztilde)
  
  l <- 1 # < Start at first adaLasso step
  pv <- 0 # < Initial p-value
  while(pv < 0.1/log(n)){ # While p-value smaller prespecified significance level
    wZe <- matrix(wZ[,which(alpha.ad[,l]!=0)], nrow=dim(Z)[1]) # Z chosen as invalid
    colnames(wZe) <- names(which(alpha.ad[,l]!=0))
    wXa <- cbind(wX, wZe)
    
    wY <- get(paste("w",deparse(substitute(Y)), sep=""))
    wDn <- get(paste("w",deparse(substitute(D)), sep=""))
    pv <- and.rub(Y=wY, D=wDn, Z=wZ, X=wXa) #df not correct like this, because actually only need to subtract nr of Ze
    print(pv)
    # < pi-value HS
    l <- l + 1
  } # end while
  
  # < Anderson-Rubin stopping algorithm
  
  a.end <<- alpha.ad[,l-1]
  l <<- l
}

# Adam function, HS-test ---- 
adamHS <- function(Y, D, tau){
  P <- dim(D)[2]
  # Vecomed
  comb <- combinations(n=length(colnames(Z)), v=colnames(Z), P, repeats.allowed=F)
  Len <- choose(length(colnames(Z)), P)
  
  bvec <- matrix(NA, nrow=Len, ncol=P)
  
  for(j in 1:Len){
    u <- comb[j,]
    Zt.end <- Z[, colnames(Z)[which(colnames(Z) %in% u)]]
    bvec[j,] <- ivreg(Y ~ D + Zt.end -1 | Z -1 )$coefficients[2:(P+1)]
  }
  plot(bvec)
  
  bmed <- apply(bvec, FUN=median, MARGIN=2)
  bmed
  #MC[i, cind:(cind+P-1)] <- bmed
  #cind <- cind + P
  
  (alpham <-solve(t(Z)%*%Z)%*%t(Z)%*%(Y - D%*%bmed))
  
  QR <- qr(Z)
  Dhat <- qr.fitted(QR, D)
  QR2 <- qr(Dhat)
  Ztilde <- qr.resid(QR2, Z)
  rm(QR, QR2, Dhat)
  
  w3 <- 1/abs(alpham)
  # Weights for adaLasso
  
  # Adaptive Lasso
  out <- glmnet(Ztilde, Y, family='gaussian', nlambda = 1000, alpha=1, standardize=F, intercept=FALSE, penalty.factor = w3)
  alpha.ad <- out$beta
  alpha.ad <- drop(alpha.ad)
  alpha.ad[alpha.ad != 0] <- 1
  alpha.ad <- unique(as.matrix(alpha.ad), MARGIN=2)
  rm(Ztilde)
  
  l <- 1 # < Start at first adaLasso step
  pv <- 0 # < Initial p-value
  while(pv < tau){ # While p-value smaller prespecified significance level
    
    Ze <- matrix(Z[,which(alpha.ad[,l]!=0)], nrow=dim(Z)[1]) # Z chosen as invali
    X <- cbind(Ze, D)
    
    res_FirstStep <- residuals(AER::ivreg(Y ~ X - 1 | Z - 1))
    
    Weight_SecondStep <- crossprod(res_FirstStep * Z);
    
    Coef_SecondStep <- solve(
      t(X) %*% Z %*% solve(Weight_SecondStep) %*%t(Z) %*% X
    ) %*% t(X) %*% Z %*% solve(Weight_SecondStep) %*% t(Z) %*% Y;
    
    res_SecondStep <- as.vector(Y - X %*% Coef_SecondStep);
    
    sd_SecondStep <- sqrt(diag(solve(
      t(X) %*% Z %*% solve(Weight_SecondStep) %*%t(Z) %*% X
    ) %*% t(X) %*% Z %*% solve(Weight_SecondStep)%*%crossprod(res_SecondStep * Z)%*%t(
      solve(
        t(X) %*% Z %*% solve(Weight_SecondStep) %*%t(Z) %*% X
      ) %*% t(X) %*% Z %*% solve(Weight_SecondStep)
    )));
    
    HansenJ_Stat <- t(res_SecondStep) %*% Z %*% solve(Weight_SecondStep) %*%
      t(Z) %*% res_SecondStep;
    
    pv <- pchisq(HansenJ_Stat, df=ncol(Z) - ncol(X), lower.tail = FALSE)
    # < pi-value HS
    l <- l + 1
  } # end while
  # < Hansen-Andrews stopping algorithm
  print(pv)
  a.end <<- alpha.ad[,l-1]
  l <<- l
}