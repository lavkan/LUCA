## ace.R (2008-03-10)

##     Ancestral Character Estimation

## Copyright 2005-2008 Hua Li

## May 2009 - exponential of a matrix is calculated directly

## Aug 2010 - constraints on rates

library("msm")

aceOptim<-function(x, phy, ip= 0, model="ER"){

    obj <- list()
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
  if (!is.factor(x)) x <- factor(x)

  nl <- nlevels(x)
  lvls <- levels(x)
  x <- as.integer(x)

  if (is.character(model)) {
      rate <- matrix(NA, nl, nl)
      if (model == "ER") np <- rate[] <- 1
      if (model == "ARD") {
          np <- nl*(nl - 1)
          rate[col(rate) != row(rate)] <- 1:np
       }
       if (model == "SYM") {
          np <- nl * (nl - 1)/2
          rate[col(rate) < row(rate)] <- 1:np
          rate <- t(rate)
          rate[col(rate) < row(rate)] <- 1:np
            }
  }  else {
            if (ncol(model) != nrow(model))
              stop("the matrix given as `model' is not square")
            if (ncol(model) != nl)
              stop("the matrix `model' must have as many rows as the number of categories in `x'")
            rate <- model
            np <- max(rate)
    }

    index.matrix <- rate
    index.matrix[cbind(1:nl, 1:nl)] <- NA
    rate[cbind(1:nl, 1:nl)] <- 0
    rate[rate == 0] <- np + 1 # to avoid 0's since we will use this an numeric indexing

    liks <- matrix(0, nb.tip + nb.node, nl)
    for (i in 1:nb.tip) liks[i, x[i]] <- 1

    phy <- reorder(phy, "pruningwise")

     Q <- matrix(0, nl, nl)
     dev <- function(p, output.liks = FALSE) {
            Q[] <- c(p, 0)[rate]
            diag(Q) <- -rowSums(Q)
            for (i  in seq(from = 1, by = 2, length.out = nb.node)) {
                j <- i + 1
                anc <- phy$edge[i, 1]
                des1 <- phy$edge[i, 2]
                des2 <- phy$edge[j, 2]
                P1 <- MatrixExp(Q * phy$edge.length[i])
                P2 <- MatrixExp(Q * phy$edge.length[j])
                liks[anc, ] <- P1 %*% liks[des1, ] * P2 %*% liks[des2, ]
            }
            if (output.liks) return(liks[-(1:nb.tip), ])
            - 2 * log(sum(liks[nb.tip + 1, ]))
        }
        #out <- nlm(function(p) dev(p), p = rep(ip, length.out = np),
        #           hessian = TRUE)
        if (model == "ER") 
        {
         out<-optim( par = rep(ip, length.out = np),  function(p) dev(p), method="L-BFGS-B", lower = 0.000001)
        }
        else
        {
        out <- constrOptim(theta = rep(ip, length.out = np),  function(p) dev(p), grad = NULL, ui = rbind(c(10,-1),c(-1,20),c(1,0),c(0,1)), ci = c(0,0,0,0))
        }
        obj$loglik <- -out$value / 2
        obj$rates <- out$par
        obj$index.matrix <- index.matrix

            lik.anc <- dev(obj$rates, TRUE)
            lik.anc <- lik.anc / rowSums(lik.anc)
            colnames(lik.anc) <- lvls
            obj$lik.anc <- lik.anc
          
        obj$root.prob <- as.numeric(obj$lik.anc[1,2])        
        obj$l_g.ratio <- obj$rates[1]/obj$rates[2]
         obj
}

