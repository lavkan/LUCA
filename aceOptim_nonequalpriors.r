## ace.R (2008-03-10)

##     Ancestral Character Estimation

## Copyright 2005-2008 Hua Li

## May 2009 - exponential of a matrix is calculated directly

## Aug 2010 - constraints on rates

library("msm")

aceOptim_nonUniformPrior<-function(x, phy, ip= 0, model="ER"){

    obj <- list()
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
  if (!is.factor(x)) x <- factor(x)

  nl <- nlevels(x)
  lvls <- levels(x)
  x <- as.integer(x)

  rate <- model
  np <- max(rate)


    index.matrix <- rate
    index.matrix[cbind(1:nl, 1:nl)] <- NA
    rate[cbind(1:nl, 1:nl)] <- 0
    rate[rate == 0] <- np + 1 # to avoid 0's since we will use this an numeric indexing

    liks <- matrix(0, nb.tip + nb.node, nl)
    for (i in 1:nb.tip) liks[i, x[i]] <- 1   ##

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
            if (output.liks) return(liks[-(1:nb.tip), ]) ##### return the whole tree likelihood
            if (dim(Q)[1] == 2){
                L0<- - 2 * log(sum(liks[nb.tip + 1, ]))
            }
            if (dim(Q)[1] == 3){
                prior<-c(0.5, 0.45, 0.05)
                L0<- - 2 * log(sum(liks[nb.tip + 1, ] * prior ))
            }
            if (dim(Q)[1] == 4){
                prior<-c(0.5, 0.45, 0.045, 0.005)
                L0<- - 2 * log(sum(liks[nb.tip + 1, ] * prior ))
            }
            L0
        }
        #out <- nlm(function(p) dev(p), p = rep(ip, length.out = np),
        #           hessian = TRUE)
        dQ = dim(Q)[1]
        out <- constrOptim(theta = rep(ip, length.out = np),  function(p) dev(p), grad = NULL, ui = rbind(diag(np), c(rep(c(-1,10),dQ-1),rep(0,np-2*(dQ-1))), c(rep(c(20,-1),dQ-1),rep(0,np-2*(dQ-1)))), ci = rep(0,length.out = np+2))
        obj$loglik <- -out$value / 2
        obj$rates <- out$par
        obj$index.matrix <- index.matrix

            lik.anc <- dev(obj$rates, TRUE)
            lik.anc <- lik.anc / rowSums(lik.anc)
            colnames(lik.anc) <- lvls
            obj$lik.anc <- lik.anc
        
        if (dQ == 2) {obj$root.prob <- as.numeric(obj$lik.anc[1,2])
                             obj$l_g.ratio <- obj$rates[2]/obj$rates[1]
                            } 
        else {obj$root.prob <- (0.45*obj$lik.anc[1,2]+0.05*obj$lik.anc[1,3])/(0.5*obj$lik.anc[1,1] + 0.45*obj$lik.anc[1,2] + 0.05*obj$lik.anc[1,3]) 
              obj$l_g.ratio <-  (obj$rates[2]+obj$rates[4])/(obj$rates[1]+obj$rates[3])
             }  
         obj
}
