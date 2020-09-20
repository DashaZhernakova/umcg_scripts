library(ggpmisc)
# idx = 2
# 
# merged_tab <- rm_na_outliers(traits_m, pheno_m, idx, method = "IQR", log_tr = F)
# colnames(merged_tab)[1] <- "phenotype"
# merged_tab <- merged_tab[(merged_tab$age < max_age) & (merged_tab$age >= min_age),]
# merged_tab <- mutate(merged_tab, ord_gender_F1M2 = ordered(gender_F1M2, levels = c('1', '2')))
# merged_tab <- mutate(merged_tab, gender_F1M2 = factor(gender_F1M2))

get_breakpoints_derivatives <- function(merged_tab, cutoff = 0.00015, min_age = 20, max_age = 80, span_val = 7){
  age_array <- seq(min_age, max_age, length = n_points)
  breakpoints <- list()
  
  d <- merged_tab[merged_tab$gender_F1M2 == 1,]
  d$phenotype <- scale(d$phenotype)
  deriv2 <- get_secod_derivative(d)
  peaks <- ggpmisc:::find_peaks(abs(deriv2), span = span_val) & abs(deriv2) > cutoff 
  breakpoints[[1]] <- age_array[peaks]
  
  d <- merged_tab[merged_tab$gender_F1M2 == 2,]
  d$phenotype <- scale(d$phenotype)
  deriv2 <- get_secod_derivative(d)
  peaks <- ggpmisc:::find_peaks(abs(deriv2), span = span_val) & abs(deriv2) > cutoff 
  breakpoints[[2]] <- age_array[peaks]
  return(breakpoints)
  #plot(y = deriv2, x = age_array[2:n_points], type = 'l')
  #cutoff <- 0.0008
  #abline(h=cutoff)
  #abline(h=-1*cutoff)
}
  

get_secod_derivative <- function(d, eps = 1e-7, n_points = 300){

  mod <- gam(phenotype ~  ba + eo + er + gr + 
               ly + mo + tr + s(age), data = d, method = "REML")
  new.x <- with(d, expand.grid(age = seq(min_age, max_age, length = n_points), 
                               ba = mean(ba), eo = mean(eo), er = mean(er), gr = mean(gr), 
                               ly = mean(ly),  mo = mean(mo), tr = mean(tr)))
  new.y <- data.frame(predict(mod, newdata = new.x, se.fit = TRUE, type = "response"))
  pdat <- data.frame(new.x, new.y)
  pdat <- rename(pdat, pred = fit, SE = se.fit)
  pdat <- mutate(pdat, lwr = pred - 1.96 * SE, upr = pred + 1.96 * SE) # calculating the 95% confidence interval
  
  
  newDF <- data.frame(pdat)
  X0 <- predict(mod, newDF, type = "lpmatrix")
  newDFeps_p <- newDF + eps
  X1 <- predict(mod, newDFeps_p, type = "lpmatrix")
  Xp <- (X1 - X0) / eps
  deriv1 <- Xp %*% coef(mod) #first derivative
  age_array <- seq(min_age, max_age, length = n_points)
  #plot(y = deriv1, x = age_array, type = 'l')
  deriv2 <- deriv1[2:length(deriv1)] - deriv1[1:length(deriv1)-1] # a very rough second derivative with epsilon = 1
  #plot(y = deriv2, x = age_array[2:n_points], type = 'l')
  deriv1 <- pdat$pred[2:nrow(pdat)] - pdat$pred[1:nrow(pdat)-1]
  deriv2 <- deriv1[2:length(deriv1)] - deriv1[1:length(deriv1)-1] # a very rough second derivative with epsilon = 1
  
  plot(y = deriv1, x = age_array[2:n_points], type = 'l')
  plot(y = deriv2, x = age_array[3:n_points], type = 'l')
   return(deriv2)
}

#mod.d <- Deriv(mod, newdata = pdat, term = "age")
# mod.dci <- confint(mod.d, term = "age")
# mod.dsig <- signifD(pdat$pred, d = mod.d[["age"]]$deriv, mod.dci[["age"]]$upper, mod.dci[["age"]]$lower)
# 
# 
# ylim <- with(pdat, range(upr, lwr,pred))
# plot(pred ~ age, data = pdat, type = "n", ylab = "pheno", ylim = ylim)
# lines(pred ~ age, data = pdat)
# lines(upr ~ age, data = pdat, lty = "dashed")
# lines(lwr ~ age, data = pdat, lty = "dashed")
# lines(unlist(mod.dsig$incr) ~ age, data = pdat, col = "red", lwd = 3)
# lines(unlist(mod.dsig$decr) ~ age, data = pdat, col = "blue", lwd = 3)






################################################
## Functions for derivatives of GAM(M) models ##
################################################
Deriv <- function(mod, n = 300, eps = 1e-7, newdata, term) {
  if(inherits(mod, "gamm"))
    mod <- mod$gam
  m.terms <- attr(terms(mod), "term.labels")
  if(missing(newdata)) {
    newD <- sapply(model.frame(mod)[, m.terms, drop = FALSE],
                   function(x) seq(min(x), max(x), length = n))
    names(newD) <- m.terms
  } else {
    newD <- newdata
  }
  newDF <- data.frame(newD) ## needs to be a data frame for predict
  X0 <- predict(mod, newDF, type = "lpmatrix")
  newDF <- newDF + eps
  X1 <- predict(mod, newDF, type = "lpmatrix")
  Xp <- (X1 - X0) / eps
  Xp.r <- NROW(Xp)
  Xp.c <- NCOL(Xp)
  ## dims of bs
  bs.dims <- sapply(mod$smooth, "[[", "bs.dim") - 1
  ## number of smooth terms
  t.labs <- attr(mod$terms, "term.labels")
  ## match the term with the the terms in the model
  if(!missing(term)) {
    want <- grep(term, t.labs)
    if(!identical(length(want), length(term)))
      stop("One or more 'term's not found in model!")
    t.labs <- t.labs[want]
  }
  nt <- length(t.labs)
  ## list to hold the derivatives
  lD <- vector(mode = "list", length = nt)
  names(lD) <- t.labs
  for(i in seq_len(nt)) {
    Xi <- Xp * 0
    want <- grep(t.labs[i], colnames(X1))
    Xi[, want] <- Xp[, want]
    df <- Xi %*% coef(mod)
    df.sd <- rowSums(Xi %*% mod$Vp * Xi)^.5
    lD[[i]] <- list(deriv = df, se.deriv = df.sd)
  }
  class(lD) <- "Deriv"
  lD$gamModel <- mod
  lD$eps <- eps
  lD$eval <- newD - eps
  lD ##return
}

confint.Deriv <- function(object, term, alpha = 0.05, ...) {
  l <- length(object) - 3
  term.labs <- names(object[seq_len(l)])
  if(missing(term)) {
    term <- term.labs
  } else { ## how many attempts to get this right!?!?
    ##term <- match(term, term.labs)
    ##term <- term[match(term, term.labs)]
    term <- term.labs[match(term, term.labs)]
  }
  if(any(miss <- is.na(term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  res <- vector(mode = "list", length = length(term))
  names(res) <- term
  residual.df <- df.residual(object$gamModel)
  tVal <- qt(1 - (alpha/2), residual.df)
  ##for(i in term.labs[term]) {
  for(i in term) {
    upr <- object[[i]]$deriv + tVal * object[[i]]$se.deriv
    lwr <- object[[i]]$deriv - tVal * object[[i]]$se.deriv
    res[[i]] <- list(upper = drop(upr), lower = drop(lwr))
  }
  res$alpha = alpha
  res
}

signifD <- function(x, d, upper, lower, eval = 0) {
  miss <- upper > eval & lower < eval
  incr <- decr <- x
  want <- d > eval
  incr[!want | miss] <- NA
  want <- d < eval
  decr[!want | miss] <- NA
  list(incr = incr, decr = decr)
}

plot.Deriv <- function(x, alpha = 0.05, polygon = TRUE,
                       sizer = FALSE, term,
                       eval = 0, lwd = 3,
                       col = "lightgrey", border = col,
                       ylab, xlab, main, ...) {
  l <- length(x) - 3
  ## get terms and check specified (if any) are in model
  term.labs <- names(x[seq_len(l)])
  if(missing(term)) {
    term <- term.labs
  } else {
    term <- term.labs[match(term, term.labs)]
  }
  if(any(miss <- is.na(term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  if(all(miss))
    stop("All terms in 'term' not found in model.")
  l <- sum(!miss)
  nplt <- n2mfrow(l)
  tVal <- qt(1 - (alpha/2), df.residual(x$gamModel))
  if(missing(ylab))
    ylab <- expression(italic(hat(f)*"'"*(x)))
  if(missing(xlab)) {
    xlab <- attr(terms(x$gamModel), "term.labels")
    names(xlab) <- xlab
  }
  if (missing(main)) {
    main <- term
    names(main) <- term
  }
  ## compute confidence interval
  CI <- confint(x, term = term)
  ## plots
  layout(matrix(seq_len(l), nrow = nplt[1], ncol = nplt[2]))
  for(i in term) {
    upr <- CI[[i]]$upper
    lwr <- CI[[i]]$lower
    ylim <- range(upr, lwr)
    plot(x$eval[,i], x[[i]]$deriv, type = "n",
         ylim = ylim, ylab = ylab, xlab = xlab[i], main = main[i], ...)
    if(isTRUE(polygon)) {
      polygon(c(x$eval[,i], rev(x$eval[,i])),
              c(upr, rev(lwr)), col = col, border = border)
    } else {
      lines(x$eval[,i], upr, lty = "dashed")
      lines(x$eval[,i], lwr, lty = "dashed")
    }
    abline(h = 0, ...)
    if(isTRUE(sizer)) {
      lines(x$eval[,i], x[[i]]$deriv, lwd = 1)
      S <- signifD(x[[i]]$deriv, x[[i]]$deriv, upr, lwr,
                   eval = eval)
      lines(x$eval[,i], S$incr, lwd = lwd, col = "blue")
      lines(x$eval[,i], S$decr, lwd = lwd, col = "red")
    } else {
      lines(x$eval[,i], x[[i]]$deriv, lwd = 2)
    }
  }
  layout(1)
  invisible(x)
}
