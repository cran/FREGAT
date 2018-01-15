# FREGAT (c) 2016 Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

pval.famBT <- function(Z) {

	w <- Z$w
	Z <- P11CholInvCo %*% Z$Z # Cor
	if (ncol(Z) > 1) {
#CC		Z <- Z %*% diag(w) # n x m
		Z <- t(t(Z) * w)
		Z <- rowMeans(Z)
	} else {
		Z <- Z * w
	}
	se.beta0 <- 1 / sum(Z * Z)
	se.beta <- se.beta0 * nullmod$total.var
	est.beta <- sum(Z * as.vector(pheno)) * se.beta0
	chi2 <- (est.beta ^ 2) / se.beta
	p <- pchisq(chi2, 1, lower.tail = F)

	return(c(p, est.beta, sqrt(se.beta)))#, m0, m1))

}

pval.MLR <- function(Z) {
	if (H2est) Z$Z <- CholInvCo %*% Z$Z

	fit <- glm(pheno ~ covariate + Z$Z - 1, family = "gaussian")
	if (stat == 'F') {
		p <- anova(fit, test = stat)[3, 6]
	} else {
		p <- anova(fit, test = stat)[3, 5]
	}
	return(p)

}

pval.famFLM <- function(Z) {

	m1 <- dim(Z$Z)[2]
	Z$Z <- t(t(Z$Z) * Z$w)
	model <- model0

	# base condition is m >= kg >= kb
	# previously ensured that
	# if g then kg > kb
	# if kg == kb then !g
	# all fourier bases are odd

	# g:
	# m >= kg => default
	# m < kg => kg <- m
		# if fourier if kg is even => kg <- kg -1
		# ensured that m >= kg 
		# if kg <= kb => MLR
		# if kg > kb => genobasis recalculated

	if (g) {
		if (m1 >= kg0) {
			genobasis <- genobasis0
			J <- J0
		} else {
			if (GVF == 'bspline') {
				order <- min(order0, m1)
				kg <- m1
			} else if (GVF == 'fourier') {
				if (m1 %% 2) kg <- m1 else kg <- (m1 - 1)
			}
			if (kg <= kb0) {  # m1 <= kg = kb, (BB, BF) -> 0B, (FF, FB) -> 0F with kb, betabasis recalculated
				test <- 'MLR'
			} else {  # m1 >= kg > kb, genobasis recalculated
				if (GVF == 'bspline') {
					genobasis <- create.bspline.basis(norder = order, nbasis = kg)
				} else if (GVF == 'fourier') {
					genobasis <- create.fourier.basis(c(0, 1), nbasis = kg)
				}
				J <- inprod(genobasis, betabasis0)
				model <- paste(GVF, kg, '-', BSF, kb0, sep = '')
			}
		}

	# !g:
	# m > kb => default
	# m <= kb => MLR

	} else {
		if (m1 > kb0) {
			betabasis <- betabasis0
		} else {
			test <- 'MLR'
		}
	}

	if (test == 'MLR') {
		test <- 'famFLM'
		return(c(pval.MLR(Z), 'MLR'))
	}

	### Calculation of matrix FI in given positions of a region
	if (g) {
		B <- eval.basis(Z$pos, genobasis)
		### Calculation of formula (1) for functional genotypes (depend on ONLY positions)
		FFF <- ginv(t(B) %*% B) %*% t(B)
		U <- Z$Z %*% t(FFF)
		UJ <- matrix( U %*% J, ncol = dim(J)[2])
	} else {
		B <- eval.basis(Z$pos, betabasis)
		UJ <- Z$Z %*% B
	}
	if (H2est) UJ <- CholInvCo %*% UJ  # making centralized and independent genotypes # Cor
	fit <- glm(pheno ~ covariate + UJ - 1, family = "gaussian")
	if (stat == 'F') {
		p <- anova(fit, test = stat)[3, 6]

	} else { p <- anova(fit, test = stat)[3, 5] }

	c(p, model)

}

pval.single <- function(Z) {
	Z <- P11CholInvCor %*% Z
	se.beta0 <- 1 / sum(Z * Z)
	se.beta <- se.beta0 * nullmod$total.var
	est.beta <- sum(Z * as.vector(pheno)) * se.beta0
	chi2 <- (est.beta ^ 2) / se.beta
	p <- pchisq(chi2, 1, lower.tail = F)
	return(c(p, est.beta, sqrt(se.beta)))
}

pval.PCA <- function(Z) {
	Gcc <- (Z$Z - rep(1, n1) %*% t(colMeans(Z$Z))) # n1 = N
	Gc <- t(Z$w * t(Gcc))
	numberPCA(yc, Gc, n.pca, var.fraction) # n.pca = N - 1
}

numberPCA <- function(yc, Gc, n, var.fraction) {
	m <- qr(Gc)$rank
	pCA <- PC(Gc, n) # n = N - 1
	COMPs0 <- as.matrix(pCA$scores[, 1:m])     # matrix of independent scores
	CPV <- pCA$importance[3, 1:m]   # Cumulative Proportion of Variance (CPV)
	comp <- which(CPV >= var.fraction)      # components for which Explained variance fraction is about 85%
	minP <- 100
	minM <- 0
	M <- min(comp)
	COMPs <- as.matrix(COMPs0[,1:M])
	m <- qr(COMPs)$rank
	GY <- as.vector(t(COMPs) %*% yc)
	CC <- t(COMPs) %*% COMPs
	if (M > 1) {
		RSS <- n - sum(GY * as.vector((CC %^% (-1)) %*% GY))
	} else { RSS <- n - GY * GY / CC }
	Fstat <- ((n - m) / m) * (n - RSS) / RSS    # F-statistic
	p <- pf(Fstat, m, n - m, lower.tail = FALSE)
	if (p < minP) minM <- M
	minP <- min(p, minP)
	c(minP, minM, CPV[minM])
}

PC <- function(X, n) {
	U <- (t(X) %*% X) / n # n = N - 1
	eX1 <- eigen(U, symmetric = TRUE)
	values <- eX1$values
	vectors <- eX1$vectors
	values[values < 0] <- 0
	sdev <- sqrt(values)
	ind.var <- sdev ^ 2
	total.var <- sum(ind.var)
	scores <- X %*% vectors
	scorenames <- paste('PC', 1:ncol(X), sep = '')
	colnames(scores) <- colnames(vectors) <- scorenames
	rownames(vectors) <- colnames(X)
	prop.var <- ind.var / total.var
	cum.var <- rep(NA, ncol(X))
	for (i in 1:ncol(X)) { cum.var[i] <- sum(prop.var[1:i]) }
	importance <- rbind(sdev, prop.var, cum.var)
	colnames(importance) <- scorenames
	rownames(importance) <- c("Standard Deviation", "Proportion of Variance", "Cumulative Proportion")
	list(values=values,vectors=vectors,scores=scores,importance=importance,sdev=sdev)
}

"%^%" <- function(U, k) {
	UUU <- eigen(U, symmetric = TRUE)  # UUU = Uvec %*% diag(Uval) %*% t(Uvec)
	Uval <- UUU$val
	Uvec <- UUU$vec
	Uvec <- Uvec[,Uval > 1e-7]
	Uval <- Uval[Uval > 1e-7]
	Uvec %*% (t(Uvec) * (Uval ^ k))
}
