# FREGAT (c) 2016 Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

'single.point' <- function (formula, phenodata, genodata, kin = NULL, nullmod, regions = NULL,
mode = 'add', ncores = 1, return.time = FALSE, impute.method = 'mean', write.file = FALSE, ...) { #, na.action = 'impute'

t0 <- proc.time()

############## CHECKS

if (missing(genodata)) stop("'genodata' not found")

ch <- check.ini(substitute(formula), phenodata, kin)
for(i in 1:length(ch)) assign(names(ch)[i], ch[[i]])

mode <- match.arg(mode, c('add', 'dom', 'rec'))
impute.method <- match.arg(impute.method, c('mean', 'blue'))

# na.action <- match.arg(na.action, c('impute', 'omit'))

k <- FALSE
gtype <- 0

if (class(genodata) %in% c('gwaa.data', 'snp.data')) {
	if (requireNamespace("GenABEL", quietly = TRUE)) {
		if (class(genodata) == 'gwaa.data') genodata <- GenABEL::gtdata(genodata)
		gtype <- 1
	} else { stop(paste("'genodata' class '", class(genodata),"' cannot be processed, 'GenABEL' package not found",sep='')) }
} else if (class(genodata) == 'character' & length(genodata) == 1) {
	#if (grepl('vcf', genodata)) stop("VCF not yet supported for single point analysis")
	tmp <- check.geno(genodata, regions, n, ...)
	for (i in 1:length(tmp)) assign(names(tmp)[i], tmp[[i]])
} else if (class(genodata) == 'snpMatrix') {
	genodata <- as(genodata, 'numeric')
} else if (!class(genodata) %in% c('matrix', 'data.frame')) {
	if (is.numeric(genodata)) {
		genodata <- as.matrix(genodata)
	} else stop("Wrong 'genodata' class")
} else if (!is.numeric(genodata[1, 1])) {
	stop("'genodata' elements should be numeric")
}

if (gtype < 2) {
	if (dim(genodata)[1] != n) stop("Dimensions of 'phenodata' and 'genodata' do not match")
	k <- dim(genodata)[2]
}

X <- check.covariates(n, formula, phenodata)

measured.ids <- as.logical(X$measured.ids)
X <- as.matrix(X$X)
n1 <- sum(measured.ids)

if (n1 != n) {
	if (gtype < 2) genodata <- genodata[measured.ids, ]
	if (H2est) kin <- kin[measured.ids, measured.ids]
}

if (gtype < 2) {
	snvnames <- colnames(genodata)
	if (is.null(snvnames)) snvnames <- paste('snv', 1:k, sep = '')
}

if (gtype < 3) {
	if (!is.null(regions)) {
		regions <- as.matrix(regions)
		rtype <- dim(regions)[2]
		if (rtype == 1) {
			if (length(regions) != k) {
				warning("Dimensions of 'regions' and 'genodata' do not match, 'regions' ignored")
				regions <- NULL
			}
		} else {
			if (rtype > 2) {
				warning("'regions' with more than two columns, only the first two will be used")
				regions <- regions[,1:2]
				rtype <- 2
			}
			snvnames.out <- unique(regions[, 1])
			if (is.null(snvnames)) stop("snv names not found in 'genodata'")
			v <- snvnames.out %in% snvnames
			k <- sum(v)
			if (k == 0) stop("No variants from 'regions' found in 'genodata'")
			if (length(snvnames.out) > k) warning(length(snvnames.out) - k, " variant(s) from 'regions' not found in 'genodata'")
			snvnames.out <- snvnames.out[v]
			sq <- match(snvnames.out, snvnames)
	#		if (gtype < 2) {
	#			genodata <- genodata[, colnames(genodata) %in% snvnames.out]
	#			gc()
	#			k <- dim(genodata)[2]
	#		}
		}
	}

	if (is.null(regions)) {
		regions <- rep('none', k)
		rtype <- 1
	}
	if (rtype == 1) {
		snvnames.out <- snvnames
		sq <- 1:k
	}
} else {
	if (is.null(regions)) stop ("'regions' should be set when using VCF file")
	regions <- as.matrix(regions)
	tmp <- check.regions(k, regions)
	for (i in 1:length(tmp)) assign(names(tmp)[i], tmp[[i]])
	k <- nreg

}
ncores <- check.cores(ncores, k)

############# NULL MODEL

t00 <- proc.time()

run.null <- check.nullmod(nullmod, X)

if (run.null) nullmod <- NullMixedModel(X[, 1], X[, -1], kin * 2, H2est = H2est, ...)

if (H2est) {
	SIGMA <- nullmod$total.var * (kin * 2 * nullmod$h2 + diag(n1) * (1 - nullmod$h2))
	SIGMAi <- chol2inv(chol(SIGMA))
	if (impute.method == 'blue') colInvOmega <- as.vector(colMeans(SIGMAi))
} else {
	if (impute.method == 'blue') impute.method <- 'mean'
}

pheno <- X[, 1] - X[, -1] %*% as.matrix(nullmod$alpha[match(colnames(X)[-1], rownames(nullmod$alpha)), 1])
X <- X[, -1]

ttt0 <- (proc.time() - t00)

########## PRELIMINARIES

if (H2est) {
	Cori <- nullmod$total.var * SIGMAi
	CholInvCor <- chol(Cori)
	covariate <- CholInvCor %*% X
	P11 <- diag(n1) - covariate %*% solve(t(covariate) %*% covariate, t(covariate))# == P11 famBT
	P11CholInvCor <- P11 %*% CholInvCor
	pheno <- CholInvCor %*% pheno  # making independent pheno
} else {
	P11CholInvCor <- diag(n1) - X %*% solve(t(X) %*% X, t(X))
}

########## ANALYSE

t00 <- proc.time()

environment(analyze.single) <- environment()
environment(pval.single) <- environment()
if (gtype == 3) {
	environment(genotypes) <- environment()
	test <- 'single'
	omit.linear.dependent <- FALSE
}

if (ncores == 1) {
	
	if (write.file != FALSE) write.table(t(c('region', 'snv', 'n', 'pvalue', 'beta', 'se.beta', 'AF')), file = write.file, col.names = F, row.names = F, quote = F)
	if (gtype == 3) {
		
		out <- matrix(NA, 0, 7)
		
		for (i in 1:nreg) {
			r <- as.character(l[i])
			Z <- genotypes()$Z
			if (is.null(Z)) {
				warning("No polymorphic variants in region ", r, ", skipped")
			} else if (length(Z) == 1) {
				if (Z == 'not in range') warning("Genotypes are not in range 0..2 in region ", r, ", skipped")
			} else {
				Z <- as.matrix(Z)
				for (i in 1:dim(Z)[2]) {
					Zi <- Z[, i]
					tmp <- c(r, colnames(Z)[i], sum(!is.na(Zi)), pval.single(Zi), mean(Zi, na.rm = T) / 2)
					if (write.file != FALSE) write.table(t(tmp), file = write.file, append = T, col.names = F, row.names = F, quote = F)
					out <- rbind(out, tmp)
				}
			}
		}
		
		out <- as.data.frame(out, stringsAsFactors = FALSE)
		out[, 3:7] <- sapply(3:7, function(x) as.numeric(out[, x]))
		colnames(out) <- c('region', 'snv', 'n', 'pvalue', 'beta', 'se.beta', 'AF')
	
	} else {
	
		out <- data.frame(snv = snvnames.out, n = NA, pvalue = NA, beta = NA, se.beta = NA, AF = NA)

		for (i in 1:k) {
			r <- snvnames.out[i]
	#		shd be: snvnames.out[i] == snvnames[s]
			s <- sq[i]
			out[i, 2:6] <- analyze.single(s)
		}

	}

} else { # gtype=3?

	nj <- floor(k / ncores)
	vj <- rep(nj, ncores)
	nadd <- k %% ncores
	if (nadd > 0) vj[1:nadd] <- nj + 1

	cl <- makeCluster(ncores)
	doParallel::registerDoParallel(cl)
	clusterExport(cl, varlist = ls(), envir = environment())

	if (gtype == 3) {

		environment(read.vcf) <- environment()
		out <- foreach::'%dopar%'(foreach::foreach(j = 1:ncores, .combine = rbind, .inorder = F), {

			i0 <- (sum(vj[1:(j - 1)]) * (j > 1) + 1)
			out <- c()

			for (i in i0:(i0 + vj[j] - 1)) {
				r <- as.character(l[i])
				Z <- genotypes()$Z
				if (!is.null(Z)) {
					if (length(Z) > 1) {
						Z <- as.matrix(Z)
						for (i in 1:dim(Z)[2]) {
							Zi <- Z[, i]
							tmp <- c(r, colnames(Z)[i], sum(!is.na(Zi)), pval.single(Zi), mean(Zi, na.rm = T) / 2)
							out <- rbind(out, tmp)
						}
					}
				}
			}
			out
		})

		stopCluster(cl)

		out <- as.data.frame(out, stringsAsFactors = FALSE)
		out[, 3:7] <- sapply(3:7, function(x) as.numeric(as.character(out[, x])))
		colnames(out) <- c('region', 'snv', 'n', 'pvalue', 'beta', 'se.beta', 'AF')

	} else {

		out <- foreach::'%dopar%'(foreach::foreach(j = 1:ncores, .combine = rbind, .inorder = F), {

	#		if (class(genodata) == 'snp.data') requireNamespace("GenABEL", quietly = TRUE)

			i0 <- (sum(vj[1:(j - 1)]) * (j > 1) + 1)
			out <- c()

			for (i in i0:(i0 + vj[j] - 1)) {
				r <- snvnames.out[i]
				s <- sq[i]
				out <- rbind(out, c(r, analyze.single(s)))
			}
			out
		})

		stopCluster(cl)

		out <- as.data.frame(out)
		out[, 2:6] <- sapply(2:6, function(x) as.numeric(as.character(out[, x])))
		colnames(out) <- c('snv', 'n', 'pvalue', 'beta', 'se.beta', 'AF')
	}
}

if (gtype < 3) {
	if (rtype == 1) {
		out <- cbind(region = as.character(regions), out)
	} else {
		out <- cbind(region = as.character(regions[, 2]), snv = as.character(regions[, 1]), out[match(regions[, 1], out[, 1]), 2:6])
	}
}
rownames(out) <- NULL

ttt <- (proc.time() - t00)

out0 <- out
out <- c()
out$results <- out0
out$nullmod <- nullmod
if (return.time) {
	out$time$null <- ttt0
	out$time$regions <- ttt
	out$time$total <- (proc.time() - t0)
}

out

}