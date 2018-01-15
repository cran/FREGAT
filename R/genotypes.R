# FREGAT (c) 2016

genotypes <- function() {#(r, reg, genodata, mode, weights, fweights, beta.par, positions, flip.genotypes, test, et cetera) {
	if (gtype < 3) {
		if (all(is.na(reg))) return(list(m0 = 0, Z = 'not found', w = NULL, pos = NULL))
	}
	w <- pos <- NULL
	if (!gtype) {
		Z <- as.matrix(genodata[, reg])
	} else if (gtype == 1) {
		Z <- as.matrix(GenABEL::as.double.snp.data(genodata[, reg]))
	} else if (gtype == 2) { # mozhno li dostat' pos is failov?
		Z <- read.plink.region(bed, snvnames, idnames, reg)
		Z <- subset.SnpMatrix(Z, measured.ids)
		Z <- SnpMatrix2numeric(Z)
	} else if (gtype == 3) {
		if (length(annoType) == 1) {
			Z <- read.vcf(genodata, geneFile, r, annoType)
			if (is.null(Z[[1]])) return(list(m0 = 0, Z = NULL, w = NULL, pos = NULL))
			Z <- t(Z[[1]])
			if (test == 'famFLM') pos <- read.vcf.pos(genodata, geneFile, r, annoType)
		} else {
			Z <- c()
			for (anno in annoType) {
				Z0 <- read.vcf(genodata, geneFile, r, anno)
				if (is.null(Z0[[1]])) next
				Z <- cbind(Z, t(Z0[[1]]))
				if (test == 'famFLM') {
					pos0 <- read.vcf.pos(genodata, geneFile, r, anno)
					pos <- c(pos, pos0)
				}
			}
		}
		if (is.null(Z)) return(list(m0 = 0, Z = NULL, w = NULL, pos = NULL))
		if (!dim(Z)[2] == 0) Z <- as.matrix(Z[measured.ids, , drop = FALSE])
		#if (weights.table) {
		
		#}
	}
	if (dim(Z)[2] == 0) return(list(m0 = 0, Z = NULL, w = NULL, pos = NULL))
	if (gtype < 3) { ## weights with vcf not covered!!
		w <- pos <- NULL
		if (test == 'famFLM') pos <- positions[reg]
		if (!is.null(weights)) w <- weights[reg]
	}

	m0 <- dim(Z)[2]

	Z[Z == -9] <- NA

	MAF <- colMeans(Z, na.rm = TRUE) / 2

	if (any(is.na(MAF))) {
		if (all(is.na(MAF))) return(list(m0 = m0, Z = NULL, w = w, pos = pos))
		index <- !is.na(MAF)
		Z <- as.matrix(Z[, index, drop = FALSE])
		MAF <- MAF[index]
		if (test == 'famFLM') pos <- as.vector(pos[index])
		if (!is.null(weights)) w <- as.vector(w[index])
	}

	index <- sapply(1:dim(Z)[2], function(x) all(Z[, x] >= 0, na.rm = T) & all(Z[, x] <= 2, na.rm = T))
	if (sum(index) == 0) return(list(m0 = m0, Z = 'not in range', w = NULL, pos = NULL))
	Z <- as.matrix(Z[, index, drop = FALSE])
	MAF <- MAF[index]
	if (test == 'famFLM') pos <- as.vector(pos[index])
	if (!is.null(weights)) w <- as.vector(w[index])

	if (any(MAF > .5) & test %in% c('famBT', 'famSKAT', 'MLR')) {
		if (flip.genotypes | test == 'MLR') {
			v <- MAF > .5
			Z[, v] <- abs(2 - Z[, v])
			MAF[v] <- 1 - MAF[v]
		} else warning("Check whether minor allele homozygotes are coded as 2 at reg ", r, ", use 'flip.genotypes = TRUE' if needed")
	}

	if (mode != 'add') {
		Z <- round(Z)
		if (mode == 'rec') Z[Z == 1] <- 0
		Z[Z == 2] <- 1
	}

	v <- sapply(1:dim(Z)[2], function(x) var(Z[, x], na.rm = TRUE))
	Z <- as.matrix(Z[, v > 0, drop = FALSE])
	if (test == 'famFLM') pos <- pos[v > 0]
	if (!is.null(weights)) { w <- w[v > 0]
	} else if (!is.null(fweights)) {
		MAF <- MAF[v > 0]
		w <- fweights(MAF)
	}
	if (dim(Z)[2] == 0) return(list(m0 = m0, Z = NULL, w = NULL, pos = NULL))

	#impute missing genotypes
	if (mode == 'add') {
		if (impute.method == 'mean') {
			for (z in 1:dim(Z)[2]) Z[is.na(Z[, z]), z] <- mean(Z[, z], na.rm = T)
		} else { # BLUE imputation
			for (z in 1:dim(Z)[2]) {
				gen <- as.vector(Z[, z]) # these are genotypes of all individuals for one marker
				BLUE <- sum(colInvOmega[!is.na(gen)] * gen[!is.na(gen)]) / sum(colInvOmega[!is.na(gen)])
				Z[is.na(gen), z] <- BLUE
			}
		}
	} else {
		if (impute.method == 'mean' & mode == 'dom') {
			for (z in 1:dim(Z)[2]) {
				q <- mean(Z[, z], na.rm = T) / 2
				Z[is.na(Z[, z]), z] <- 2 * q - q ^ 2
			}
		} else if (impute.method == 'mean' & mode == 'rec') {
			for (z in 1:dim(Z)[2]) {
				q <- mean(Z[, z], na.rm = T) / 2
				Z[is.na(Z[, z]), z] <- q ^ 2
			}
		} else if (impute.method == 'blue' & mode == 'dom') {
			for (z in 1:dim(Z)[2]) {
				gen <- as.vector(Z[, z])
				q <- sum(colInvOmega[!is.na(gen)] * gen[!is.na(gen)]) / sum(colInvOmega[!is.na(gen)]) / 2
				Z[is.na(gen), z] <- 2 * q - q ^ 2
			}
		} else if (impute.method == 'blue' & mode == 'rec') {
			for (z in 1:dim(Z)[2]) {
				gen <- as.vector(Z[, z])
				q <- sum(colInvOmega[!is.na(gen)] * gen[!is.na(gen)]) / sum(colInvOmega[!is.na(gen)]) / 2
				Z[is.na(gen), z] <- q ^ 2
			}
		}
	}

	if (omit.linear.dependent) {
		# detect and eliminate linear-dependent genetic variants
		dqr <- qr(Z) # QR-decomposition
		index <- dqr$pivot[1:dqr$rank] # indexes of dependent genetic variants
		Z <- as.matrix(Z[, index]) # elimination
		if (test == 'famFLM') pos <- as.vector(pos[index])
		if (!is.null(w)) w <- as.vector(w[index])
	}
	if (test == 'famBT') {
		return(list(m0 = m0, Z = Z, w = w, pos = pos))
	}

	if (test == 'famFLM') {
		if (any(is.na(pos))) {
			pos <- NULL
		} else {
			v <- order(pos)
			Z <- as.matrix(Z[, v])
			w <- as.vector(w[v])
			pos <- pos[v]
			if (sum(duplicated(pos)) > 0) {
				for (i in 1:length(pos)) {
					pos[duplicated(pos)] <- pos[duplicated(pos)] + 1
					if (sum(duplicated(pos)) == 0) break()
				}
			}
		if (max(pos) > min(pos)) pos <- (pos - min(pos)) / (max(pos) - min(pos))
		if (max(pos) == min(pos)) pos <- .5
		}
	if (flip.genotypes & dim(Z)[2] > 1) Z <- flipper(Z)
	}
	return(list(m0 = m0, Z = Z, w = w, pos = pos))
}

flipper <- function(Z) {
	for (j in 1:(dim(Z)[2] - 1)) {
		A_B <- mean((Z[,j+1]-1)*(Z[,j]-1))
		if (A_B < 0) { 
		Z[, j+1] <- (2 - Z[, j + 1])
		}
	}
	return(Z)
}



read.vcf <- function(genodata, geneFile, r, annoType) {
	if (getRversion() >= "3.2.0") {
		invisible(capture.output(invisible(capture.output(Z <- seqminer::readVCFToMatrixByGene(genodata, geneFile, r, annoType), type = 'output')), type = 'message'))
	} else { Z <- seqminer::readVCFToMatrixByGene(genodata, geneFile, r, annoType) }
	return(Z)
}

read.vcf.pos <- function(genodata, geneFile, r, annoType) {
	if (getRversion() >= "3.2.0") {
		invisible(capture.output(invisible(capture.output(pos <- seqminer::readVCFToListByGene(genodata, geneFile, r, annoType, 'POS', '', '')$POS, type = 'output')),  type = 'message'))
	} else { pos <- seqminer::readVCFToListByGene(genodata, geneFile, r, annoType, 'POS', '', '')$POS }
	return(pos)
}