# FREGAT (c) 2017 Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

analyze.region <- function() {

	Z <- genotypes()
	m0 <- Z$m0
	if (is.null(Z$Z)) {
		warning("No polymorphic variants in region ", r, ", skipped")
		return(c(r, m0, 0, rep(NA, lgt - 3)))
	}
	if (length(Z$Z) == 1) {
		if (Z$Z == 'not in range') warning("Genotypes are not in range 0..2 in region ", r, ", skipped")
		if (Z$Z == 'not found') warning("Cannot locate region ", r, ", skipped")
		return(c(r, m0, NA, rep(NA, lgt - 3)))
	}
	m1 <- dim(Z$Z)[2]
	if (test == 'famFLM' & is.null(Z$pos)) {
		warning("Some positions are not assigned in region ", r, ", skipped")
		return(c(r, m0, m1, rep(NA, lgt - 3)))
	}
	if (test == 'famSKAT') {
		if (rho & qr(Z$Z)$rank > 1) return(c(r, m0, m1, pval.famSKATO(Z)))
	}
	return(c(r, m0, m1, pval.region(Z)))

}


analyze.single <- function(snv) {

	if (!gtype) {
		Z <- as.matrix(genodata[, snv])
	} else if (gtype == 1) {
		Z <- as.matrix(GenABEL::as.double.snp.data(genodata[, snv]))
	} else if (gtype == 2) {
		Z <- read.plink.region(bed, snvnames, idnames, snv)
		Z <- subset.SnpMatrix(Z, measured.ids)
		Z <- SnpMatrix2numeric(Z)
	} else if (gtype == 3) {
		if (length(annoType) == 1) {
			Z <- read.vcf.single(genodata, geneFile, reg, annoType)
			if (is.null(Z[[1]])) return(list(m0 = 0, Z = NULL, w = NULL, pos = NULL))
			Z <- t(Z[[1]])
		} else { 
			Z <- c()
			for (anno in annoType) {
				Z0 <- read.vcf.single(genodata, geneFile, reg, anno)
				if (is.null(Z0[[1]])) next
				Z <- cbind(Z, t(Z0[[1]]))
			}
		}
		if (is.null(Z)) return(list(m0 = 0, Z = NULL, w = NULL, pos = NULL))
		if (!dim(Z)[2] == 0) Z <- as.matrix(Z[measured.ids, ])
	}

	n11 <- sum(!is.na(Z))
	if (n11 == 0) {
		warning("Zero call rate variant ", r, ", skipped")
		return(c(0, NA, NA, NA, NA))
	}

	Z[Z == -9] <- NA

	MAF <- mean(Z, na.rm = TRUE) / 2
	if (MAF > .5) warning("MAF > 0.5 at variant ", r)

	if (mode != 'add') {
		Z <- round(Z)
		index <- Z >= 0 & Z <= 2
		if (sum(index, na.rm = T) == 0) {
			warning("Genotypes are not in range 0..2 at variant ", r, ', skipped')
			return(c(n11, NA, NA, NA, NA))
		}
		Z[!index] <- NA 
		if (mode == 'rec') Z[Z == 1] <- 0
		Z[Z == 2] <- 1
	}
	
	#impute missing genotypes
	if (mode == 'add') {
		if (impute.method == 'mean') {
			Z[is.na(Z)] <- mean(Z, na.rm = T)
		} else { # BLUE imputation
			Z[is.na(Z)] <- sum(colInvOmega[!is.na(Z)] * Z[!is.na(Z)]) / sum(colInvOmega[!is.na(Z)])
		}
	} else {
		if (impute.method == 'mean') {
			q <- mean(Z, na.rm = T) / 2
		} else {
			q <- sum(colInvOmega[!is.na(Z)] * Z[!is.na(Z)]) / sum(colInvOmega[!is.na(Z)]) / 2
		}
		if (mode == 'dom') {
			Z[is.na(Z)] <- 2 * q - q ^ 2
		} else {
			Z[is.na(Z)] <- q ^ 2
		}
	}

	n11 <- n1

	v <- var(Z)
	if (v == 0) {
		warning("Nonpolymorphic variant ", r, ", skipped")
		return(c(n11, NA, NA, NA, 0))
	}

	return(c(n11, pval.single(Z), MAF))
}

