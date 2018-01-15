\name{PCA}
\alias{PCA}
\title{Principal Components Analysis}
\description{
Test for association between a trait and principal components of genotypes within a region, on summary statistics
}
\usage{
PCA(formula, phenodata, genodata, kin = NULL, nullmod, regions = NULL,
sliding.window = c(20, 10), mode = "add", ncores = 1, return.time = FALSE,
beta.par = c(1, 1), weights = NULL, var.fraction = .85, impute.method = "mean",
write.file = FALSE, ...)
}

\arguments{
	\item{formula}{referring to the column(s) in \code{phenodata} to be analyzed as outcome and,
	if needed, covariates.}

	\item{phenodata}{a data frame containing columns mentioned in \code{formula}: trait to analyze and,
	if needed, covariates. Individuals not measured for trait or covariates will be omitted.}

	\item{genodata}{an object with genotypes to analyze. Several formats are allowed:\cr
	- a data frame or matrix (with individuals in the rows and genetic variants in the columns)
	containing genotypes coded as AA = 0, Aa = 1 and aa = 2, where a is a minor allele.\cr
	- for PLINK binary data format, a character string indicating a *.bed file name (*.bim and *.fam
	files should have the same prefix). This will make use of \code{read.plink()} function.\cr
	- for VCF format, a character string indicating a *vcf.gz file name. This will require
	\code{seqminer} R-package to be installed. Its \code{readVCFToMatrixByGene()} function will be
	used to read VCF file gene-wise. The function also requires a geneFile, a text file listing all
	genes in refFlat format (see Examples below). VCF file should be bgzipped and indexed by Tabix.\cr
	- an object of \code{gwaa.data} or \code{snp.data} class (this will require
	\code{GenABEL} R-package to be installed).}

	\item{kin}{a square symmetric matrix giving the pairwise kinship coefficients between analyzed
	individuals. Under default \code{kin = NULL} all individuals will be considered as unrelated.}

	\item{nullmod}{an object containing parameter estimates under the null model. Setting \code{nullmod}
	allows to avoid re-estimation of the null model that does not depend on genotypes and can be
	calculated once for a trait. If not set, the null model parameters will be estimated within the function.
	The \code{nullmod} object in proper format can be obtained by \code{null.model()} function
	or any analysis function in \code{FREGAT}.}

	\item{regions}{an object assigning regions to be analyzed. This can be:\cr
	- a vector of length equal to the number of genetic variants assigning the region
	for each variant (see Examples).\cr
	- a data frame / matrix with names of genetic variants in the first column and names of regions
	in the second column (this format allows overlapping regions).\cr
	- for VCF format, a character vector with names of genes to analyze.\cr
	If NULL, \code{sliding.window} parameters will be used.}

	\item{sliding.window}{the sliding window size and step. Has no effect if \code{regions} is defined.}

	\item{mode}{the mode of inheritance: "add", "dom" or "rec" for additive, dominant or recessive mode,
	respectively. For dominant (recessive) mode genotypes will be recoded as AA = 0, Aa = 1 and aa = 1
	(AA = 0, Aa = 0 and aa = 1), where a is a minor allele. Default mode is additive.}

	\item{ncores}{number of CPUs for parallel calculations. Default = 1.}

	\item{return.time}{a logical value indicating whether the running time should be returned.}

	\item{beta.par}{two positive numeric shape parameters in the beta distribution to assign weights 
	for each genetic variant as a function of MAF in the default weights function (see Details). Default = c(1, 1).}

	\item{weights}{a numeric vector or a function of minor allele frequency (MAF) to assign weights
	for each genetic variant in the weighted kernels. Has no effect if one of unweighted kernels
	was chosen. If NULL, the weights will be calculated using the beta distribution (see Details).}

	\item{var.fraction}{minimal proportion of genetic variance within region that should be explained by principal
	components used (see Details for more info).}

	\item{impute.method}{a method for imputation of missing genotypes. It can be either "mean" (default)
	or "blue". If "mean" the genotypes will be imputed by the simple mean values. If "blue"
	the best linear unbiased estimates (BLUEs) of mean genotypes will be calculated
	taking into account the relationships between individuals [McPeek, et al., 2004,
	DOI: 10.1111/j.0006-341X.2004.00180.x] and used for imputation.}

	\item{write.file}{output file name. If specified, output (as it proceeds) will be written 
	to the file.}

	\item{...}{other arguments that could be passed to \code{null()}, \code{read.plink()}\cr
	and \code{readVCFToMatrixByGene()}.}

}
\details{
	PCA test is a useful tool to detect association between genetic variants of a region and a trait
	when genetic variants are strongly correlated. PCA test is based on the spectral decomposition of
	genotype matrix. The number of top principal components will be chosen
	in such a way that >= \code{var.fraction} of region variance can be explained by these PCs.
	By default,	\code{var.fraction} = 0.85, i.e. >= 85\% of region variance will be explained by PCs.
	If \code{var.fraction} = 1 then the results of PCA test and MLR-based test are identical.
	
	\code{beta.par = c(a, b)} can be used to set weights for genetic variants.
	Given the shape parameters of the beta function, \code{beta.par = c(a, b)}, 
	the weights are defined using probability density function of the beta distribution:\cr
	\cr
	\eqn{W_{i}=(B(a,b))^{^{-1}}MAF_{i}^{a-1}(1-MAF_{i})^{b-1} },\cr
	\cr
	where \eqn{MAF_{i}} is a minor allelic frequency for the \eqn{i^{th}} genetic variant in the region,
	which is estimated from genotypes, and \eqn{B(a,b)} is the beta function. This way of defining weights
	is the same as in original SKAT (see [Wu, et al., 2011] for details).
}
\value{
	
	A list with values:

	\item{results}{a data frame containing P values, numbers of variants
	and filtered variants for each of analyzed regions.
	It also contains the number of the principal components
	used for each region and the proportion of genetic variance
	they make up.}

	\item{nullmod}{an object containing the estimates of the null model parameters:
	heritability (h2), total variance (total.var), estimates of fixed effects of covariates (alpha),
	the gradient (df), and the total log-likelihood (logLH).}

	\item{sample.size}{the sample size after omitting NAs.}

	\item{time}{If \code{return.time = TRUE} a list with running times for null model, regional analysis and total analysis is returned. See \code{proc.time()} for output format.}

}
\references{
	Jolliffe, I.T. A note on the use of principal components in regression. J R Stat Soc Ser C 31, 300-303 (1982).
	}
\examples{

data(example.data)

## Run PCA with sliding window (default):
out <- PCA(trait ~ age + sex, phenodata, genodata, kin)

## Run PCA with regions defined in snpdata$gene and with
## null model parameters obtained in the first run:
out <- PCA(trait ~ age + sex, phenodata, genodata, kin,
	out$nullmod, regions = snpdata$gene)

## Run PCA parallelized on two cores (this will require
## 'foreach' and 'doParallel' R-packages installed and
## cores available):
out <- PCA(trait ~ age + sex, phenodata, genodata, kin,
	out$nullmod, ncores = 2)

## Run PCA with genotypes in VCF format:
VCFfileName <- system.file(
	"testfiles/1000g.phase1.20110521.CFH.var.anno.vcf.gz",
	package = "FREGAT")
geneFile <- system.file("testfiles/refFlat_hg19_6col.txt.gz",
	package = "FREGAT")
phe <- data.frame(trait = rnorm(85))
out <- PCA(trait, phe, VCFfileName, geneFile = geneFile,
	reg = "CFH", annoType = "Nonsynonymous")

## Run PCA with genotypes in PLINK binary data format:
bedFile <- system.file("testfiles/sample.bed",
	package = "FREGAT")
phe <- data.frame(trait = rnorm(120))
out <- PCA(trait, phe, bedFile)

}
