\name{famFLM}
\alias{famFLM}
\title{family Functional Linear Model}
\description{
A region-based association test for familial or population data under functional linear models (functional data analysis approach)
}
\usage{
famFLM(formula, phenodata, genodata, kin = NULL, nullmod,
regions = NULL, sliding.window = c(20, 10), mode = "add",
ncores = 1, return.time = FALSE, beta.par = c(1, 1),
weights = NULL, positions = NULL, GVF = FALSE,
BSF = "fourier", kg = 30, kb = 25, order = 4, stat = "F",
flip.genotypes = FALSE, impute.method = 'mean',
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
	for each genetic variant as a function of MAF (see Details). Default = c(1, 1) corresponds
	to standard unweighted FLM. Has no effect if \code{weights} are defined.}

	\item{weights}{a numeric vector or a function of minor allele frequency (MAF) to assign weights
	for each genetic variant in the weighted kernels. Has no effect if one of unweighted kernels
	was chosen. If NULL, the weights will be calculated using the beta distribution (see Details).}

	\item{positions}{a vector of physical positions for genetic variants in \code{genodata}.
	Not used when VCF file supplied.}

	\item{GVF}{a basis function type for Genetic Variant Functions. Can be set to
	"bspline" (B-spline basis) or "fourier" (Fourier basis). The default \code{GVF = FALSE}
	assumes beta-smooth only. If \code{GVF = TRUE} the B-spline basis will be used.}

	\item{BSF}{a basis function type for beta-smooth. Can be set to "bspline" (B-spline basis) or
	"fourier" (Fourier basis, default).}

	\item{kg}{the number of basis functions to be used for \code{GVF} (default = 30, has no effect
	under \code{GVF = FALSE}).}

	\item{kb}{the number of basis functions to be used for \code{BSF} (default = 25).}

	\item{order}{a polynomial order to be used in "bspline". Default = 4 corresponds to the cubic B-splines.
	as no effect if only Fourier bases are used.}

	\item{stat}{the statistic to be used to calculate the P values. One of "F" (default), "Chisq", "LRT".}

	\item{flip.genotypes}{a logical value indicating whether the genotypes of some genetic variants should be
	flipped (relabeled) for their better functional representation [Vsevolozhskaya, et al., 2014]. Default = FALSE.}

	\item{impute.method}{a method for imputation of missing genotypes. It can be either "mean" (default)
	or "blue". If "mean" the genotypes will be imputed by the simple mean values. If "blue"
	the best linear unbiased estimates (BLUEs) of mean genotypes will be calculated
	taking into account the relationships between individuals [McPeek, et al., 2004,
	DOI: 10.1111/j.0006-341X.2004.00180.x] and used for imputation.}

	\item{write.file}{output file name to write results as they come (sequential mode only).}

	\item{...}{other arguments that could be passed to \code{null()}, \code{read.plink()}\cr
	and \code{readVCFToMatrixByGene()}.}
}
\details{
	The test assumes that the effects of multiple genetic variants
	(and also their genotypes if GVFs are used) can be described
	as a continuous function, which can be modelled through B-spline
	or Fourier basis functions. When the number of basis functions
	(set by \eqn{Kg} and \eqn{Kb}) is less than the number of variants
	within the region, the famFLM test may have an advantage of using
	less degrees of freedom [Svishcheva, et al., 2015].\cr

	Several restrictions exist in combining B-spline or Fourier bases
	for construction of GVFs and BSF [Svishcheva, et al., 2015], and
	the famFLM function takes them into account. Namely:\cr

	1) \eqn{m \geq Kg \geq Kb}, where \eqn{m} is the number of polymorphic
	genetic variants within a region.\cr

	2) Under \eqn{Kg = Kb}, B-B and B-F models are equivalent
	to 0-B model, and F-F and F-B models are equivalent to 0-F model.
	0-B and 0-F models will be used for these cases, respectively.\cr

	3) Under \eqn{m = Kb}, 0-B and 0-F models are equivalent to a
	standard multiple linear regression, and it will be used for these cases.\cr

	4) When Fourier basis is used, the number of basis functions should be
	an odd integer. Even values will be changed accordingly.\cr

	Because of these restrictions, the model in effect may not always
	be the same as it has been set. The ultimate model name is returned in
	results in the "model" column (see below).
	
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
	and informative polymorphic variants for each of analyzed regions.
	It also contains the names of the functional models used for each region
	(it may not always coincide with what was set, because of restrictions described
	in Details section). The first part of the name relates to the functional basis
	of GVFs and the second one to that of BSF, e.g. "F30-B25" means that 30 Fourier
	basis functions were used for construction of GVFs and 25 B-spline basis functions
	were used for construction of BSF. "0-F25" means that genotypes were not smoothed
	and 25 Fourier basis functions were used for beta-smooth. "MLR" means that standard
	multiple linear regression was applied.}

	\item{nullmod}{an object containing the estimates of the null model parameters: heritability (h2),
	total variance (total.var), estimates of fixed effects of covariates (alpha), the gradient (df), and the total log-likelihood (logLH).}

	\item{sample.size}{the sample size after omitting NAs.}

	\item{time}{If \code{return.time = TRUE} a list with running times for null model, regional analysis and total analysis is returned. See \code{proc.time()} for output format.}
}
\references{
	Svishcheva G.R., Belonogova N.M. and Axenovich T.I. (2015) Region-based association test for familial data under functional linear models. PLoS ONE 10(6): e0128999.\cr
	Vsevolozhskaya O.A., et al. (2014) Functional Analysis of Variance for Association Studies. PLoS ONE 9(9): e105074.\cr
	Wu M.C., et al. (2011) Rare-variant association testing for sequencing data with the sequence kernel association test. Am. J. Hum. Genet., Vol. 89, P. 82-93.
	}
\examples{

data(example.data)

## Run famFLM with sliding window (default):
out <- famFLM(trait ~ age + sex, phenodata, genodata, kin,
	positions = snpdata$position)

## Run famFLM with regions defined in snpdata$gene and with
## null model parameters obtained in the first run:
out <- famFLM(trait ~ age + sex, phenodata, genodata, kin,
	out$nullmod, positions = snpdata$position,
	regions = snpdata$gene)

## Run famFLM parallelized on two cores (this will require
## 'foreach' and 'doParallel' R-packages installed and
## cores available):
out <- famFLM(trait ~ age + sex, phenodata, genodata, kin,
	out$nullmod, positions = snpdata$position, ncores = 2)

## Run MLR with genotypes in VCF format:
VCFfileName <- system.file(
	"testfiles/1000g.phase1.20110521.CFH.var.anno.vcf.gz",
	package = "FREGAT")
geneFile <- system.file("testfiles/refFlat_hg19_6col.txt.gz",
	package = "FREGAT")
phe <- data.frame(trait = rnorm(85))
out <- famFLM(trait, phe, VCFfileName, geneFile = geneFile,
	reg = "CFH", annoType = "Nonsynonymous",
	flip.genotypes = TRUE)

## Run famFLM with genotypes in PLINK binary data format:
bedFile <- system.file("testfiles/sample.bed",
	package = "FREGAT")
data <- read.plink(bedFile)
phe <- data.frame(trait = rnorm(120))
out <- famFLM(trait, phe, bedFile, positions = data$map$position)
}
