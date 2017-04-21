\name{FFBSKAT}
\alias{FFBSKAT}
\title{Fast Family-Based SKAT}
\description{
A fast regional association analysis in related or population  samples
}
\usage{
FFBSKAT(formula, phenodata, genodata, kin = NULL, nullmod,
regions = NULL, sliding.window = c(20, 10), mode = "add",
ncores = 1, return.time = FALSE, kernel = "linear.weighted",
beta.par = c(1, 25), weights = NULL, method = "kuonen",
acc = 1e-8, lim = 1e+6, return.variance.explained = FALSE,
reml = TRUE, flip.genotypes = FALSE, impute.method = 'mean',
rho = FALSE, write.file = FALSE, ...)
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

	\item{kernel}{one of "linear.weighted" (default), "quadratic", "IBS", "IBS.weighted", "2wayIX".
	See Details for "linear.weighted" kernel description and [Wu, et al., 2011] for other kernel types.
	"2wayIX" kernel considers SNP-SNP interaction terms along with main effects. For
	"linear.weighted" and "IBS.weighted" kernels, weights can be varied by defining
	\code{weights} or \code{beta.par}.}

	\item{beta.par}{two positive numeric shape parameters in the beta distribution to assign weights 
	for each SNP in weighted kernels (see Details). Default = c(1, 25) is recommended for analysis
	of rare variants. Has no effect for unweighted kernels or if \code{weights} are defined.}

	\item{weights}{a numeric vector or a function of minor allele frequency (MAF) to assign weights
	for each genetic variant in the weighted kernels. Has no effect if one of unweighted kernels
	was chosen. If NULL, the weights will be calculated using the beta distribution (see Details).}

	\item{method}{either "kuonen" or "davies". Method for computing the P value (see Details).
	Default = "kuonen".}

	\item{acc}{accuracy parameter for "davies" method.}

	\item{lim}{limit parameter for "davies" method.}

	\item{return.variance.explained}{a logical value indicating whether the (marginal) variance explained by
	each region should be returned. Default = FALSE for faster performance.}

	\item{reml}{a logical value indicating whether the restricted maximum likelihood should be used to estimate
	the variance explained by each region. Default = TRUE for faster performance.\cr
	Has no effect if \code{return.variance.explained = FALSE}.}

	\item{flip.genotypes}{a logical value indicating whether the genotypes of some genetic variants should be
	flipped (relabeled) to ensure that all minor allele frequencies (MAFs) < 0.5. Default = FALSE, with warning of any MAF > 0.5.}

	\item{impute.method}{a method for imputation of missing genotypes. It can be either "mean" (default)
	or "blue". If "mean" the genotypes will be imputed by the simple mean values. If "blue"
	the best linear unbiased estimates (BLUEs) of mean genotypes will be calculated
	taking into account the relationships between individuals [McPeek, et al., 2004,
	DOI: 10.1111/j.0006-341X.2004.00180.x] and used for imputation.}

	\item{rho}{If TRUE the optimal test is used [Lee, et al., 2012]. \code{rho} can be a vector of grid
	values from 0 to 1. The default grid is (0 : 10) / 10.}

	\item{write.file}{output file name to write results as they come (sequential mode only).}

	\item{...}{other arguments that could be passed to \code{null()}, \code{read.plink()}\cr
	and \code{readVCFToMatrixByGene()}.}
}
\details{
	By default, FFBSKAT uses the linear weighted kernel function to set the inter-individual
	similarity matrix \eqn{K = GWWG^T}, where \eqn{\mathit{G}} is the \eqn{\mathit{n\times p}}
	genotype matrix for \eqn{\mathit{n}} individuals and \eqn{\mathit{p}} genetic variants in the region, 
	and \eqn{\mathit{W}} is the \eqn{\mathit{p\times p}} diagonal weight matrix. Given the shape parameters
	of the beta function, \code{beta.par = c(a, b)}, 
	the weights are defined using probability density function of the beta distribution:\cr
	\cr
	\eqn{W_{i}=(B(a,b))^{^{-1}}MAF_{i}^{a-1}(1-MAF_{i})^{b-1} },\cr
	\cr
	where \eqn{MAF_{i}} is a minor allelic frequency for the \eqn{i^{th}} genetic variant in the region,
	which is estimated from genotypes, and \eqn{B(a,b)} is the beta function. This way of defining weights
	is the same as in original SKAT (see [Wu, et al., 2011] for details). \code{beta.par = c(1, 1)} corresponds
	to the unweighted SKAT.
	The formula: \cr
	\cr
	\eqn{Q=0.5\tilde{y}^{T}\Omega^{-1}K\Omega^{-1}\tilde{y}}\cr
	\cr
	is used to calculate score statistic, where \eqn{\tilde{y}} and \eqn{\Omega} are environmental
	residuals and covariance matrix obtained under the null hypothesis, respectively. Depending on
	the method option chosen, either Kuonen or Davies method is used to calculate P values from the
	score statistic Q. Both an Applied Statistics algorithm that inverts the characteristic function
	of the mixture chisq [Davies, 1980] and a saddlepoint approximation [Kuonen, 1999] are nearly exact,
	with the latter usually being a bit faster.
	For other kernel types, see [Wu, et al., 2011].
}
\value{
	A list with values:

	\item{results}{a data frame containing P values, numbers of variants and polymorphic variants for each of analyzed regions.\cr
	If \code{return.variance.explained = TRUE} it contains also the column with marginal amounts of variance explained by each region. If \code{reml = FALSE} the new estimates of heritability (h2) and total variance with corresponding total log-likelihood are also returned.}

	\item{nullmod}{an object containing the estimates of the null model parameters: heritability (h2), total variance (total.var),
	estimates of fixed effects of covariates (alpha), the gradient (df), and the total log-likelihood (logLH).}

	\item{sample.size}{the sample size after omitting NAs.}

	\item{time}{If \code{return.time = TRUE} a list with running times for null model, regional analysis and total analysis is returned. See \code{proc.time()} for output format.}
}
\references{
	Svishcheva G.R., Belonogova N.M. and Axenovich T.I. (2014) FFBSKAT: Fast Family-Based Sequence Kernel Association Test. PLoS ONE 9(6): e99407. doi:10.1371/journal.pone.0099407\cr
	Davies R.B. (1980) Algorithm AS 155: The Distribution of a Linear Combination of chi-2 Random Variables, Journal of the Royal Statistical Society. Series C (Applied Statistics), Vol. 29, N 3, P. 323-333.\cr
	Kuonen D. (1999) Saddlepoint Approximations for Distributions of Quadratic Forms in Normal Variables. Biometrika, Vol. 86, No. 4, P. 929-935.\cr
	Wu M.C., et al. (2011) Rare-variant association testing for sequencing data with the sequence kernel association test. Am. J. Hum. Genet., Vol. 89, P. 82-93.\cr
	Lee S., et al. (2012) Optimal unified approach for rare variant association testing with application to small sample case-control whole-exome sequencing studies. American Journal of Human Genetics, 91, 224-237.
	}
\examples{

data(example.data)  

## Run FFBSKAT with sliding window (default):
out <- FFBSKAT(trait ~ age + sex, phenodata, genodata, kin)

## Run FFBSKAT with regions defined in snpdata$gene and with
## null model obtained in first run:
out <- FFBSKAT(trait ~ age + sex, phenodata, genodata, kin,
out$nullmod, regions = snpdata$gene)

## Run FFBSKAT parallelized on two cores (this will require
## 'foreach' and 'doParallel' R-packages installed and
## cores available):
out <- FFBSKAT(trait ~ age + sex, phenodata, genodata, kin,
	out$nullmod, ncores = 2)

## Run FFBSKAT with genotypes in VCF format:
VCFfileName <- system.file(
	"testfiles/1000g.phase1.20110521.CFH.var.anno.vcf.gz",
	package = "FREGAT")
geneFile <- system.file("testfiles/refFlat_hg19_6col.txt.gz",
	package = "FREGAT")
phe <- data.frame(trait = rnorm(85))
out <- FFBSKAT(trait, phe, VCFfileName, geneFile = geneFile,
	reg = "CFH", annoType = "Nonsynonymous",
	flip.genotypes = TRUE)

## Run FFBSKAT with genotypes in PLINK binary data format:
bedFile <- system.file("testfiles/sample.bed",
	package = "FREGAT")
phe <- data.frame(trait = rnorm(120))
out <- FFBSKAT(trait, phe, bedFile)

}
