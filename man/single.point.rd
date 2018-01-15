\name{single.point}
\alias{single.point}
\title{A single-point association test}
\description{
A single-point association test for familial or population data
}
\usage{
single.point(formula, phenodata, genodata, kin = NULL, nullmod,
regions = NULL, mode = "add", ncores = 1, return.time = FALSE,
impute.method = 'mean', write.file = FALSE, ...)
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
	- an object of \code{gwaa.data} or \code{snp.data} class (this will require
	\code{GenABEL} R-package to be installed).}

	\item{kin}{a square symmetric matrix giving the pairwise kinship coefficients between analyzed
	individuals. Under default \code{kin = NULL} all individuals will be considered as unrelated.}

	\item{nullmod}{an object containing parameter estimates under the null model. Setting \code{nullmod}
	allows to avoid re-estimation of the null model that does not depend on genotypes and can be
	calculated once for a trait. If not set, the null model parameters will be estimated within the function.
	The \code{nullmod} object in proper format can be obtained by \code{null.model()} function
	or any analysis function in \code{FREGAT}.}

	\item{regions}{an object assigning regions for genetic variants in \code{genodata} can be either\cr
	- a vector of length equal to the number of genetic variants assigning the region
	for each variant (see Examples).\cr
	- a data frame / matrix with names of genetic variants in the first column and names of regions
	in the second column (this format allows overlapping regions).\cr
	If NULL, all variants in \code{genodata} will be analysed.}

	\item{mode}{the mode of inheritance: "add", "dom" or "rec" for additive, dominant or recessive mode,
	respectively. For dominant (recessive) mode genotypes will be recoded as AA = 0, Aa = 1 and aa = 1
	(AA = 0, Aa = 0 and aa = 1), where a is a minor allele. Default mode is additive.}

	\item{ncores}{number of CPUs for parallel calculations. Default = 1.}

	\item{return.time}{a logical value indicating whether the running time should be returned.}

	\item{impute.method}{a method for imputation of missing genotypes. It can be either "mean" (default)
	or "blue". If "mean" the genotypes will be imputed by the simple mean values. If "blue"
	the best linear unbiased estimates (BLUEs) of mean genotypes will be calculated
	taking into account the relationships between individuals [McPeek, et al., 2004,
	DOI: 10.1111/j.0006-341X.2004.00180.x] and used for imputation.}

	\item{write.file}{output file name. If specified, output (as it proceeds) will be written 
	to the file.}

	\item{...}{other arguments that could be passed to \code{null()} and \code{read.plink()}.}
}
\value{
	A list with values:

	\item{results}{a data frame containing P values, estimates of betas and their s.e., MAFs and sample size for each of analyzed variants.}

	\item{nullmod}{an object containing the estimates of the null model parameters: heritability (h2),
	total variance (total.var), estimates of fixed effects of covariates (alpha), the gradient (df), and the total log-likelihood (logLH).}

	\item{time}{If \code{return.time = TRUE} a list with running times for null model, regional analysis and total analysis is returned. See \code{proc.time()} for output format.}
}
\examples{

data(example.data)

## Run single.point analysis for all variants in genodata (default):
out <- single.point(trait ~ age + sex, phenodata, genodata, kin)

## Run single.point analysis for regions defined in snpdata$gene and with
## null model parameters obtained in the first run:
out <- single.point(trait ~ age + sex, phenodata, genodata, kin,
	out$nullmod, regions = snpdata$gene)

## Run single.point analysis parallelized on two cores (this will require
## 'foreach' and 'doParallel' R-packages installed and
## cores available):
out <- single.point(trait ~ age + sex, phenodata, genodata, kin,
	out$nullmod, ncores = 2)

## Run single.point analysis with genotypes in PLINK binary data format:
bedFile <- system.file("testfiles/sample.bed",
	package = "FREGAT")
phe <- data.frame(trait = rnorm(120))
out <- single.point(trait, phe, bedFile)
}
