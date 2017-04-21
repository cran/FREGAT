\name{null}
\alias{null}
\title{Fitting the null model}
\description{
Gives estimation of model parameters under the null hypothesis
}
\usage{
null(formula, phenodata, kin = NULL, opt.method = 'optimize',
ih2 = 0.3, eps = 1.e-04)
}

\arguments{
	\item{formula}{referring to the column(s) in \code{phenodata} to be analyzed as outcome and,
	if needed, covariates.}

	\item{phenodata}{a data frame containing columns mentioned in \code{formula}: trait to analyze and,
	if needed, covariates. Individuals not measured for trait or covariates will be omitted.}

	\item{kin}{a square symmetric matrix giving the pairwise kinship coefficients between analyzed
	individuals. Under default \code{kin = NULL} all individuals will be considered as unrelated.}

	\item{opt.method}{optimization method, one of "optimize" (default), "optim" and "nlminb". Corresponding
	R functions will be used in optimization.}

	\item{ih2}{initial value for h2. Default = 0.3.}

	\item{eps}{epsilon (precision) value for optimization. Default = 1e-04.}

}
\details{
	The function performs one-dimensional optimization for h2 with analytical calculations for the other parameters.
}
\value{
	A list with values:

	\item{h2}{estimate of the heritability}

	\item{total.var}{estimate of the total variance}

	\item{alpha}{estimates of fixed effects of covariates}

	\item{df}{the gradient}

	\item{logLH}{the total log-likelihood}

	\item{p.normality}{p value of Shapiro-Wilk normality test for the null model residuals.}
}
\examples{

data(example.data)

## Run the null model:
nullmod <- null(trait ~ age + sex, phenodata, kin)

## SKAT with the null model object obtained in the first run:
out <- FFBSKAT(trait ~ age + sex, phenodata, genodata, kin, nullmod)

}
