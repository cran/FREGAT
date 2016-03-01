\name{example.data}
\alias{example.data}
\alias{genodata}
\alias{phenodata}
\alias{snpdata}
\alias{kin}
\title{A small example data set}
\description{
	\code{genodata}
	A matrix containing genotypes of 50 genetic variants (given in columns) in 66 individuals
	(given in rows). Three genotypes are coded as 0, 1 and 2.

	\code{phenodata}
	A data frame containing "trait", "sex" and "age" columns: a quantitative trait to be analyzed
	and its covariates.

	\code{snpdata}
	A data frame with descriptive information on 50 genetic variants in \code{genodata}.
	The important column is "gene": it assigns each variant to a certain gene region.

	\code{kin}
	A kinship matrix for the 66 individuals.
}
\usage{
	data(example.data)
}
