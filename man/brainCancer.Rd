\name{brainCancer}
\alias{brainCancer}
\docType{data}
\title{The brain cancer data set}
\description{2 sets containing the gene expression profiles of 55 and 65 brain cancer patients respectively.
}

\usage{
data(brainCancer)
}

\format{
\code{brainCancer} is a list of 2 components: train and test. "train" is a numeric matrix with 55 samples (rows) across 5000 genes (columns). "test" is a numeric matrix with 65 samples (rows) across the same 5000 genes (columns).     
}

\source{Horvath S, Zhang B, Carlson M, Lu K, Zhu S, Felciano R, Laurance M, Zhao W, Shu Q, Lee Y, Scheck A,
Liau L, Wu H, Geschwind D, Febbo P, Kornblum H, TF C, Nelson S, Mischel P: Analysis of Oncogenic Signaling
Networks in Glioblastoma Identifies ASPM as a Novel Molecular Target. Proc Natl Acad Sci U S A 2006,
103(46):17402-7.} 
  
\references{
Lin Song, Peter Langfelder, Steve Horvath: Random generalized linear model: a highly accurate and interpretable ensemble predictor. BMC Bioinformatics (2013)
}

\author{Lin Song, Steve Horvath}

\examples{
data(brainCancer)
}
