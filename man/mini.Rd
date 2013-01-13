\name{mini}
\alias{mini}
\docType{data}
\title{An example data set derived from the brain cancer data set}
\description{This example contains one training set, one test set, a corresponding binary outcome and a corresponding continuous outcome. Outcomes are gene traits derived from the brain cancer data set.
}

\usage{
data(mini)
}


\format{
\code{mini} is a list of 6 components: x, xtest, yB, yBtest, yC and yCtest. "x" is a numeric matrix with 55 samples (rows) across 4999 genes (columns). "xtest" is a numeric matrix with 65 samples (rows) across the same 4999 genes (columns). They are subsets of the original \code{brainCancer} data. One gene is left out. "yC" and "yCtest" are continuous outcomes equal to the expression values of the left out gene in the training and test set respectively. "yB" and "yBtest" are binary outcomes dichotomized at the median based on "yC" and "yCtest".      
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
data(mini)
}
