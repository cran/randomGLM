\name{thinRandomGLM}
\alias{thinRandomGLM}
\title{Random generalized linear model predictor thinning}
\description{ This function allows the user to define a "thinned" version of a random  generalized linear
model predictor by focusing on those features that occur relatively frequently. 
}

\usage{
thinRandomGLM(rGLM, threshold)
}

\arguments{
  \item{rGLM}{a \code{randomGLM} object such as one returned by \code{\link{randomGLM}}. }

  \item{threshold}{ integer specifying the minimum of times a feature was selected across the bags in
\code{rGLM} for the feature to be kept. Note that only features selected \code{threshold +1} times and more
are retained. For the purposes of this count, appearances in interactions are not
counted. Features that appear 
\code{threshold} times or fewer are removed from the underlying regression models when the models are re-fit.}

}
\details{
The function "thins out" (reduces) a previously-constructed random generalized linear model predictor by
removing rarely selected features and refitting each (generalized) linear model (GLM). 
Each GLM (per bag) is refit using only those
features that occur more than \code{threshold} times across the \code{nBags} number of bags. The
occurrence count excludes interactions (in other words, the threshold will be applied to the first row of
\code{timesSelectedByForwardRegression}).
}

\value{
The function returns a valid \code{randomGLM} object (see \code{\link{randomGLM}} for details) that can be
used as input to the predict() method (see \code{\link{predict.randomGLM}}). The returned object contains a
copy of  the input \code{rGLM} in which the following components were modified:

 \item{predictedOOB}{the updated continuous prediction (if \code{classify} is \code{FALSE}) 
or predicted classification (if \code{classify} is \code{TRUE}) of the input data based on out-of-bag 
samples.
}
  \item{predictedOOB.response}{In case of a binary outcome, the updated predicted probability of each
outcome
specified by \code{y} based on out-of-bag samples. In case of a continuous outcome, this is the predicted
value based on out-of-bag samples (i.e., a copy of \code{predictedOOB}).}

  \item{featuresInForwardRegression}{features selected by forward selection in each bag. A list with one
component per bag. Each component
is a matrix with \code{maxInteractionOrder} rows.
Each column represents one
interaction obtained by multiplying the features indicated by the entries in each column (0 means no
feature, i.e. a lower order interaction).  }

  \item{coefOfForwardRegression}{coefficients of forward regression. A list with one
component per bag. Each component is a vector giving the coefficients of the model determined by forward
selection in the corresponding bag. The order of the coefficients is the same as the order of the terms in
the corresponding component of \code{featuresInForwardRegression}. }

 \item{interceptOfForwardRegression}{a vector with one component per bag giving the intercept of the
regression model in each bag.}

  \item{timesSelectedByForwardRegression}{a matrix of \code{maxInteractionOrder} rows and number of features
columns. Each entry gives the number of times the corresponding feature appeared in a predictor model at the
corresponding order of interactions. Interactions where a single feature enters more than once (e.g., a
quadratic interaction of the feature with itself) are counted once.
}
  \item{models}{the "thinned" regression models for each bag. }

}


  
\references{
Lin Song, Peter Langfelder, Steve Horvath: Random generalized linear model: a highly accurate and interpretable ensemble predictor. BMC Bioinformatics (2013)
}

\author{Lin Song, Steve Horvath, Peter Langfelder}

\examples{

## binary outcome prediction
# data generation
data(iris)
# Restrict data to first 100 observations
iris=iris[1:100,]
# Turn Species into a factor
iris$Species = as.factor(as.character(iris$Species))
# Select a training and a test subset of the 100 observations
set.seed(1)
indx = sample(100, 67, replace=FALSE)
xyTrain = iris[indx,]
xyTest = iris[-indx,]
xTrain = xyTrain[, -5]
yTrain = xyTrain[, 5]

xTest = xyTest[, -5]
yTest = xyTest[, 5]

# predict with a small number of bags - normally nBags should be at least 100.
RGLM = randomGLM(xTrain, yTrain, nCandidateCovariates=ncol(xTrain), nBags=30, keepModels = TRUE, nThreads = 1)
table(RGLM$timesSelectedByForwardRegression[1, ])
# 0  7 23 
# 2  1  1 

thinnedRGLM = thinRandomGLM(RGLM, threshold=7)
predicted = predict(thinnedRGLM, newdata = xTest, type="class")
predicted = predict(RGLM, newdata = xTest, type="class")

}
\keyword{misc}
