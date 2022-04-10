\name{predict.randomGLM}
\alias{predict.randomGLM}
\title{Prediction from a random generalized linear model predictor }
\description{ Implements a predict method on a previously-constructed random  generalized linear model predictor and new data.
}

\usage{
\method{predict}{randomGLM}(object, newdata, type=c("response", "class"), 
                 thresholdClassProb = object$details$thresholdClassProb, ...)
}

\arguments{
  \item{object}{a \code{randomGLM} object such as one returned by \code{\link{randomGLM}}. }

  \item{newdata}{specification of test data for which to calculate the prediction.}

  \item{type}{type of prediction required. Type "response" gives the fitted probabilities for classification, the fitted values for regression. Type "class" applies only to classification, and produces the predicted class labels.}

  \item{thresholdClassProb}{the threshold of predictive probabilities to arrive at classification. Takes values between 0 and 1. Only used for binary outcomes.}

  \item{...}{other arguments that may be passed to and from methods. Currently unused.}
}

\details{
The function calculates prediction on new test data. It only works if \code{object} contains the regression
models that were used to construct the predictor (see argument \code{keepModels} of the function
\code{\link{randomGLM}}). 

If the predictor was trained on a multi-class response, the prediction is applied to each of the
representing binary variables (see \code{\link{randomGLM}} for details).
}

\value{
For continuous prediction, the predicted values. For classification of binary response, predicted class when
\code{type="class"}; or a two-column matrix giving the class probabilities if \code{type="response"}.

If the predictor was trained on a multi-class response, the returned value is a matrix of "cbind"-ed results
for the representing individual binary variables (see \code{\link{randomGLM}} for details).
}

\references{
Lin Song, Peter Langfelder, Steve Horvath: Random generalized linear model: a highly accurate and interpretable ensemble predictor. BMC Bioinformatics (2013)
}

\author{Lin Song, Steve Horvath and Peter Langfelder.
}

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

# predict with a small number of bags 
# - normally nBags should be at least 100.
RGLM = randomGLM(
   xTrain, yTrain, 
   nCandidateCovariates=ncol(xTrain), 
   nBags=30, 
   keepModels = TRUE, nThreads = 1)

predicted = predict(RGLM, newdata = xTest, type="class")
table(predicted, yTest)

## continuous outcome prediction

x=matrix(rnorm(100*20),100,20)
y=rnorm(100)

xTrain = x[1:50,]
yTrain = y[1:50]
xTest = x[51:100,]
yTest = y[51:100]

RGLM = randomGLM(
   xTrain, yTrain, 
   classify=FALSE, 
   nCandidateCovariates=ncol(xTrain), 
   nBags=10, 
   keepModels = TRUE, nThreads = 1)

predicted = predict(RGLM, newdata = xTest)
}
\keyword{misc}
