\name{randomGLM}
\alias{randomGLM}
\title{Random generalized linear model predictor}
\description{
Ensemble predictor comprised of individual generalized linear model predictors. 
}

\usage{
randomGLM(
  # Input data
  x, y, xtest = NULL,
  weights = NULL,

  # Which columns in x are categorical?
  categoricalColumns = NULL,
  maxCategoricalLevels = 2,

  # Include interactions?
  maxInteractionOrder = 1,
  includeSelfinteractions = TRUE,

  # Prediction type: type can be used to set 
  # the prediction type in a simplified way...
  type = c("auto", "linear", "binary", "count", "general", "survival"),

  # classify is retained mostly for backwards compatibility
  classify = switch(type, 
    auto = !is.Surv(y) & (is.factor(y) | length(unique(y)) < 4),
    linear = FALSE,
    binary = TRUE ,
    count = FALSE,
    general = FALSE,
    survival = FALSE),

  # family can be used to fine-tune the underlying regression model
  family = switch(type, 
    auto = NULL,
    linear = gaussian(link="identity"),
    binary = binomial(link=logit),
    count = poisson(link = "log"),
    general = NULL,
    survival = NULL),

  # Multi-level classification options - only apply to classification 
  # with multi-level response
  multiClass.global = TRUE,
  multiClass.pairwise = FALSE,
  multiClass.minObs = 1,
  multiClass.ignoreLevels = NULL,

  # Sampling options
  nBags = 100,
  replace = TRUE,
  sampleBaggingWeights = NULL,
  nObsInBag = if (replace) nrow(x) else as.integer(0.632 * nrow(x)),
  nFeaturesInBag = ceiling(ifelse(ncol(x)<=10, ncol(x),
        ifelse(ncol(x)<=300, (1.0276-0.00276*ncol(x))*ncol(x), ncol(x)/5))),
  minInBagObs = min( max( nrow(x)/2, 5), 2*nrow(x)/3),
  maxBagAttempts = 100*nBags,
  replaceBadBagFeatures = TRUE,

  # Individual ensemble member predictor options
  nCandidateCovariates=50,
  corFncForCandidateCovariates= cor,
  corOptionsForCandidateCovariates = list(method = "pearson", use="p"),
  mandatoryCovariates = NULL,
  interactionsMandatory = FALSE,
  keepModels = is.null(xtest),

  # Miscellaneous options
  thresholdClassProb = 0.5,
  interactionSeparatorForCoefNames = ".times.",
  randomSeed = 12345,
  nThreads = NULL,
  verbose =0 )
}

\arguments{
  \item{x}{a matrix whose rows correspond to observations (samples) and whose columns correspond to features (also known as covariates or variables). Thus, \code{x} contains the training data sets.}

  \item{y}{outcome variable corresponding to the rows of \code{x}: at this point, one can either use a binary class outcome (factor variable) or a quantitative outcome (numeric variable). }

  \item{xtest}{an optional matrix of a second data set (referred to as test data set while the data in
\code{x} are interpreted as training data). The number of rows can (and typically will) be different from
the number of rows in \code{x}.} 

  \item{weights}{optional specifications of sample weights for the regression models. Not to be confused
with \code{sampleBaggingWeights} below.}

  \item{categoricalColumns}{optional specifications of columns that are to be treated as categorical. If
not given, columns with at most \code{maxCategoricalLevels} (see below) unique values will be considered
categorical.}

  \item{maxCategoricalLevels}{columns with no more than this number of unique values will be considered
categorical.}

  \item{maxInteractionOrder}{ integer specifying the maximum interaction level. The default is to have no
interactions; numbers higher than 1 specify interactions up to that order. For example, 3 means quadratic
and cubic interactions will be included. Warning: higher order interactions greatly increase the computation
time. We see no benefit of using maxInteractionOrder>2.} 

  \item{includeSelfinteractions}{logical: should self-interactions be included?}

  \item{type}{character string specifying the type of the response variable. Recognized values are (unique
abbreviations of) 
\code{"auto"}, \code{"linear"}, \code{"binary"}, \code{"count"}, \code{"general"}, and \code{"survival"}. See
Details for what the individual types mean.
}

  \item{classify}{logical indicating whether the response is a categorical variable. This argument is
present mainly for backwards compatibility; please use \code{type} above to specify the type of the response
variable. If \code{TRUE} the response \code{y} will be interpreted as a binary variable and
logistic regression will be used. If \code{FALSE} the response \code{y} will be interpreted as a
quantitative numeric variable and a least squares regression model will be used to arrive at base learners.
Multi-level classification is split into a series of binary classification problems according to the
\code{multiClass...} arguments described below.}

  \item{family}{Specification of family (see \code{\link[stats]{family}}) for general linear model fitting
(see \code{\link[stats]{glm}}). Default values are provided for most of the specific types but can be
overriden here (for example, if a different link function is desired). There is no default value for
\code{type = "general"} and the user must specify a valid family. In contrast, this argument must be
\code{NULL} when \code{type = "survival"}. }

  \item{multiClass.global}{for multi-level classification, this logical argument controls whether binary
variables of the type "level vs. all others" are included in the series of binary variables to which
classification is applied. }

  \item{multiClass.pairwise}{for multi-level classification, this logical argument controls whether binary 
variables of the type "level A vs. level B" are included in the series of binary variables to which 
classification is applied. }

  \item{multiClass.minObs}{an integer specifying the minimum number of observations for each level for the
level to be considered when creating "level vs. all" and "level vs. level" binary variables. }

  \item{multiClass.ignoreLevels}{optional specifications of the values (levels) of the input response
\code{y} that are to be ignored when constructing level vs. all and level vs. level binary responses. Note
that observation with these values will be included in the "all" but will not have their own "level vs. all"
variables.}

  \item{nBags}{number of bags (bootstrap samples) for defining the ensemble predictor, i.e. this also
corresponds to the number of individual GLMs.} 

  \item{replace}{logical. If \code{TRUE} then each bootstrap sample (bag) is defined by sampling with
replacement. Otherwise, sampling is carried out without replacement. We recommend to choose \code{TRUE}. } 
 
  \item{sampleBaggingWeights}{weights assigned to each observations (sample) during bootstrap sampling. Default \code{NULL} corresponds to equal weights.}

  \item{nObsInBag}{number of observations selected for each bag. Typically, a bootstrap sample (bag) has the
same number of observations as in the original data set (i.e. the rows of \code{x}).}  

  \item{nFeaturesInBag}{number of features randomly selected for each bag. Features are randomly selected
without replacement. If there are no interaction terms, then this number should be smaller than or equal to
the number of rows of \code{x}.}  

  \item{minInBagObs}{minimum number of unique observations that constitute a valid bag. If the sampling
produces a bag with fewer than this number of unique observations, the bag is discarded and re-sampled again
until the number of unique observations is at least \code{minInBagObs}. This helps prevent too few unique
observations in a bag which would lead to problems with model selection. }

  \item{maxBagAttempts}{Maximum number of bagging attempts.}
 
  \item{replaceBadBagFeatures}{If a feature in a bag contains missing data, should it be replaced?}

  \item{nCandidateCovariates}{Positive integer. The number of features that are being considered for forward
selection in each GLM (and in each bag). For each bag, the covariates are being chosen according their
highest absolute correlation with the outcome. In case of a binary outcome, it is first turned into a binary
numeric variable.  } 

  \item{corFncForCandidateCovariates}{the correlation function used to select candidate covariates. Choices
include \code{\link[stats]{cor}} or biweight midcorrelation, \code{bicor}, implemented in the
package WGCNA. The biweight mid-correlation is a robust alternative to the Pearson correlation.}

  \item{corOptionsForCandidateCovariates}{list of arguments to the correlation function. Note that robust correlations are sometimes problematic for binary class outcomes. When using the robust
correlation \code{bicor}, use the argument \code{"robustY=FALSE"}.}

  \item{mandatoryCovariates}{indices of features that are included as mandatory covariates in each GLM model. The default is no mandatory features. This allows the user to "force" variables into each GLM.} 

  \item{interactionsMandatory}{logical: should interactions of mandatory covariates be mandatory as well? Interactions are only included up to the level specified in \code{maxInteractionOrder}. }

  \item{keepModels}{logical: should the regression models for each bag be kept? The models are necessary for future predictions using the \code{predict} function, predict() generic.}

  \item{thresholdClassProb}{number in the interval [0,1]. Recommended value 0.5. This parameter is only
relevant in case of a binary outcome, i.e. for a logistic regression model. Then this threshold will be
applied to the predictive class probabilities to arrive at binary outcome (class outcome). }

 \item{interactionSeparatorForCoefNames}{a character string that will be used to separate feature names when
forming names of interaction terms. This is only used when interactions are actually taken into account (see
\code{maxInteractionLevel} above) and only affects coefficient names in models and columns names in returned
\code{featuresInForwardRegression} (see output value below). We recommend setting it so the interaction
separator does not conflict with any feature name since this may improve interpretability of the results.}


  \item{randomSeed}{NULL or integer. The seed for the random number generator. If NULL, the seed will not be
set. If non-NULL and the random generator has been initialized prior to the function call, the latter's
state is saved and restored upon exit.} 

  \item{nThreads}{number of threads (worker processes) to perform the calculation. If not given, will be
determined automatically as the number of available cores if the latter is 3 or less, 
and number of cores minus 1 if the number of available cores is 4 or more. Invalid entries (missing value,
zero or negative values etc.) are changed to 1, with a warning.}

  \item{verbose}{value 0 or 1 which determines the level of verbosity. Zero means silent, 1 reports the bag number the function is working on. At this point verbose output only works if \code{nThreads}=1}
}

\details{

The function \code{randomGLM} can be used to predict a variety of different types of outcomes. The outcome
type is specified by the argument \code{type} as follows:

If \code{type = "auto"}, the function will attempt to determine the response type automatically. If the response
is not a \code{\link[survival]{Surv}} object and is a factor or has 3 or fewer unique values, 
it is assumed to be categorical. If the number of unique values is 2, the function uses logistic regression;
for categorical responses with more 3 or more possible values, see below.

If \code{type = "linear"}, the responses is assumed to be numeric with Gaussian errors, and the function
will use linear regression.

If \code{type = "binary"}, the response is assumed to be categorical (at this point not necessarily binary
but that may change in the future). If the response has 2 levels, logistic regression (binomial family with
the logit link) is used. If the response has more than 2 levels, see below.

If \code{type = "count"}, the response is assumed to represent counts with Poisson-distributed errors, and
Poisson regression (poisson family with the logarithmic link) is used.

If \code{type = "general"}, the function does not make assumption about the response type and the user must
specify an appropriate family (see \code{\link{stats}{family}}).

If \code{type = "survival"}, the function assumes the response is a censored time, that is a
\code{\link[survival]{Surv}} object. In this case the argument \code{family} must be \code{NULL} (the
default) and the function uses Cox proportional hazard regression implemented in function
\code{\link[survival]{coxph}}

The function proceeds along the following steps:

Step 1 (bagging): \code{nBags} bootstrapped data sets are being generated based on random sampling from the
original training data set (\code{x},\code{y}). If a bag contains less than \code{minInBagObs} unique
observations or it contains all observations, it is discarded and re-sampled again.

Step 2 (random subspace): For each bag, \code{nFeaturesInBag} features are randomly selected (without
replacement) from the columns of \code{x}. Optionally, interaction terms between the selected features can
be formed (see the argument \code{maxInteractionOrder}). 

Step 3 (feature ranking): In each bag, features are ranked according to their correlation with the outcome
measure. Next the top \code{nCandidateCovariates} are being considered for forward selection in each GLM
(and in each bag).  

Step 4 (forward selection): Forward variable selection is employed to define a multivariate GLM model of the
outcome in each bag. 

Step 5 (aggregating the predictions): Prediction from each bag are aggregated. In case, of a quantitative
outcome, the predictions are simply averaged across the bags. 

Generally, \code{nCandidateCovariates}>100 is not recommended, because the forward
selection process is 
time-consuming. If arguments \code{"nBags=1, replace=FALSE, nObsInBag=nrow(x)"} are used, 
the function becomes a forward selection GLM predictor without bagging. 

Classification of multi-level categorical responses is performed indirectly by turning the single
multi-class response into a set of binary variables. The set can include two types of binary variables:
Level vs. all others (this binary variable is 1 when the original response equals the level and zero
otherwise), and level A vs. level B (this binary variable is 0 when the response equals level A, 1 when the
response equals level B, and NA otherwise). 
For example, if the input response \code{y} contains observations with values (levels) "A", "B", 
"C", the binary variables
will have names "all.vs.A" (1 means "A", 0 means all others), "all.vs.B",
"all.vs.C", and optionally also "A.vs.B" (0 means "A", 1 means "B", NA means neither "A" nor "B"), "A.vs.C",  
and "B.vs.C". 
Note that using pairwise level vs. level binary variables be
very time-consuming since the number of such binary variables grows quadratically with the number of levels
in the response. The user has the option to limit which levels of the original response will have their
"own" binary variables, by setting the minimum observations a level must have to qualify for its own binary
variable, and by explicitly enumerating levels that should not have their own binary variables. Note that
such "ignored" levels are still included on the "all" side of "level vs. all" binary variables. 

At this time the predictor does not attempt to summarize the binary variable classifications into a single
multi-level classification. 

Training this predictor on data with fewer than 8 observations is not recommended (and the
function will warn about it). Due to the bagging
step, the number of unique observations in each bag is less than the number of observations in the input
data; the low number of unique observations can (and often will) lead to an essentially perfect fit which
makes it impossible to perfrom meaningful stepwise model selection.  

Feature names: In general, the column names of input \code{x} are assumed to be the feature names. If
\code{x} has no column names (i.e., \code{colnames(x)} is \code{NULL}), stadard column names of the form
\code{"F01", "F02", ...} are used. If \code{x} has non-NULL column names, they are turned into valid and
unique names using the function \code{\link{make.names}}. If the function \code{\link{make.names}} returns
names that are not the same as the column names of \code{x}, the component \code{featureNamesChanged} will
be \code{TRUE} and the component \code{nameTranslationTable} contains the information about input and actual
used feature names. The feature names are used as predictor names in the individual models in each bag.

}

\value{
  The function returns an object of class \code{randomGLM}. For continuous prediction or two-level
classification, this is a list with the following components:
  \item{predictedOOB}{the continuous prediction (if \code{classify} is \code{FALSE}) 
or predicted classification (if \code{classify} is \code{TRUE}) of the input data based on out-of-bag
samples.
}
  \item{predictedOOB.response}{In case of a binary outcome, this is the predicted probability of each outcome
specified by \code{y} based on out-of-bag samples. In case of a continuous outcome, this is the predicted
value based on out-of-bag samples (i.e., a copy of \code{predictedOOB}).}

  \item{predictedTest.cont}{if test set is given, the predicted probability of each outcome specified by
\code{y} for test data for binary outcomes. In case of a continuous outcome, this is the test set predicted
value. }

  \item{predictedTest}{if test set is given, the predicted classification for test data. Only for binary
outcomes.} 

  \item{candidateFeatures}{candidate features in each bag. A list with one component per bag. Each component
is a matrix with \code{maxInteractionOrder} rows and \code{nCandidateCovariates} columns. 
Each column represents one
interaction obtained by multiplying the features indicated by the entries in each column (0 means no
feature, i.e. a lower order interaction). }

  \item{featuresInForwardRegression}{features selected by forward selection in each bag. A list with one
component per bag. Each component
is a matrix with \code{maxInteractionOrder} rows.
Each column represents one
interaction obtained by multiplying the features indicated by the entries in each column (0 means no 
feature, i.e. a lower order interaction). The column names contain human-readable names for the terms. }

  \item{coefOfForwardRegression}{coefficients of forward regression. A list with one
component per bag. Each component is a vector giving the coefficients of the model determined by forward
selection in the corresponding bag. The order of the coefficients is the same as the order of the terms in
the corresponding component of \code{featuresInForwardRegression}. }

 \item{interceptOfForwardRegression}{a vector with one component per bag giving the intercept of the
regression model in each bag.}


  \item{bagObsIndx}{a matrix with \code{nObsInBag} rows and \code{nBags} columns, giving the indices of
observations selected for each bag.} 

  \item{timesSelectedByForwardRegression}{a matrix of \code{maxInteractionOrder} rows and number of features
columns. Each entry gives the number of times the corresponding feature appeared in a predictor model at the
corresponding order of interactions. Interactions where a single feature enters more than once (e.g., a
quadratic interaction of the feature with itself) are counted once. 
}
  \item{models}{the regression models for each bag. Predictor features in each bag model are named using
their}

 \item{featureNamesChanged}{logical indicating whether feature names were copied verbatim from column names
of \code{x} (\code{FALSE}) or whether they had to be changed to make them valid and unique names
(\code{TRUE}).}

  \item{nameTranslationTable}{ only present if above \code{featureNamesChanged} is \code{TRUE}. 
A data frame with three columns and one row per input feature (column of input \code{x}) giving the  
feature number, original feature name, and modified feature name that is used for model fitting.}

In addition, the output value contains a copy of several input arguments. These are included to facilitate
prediction using the \code{predict} method. These returned values should be considered undocumented and may
change in the future.

In the multi-level classification classification case, the returned list (still considered a valid
\code{randomGLM} object) contains the following components:

  \item{binaryPredictors}{a list with one component per binary variable, containing the \code{randomGLM} 
predictor trained on that binary variable as the response. The list is named by the corresponding binary
variable. For example, if the input response \code{y} contains observations with values (levels) "A", "B",
"C", the binary variables (and components of this list) 
will have names "all.vs.A" (1 means "A", 0 means all others), "all.vs.B",
"all.vs.C", and optionally also "A.vs.B" (0 means "A", 1 means "B", NA means neither "A" nor "B"), "A.vs.C",
and "B.vs.C".  }

  \item{predictedOOB}{a matrix in which columns correspond to the binary variables and rows to samples,
containing the predicted binary classification for each binary variable. Columns names and meaning of 0 and
1 are described above. }

  \item{predictedOOB.response}{a matrix with two columns per binary variable, giving the class probabilities
for each of the two classes in each binary variables. Column names contain the variable and class names.}

  \item{levelMatrix}{a character matrix with two rows and one column per binary variable, giving the level
 corresponding to value 0 (row 1) and level corresponding to value 1 (row 2). This encodes the same
information as the names of the \code{binaryPredictors} list but in a more programmer-friendly way.}

  If input \code{xTest} is non-NULL, the components \code{predictedTest} and \code{predictedTest.response}
contain test set predictions analogous to \code{predictedOOB} and \code{predictedOOB.response}.
}
  
\references{
Lin Song, Peter Langfelder, Steve Horvath: Random generalized linear model: a highly accurate and interpretable ensemble predictor. BMC Bioinformatics 2013 Jan 16;14:5. doi: 10.1186/1471-2105-14-5. 
}

\author{Lin Song, Steve Horvath, Peter Langfelder. 
The function makes use of the \code{glm} function and other standard R functions. }

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
   xTest, 
   nCandidateCovariates=ncol(xTrain), 
   nBags=30, nThreads = 1)

yPredicted = RGLM$predictedTest
table(yPredicted, yTest)


## continuous outcome prediction

x=matrix(rnorm(100*20),100,20)
y=rnorm(100)

xTrain = x[1:50,]
yTrain = y[1:50]
xTest = x[51:100,]
yTest = y[51:100]

RGLM = randomGLM(
  xTrain, yTrain, 
  xTest, classify=FALSE, 
  nCandidateCovariates=ncol(xTrain), 
  nBags=10, 
  keepModels = TRUE, nThreads = 1)
}
\keyword{misc}

