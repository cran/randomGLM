# Peter Langfelder's additions and changes

# 1.06:

#  . Add the capability to work with factors. There are a few issues with this:

#    1. Sampling: all levels of a factor should be present in every bag. If they are not, it will lead to
#    problems when the factor has more than 2 levels: n-1 levels are present in the bag so the factor
#    variable is retained for the models, and the predict() function will throw an error.

#    2. Generating interactions: I should be able to use the function model.matrix to generate the
#    interactions. It works automatically with factors, uses formula as input, and basically generates the
#    same thing that my .generateInteractions will do with numeric variables. The workflow for generating
#    interactions should be this: prepare the formula as a character string and call terms() on it. Part of
#    the terms() value looks a lot like my interaction matrix. Save the result as part of the output for
#    each bag; call model.matrix on it to get the actual interaction variables. These can be used to do
#    variable selection and final model fitting.
#
#    Alternatively, can do this on my own. First, binarize all categorical variables into values 0, 1. Then
#    generate interactions as before except binary variables can only enter with powers 0 or 1. This will
#    entail keeping track of which variables are binary, and changing the interaction matrix generator
#    appropriately. Further, for prediction on test data, need to keep information on the original levels
#    of the categorical predictors so the test set predictors can be appropriately binarized.

#    Note that there's no need to generate interactions of binary variables in a special way; a simple
#    product of binary indicators is enough. This is obvious for 2 binary variables a, b: for example, 
#    a & !b equals a-a*b. By induction this is also true for any number of binary variables: a combination
#    of n variables can be expressed as a combination of 2 variables: the combination of n-1 variables and
#    the last variable, so if I can express combinations of n-1 variables, I can also express combinations
#    of n variables. 


#  . To accomodate categorical variables:

#    Will keep track of which variables were categorical, what their levels were, and which
#    prediction variable they correspond to.

#    Put the binarized categorical variables into a separate xCat variable for the bagged prediction? They
#    will have to be handled with more care for bagging and during prediction anyway. 

#  . convert code to x being a data frame. This is necessary to correctly handle possible factors in the
#    input x and xTest.

#  . Should also add an option to in-/exclude self-interactions since such terms are non-standard in normal
#    statistical model fitting

# Main change: unifying all levels of interactions in a single function. 

# Introducing the concept of an interaction matrix that specifies how the input
# features get multiplied into the predictors used in the underlying models. 
# The matrix has maxInteractionOrder rows and as many columns as needed. Each
# column corresponds to one term in for the (g)lm models. The entries in a
# column give indices of the input features that are multiplied in the term; 0
# index means "none". Thus, if one has 3 input features x1, x2, x3, an interaction
# matrix with up to 3-way interactions may look like (first line is column name)

# c1 c2 c3 c4 c5 c6 c7
#  1  2  3  1  1  1  2 
#  0  0  0  1  2  3  2
#  0  0  0  0  2  0  0

# This means that term 1 (column 1) is x1, term 2 is x2, term 3 is x3, term 4 is
# x1*x1, term5 is x1*x2*x2, term 6 in x1*x3, term 7 is x2*x2. No particular order
# is assumed in any of the rows or columns. The function that generates this
# matrix, named .interactionMatrix, 
# also assigns unique names to each column which is helpful when identifying
# individual terms and extracting their content.
#
# The advantage of this representation is that it makes it very convenient to
# actually generate the interaction terms from given numeric data. The generation
# of the interactions is implemented in function .generateInteractions.
#
# Furthermore, given an interaction matrix, it is easy to count the occurences of
# each input feature. This is implemented in function .countsInInteractionMatrix
# which counts the number of times each input feature appears in each level of
# interactions in a given interaction matrix.


# . Removed scaling of predictors

## change to 3.0
# 1. function name change to randomGLM
# 2. parameter change: corFncForCandidateCovariates, corOptionsForCandidateCovariates
# 3. change all "gene" to "feature"
# 4. add interceptOfForwardRegression

# For now, include required packages

#library(gtools)
#library(MASS)

#=====================================================
#
# Helper functions
#
#=====================================================

.spaste = function(...)
{
  paste(..., sep = "");
}

# prettier print than the standard "print" function.

.cat.nl = function(...)
{
  cat(.spaste(..., "\n"));
}

.prependZeros = function(x, len = max(nchar(x)))
{
  lengths = nchar(x);
  if (len < max(lengths)) stop("Some entries of 'x' are too long.");
  out = as.character(x);
  n = length(x);
  for (i in 1:n) if (lengths[i] < len)
    out[i] = .spaste( paste(rep("0", len-lengths[i]), collapse = ""),
                     x[i]);

  out;
}

.is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  { abs(x - round(x)) < tol }


  
#=====================================================
#
# Generator of interactions
#
#=====================================================

# My own function for generating combinations since the one in gtools fails with large n's.
# Combinations with repeats

.combinationsWithReplacement = function(n, order)
{
  if (order==1) return(matrix(c(1:n), 1, n));

  # If order is not 1, calculate result using recursion
  nOut = choose(n+order-1, order);
  out = matrix(0, order, nOut);
  sub = .combinationsWithReplacement(n, order-1)
  index = 1;
  for (i in 1:n)
  {
    n1 = ncol(sub);
    out[, index:(index + n1 - 1)] = rbind( rep(i, ncol(sub)), sub);
    index = index + n1;
    sub = sub[, colSums(sub==i) == 0, drop = FALSE];
  }
  out;
}

# Assume that order is at most n. If not, order will be adjusted
.combinationsWithoutReplacement = function(n, order)
{
  if (order > n) order = n;
  if (order==1) return(matrix(c(1:n), 1, n));

  # If order is not 1, calculate result using recursion
  nOut = choose(n, order);
  out = matrix(0, order, nOut);
  sub = .combinationsWithoutReplacement(n-1, order-1) + 1  
  index = 1;
  for (i in 1:(n-1))
  {
    n1 = ncol(sub);
    if (n1 > 0)
    {
      out[, index:(index + n1 - 1)] = rbind( rep(i, ncol(sub)), sub);
      index = index + n1;
      sub = sub[, colSums(sub==i+1) == 0, drop = FALSE];
    }
  }
  out;
}

.combinations = function(nWith, nWithout, order)
{
  nAll = nWith + nWithout;
  if (order==1) return(matrix(c(1:nAll), 1, nAll));
  if (nWith==0) return (.combinationsWithoutReplacement(nWithout, order));

  nOut = 0;
  for (o in 0:max(order, nWithout))
    nOut = nOut + choose(nWith+o-1, o) * choose(nWithout, order-o)
  out = matrix(0, order, nOut);
  sub = .combinations(nWith, nWithout, order-1);
  index = 1;
  for (i in 1:(nAll- (nWithout>0)))
  {
    wr = i<=nWith;
    n1 = ncol(sub);
    if (n1 > 0)
    {
      out[, index:(index + n1 - 1)] = rbind( rep(i, ncol(sub)), sub);
      index = index + n1;
      if (wr) {
        sub = sub[, colSums(sub==i) == 0, drop = FALSE];
        if (i==nWith & nWithout > 0)
           sub = sub[, colSums(sub==i+1) == 0, drop = FALSE];
      } else
        sub = sub[, colSums(sub==i+1) == 0, drop = FALSE];
    }
  }
  out;
}


 
# Generates an interaction matrix of all interactions up to specified maxOrder of
# variables indexed 1:nWith (with replacement) and (nWith+1):(nWithout) (without replacement). 
# Optionally also sets column names with some flexibility.
# Empty slots in the matrix are marked with a zero.

.interactionMatrix = function(nWith, nWithout, maxOrder, 
                 index.without,
                 setColNames = TRUE,
                 originalNames = c(1:(nWith+nWithout)),
                 featureSeparator = ".")
{
  # if (maxOrder==1) return(matrix(c(1:(nWith+nWithout)), 1, nWith + nWithout));

  combs = .combinations(nWith + 1, nWithout, maxOrder)[, -1, drop = FALSE]; 
            # drop the 1 1 1 1... column (intercept)
  if (length(index.without) > 0)
  {
    generatedOrder = c(1:(nWith + nWithout +1));
    targetOrder = c(generatedOrder[-(index.without+1)], generatedOrder[(index.without+1)]);
    out = targetOrder[combs];
    dim(out) = dim(combs);  
  }  else
    out = combs;

  # Perform a few tricks to get the order of the variables nice.

  # Move the number 1 variable to last

  last = nWith + nWithout + 2;
  out[combs==1] = last;
  out = out-1;

  # Flip the row order to have the empty slots on the bottom and sort each column of the matrix.
  out = -t(matsort(-t(out [ c(maxOrder:1), , drop = FALSE])));

  nameMatrix = array(c(originalNames, "")[out], dim = dim(out));
  if (setColNames) 
  {
    colnames = apply(nameMatrix, 2, function(x) {paste(x[x!=""], collapse = featureSeparator)});
    colnames(out) = colnames;
  }
  rownames(out) = .spaste("Feature.", c(1:maxOrder));

  # Get the order of columns
  orderNames = 1:(nWith + nWithout);
  orderNameMatrix = array(c(orderNames, 0)[out], dim = dim(out))
  interactionOrder = colSums(out!=last-1);
  orderDF = as.data.frame(cbind(interactionOrder, t(orderNameMatrix)));
  order = do.call("order", orderDF);

  out2 = out[, order, drop = FALSE];
  out2[out2==last-1] = 0;
  out2;
}


# utility function used by .generateInteractions below.
.generateInteractions.1row = function(x, interactionMatrix)
{
  mat = x[interactionMatrix];
  dim(mat) = dim(interactionMatrix);
  apply(mat, 2, prod);
}

# Generates interactions from a given matrix of numeric data for "original" or
# "input" features and an interaction matrix. 
# the trick to generate everything in one step is to add a dummy "feature" equal to
# 1 in all rows to the actual features and use it as "feature with index 0". 
# The function optionally transfers the column names from the interaction matrix to
# the result, which makes the columns names of the result unique and somewhat
# descriptive.

# Since the function is always called with a valid interaction matrix, remove the code that generates it if
# it's empty.

.generateInteractions = function(x, #maxOrder, binaryIndex, includeSelfInteractions,
                                 interactionMatrix, x1 = NULL,
                                 setColNames = FALSE)
                                 #originalNames= c(1:ncol(x)))
{
  nAll = ncol(x);
  #if (is.null(interactionMatrix)) 
  #{
  #  if (includeSelfInteractions)
  #  {
  #    nWithout = length(binaryIndex)
  #  } else {
  #    nWithout = nAll;
  #    binaryIndex = 1:nAll;
  #  }
  #  nWith = nAll - nWithout;
  #  interactionMatrix = .interactionMatrix(nWith, nWithout, maxOrder, index.without = binaryIndex, 
  #                                         setColNames = setColNames, 
  #                                         originalNames = originalNames);
  #  
  #  
  #}
  maxOrder= nrow(interactionMatrix);
  if (maxOrder==1)
  { 
    # If maxOrder is 1, the result is simply x but the order may be changed.
    x.int = x[, as.vector(interactionMatrix), drop = FALSE];
  } else {
    # The non-trivial case
    if (is.null(x1)) x1 = cbind(x, rep(1, nrow(x)));
    interactionMatrix[ interactionMatrix==0 ] = nAll+1;

    # If the number of variables (columns in interactionMatrix) is 1, the result needs no transposing;
    # otherwise it needs to be transposed.
    if (ncol(interactionMatrix)==1)
    {
      x.int = as.matrix(apply(x1, 1, .generateInteractions.1row, interactionMatrix))
    } else 
      x.int = t(apply(x1, 1, .generateInteractions.1row, interactionMatrix));
  }
  if (setColNames) colnames(x.int) = colnames(interactionMatrix);
  x.int
}

# This function counts the number of times each input feature appears in each level
# of interactions in a given interaction matrix.

.countsInInteractionMatrix = function(im, nFeatures)
{
  maxLevel = nrow(im);
  counts = matrix(0, maxLevel, nFeatures);
  level = maxLevel - colSums(im==0);
  for (l in 1:maxLevel)
  {
    mat1 = im[ 1:l, level==l, drop = FALSE];
    mat1.unique = sapply(as.data.frame(mat1), unique);
    counts1 = table(unlist(mat1.unique));
    where = as.numeric(names(counts1));
    counts[ l, where] = as.numeric(counts1);
  }
  rownames(counts) = .spaste("Level.", c(1:maxLevel));
  counts;
}


# Translate ordinal data using a dictionary

.translate = function(data, dictionary)
{
  translated = dictionary[ match(data, dictionary[, 1]), 2];
  attributes(translated) = attributes(data);
  translated;
}

#=======================================================================================================
#
# binarizeCategoricalVar
#
#=======================================================================================================

# Assumes x is a vector but can easily be modified to also work with matrices.

.binarizeCategoricalVar = function(x, minCount = 3, val1 = 0, val2 = 1, nameSep = ".vs.", namePrefix = "",
                                  nameForAll = "all",
                                  ignore = NULL, includePairwise = TRUE,
                                  includeLevelVsAll = FALSE, levels = NULL, levelOrder = NULL,
                                  dropFirstLevelVsAll = FALSE,
                                  dropUninformative = TRUE)

{
  if (is.null(levels))
  {
    tab = table(x);
    levels0 = names(tab);
    tab = tab[ tab >= minCount & !(levels0 %in% ignore) ];
    levels = names(tab);
  }
  if (!is.null(levelOrder))
  {
    order = match(levelOrder, levels);
    order = order[is.finite(order)];
    levels0 = levels[order];
    levels1 = levels[ !levels %in% levels0];
    levels = c(levels0, levels1);
  }
  nSamples = length(x);
  nLevels = length(levels)
  nBinaryVars = includePairwise * nLevels * (nLevels - 1)/2 +
                      includeLevelVsAll * (nLevels - dropFirstLevelVsAll)
  if (nBinaryVars==0)
  {
    if (dropUninformative)
    {
      return(NULL)
    } else {
      out = as.matrix(rep(val2, nSamples));
      colnames(out) = levels[1];
      return(out);
    }
  }
  out = matrix(NA, nSamples, nBinaryVars)
  levelTable = matrix("", 2, nBinaryVars);
  ind = 1;
  names = rep("", nBinaryVars);
  if (includePairwise)
  {
    for (v1 in 1:(nLevels-1)) for (v2 in (v1+1):nLevels)
    {
       out[ x==levels[v1], ind] = val1;
       out[ x==levels[v2], ind] = val2;
       names[ind] = .spaste(namePrefix, levels[v1], nameSep, levels[v2]);
       levelTable[, ind] = levels[ c(v1, v2)];
       ind = ind + 1;
    }
  }
  if (includeLevelVsAll)
    for (v1 in (1 + as.numeric(dropFirstLevelVsAll)):nLevels)
    {
      out[, ind] = c(val1, val2) [ as.numeric(x==levels[v1])+1 ];
      names[ind] = .spaste(namePrefix, nameForAll, nameSep, levels[v1]);
      levelTable[, ind] = c(nameForAll, levels[v1]);
      ind = ind+1;
    }
  colnames(out) = names;
  colnames(levelTable) = names;
  rownames(levelTable) = .spaste("Value.", c(val1, val2));
  attr(out, "includedLevels") = levelTable;
  out;
}

# Split factors into independent binary variables. Leave out the last level.
# It is important to always leave out the same level, i.e. not make the left-out level dynamic.

.binarizeCategoricalColumns = function(df, categoricalIndicator, levels)
{
  if (sum(categoricalIndicator)==0) 
  {
    attr(df, "originalNames") = names(df);
    attr(df, "indexOfOrigin") = c(1:ncol(df));
    attr(df, "isCategorical") = rep(FALSE, ncol(df));
    return(df);
  }

  if (length(levels)>0 && (length(levels)!=sum(categoricalIndicator)))
    stop(".binarizeCategoricalColumns: Internal error: length of 'levels' is not the same as",
         " length of 'categoricalFactors'");

  if (is.null(levels))
  {
    binFactors.list = mapply(.binarizeCategoricalVar, x = df[categoricalIndicator], 
                      nameForAll = names(df)[categoricalIndicator], 
                      MoreArgs = list(minCount = 0, includePairwise = FALSE, includeLevelVsAll = TRUE, 
                      nameSep = ".", dropFirstLevelVsAll = TRUE,
                      dropUninformative = FALSE), SIMPLIFY = FALSE);
  } else
    binFactors.list = mapply(.binarizeCategoricalVar, x = df[categoricalIndicator], 
                      nameForAll = names(df)[categoricalIndicator], levels = levels,
                      MoreArgs = list(minCount = 0, includePairwise = FALSE, includeLevelVsAll = TRUE,
                      nameSep = ".", dropFirstLevelVsAll = TRUE, dropUninformative = FALSE), SIMPLIFY = FALSE);
  
  factorNames = mapply(rep, names(df)[categoricalIndicator], sapply(binFactors.list, ncol), SIMPLIFY = FALSE)

  if (any(!categoricalIndicator)) 
  {
    out.0 = cbind(df[ , !categoricalIndicator], do.call(cbind, binFactors.list));
    names.0 = c(names(df)[!categoricalIndicator], unlist(factorNames));
    order = unlist(lapply(names(df), function(x, l) {which(x==l)}, names.0));
    df.out = as.data.frame(out.0[, order]);
  } else {
    df.out = as.data.frame(do.call(cbind, binFactors.list));
    names.0 = unlist(factorNames);
    order = unlist(lapply(names(df), function(x, l) {which(x==l)}, names.0));
  }

  names(df.out) = make.names(names(df.out));

  names.0.ord =  names.0[order];
  attr(df.out, "originalNames") = names.0.ord;
  attr(df.out, "indexOfOrigin") = match(names.0.ord, names(df));
  attr(df.out, "isCategorical") = names(df.out)!=names.0.ord

  df.out;

}

.nonMissingLevels = function(x)
{
  u = unique(x);
  sort(u[!is.na(u)]);
}

.nNonMissingLevels = function(x)
{
  u = unique(x);
  length(u) - sum(is.na(u));
}
.categoricalColumns = function(x, maxLevels)
{
  nLevels = sapply(x, .nNonMissingLevels);
  fac = sapply(x, is.factor);
  
  nLevels <= maxLevels | fac;
}


#==========================================================================================
#
# .createModel and .predictFromModel
#
#==========================================================================================
# Depending on type and family, get the appropriate model for the formula and data.

.createModel = function(formula, modelData, weights, type, family)
{
    switch(type, linear = lm(as.formula(formula), data = modelData, weights = weights, model = FALSE),
                 count = ,
                 binary = glm(as.formula(formula), data = modelData, weights = weights, 
                              family = family, model = FALSE),
                 general = if (substring(family$family, 1, nchar("Negative Binomial"))=="Negative Binomial")
                               glm.nb(as.formula(formula), data = modelData, weights = weights,
                                      link = match.fun(family$link), model = FALSE) else
                               glm(as.formula(formula), data = modelData, weights = weights, 
                                   family = family, model = FALSE),
                 survival = coxph(as.formula(formula), data = modelData, weights, model = FALSE))
}

.predictFromModel = function(model, newdata, modelType)
{
  switch(modelType, linear = predict(model, newdata = newdata),
                    count = ,
                    binary = ,
                    general = predict(model, newdata = newdata, type = "response"),
                    survival = predict(model, newdata, type = "lp") )
}

.keepCoeffs = function(coeffs, modelType)
{
  n = nrow(coeffs);
  drop1 = rownames(coeffs)[1] == "(Intercept)"
  if (drop1) return(2:n)
  c(1:n);
}

.haveIntercept = function(coeffs, modelType)
{
  return(rownames(coeffs)[1] == "(Intercept)")
}

#===========================================================================================
#
# .significance: calculate significance of features for outcome.
#
#===========================================================================================

.significance = function(x, y, 
                         type, family,
                         corFncForCandidateCovariates, corOptionsForCandidateCovariates)
{
  if (type=="survival")
  {
    y.cor = residuals(coxph(y~ 1, model = TRUE), type = "deviance");
  } else 
    y.cor = y;

  corOptionsForCandidateCovariates$x = x;
  corOptionsForCandidateCovariates$y = y.cor;

  as.vector(do.call(corFncForCandidateCovariates, corOptionsForCandidateCovariates));
}

#===========================================================================================
#
# Identify bad features: features with zero variance or with missing data.
# Requires the function colSds from matrixStats.
#
#===========================================================================================

.goodFeatures = function(x, returnIndex = FALSE)
{
  indxMiss = colSums(is.na(x))>0
  indxVar = colSds(as.matrix(x), na.rm = TRUE) ==0
  out = !(indxMiss | indxVar);
  if (returnIndex) which(out) else out;
}


#=====================================================
#
# forwardSelection
#
#=====================================================

# Note: yBag is assumed to be a matrix.

.forwardSelection = function(xBag, yBag, xTestBag, 
                             weights.Bag,
                             classify, 
                             binaryIndicator,  # Need indicator instead of an index 
                             maxInteractionOrder,
                             includeSelfinteractions,
                             nCandidateCovariates,  
                             corFncForCandidateCovariates, corOptionsForCandidateCovariates, 
                             NmandatoryCovariates,
                             interactionsMandatory,
                             keepModel,
                             interactionSeparatorForCoefNames,
                             type,
                             family,
                             responseName)
{  
  # remove features with missing values or variance=0
  keepFeatures = .goodFeatures(xBag, returnIndex = FALSE);
  if (!any(keepFeatures)) 
    stop("All predictor features have been removed from a bag due to the presence of missing data or",
         "\ndue to the feature being constant. Please remove constant features and ",
         "\neither impute missing data or remove features and/or samples with missing data.");

  xBag = xBag[, keepFeatures, drop=FALSE]
  xTestBag = xTestBag[, keepFeatures, drop=FALSE]
  binaryIndicator = binaryIndicator[keepFeatures];

  nFeatures = ncol(xBag);

  if (NmandatoryCovariates>0)
  {
    mandatCovars = c(1:NmandatoryCovariates)
    # if removed features include mandatory cov, then NmandatoryCovariates should be decreased
    mandatCovars = mandatCovars[ keepFeatures[mandatCovars]];
    NmandatoryCovariates = length(mandatCovars);
  } else
    mandatCovars = numeric(0);

  nonMandatCovars = setdiff( c(1:nFeatures), mandatCovars);

  # Add interaction terms:
  # generate the interaction matrix

  if (includeSelfinteractions)
  { 
    nWithout = sum(binaryIndicator);
  } else {
    nWithout = nFeatures;
    binaryIndicator = rep(TRUE, nFeatures);
  }

  nWith = nFeatures - nWithout;

  interactionMatrix = .interactionMatrix(nWith, nWithout, maxOrder = maxInteractionOrder, 
                                    index.without = which(binaryIndicator), 
                                    originalNames = colnames(xBag),
                                    setColNames = TRUE, featureSeparator = interactionSeparatorForCoefNames);
  # Identify mandatory interactions. If interactions of mandatory covariates are also mandatory, add all
  # interactions where at least one term is mandatory (which will be many, so must be used with
  # caution)
  if (interactionsMandatory)
  {
    #mandatoryInteractions = apply( interactionMatrix, 2, function(x) { any(x %in% mandatCovars) } );
    inMandat = interactionMatrix %in% mandatCovars;
    dim(inMandat) = dim(interactionMatrix);
    mandatoryInteractions = colSums(inMandat) > 0;
  } else 
    mandatoryInteractions = mandatCovars;

  nMandatoryInteractions = length(mandatoryInteractions);
  if (nMandatoryInteractions > nCandidateCovariates)
     stop("Number of mandatory interactions is larger than number of candidate covariates.");

  nInteractions = ncol(interactionMatrix);
  nonMandatInteractions = setdiff( c(1:nInteractions), mandatoryInteractions);

  x.int = .generateInteractions(xBag, interactionMatrix = interactionMatrix, setColNames = TRUE);
  xTest.int = .generateInteractions(xTestBag, interactionMatrix = interactionMatrix, setColNames = TRUE);


  # calculate feature significance  
  # corOptionsForCandidateCovariates$x = x.int
  # corOptionsForCandidateCovariates$y = yBag
  absGS = abs(.significance(x.int, yBag, type = type, family = family,
                            corFncForCandidateCovariates = corFncForCandidateCovariates,
                            corOptionsForCandidateCovariates = corOptionsForCandidateCovariates));
                            
  ## nCandidateCovariates could be smaller than indicated due to missing data in x.
  nCandidateCovariates = min(nCandidateCovariates, ncol(x.int))

  ## get indices of candidate cov
  rank = rank(-absGS[nonMandatInteractions], ties.method="f")
  indx = c(mandatoryInteractions, nonMandatInteractions[rank<=(nCandidateCovariates-nMandatoryInteractions)]);

  x.int = x.int[, indx, drop=FALSE]
  xTest.int = xTest.int[, indx, drop=FALSE]
  absGS = absGS[indx]

  # output candidate cov
  candidateFeatures = interactionMatrix[, indx, drop = FALSE];

  # index of most significant feature, used in initial model.
  featureMax = which.max(absGS)

  ## define initial model and full model for binary and continuous outcome. Mandatory covariates must show
  # up in final model, so put them in initial model.

  # the data frame modelData will contain the predictors and the response. Note that the apparently unusual
  # way of how data.frame handles arguments that are itself matrices comes in handy here.
  # In other words, even if yBag is a matrix, it is represented by a single "column" in modelData.
  rownames(x.int) = make.names(rownames(x.int), unique = TRUE);
  modelData = data.frame(x.int, response = yBag);
  names(modelData)[ ncol(modelData) ] = responseName;
  predictorNames = colnames(x.int);

  # Need to remove large variables from the current environment because a copy of this environment is kept
  # in the formula and models below. 
  rm(xBag, xTestBag, x.int, corOptionsForCandidateCovariates, interactionMatrix);

  initialFormula = paste(responseName, "~", 
               paste(predictorNames[ if (nMandatoryInteractions>0) mandatoryInteractions else featureMax ], 
                     collapse = " + "));
  initialModel = .createModel(initialFormula, modelData, weights = weights.Bag, type = type, family = family);

  upperFormula = paste(responseName, "~ .")
  upperModel = .createModel(upperFormula, modelData, weights = weights.Bag, type = type, family = family);

  ## forward model selection
  model = stepAIC(initialModel,
		  scope = list(upper = upperModel), 
		  direction="forward", 
		  trace=FALSE)

  ## output selected feature (by their names) and their coefficients, is there any smarter way to fish out
  # which features are selected into model? 

  coeffs = summary(model)$coefficients;
  keepCoeffs = .keepCoeffs(coeffs, type);
  haveIntercept = .haveIntercept(coeffs, type)
  selected = rownames(coeffs)[keepCoeffs]
  featuresInForwardRegression = candidateFeatures[, match(selected, colnames(candidateFeatures)), 
                                                    drop = FALSE];

  # output candidate covariates, taking into account the fact that some features may have been dropped.
  dict1 = cbind(0:nFeatures, c(0, which(keepFeatures)))
  candidateFeatures = .translate(candidateFeatures, dict1);
  featuresInForwardRegression = .translate(featuresInForwardRegression, dict1);

  coefOfForwardRegression = coeffs[keepCoeffs,1]
  interceptOfForwardRegression = if (haveIntercept) coeffs[1,1] else NA; 

  ## outHat is piHat for binary outcome and yHat for quantitative outcome
  outHat = .predictFromModel(model, newdata = as.data.frame(xTest.int), modelType = type)

  out = list(predicted = outHat, 
             retainedFeatures = which(keepFeatures),
             candidateFeatures = candidateFeatures, 
             featuresInForwardRegression = featuresInForwardRegression, 
             coefOfForwardRegression = coefOfForwardRegression, 
             interceptOfForwardRegression = interceptOfForwardRegression,
             model = if (keepModel) model else NULL )

  # The environment of this function is kept in the model returned by stepAIC. Thus, delete everything but
  # the output value so the environment doesn't take up too much memory.
  varList = ls(all.names = TRUE);
  rm(list = setdiff( varList, c("out")));
  out;
}


#=====================================================
#
# Handling of multi-threaded calculations
#
#=====================================================

.disableThreads = function(clusterInfo = list())
{
  if (!is.null(clusterInfo$cluster))
      try( stopCluster(cluster), silent = TRUE);
}

.enableThreads = function(nThreads, verbose)
{
  if (is.null(nThreads)) nThreads = max(ceiling(detectCores()*3/4), detectCores()-1);
  if (is.na(nThreads)) nThreads = 1;

  if (nThreads < 1) 
  {
    warning("In function randomGLM: 'nThreads' is below 1. Will use serial execution.");
    nThreads = 1;
  }

  if (nThreads > 1)
  {
    if (verbose > 1) .cat.nl("Will use parallel calculation with ", nThreads, " workers.");
    if (.Platform$OS.type=="windows")
    {
      # On Windows: create a cluster manually
      # and  export the parent evinronment of this function for randomGLM to work as well, plus
      # the environment of packages gtools and MASS that are needed. 
      cluster = makePSOCKcluster(nThreads, outfile = "");
      #assign(".randomGLMparallelCluster", cluster, pos = ".GlobalEnv");
      clusterExport(cluster, varlist = ls(envir = parent.env(environment()), all.names = TRUE), 
                             envir = parent.env(environment()))
      clusterCall(cluster, library, package = "MASS", character.only = TRUE);
      clusterCall(cluster, library, package = "gtools", character.only = TRUE);
      registerDoParallel(cluster);
      out = list(nThreads = nThreads, cluster = cluster);
    } else {
      # On linux, simply register a parallel backend with nThreads workers
      registerDoParallel(nThreads);
      out = list(nThreads = nThreads, cluster = NULL);
    }
  } else {
    if (verbose > 1) .cat.nl("Will use serial calculation with a single worker process.");
    registerDoSEQ();
    out = list(nThreads = nThreads, cluster = NULL);
  }
  out; 
}

#=================================================================================================
#
# main user-level function randomGLM
#
#=================================================================================================

randomGLM = function(
  # Input data
  x, y, xtest = NULL, 

  weights = NULL,

  # Which columns in x are categorical?
  categoricalColumns = NULL,
  maxCategoricalLevels = 2,

  # Include interactions?
  maxInteractionOrder = 1,
  includeSelfinteractions = TRUE,

  # Prediction type: type can be used to set the prediction type in a simplified way...
  type = c("auto", "linear", "binary", "count", "general", "survival"),

  # classify is retained mostly for backwards compatibility
  classify = switch(type, auto = !is.Surv(y) & (is.factor(y) | length(unique(y)) < 4),
                          linear = FALSE,
                          binary = TRUE ,
                          count = FALSE,
                          general = FALSE,
                          survival = FALSE),

  # family can be used to fine-tune the underlying regression model
  family = switch(type, auto = NULL,
                        linear = gaussian(link="identity"),
                        binary = binomial(link=logit),
                        count = poisson(link = "log"),
                        general = NULL,
                        survival = NULL),

  # Multi-level classification options - only apply to classification with multi-level response
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
{

  # save original data 
  ySaved = y;
  xSaved = x;

  type = match.arg(type);
  
  # if y is binary, extract y levels
  if (classify)
  {
    if (! type %in% c("auto", "binary"))
      stop("Inconsistent 'type' and 'classify': 'classify' can be TRUE only\n", 
           "   when 'type' is \"auto\" or \"binary\"");

    originalYLevels = sort(unique(y));

    # If y has more than 2 levels, do classification on binarized variables. 
    if (length(originalYLevels)>2) 
    {
      if (is.na(multiClass.minObs)) multiClass.minObs = 0;
      if (multiClass.minObs < 1) 
      {
         .cat.nl("Warning: invalid input of 'multiClass.nimObs' changed to 1.");
         multiClass.minObs = 1;
      }
      if (length(originalYLevels) > length(y)/multiClass.minObs | length(originalYLevels) == length(y))
      {
        stop("The response 'y' has too many levels for classification.\n", 
             "   Perhaps you should set 'classify = FALSE' or an appropriate 'type'?");
      } else {
        .cat.nl("randomGLM: transforming multi-level response to a series of binary variables.");
      }
   
      yBin = .binarizeCategoricalVar(as.character(y),
                                  minCount = multiClass.minObs, 
                                  val1 = 0, val2 = 1, nameSep = ".vs.", namePrefix = "",
                                  ignore = multiClass.ignoreLevels, 
                                  includePairwise = multiClass.pairwise,
                                  includeLevelVsAll = multiClass.global, 
                                  levelOrder = NULL);
      nY = ncol(yBin);
      yBinNames = colnames(yBin);
      yBinLevels = attr(yBin, "includedLevels");
      
      # Apply randomGLM recursively to each column of yBin.

      out = list(binaryPredictors = list());
      for (iy in 1:nY)
      {
        if (verbose > 0)
          .cat.nl("..Working on binary variable ", yBinNames[iy], " (", iy, " of ", nY, ")");
        out$binaryPredictors[[iy]] = randomGLM(x = x, y = yBin[, iy], xtest = xtest,
                              weights = weights,
                              categoricalColumns = categoricalColumns,
                              maxCategoricalLevels = maxCategoricalLevels,
                              maxInteractionOrder = maxInteractionOrder,
                              includeSelfinteractions = includeSelfinteractions,
                              type = type,
                              classify = classify,
                              family = family,
                              nBags = nBags,
                              replace = replace,
                              sampleBaggingWeights = sampleBaggingWeights,
                              nObsInBag = nObsInBag,
                              nFeaturesInBag = nFeaturesInBag,
                              minInBagObs = minInBagObs,
                              maxBagAttempts = maxBagAttempts,
                              replaceBadBagFeatures = replaceBadBagFeatures,
                              nCandidateCovariates = nCandidateCovariates,
                              corFncForCandidateCovariates = corFncForCandidateCovariates,
                              corOptionsForCandidateCovariates = corOptionsForCandidateCovariates,
                              mandatoryCovariates = mandatoryCovariates,
                              interactionsMandatory = interactionsMandatory,
                              keepModels = keepModels,
                              thresholdClassProb = thresholdClassProb,
                              interactionSeparatorForCoefNames = interactionSeparatorForCoefNames,
                              randomSeed = randomSeed,
                              nThreads = nThreads,
                              verbose = verbose - 1);

      }

      names(out$binaryPredictors) = yBinNames;

      out$predictedOOB = as.matrix(sapply(out$binaryPredictors, getElement, "predictedOOB"));
      colnames(out$predictedOOB) = .spaste("PredictionFor.", yBinNames);

      
      out$predictedOOB.response = do.call(cbind, 
                         lapply(out$binaryPredictors, getElement, "predictedOOB.response"));
      responseNames = .spaste(rep(yBinNames, rep(2, nY)), ".ProbabilityOfClass.", as.vector(yBinLevels));
      colnames(out$predictedOOB.response) = responseNames;

      rownames(out$predictedOOB.response) = rownames(out$predictedOOB) = rownames(x);
      if (!is.null(xtest))
      {
        out$predictedTest = sapply(out$binaryPredictors, getElement, "predictedTest");
        colnames(out$predictedTest) = yBinNames;

        out$predictedTest.response = do.call(cbind, 
                         lapply(out$binaryPredictors, getElement, "predictedTest.response"));
        colnames(out$predictedOOB.response) = responseNames;
        rownames(out$predictedTest) = rownames(out$predictedTest.response) = rownames(xtest);
      }

      out$levelMatrix = yBinLevels;
      out$thresholdClassProb = thresholdClassProb;

      class(out) = c("randomGLM", class(out));
     
      return(out);
    }

    y = as.numeric(as.factor(y))-1;
    # The next 3 lines are not needed since the results are always c(0, 1), 0, and 1.
    numYLevels = sort(unique(y));
    minY = min(y, na.rm = TRUE);
    maxY = max(y, na.rm = TRUE);
  } 

  # If type is "auto", guess the actual type. The value "auto" is invalid in internal code.
  if (type=="auto")
  {
    if (is.Surv(y)) 
    {
       type = "survival"
    } else if (classify) {
       type = "binary"
    } else if (all(.is.wholenumber(y))) {
       type = "count"
    } else type = "linear"
    # The default family for type="auto" is NULL; get the appropriate family here now.
    if (is.null(family)) 
       family = switch(type, linear = gaussian(link="identity"),
                             binary = binomial(link=logit),
                             count = poisson(link = "log"),
                             survival = NULL);

    if (verbose > 2) .cat.nl("Determined type ", type);
  } else {
    if (classify && type!="binary")
      stop("Inconsistent input: when 'classify' is TRUE, 'type' must be 'binary' or 'auto'.");
  }

  if (type!="survival" && is.null(family))
    stop("'family' could not be set automatically (most likely because the input 'type' is 'general').\n",
         "  Please set an appropriate 'family' or use a 'type' that automatically selects 
            an appropriate 'family'.");

  # Check that the response y is numeric, unless type is 'survival'

  responseTypeOK = switch(type, linear = is.numeric(y),
                          binary = is.numeric(y),
                          count = is.numeric(y) && all(.is.wholenumber(y)),
                          general = is.numeric(y),
                          survival = is.Surv(y));

  if (!responseTypeOK)
    stop(.spaste("The response 'y' is of incorrect type. For 'type' \"", type, "\", the response\n",
                "must be ", switch(type, linear = "numeric",
                          binary = "vector with 2 values",
                          count = "numeric with integer values",
                          general = "numeric",
                          survival = "a Survival object"), "."));

  corFncForCandidateCovariates = match.fun(corFncForCandidateCovariates);

  featureNames.original = colnames(x);
  if (is.null(dim(x)))  
  {
     x = data.frame(x = x);
  } else { 
     x = as.data.frame(x);
  }

  #keepCols = sapply(x, function(.x) !all(.x==.x[1]));
  #if (!all(keepCols)) 
  #x = x[, keepCols];
  #featureNames.original = featureNames.original[keepCols];

  # PL: Why this line? Because this drops all col- and row-names? 
  #    Edit: it also converts factors to numbers, so it makes more sense than as.matrix(x)
  # x = matrix(as.numeric(x), nrow(x), ncol(x))

  nFeatures.original = ncol(x);
  if (is.null(featureNames.original)) 
  {
     featureNames.original = featureNames = .spaste("F", .prependZeros(c(1:nFeatures.original)));
     colnames(x) = featureNames.original;
     namesChanged = FALSE;
     nameTranslationTable = data.frame(Column = c(1:nFeatures.original), 
                                       OriginalName = rep(NA, nFeatures.original),
                                       CoefficientName = featureNames);
  } else {
     featureNames = make.names(featureNames.original, unique = TRUE);
     if (isTRUE(all.equal(featureNames, featureNames.original)))
     {
       namesChanged = FALSE;
       nameTranslationTable = data.frame(Column = c(1:nFeatures.original), 
                                         OriginalName = featureNames,
                                         CoefficientName = featureNames);
     } else {
       namesChanged = TRUE;
       nameTranslationTable = data.frame(Column = c(1:nFeatures.original), OriginalName = featureNames.original,
                                    CoefficientName = featureNames);
       colnames(x) = featureNames;
     }
  }

  doTest = !is.null(xtest);
  if (doTest)
  { 
    xtestSaved = xtest;

    if (is.null(dim(xtest))) 
    {
       xtest = data.frame(x = xtest);
    } else
       xtest = as.data.frame(xtest);

    if (ncol(x)!=ncol(xtest))
      stop("Number of learning and testing predictors (columns of x, xtest) must equal.");

    if (!is.null(colnames(xtestSaved)))
    {
      if (!isTRUE(all.equal(colnames(xtestSaved), featureNames.original)))
        stop("Column names of 'x' and 'xtest' disagree.");
    } 
    colnames(xtest) = colnames(x)
    nTestSamples = nrow(xtest);
    # matrix for test set predicted values across bags 
    predictedTestMat =  matrix(NA, nTestSamples, nBags);
  } 

  # Handling of factors and other categorical features.

  # if categoricalColumns is not logical, turn it into logical indicator.
  if (!is.logical(categoricalColumns))
    categoricalColumns = c(1:nFeatures.original) %in% categoricalColumns;

  # Add columns deemed categorical to the indicator
  isCategorical.all = categoricalColumns | .categoricalColumns(x, maxLevels = maxCategoricalLevels);
  isFactor.original = sapply(x, is.factor);

  if (sum(isCategorical.all) > 0)
  {
    xLevels = lapply(x[isCategorical.all], .nonMissingLevels);

    # If test data were suplied, check whether the factor positions and levels agree between training and
    # test data 

    # Binarize the factors and adjust mandatory covariates, feature names, number of features etc.
    x = .binarizeCategoricalColumns(x, isCategorical.all, xLevels);
    bin2original = attr(x, "indexOfOrigin");
    
    if (!is.null(mandatoryCovariates)) 
    {
      mandatoryCovariates.bin = which(bin2original %in% mandatoryCovariates);
    } else
      mandatoryCovariates.bin = NULL;
    featureNames.bin = colnames(x)
    nFeatures.bin = ncol(x);
    isCategorical.bin = attr(x, "isCategorical");

    namesChanged = TRUE;
    nameTranslationTable = cbind(nameTranslationTable[bin2original, ], 
                                 BinarizedCoefficientName = featureNames.bin);
    if (doTest) 
    {
      isCategorical.all.test = categoricalColumns | 
                    .categoricalColumns(xtest, maxLevels = maxCategoricalLevels);
      if (any(isCategorical.all & !isCategorical.all.test))
        stop("All columns that are categorical in training data must also be categorical in test data."); 
      isFactor.test = sapply(xtest, is.factor);
      if (!isTRUE(all.equal(isFactor.original, isFactor.test))) 
         stop("Factor variables in training data 'x' must\n", 
              "also be factor variables in test data 'xtest' and vice-versa.")

      levels.test = lapply(xtest[isCategorical.all], .nonMissingLevels);

      consistent = mapply(function(x1, x2) {all(x2 %in% x1)}, xLevels, levels.test);
      if (!all(consistent))
      {
        stop(.spaste("The following categorical variables have test levels that", 
                    " do not appear in training data:\n",
             paste(colnames(x)[isCategorical.all][!consistent], collapse = ", ")));
      }
      xtest = .binarizeCategoricalColumns(xtest, isCategorical.all, xLevels)
    }
  } else {
    isCategorical.bin = rep(FALSE, nFeatures.original);
    nFeatures.bin = nFeatures.original
    featureNames.bin = featureNames.original;
    mandatoryCovariates.bin = mandatoryCovariates;
    xLevels = NULL;
  }

  nFeatures.eff = nFeatures.bin;
  featureNames.eff = featureNames.bin;
  isCategorical.eff = isCategorical.bin
  mandatoryCovariates.eff = mandatoryCovariates.bin;

  # Choose a response name that does not conflict with any of the column names in x.
  responseName = make.unique(c(featureNames.eff, "yResponse"))[nFeatures.eff + 1];

  # Important: turn the response into a matrix. This allows me to use one code for simple responses as well
  # as for survival (where the response is a matrix with 2 columns).

  if (type!="survival") y = as.matrix(y);

  if (nrow(y) != nrow(x))
      stop("x and y must have the same number of observations.")

  nSamples = nrow(y);

  if (!is.null(weights) && length(weights)!=nSamples) 
    stop("If 'weights' are given, they must have the same length (number of samples) as 'y' and 'x'.");

  if (nSamples < 8) 
  {
    .cat.nl("*****************************************************************\n",
            "* Warning in randomGLM: there are 7 or fewer observations.\n",
            "*   This may be too few to perform meaningful model selection\n", 
            "*   on in-bag (i.e., even fewer) samples.\n",
            "*   Model selection algorithm will likely output additional warnings.\n",
            "*   The resulting predictor should be used with caution.\n",
            "*****************************************************************");
  }

  nonMandatCovars.eff = setdiff( c(1:nFeatures.eff), mandatoryCovariates.eff);
  nMandatoryCovariates.eff = length(mandatoryCovariates.eff);


  # nFeaturesInBag shouldn't be greater than total number of features.
  if (nFeatures.eff<nFeaturesInBag)
  {
    nFeaturesInBag = ncol(x)
    .cat.nl("Warning in randomGLM: nFeaturesInBag is larger than the effective number of features.\n",
            "   Will use nFeaturesInBag equal to the number of features."); 
  }

  # nCandidateCovariates shouldn't be greater than nFeaturesInBag
  if (nCandidateCovariates>nFeaturesInBag)
  {
     .cat.nl("Warning in randomGLM: nCandidateCovariates is larger than nFeaturesInBag.\n",
             "  Will use nCandidateCovariates=nFeaturesInBag");
     nCandidateCovariates = nFeaturesInBag;
  }

  mandatoryCovarsGiven = length(mandatoryCovariates.eff)>0
  if (mandatoryCovarsGiven  & nCandidateCovariates <= nMandatoryCovariates.eff)
  {
    stop("Error: number of mandatoryCovariates >= nCandidateCovariates")
  }

  if (thresholdClassProb<0 | thresholdClassProb>1)
    stop("'thresholdClassProb' must be between 0  and 1.")

  # Check that the given family makes sense.

  if (type!="survival")
  {
    if (inherits(family, "character")) family = match.fun(family) ();
    if (!inherits(family, "family"))
       stop("Argument 'family' must be a valid 'family' object (see help(family))\n",
            "   or a character string specifying a valid family.");
  }

  # set seed
  if (!is.null(randomSeed))
  {
    if (exists(".Random.seed"))
    {
      saved.seed = .Random.seed;
      on.exit(.Random.seed <<- saved.seed);
    } 
    set.seed(randomSeed);
  }

  # matrix for predicted values in the training set 
  #predictedMat = matrix(NA, nSamples, nBags);
  # matrix for other outputs
  #featuresInForwardRegression = candidateFeatures = coefOfForwardRegression = list();
  #interceptOfForwardRegression = rep(NA, nBags)
  #bagObsIndx = matrix(NA, nBags, nObsInBag)
  #models = list();

  if (is.null(nThreads) || (nThreads > 1)) clusterInfo = .enableThreads(nThreads, verbose) 
     else clusterInfo = list(nThreads = nThreads);
  on.exit(.disableThreads(clusterInfo), add = TRUE)
  nThreads = clusterInfo$nThreads

  combinePredictors = function(...)
  {
    preds = list(...);
    out = list();

    out$predictedMat = sapply(preds, getElement, "predicted");
    if (doTest) 
    {
      out$predictedTestMat = sapply(preds, getElement, "predictedTest");
      #print(dim(out$predictedTestMat));
      # If the test set consists of 1 sample, this will be a vector that needs to be turned into a
      # single-row matrix
      if (is.null(dim(out$predictedTestMat)))
      {
        if (length(preds) > 1)
        {
           out$predictedTestMat = t(as.matrix(out$predictedTestMat));
        } else
           out$predictedTestMat = as.matrix(out$predictedTestMat);
      }
    }

    out$candidateFeatures = lapply(preds, getElement, "candidateFeatures");
    out$featuresInForwardRegression = lapply(preds, getElement, "featuresInForwardRegression");
    out$coefOfForwardRegression = lapply(preds, getElement, "coefOfForwardRegression");
    out$interceptOfForwardRegression = sapply(preds, getElement, "interceptOfForwardRegression");
    out$models = lapply(preds, getElement, "model");
    out;
  }

  # Prepare out of bag samples and in-bag features. This obviates the problems associated with splitting
  # random number generation across worker processes.

  # Since categorical predictors are binarized, they don't need any special attention - if they don't vary,
  # they will be removed from the bag later.

  bagFeatures = matrix(NA, nFeaturesInBag, nBags);
  bagObsIndx = matrix(NA, nObsInBag, nBags);
  nBagAttempts = 0;
  for (bag in 1:nBags)
  {
    yBagVar = 0;
    while (yBagVar==0)
    {
      nBagAttempts = nBagAttempts + 1;
      if (nBagAttempts > maxBagAttempts)
        stop("Bagging step failed. This could mean that factor variables have\n",
             "too many levels with too few samples at each level, or\n",
             "the response is (nearly) constant.");

      # sample indices for each bag
      bagSamples = sample(nSamples, nObsInBag, replace = replace, prob=sampleBaggingWeights);
      yBag = y[bagSamples, , drop = FALSE];
      yBagVar = var(yBag[, 1], na.rm=TRUE)
      # If there are no out-of-bag samples, force re-sampling as well
      # If the number of in-bag samples is less minInBagObs, re-sample again
      nUniqueInBag = length(unique(bagSamples));
      if (nUniqueInBag==nSamples | nUniqueInBag < minInBagObs) yBagVar = 0

      if (yBagVar == 0) next;

      if (replaceBadBagFeatures)
      {
        goodIndex = .goodFeatures( x[bagSamples, ], returnIndex = FALSE);
        mandatoryCovariates.bag = mandatoryCovariates.eff[goodIndex];
        nMandatoryCovariates.bag = length(mandatoryCovariates.bag);
      } else {
        mandatoryCovariates.bag = mandatoryCovariates.eff
        nMandatoryCovariates.bag = nMandatoryCovariates.eff
        goodIndex = c(1:nFeatures.eff);
      }

      featurePool = nonMandatCovars.eff[goodIndex];
      # If no features are "good", try bagging again.
      if (length(featurePool) == 0) { yBagVar = 0; next;}

      # If there aren't enough features to fill the bag completely, issue a warning.
      if (length(featurePool) < nFeaturesInBag - nMandatoryCovariates.bag)
        warning(.spaste("Bagging bag ", bag, ": number of valid bag features is not sufficient \n",
                     "  for sampling. This may impact the performance of the predictor."));

      bagFeatures.1 = c(mandatoryCovariates.bag,
                  sample(featurePool, min(length(featurePool), nFeaturesInBag - nMandatoryCovariates.bag)))

    }
    bagObsIndx[, bag] = bagSamples;
    bagFeatures[1:length(bagFeatures.1), bag] = bagFeatures.1;
  }

  # For convenience: the iteration over each bag is put in a macro-like function
  singleBagIteration = function(bag, verbose)
  {
    if (verbose>0) {.cat.nl("..bag ", bag)}
    #mem.last = NA;
    #if (verbose > 5) {m = gc()[2,2]; .cat.nl("  Step 1: ", m, ", diff: ", m - mem.last); mem.last = m}
    out = list();
    bagSamples = bagObsIndx[, bag];
    # keep track of in bag and oob samples
    oob = c(1:nSamples)[-unique(bagSamples)];
    nOOB = length(oob);
    features = bagFeatures[, bag]
    features = features[is.finite(features)];
    nFeatures.bag = length(features);

    xBag = x[bagSamples, features, drop = FALSE];
    yBag = y[bagSamples, , drop = FALSE];

    isCategorical.bag = isCategorical.eff[ features ];

    if (doTest)
    {
      xTestBag = rbind(x[oob, features, drop = FALSE], xtest[, features, drop = FALSE]);
    } else {
      xTestBag = x[oob, features, drop = FALSE];
    }
    # Here I only need to pass the number of mandatory covariates to function forwardSelection, because
    # they're saved at the beginning of features. 

    pr = .forwardSelection(xBag, yBag, xTestBag, 
                          weights.Bag = if (is.null(weights)) NULL else weights[ bagSamples],
                          classify=classify, 
                          binaryIndicator = isCategorical.bag,
                          maxInteractionOrder = maxInteractionOrder,
                          includeSelfinteractions = includeSelfinteractions,
                          nCandidateCovariates = nCandidateCovariates,
                          corFncForCandidateCovariates = corFncForCandidateCovariates, 
                          corOptionsForCandidateCovariates = corOptionsForCandidateCovariates, 
                          NmandatoryCovariates = nMandatoryCovariates.eff,
                          interactionsMandatory = interactionsMandatory,
                          keepModel = keepModels,
                          interactionSeparatorForCoefNames = interactionSeparatorForCoefNames,
                          type = type,
                          family = family,
                          responseName = responseName);

    #if (verbose > 5) {m = gc()[2,2]; .cat.nl("  Step 3: ", m, ", diff: ", m - mem.last); mem.last = m}
    # get output
    out$predicted = rep(NA, nSamples);
    out$predicted[oob] = pr$predicted[1:nOOB];
    if (doTest) {
      out$predictedTest = pr$predicted[(nOOB+1):(nOOB + nTestSamples)];
    }
    # to extract selected feature indices from feature name "feature?", use substring
    dictionary = rbind(c(0,0), cbind(1:nFeatures.bag, features));
    out$candidateFeatures = .translate(pr$candidateFeatures, dictionary);
    out$featuresInForwardRegression = .translate(pr$featuresInForwardRegression, dictionary);
    out$coefOfForwardRegression = pr$coefOfForwardRegression;
    out$interceptOfForwardRegression = pr$interceptOfForwardRegression

    if (keepModels) out$model = pr$model;
    # For large problems: this may be necessary to prevent exhausting the entire memory of the system.
    rm(xBag, pr, xTestBag, features, bagSamples, oob, yBag); 
    #if (verbose > 5) {m = gc()[2,2]; .cat.nl("  Step 4: ", m, ", diff: ", m - mem.last); mem.last = m}

    # Result value for each iteration.
    out
  }

  # loop over bags. Try two different version of the same code, one for parallel and one for serial
  # execution.
  if (nThreads > 1)
  {
    ensemble = foreach (bag = 1:nBags, .combine = combinePredictors, .multicombine = TRUE, 
                                       .maxcombine = nBags) %dopar%
      singleBagIteration(bag, verbose = verbose)
  } else {
    bagRes = list();
    for (bag in 1:nBags)
    {
      bagRes[[bag]] = singleBagIteration(bag, verbose = verbose);
      #if (verbose > 5) {tmp = gc(); .cat.nl("  In main loop: ", tmp[2,2]);}
      #.cat.nl("Size of information for each bag:");
      #print(object.size(bagRes[[bag]]), units = "auto");
    }
    ensemble = do.call(combinePredictors, bagRes);
  }

  featuresInForwardRegression.all = do.call(cbind, ensemble$featuresInForwardRegression);
  timesSelectedByForwardRegression = .countsInInteractionMatrix(featuresInForwardRegression.all, nFeatures.eff)

  colnames(bagObsIndx) = names(ensemble$candidateFeatures) = 
     names(ensemble$featuresInForwardRegression) = names(ensemble$coefOfForwardRegression) = 
     names(ensemble$interceptOfForwardRegression) = .spaste("Bag",1:nBags)

  if (ncol(x)==ncol(xSaved) && !is.null(colnames(xSaved)))
     colnames(timesSelectedByForwardRegression) = colnames(xSaved);

  # average predictive prob over bags, but if all bags give NA, the average should also give NA.
  predictedOOB.response1 = rowMeans(ensemble$predictedMat, na.rm = TRUE);
  predictedOOB.response1[rowSums(!is.na(ensemble$predictedMat))== 0] = NA;
  
  # recover original sample names
  if (!is.null(rownames(xSaved))) {
    names(predictedOOB.response1) = rownames(xSaved)
  }

  # prepare basic output
  out = list(predictedOOB.response = predictedOOB.response1,
             predictedOOB = predictedOOB.response1,
             candidateFeatures= ensemble$candidateFeatures,
             featuresInForwardRegression = ensemble$featuresInForwardRegression,
             coefOfForwardRegression = ensemble$coefOfForwardRegression,
             interceptOfForwardRegression = ensemble$interceptOfForwardRegression,
             bagObsIndx = bagObsIndx,
             timesSelectedByForwardRegression = timesSelectedByForwardRegression,
             models = if (keepModels) ensemble$models else NULL,
             featureNamesChanged = namesChanged,
             nameTranslationTable = nameTranslationTable,
             responseName = responseName,
             details = list(
                 type = type,
                 classify = classify,
                 family = family, 
                 nFeatures.original = nFeatures.original,
                 nFeatures.bin = nFeatures.bin,
                 nFeatures.eff = nFeatures.eff,
                 maxInteractionOrder = maxInteractionOrder,
                 yLevels = if (classify) originalYLevels else NULL,
                 maxCategoricalLevels = maxCategoricalLevels,
                 xLevels = xLevels,
                 categoricalColumns = categoricalColumns,
                 isCategorical.all = isCategorical.all,
                 isCategorical.bin = isCategorical.bin,
                 isCategorical.eff = isCategorical.eff,
                 mandatoryCovariates.bin = mandatoryCovariates.bin,
                 mandatoryCovariates.eff = mandatoryCovariates.eff,
                 isFactor.original = isFactor.original,
	         x.original = xSaved,
	         y.original = ySaved,
                 x = x,
                 y = y,
                 weights = weights,
                 thresholdClassProb = thresholdClassProb))
  # add test set output
  if (doTest) {
    predictedTest.response1 = rowMeans(ensemble$predictedTestMat, na.rm = TRUE);
    predictedTest.response1[rowSums(!is.na(ensemble$predictedTestMat))== 0] = NA;
    if (!is.null(rownames(xtestSaved))) {
       names(predictedTest.response1) = rownames(xtestSaved)
    }

    out$predictedTest.response = predictedTest.response1
    out$predictedTest = predictedTest.response1
  }

  # add output for binary outcomes
  if (classify)
  {
    predictedOOB = ifelse(predictedOOB.response1>thresholdClassProb, 1, 0)
    predictedOOB = originalYLevels[predictedOOB+1]

    predictedOOB.response = cbind(1-predictedOOB.response1, predictedOOB.response1)
    colnames(predictedOOB.response) = as.character(originalYLevels)

    out$predictedOOB = predictedOOB
    out$predictedOOB.response = predictedOOB.response
	
    if (doTest) {
      predictedTest = ifelse(predictedTest.response1>thresholdClassProb, 1, 0)
      predictedTest = originalYLevels[predictedTest+1]

      predictedTest.response = cbind(1-predictedTest.response1, predictedTest.response1)
      colnames(predictedTest.response) = as.character(originalYLevels)

      out$predictedTest.response = predictedTest.response;
      out$predictedTest = predictedTest;
    }
  } 
  class(out) = c("randomGLM", class(out));
  out
}

#============================================================================
#
# predict.randomGLM
#
#============================================================================

# Internal prediction function. This function will also be called from thin.randomGLM to generate prediction
# from the thinned predictor.

.predict.internal = function(object, newdata, type, thresholdClassProb,
                             returnBothTypes = FALSE)
{
  if (!is.null(newdata))
  {

    # Check that the structure of newdata is compatible with object.
    isFactor.new = sapply(as.data.frame(newdata), is.factor);
    if (!isTRUE(all.equal(as.vector(isFactor.new), as.vector(object$details$isFactor.original))))
      stop("is.factor on 'newdata' does not agree with factor indicator in original 'object'. \n",
           "  Features that were factors in the training data must\n",
           "  also be factors in test data and vice-versa.");

    # Check that names agree. It is not necessary to change the names since 
    # names on newdata bags are set from the interaction matrix. 
    if (!is.null(colnames(newdata)) && 
          (!all.equal(colnames(newdata), colnames(object$details$x.original))))
      stop("Column names of 'newdata' differ from column names of training data.\n", 
           "  Hint: If 'newdata' has column names, they must be the same\n",
           "  as the column names of the training data. \n",
           "  Mismatched column names often indicate mismatched variables which will\n",
           "  lead to wrong predictions or hard-to-understand errors.");

    if (any(object$details$isCategorical.all))
    {
      # Check that all levels of categorical variables present in newdata
      # were also present in the training data.

      levels.new = lapply(newdata[object$details$isCategorical.all], .nonMissingLevels);
      consistent = mapply(function(x1, x2) {all(x2 %in% x1)}, object$details$xLevels, levels.new);
      if (!all(consistent))
      {
        stop(.spaste("The following categorical variables have test levels that",
                    " do not appear in training data:\n",
             paste(colnames(object$details$x.original)[object$details$isCategorical.all][!consistent], 
                   collapse = ", ")));
      }

      # Binarize categorical columns of 'newdata'.
      newdata = .binarizeCategoricalColumns(newdata, object$details$isCategorical.all, 
                                                     object$details$xLevels);
    }

    nSamples = nrow(newdata)
    newdata.1 = cbind(newdata, rep(1, nSamples))
    # colnames(newdata) = make.names(colnames(newdata), unique = TRUE);
  } else {
    nSamples = nrow(object$details$y);
    x.1 = cbind(object$details$x, rep(1, nSamples));
  }

  nBags = length(object$models)

  predictedMat = matrix(NA, nSamples, nBags)

  for (b in 1:nBags) if (inherits(object$models[[b]], "lm"))
  {
    bagIM = object$featuresInForwardRegression[[b]];
    if (!is.null(newdata))
    {
      bagNewData = .generateInteractions(x = newdata, x1 = newdata.1, interactionMatrix = bagIM,
                                         setColNames = TRUE)
      
      predictedMat[, b] = .predictFromModel(object$models[[b]], newdata = as.data.frame(bagNewData), 
                                  modelType = object$details$type);
    } else {
      oob = c(1:nSamples)[-unique(object$bagObsIndx[, b])];
      bagNewData = .generateInteractions(x = object$details$x[oob, ], x1 = x.1[oob, ], interactionMatrix = bagIM,
                                         setColNames = TRUE)
      predictedMat[oob, b] = .predictFromModel(object$models[[b]], newdata = as.data.frame(bagNewData), 
                                  modelType = object$details$type);
    }
  }

  predicted.response = rowMeans(predictedMat, na.rm = TRUE)
  predicted.response[rowSums(!is.na(predictedMat))== 0] = NA

  names(predicted.response) = if (is.null(newdata)) {
              if (is.null(rownames(object$details$y.original))) rownames(object$details$x.original) else 
                                                        rownames(object$details$y.original) 
                                 } else rownames(newdata)

  if (type=="response" | returnBothTypes)
  {
    if (object$details$classify)
    {
      out.response = cbind(1-predicted.response, predicted.response)
      colnames(out.response) = as.character(object$details$yLevels)
    }else {
      out.response = predicted.response
    }
  }

  if (type=="class" | returnBothTypes)
  {
    # Note: type == "class" only makes sense for classification. 
    # For continuous prediction, put the continuous prediction here as well.
    if (object$details$classify)
    {
      predicted.round = ifelse(predicted.response>thresholdClassProb, 1, 0)
      out.class = object$details$yLevels[predicted.round+1]
    } else
      out.class = predicted.response;
  }

  if (returnBothTypes)
  {
    return(list(response = out.response, class = out.class));
  } else if (type=="class")
    return(out.class);

  out.response;
}
  
#===============================================================================================
#
# user-level predict() function
#
#===============================================================================================
  

predict.randomGLM = function(object, newdata, type=c("response", "class"), 
                             thresholdClassProb = object$details$thresholdClassProb, ...)
{
  type = match.arg(type)

  if (!is.null(object$binaryPredictors))
  {
    predictions = do.call(cbind, lapply(object$binaryPredictors, predict.randomGLM, 
                            newdata = newdata, type = type, thresholdClassProb = thresholdClassProb, ...))

    if (type=="response")
    {
       colnames(predictions) = colnames(object$predictedOOB.response);
    } else 
       colnames(predictions) = colnames(object$predictedOOB);
    return(predictions);
  }
    
  if (is.null(object$models))
    stop("The 'object' object must contain the undelying models for prediction.\n",
         "   Please re-run the randomGLM function with argument 'keepModels = TRUE' and try again.");

  # If new data is not given, return already calculated prediction. We would need the out-of-bag data for a
  # re-prediction and we don't have them.

  if (missing(newdata))
  {
    stop("valid 'newdata' must be given.")
  }
    
  if (ncol(newdata)!=object$details$nFeatures.original)
    stop("Number of columns in 'newdata' differs from the number of features\n",
         "     in the original training data.");


  if (type=="class" & !object$details$classify)
    stop("type='class' is only valid in classification.")

  if (thresholdClassProb<0 | thresholdClassProb>1)
    stop("Error: thresholdClassProb takes values between 0  and 1.")

  # Call the internal prediction function and return its result.
  .predict.internal(object = object, newdata = newdata, type = type, 
                    thresholdClassProb = thresholdClassProb, returnBothTypes = FALSE);
}

#============================================================================
#
# thin.randomGLM
#
#============================================================================

thinRandomGLM = function(rGLM, threshold)
{

  # Check if rGLM corresponds to multi-level repsonse.
  if (!is.null(rGLM$binaryPredictors))
  {
    out = rGLM;
    out$binaryPredictors = lapply(rGLM$binaryPredictors, thinRandomGLM, 
                               threshold = threshold);

    # Create the main prediction
    yBinNames = names(rGLM$binaryPredictors);
    yBinLevels = rGLM$levelMatrix
    nY = length(rGLM$binaryPredictors);

    names(out$binaryPredictors) = yBinNames;

    out$predictedOOB = as.matrix(sapply(out$binaryPredictors, getElement, "predictedOOB"));
    colnames(out$predictedOOB) = .spaste("PredictionFor.", yBinNames);

    out$predictedOOB.response = do.call(cbind,
                       lapply(out$binaryPredictors, getElement, "predictedOOB.response"));
    responseNames = .spaste(rep(yBinNames, rep(2, nY)), ".ProbabilityOfClass.", as.vector(yBinLevels));
    colnames(out$predictedOOB.response) = responseNames;

    rownames(out$predictedOOB.response) = rownames(out$predictedOOB) = rownames(rGLM$details$x.original);

    return(out);
  }

  nBags = length(rGLM$featuresInForwardRegression)
  if (nBags != length(rGLM$models))
  if (threshold<0)
    stop("'threshold' must be positive.")

  x = rGLM$details$x
  y = rGLM$details$y

  times = rGLM$timesSelectedByForwardRegression[1, ];
  if (threshold >= max(times))
    stop("Specified threshold removes all features in predictor.")

  ## so far, threshold only applies to no interaction level.
  keepF = which(times > threshold)
  screenF = function(input, keepF) { all(is.element(input, keepF)) }

  keepS = rGLM$bagObsIndx
  
  models = featuresInForwardRegression = coefOfForwardRegression = list();
  interceptOfForwardRegression = rep(NA, nBags);
  for (b in 1:nBags)
  {
    yBag = y[ keepS[,b], , drop = FALSE]
    xBag = x[ keepS[,b], , drop = FALSE]
    if (is.null(rGLM$details$weights)) 
    {
      wBag = NULL
    } else
      wBag = rGLM$details$weights[keepS[, b]];

    bagIM = rGLM$featuresInForwardRegression[[b]]
    keepIM = apply(bagIM, 2, screenF, keepF)
    if (sum(keepIM)==0)
    {
      featuresInForwardRegression[[b]] = models[[b]] = coefOfForwardRegression[[b]] = NA;
    } else {
      bagIM = bagIM[, keepIM, drop=FALSE]
      featuresInForwardRegression[[b]] = bagIM

      xBag.int = .generateInteractions(x = xBag, interactionMatrix = bagIM, setColNames = TRUE)

      modelData = cbind(data.frame(xBag.int), data.frame(yBag = yBag))
      colnames(modelData)[ncol(modelData)] = rGLM$responseName;
      models[[b]] = .createModel( formula = paste(rGLM$responseName, "~ ."),
                                  modelData = modelData,
                                  weights = wBag,
                                  type = rGLM$details$type,
                                  family = rGLM$details$family);

      coeffs = summary(models[[b]])$coefficients;
      keepCoeffs = .keepCoeffs(coeffs, rGLM$details$type);
      haveIntercept = .haveIntercept(coeffs, rGLM$details$type)
      coefOfForwardRegression[[b]] = coeffs[keepCoeffs,1]
      interceptOfForwardRegression[b] = if (haveIntercept) coeffs[1,1] else NA;

    }
  }
  
  featuresInForwardRegression.all = do.call(cbind, featuresInForwardRegression)
  timesSelectedByForwardRegression = .countsInInteractionMatrix(featuresInForwardRegression.all, 
                                                                rGLM$details$nFeatures.eff)
  if (!is.null(colnames(x))) {
    colnames(timesSelectedByForwardRegression) = colnames(x)
  }

  names(models) = names(featuresInForwardRegression) = names(coefOfForwardRegression) = .spaste("Bag",1:nBags)
  # Copy most of the information from the original rGLM since that information is left intact.
  out = rGLM;
  out$models = models;
  out$featuresInForwardRegression = featuresInForwardRegression;
  out$timesSelectedByForwardRegression = timesSelectedByForwardRegression;
  out$coefOfForwardRegression = coefOfForwardRegression;
  out$interceptOfForwardRegression = interceptOfForwardRegression;

  # Get predictions from the new models, both response and class.

  prediction = .predict.internal(out, newdata = NULL, type = "response",
                                 thresholdClassProb = rGLM$details$thresholdClassProb,
                                 returnBothTypes = TRUE);


  # Change the elements that need to be changed.
  
  out$predictedOOB.response = prediction$response;
  out$predictedOOB = prediction$class;

  if (!inherits(out, "randomGLM")) class(out) = c("randomGLM", class(out))
  out
}


