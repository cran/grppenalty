\name{cv.grppenalty}
\Rdversion{2.0}
\alias{cv.grppenalty} 
\title{Tuning parameter selection for the concave 1-norm and 2-norm group penalties}
\description{Tuning parameter selection for the concave 1-norm and
  2-norm group penalties by k-fold cross validation. }
\usage{
cv.grppenalty(y, x, index, family = "gaussian", type = "l1",
penalty = "mcp", kappa = 1/2.7, nfold = 5,
regular.only = TRUE, nlambda = 100, lambda.min = 0.01,
epsilon = 1e-3, maxit = 1e+3, seed = 1000)
}
\arguments{
  \item{y}{outcome of interest. A vector of continuous response in
    linear models or a vector of 0 or 1 in logistic models.} 
  \item{x}{the design matrix of penalized variables. By default, an intercept
    vector will be added when fitting the model.}
  \item{index}{group index of penalized variables.}
    \item{family}{a character indicating the distribution of
    outcome. Either "gaussian" or "binomial" can be specified.}
  \item{type}{a character specifying the type of grouped penalty. Either
    "l1" or "l2" can be specified, with "l1" being the default. See
    following details for more information. } 
  \item{penalty}{a character specifying the penalty. One of "mcp" or
    "scad" should be specified, with "mcp" being the default. }
  \item{kappa}{the regularization parameter kappa, either one value or
    an increasing vector of values can be specified. The value of kappa
    should be in the range of [0,1). }
  \item{nfold}{the k value for the k-fold cross validation.}
  \item{regular.only}{when selecting the tuning parameter, should the
    selection process be limited to the regular solution only or not? The
    regular solution refers to the models with df<n, with df the degree
    of freedom/non-zero coefficients, and n the sample size.}
  \item{nlambda}{a integer value specifying the number of grids along the
    penalty parameter lambda.}
  \item{lambda.min}{a value specifying how to determine the minimal value of
    penalty parameter lambda. We define lambda_min=lambda_max*lambda.min.
    We suggest lambda.min is 0.0001 if n>p; 0.01 otherwise.}
  \item{epsilon}{a value specifying the converge criterion of algorithm.}
  \item{maxit}{an integer value specifying the maximum number of iterations for
    each coordinate.}
  \item{seed}{a integer which is used as the random seed for generating
    cross-validation index.}
}
\details{
  The package implements the concave 1-norm and 2-norm group penalties in
  linear and logistic regression models. The concave 1-norm group
  penalty is defined as rho(|beta|_1;d*lambda,kappa) with |beta|_1 being
  the L1 norm of the coefficients and d being the group size. The
  concave 2-norm group penalty is defined as
  rho(|beta|_2;sqrt(d)*lambda,kappa) with |beta|_2 being the L2 norm of
  the coefficients. Here rho() is the concave function, in current
  implementation, we only consider the smoothly clipped  absolute deviation
  (SCAD) penalty  and minimum concave penalty (MCP).

  The concave 1-norm group penalties, i.e. 1-norm gSCAD or gMCP, perform
  variable selection at group and individual levels under proper tuning
  parameters. The concave 2-norm group penalties, i.e. 2-norm gSCAD or
  gMCP selects variable at group level, i.e. the variables in the same
  group are dropped or selected at the same time. One advantage of of
  the 1-norm group penalty is that it is robust to mis-specified group
  information. The 2-norm group penalty is, however, affected by the
  mis-specified group information. The concave 2-norm group penalty
  includes group Lasso as a special case when the regularization
  parameter kappa=0. Hence, setting kappa=0 in the 2-norm group penalty
  returns the group Lasso solutions.
  
  We use the coordinate descent algorithm (CDA) to compute the solution for
  both the 1-norm and 2-norm group penalties. The solution path is
  computed along kappa. That is we use the solution at kappa=0 to
  initiate the computation for a given penalty parameter lambda. In
  general, we suggest treating both the regularization parameter kappa
  and penalty parameter lambda as tuning parameters and use data-driven
  approach to select optimal kappa and lambda. However, this practice requires
  heavy computation, thus, a  particular kappa (1/2.7 for gMCP and 1/3.7
  for gSCAD) is recommended to reduce the computational time. 

  For tuning parameter selection, we implement the cross-validation
  approach for both linear and logistic models. In linear model, we use
  the predictive mean square error (PMSE) as the index quantity. The
  tuning parameter(s) corresponding to the solution with the minimum pmse
  is selected. In logistic model, the k-fold cross-validated area under ROC curve
  (CV-AUC) is used as the index quantity. The tuning parameter(s)
  corresponding to the solution with the maximum CV-AUC is selected.
}
\value{
  A list of four elements is returned.
  \item{selected}{a list of 4 elements corresponding to the selected
  tuning parameter, with the 1st is the CV-PMSE for linear model and CV-AUC
  for logistic model, the 2nd is the regression coefficients, the 3rd is
  kappa value and the 4th is the lambda value. }
  \item{cv.values}{a matrix corresponding to the predictive
    performance in the CV process, mainly for plotting purpose.}
  \item{kappas}{a matrix of regularization parameter kappa corresponding
  to the cv.values.}
  \item{lambdas}{a matrix of penalty parameter lambdas corresponding to
  the cv.values.}
}
\author{Dingfeng Jiang}
\references{Jiang, D., Huang, J., Zhang, Y. (2011). The cross-validated
  AUC for MCP-Logistic regression with high-dimensional
  data. \emph{Statistical Methods in Medical Research}, online first.
  
  Yuan, M., Lin, Y. (2006). Model selection and estimation in regression
  with grouped variables. \emph{Journal of Royal Statistical Society
    Series B}, 68 (1): 49 - 67.
  
  Meier, L., van de Geer, S., B\''{u}hlmann, P., (2008). The group lasso
  for logistic regression. \emph{Journal of Royal Statistical Society
    Series B}, 70 (1): 53 - 71
}
\seealso{\code{grppenalty}}
\examples{
set.seed(10000)
n=100
ybi=rbinom(n,1,0.4)
yga=rnorm(n)
p=20
x=matrix(rnorm(n*p),n,p)
index=rep(1:10, each =2)
out=cv.grppenalty(yga, x, index, "gaussian", "l1", "mcp",  1/2.7)
## out=cv.grppenalty(yga, x, index, "gaussian", "l2", "mcp", 1/2.7)
## out=cv.grppenalty(yga, x, index, "gaussian", "l1", "scad",  1/2.7)
## out=cv.grppenalty(yga, x, index, "gaussian", "l2", "scad",  1/2.7)
## multiple kappas
## out=cv.grppenalty(yga, x, index, "gaussian", "l1", "mcp",  c(0,0.1,1/2.7))

## out=cv.grppenalty(ybi, x, index, "binomial", "l1", "mcp", 1/2.7)
## out=cv.grppenalty(ybi, x, index, "binomial", "l2", "mcp",  1/2.7)
## out=cv.grppenalty(ybi, x, index, "binomial", "l1", "scad",  1/2.7)
## out=cv.grppenalty(ybi, x, index, "binomial", "l2", "scad",  1/2.7)
}
\keyword{models}
\keyword{regression}
