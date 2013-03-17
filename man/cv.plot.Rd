\name{cv.plot}
\Rdversion{2.0}
\alias{cv.plot}
\title{Plot the cross validation performance}
\description{Plot the cross validation performance, for linear model
  CV-PMSE will be plotted, for logistic model CV-AUC will be plotted. }
\usage{
cv.plot(cv.out)
}
\arguments{
  \item{cv.out}{the object from the cv.grppenalty function}
}
\details{
  The cv.plot shows the cross validation performance relative to the
  kappa and lambda. This is to visualize the overall cross validation
  process.
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
index=rep(1:10, each = 2)
## cv.out=cv.grppenalty(yga, x, index, "gaussian", "l1", "mcp",  1/2.7)
## cv.plot(cv.out)
## multiple kappas
cv.out=cv.grppenalty(yga, x, index, "gaussian", "l1", "mcp",  c(0,0.1,1/2.7))
cv.plot(cv.out)
}
\keyword{models}
\keyword{regression}
