\name{path.plot}
\Rdversion{2.0}
\alias{path.plot}
\title{Plot the solution path for the concave 1-norm and 2-norm group penalties}
\description{Plot the coefficient profiles of each penalized variable for the
  concave 1-norm and 2-norm group penalties. }
\usage{
path.plot(out)
}
\arguments{
  \item{out}{the object return from the grppenalty function}
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

  The solution path plot shows the coefficients profile of variables in
  x. It illustrates the effect of penalization. Depending on the kappa
  value specified, a total of nkappa figures will be plotted.
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
out=grppenalty(yga, x, index, "gaussian", "l1", "mcp",  c(0,1/2.7))
path.plot(out)

}
\keyword{models}
\keyword{regression}
