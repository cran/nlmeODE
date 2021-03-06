\name{nlmeODE}
\alias{nlmeODE}
\title{Non-linear mixed-effects modelling in nlme using differential equations}
\description{This package combines the \code{deSolve} and \code{nlme} packages for mixed-effects modelling.}
\usage{nlmeODE(model,data,LogParms,JAC,SEQ,rtol,atol,tcrit,hmin,hmax)}

\arguments{
\item{model}{A list including the following elements: \code{DiffEq} (A list of formulas containing the ODE's for the system in the same order as the compartment numbers, i.e. formula i is the ODE for compartment i), \code{ObsEq} (A formula specifying which state is observed along with possible scaling parameters),  \code{Parms} (A vector with the names of the parameters used in \code{DiffEq} followed by the parameters in \code{ObsEq} and the initial state parameters), \code{States} (A vector with the names of the states in \code{DiffEq}), \code{Init} (A logical vector with the same length as \code{States} specifying whether initial state estimates should be obtained for the particular state when TRUE)}
\item{data}{\code{groupedData} object with a formula specifying which columns are the dependent and independent variables, and grouping factor. Optional columns in the \code{groupedData} object are the dose \code{Dose}, dosing compartment \code{Cmt}, rate of infusion \code{Rate}, covariate measurements, etc.}
\item{LogParms}{If TRUE, the parameters are reparameterized in terms of the logarithm of the parameters. Default is TRUE.}
\item{JAC}{If TRUE, the Jacobian of the system of ODE's is computed and passed to the ODE solver. In some circumstances, supplying the Jacobian can speed up the computations if the system is stiff. Default is FALSE.}
\item{SEQ}{A logical value, that when TRUE, adds a gradient attribute to the returned value calculated by simultaneous solution of the sensitivity equations associated with the system of ODE's. Default is FALSE.}
\item{rtol}{Relative error tolerance for lsoda. Default is 1E-4.}
\item{atol}{Absolute error tolerance for lsoda. Default is 1E-4.}
\item{tcrit}{Time beyond which the integration should not proceed. Default is NULL.}
\item{hmin}{Minimum value of the integration stepsize. Default is 0.}
\item{hmax}{Maximum value of the integration stepsize. Default is Inf.}
}
\value{A function compatible with \code{nlme} estimation for systems of ordinary differential equations (ODE's).}
\note{For examples of commonly used PK/PD models see \code{\link{PKPDmodels}}.}

\author{Christoffer W. Tornoe <\email{ctornoe@gmail.com}>}
\seealso{\code{\link[nlme]{nlme}}, \code{\link[deSolve]{lsoda}}}

\references{Tornoe, C. W. et al. (2004a) "Non-linear mixed-effects pharmacokinetic/pharmacodynamic modelling in NLME using differential equations", Computer Methods and Programs in Biomedicine, 76(1), 31-40 and Tornoe, C. W. et al. (2004b) "Pharmacokinetic/Pharmacodynamic Modelling of GnRH Antagonist Degarelix: A Comparison of the Non-linear Mixed-Effects Programs NONMEM and NLME",Journal of Pharmacokinetics and Pharmacodynamics, 31(6), 441-461.}

\keyword{models}
