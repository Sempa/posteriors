\name{pr_t_fun_loglog_cubic}
\alias{pr_t_fun_loglog_cubic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A cloglog function to generate P_r(t) function.
%%  ~~function to do ... ~~
}
\description{
Pr_t function used to generate Pr_t values using a logit link function. This function receives a dataset, ODn threshold or ODn and interval between last negative to first positive.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
func_logit_cubic(data_set, ODn_step, t_since_inf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_set}{
  Patient dataset to be evaluated
%%     ~~Describe \code{x} here~~
}
  \item{ODn_step}{
  A vector of ODn thresholds to be evaluated
%%     ~~Describe \code{x} here~~
}
\item{t_since_inf}{
  Timepoints within the intertest interval at which the function should be evaluated.
%%     ~~Describe \code{x} here~~
}
}
\details{
This function generates P_r(t) data at each ODn threshold for every 1000 daily step size time interval. The ODn step size is 0.01 from 0.01 to 5.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Joseph B. Sempa
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

dt <- data.frame(pr_t = pr_t_fun_loglog_cubic(
data_set = data_generate_pr_t,
ODn_step = seq(0.01, 5, 0.01),
t_since_inf = seq(0, 1000, 1))
)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
