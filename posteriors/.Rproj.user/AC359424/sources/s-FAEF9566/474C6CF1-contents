\name{likelihood_param_quad_function}
\alias{likelihood_param_quad_function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GLM parameter generation function.
%%  ~~function to do ... ~~
}
\description{
Generates the glm model parameters used in the computation of likelihood of the target ODn.using a quadratic function in time since HIV infection
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
likelihood_param_quad_function(dat, target_ODn, target_ODnstepsize, t_since_ln)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
  A dataset of Pr_t values by ODn threshold and per day over 1000 daystime that's brought in as a matrix.
%%     ~~Describe \code{x} here~~
}
\item{target_ODn}{
Patient ODn as measurement, after the first HIV positive test.
%%     ~~Describe \code{x} here~~
}
\item{target_ODnstepsize}{
Resolution around the target ODn. This depends on the dynamic range of the assay. We recommend that for assays like Sedia LAg, whose highest values are around 7, we use ODn step sizes of 0.1.
%%     ~~Describe \code{x} here~~
}
\item{t_since_ln}{
Timepoints within the intertest interval at which the function should be evaluated. smaller timeframes amy take longer time to run depending on your computer specifications. The P_r(t) data hass descrete time steps.
%%     ~~Describe \code{x} here~~
}
}
\details{
We evaluate the gardient around 6 ODn values (including the target ODn). The 5 ODn varles are chosen as the closest values, based on a specific target ODn step size, around the target ODn.
The function uses the P_r(t) data generated for the assay (loaded as part of the package) in a quadratic GLM function for each of the 6 chosen target ODn (including the ODn in question) and ouputs slope parameters which are to be used later in the likelihood_fun function.
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
