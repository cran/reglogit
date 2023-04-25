\name{pima}
\alias{pima}
\docType{data}
\title{
Pima Indian Data
}
\description{
A population of women who were at least 21 years old,
 of Pima Indian heritage and living near Phoenix, Arizona, 
was tested for diabetes according to World Health Organization 
criteria. The data were collected by the US National Institute 
of Diabetes and Digestive and Kidney Diseases. 
}
\usage{data(pima)}
\format{
  A data frame with 768 observations on the following 9 variables.
  \describe{
    \item{\code{npreg}}{number of pregnancies}
    \item{\code{glu}}{plasma glucose concentration in an oral glucose
    tolerance test}
    \item{\code{bp}}{diastolic blood pressure (mm Hg)}
    \item{\code{skin}}{tricepts skin fold thickness (mm)}
    \item{\code{serum}}{2-hour serum insulin (mu U/ml)}
    \item{\code{bmi}}{mody mass index (weight in kg/(height in m)^2)}
    \item{\code{ped}}{diabetes pedigree function}
    \item{\code{age}}{age in years}
    \item{\code{y}}{classification label: \code{1} for diabetic}
  }
}
\source{
Smith, J. W., Everhart, J. E., Dickson, W. C., Knowler, W. C. and
Johannes, R. S. (1988) 
\emph{Using the ADAP learning algorithm to forecast the onset of
 diabetes mellitus. In Proceedings of the Symposium on Computer 
Applications in Medical Care (Washington, 1988)}, ed. 
R. A. Greenes, pp. 261-265. Los Alamitos, CA: IEEE Computer Society Press.

Ripley, B.D. (1996) \emph{Pattern Recognition and Neural Networks.}
 Cambridge: Cambridge University Press.

UCI Machine Learning Repository\cr
\url{http://archive.ics.uci.edu/ml/datasets/Pima+Indians+Diabetes}
}

\examples{
data(pima)
## see reglogit documentation for an example using this data
}
\keyword{datasets}
