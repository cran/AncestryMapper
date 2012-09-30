\name{AncestryMapper}
\alias{AncestryMapper}
\docType{package}
\title{Ancestry Mapper
}
\description{
Package computes the genetic distance, defined as the euclidean distance amongst a sample set of individuals and the 51 HGDP reference panels. The package allows for the visualisation of the relationship of sample individuals to the 51 HGDP reference populations, thus permitting the user to assess the relationship of individuals to different geographic groupings.
}
\details{
\tabular{ll}{
Package: \tab AncestryMapper\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-09-27\cr
License: \tab None\cr
LazyLoad: \tab yes\cr
}

}
\author{
Tiago R. Magalh\~aes, Darren J. Fitzpatrick

Maintainer: Tiago R. Magalh\~aes  <tiago17@gmail.com>
}
\references{
Magalh\~aes TR et al. (2012) PLoS ONE 7(11): e49438. doi:10.1371/journal.pone.0049438 
}
\keyword{ package }
\keyword{ genetic distance }
\keyword{ ancestry mapping }
\keyword{ population genetics }
\seealso{
image dist
}
\examples{
\dontrun{
library(AncestryMapper)
HGDP_References <- system.file('extdata',
                               'HGDP_References.txt',
                                package='AncestryMapper')

HGDP_500SNPs <- system.file('extdata',
                            'HGDP_500SNPs.ped',
                             package='AncestryMapper')

HGDP_Phenotypes <- system.file('extdata',
                               'HGDP_Phenotypes.txt',
                                package='AncestryMapper')

genetic.distance <- calculateAMids(pedtxtFile=HGDP_500SNPs,
                                   fileReferences=HGDP_References)
plotAMids(AMids=genetic.distance, phenoFile=HGDP_Phenotypes)}
}
