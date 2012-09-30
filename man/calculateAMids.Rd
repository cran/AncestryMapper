\name{calculateAMids}
\alias{calculateAMids}
\title{Calculate genetic distances}
\description{Calculates genetics distance between a sample population and HGDP references.}
\usage{
calculateAMids(pedtxtFile, fileReferences)
}
\arguments{
  \item{pedtxtFile}{Name of the ped file to be used. The ped file should include all 51 HGDP references and the individuals for which the user wishes to calculate the genetic distance.
}
  \item{fileReferences}{
File detailing the 51 references and their order of appearance in the plot. This file, HGDP references is provided with the package.
}
}
\details{
Program computes the genetic distances amongst individuals and the 51 HGDP references.
}
\value{
calculateAMids returns a dataframe containing the genetic distance of each sample individual to the 51 HGDP references. This dataframe is subsequently used as input to the plotAMids function.
}
\references{
Magalh\~aes TR et al. (2012) PLoS ONE 7(11): e49438. doi:10.1371/journal.pone.0049438
}
\author{
Tiago R. Magalh\~aes, Darren J. Fitzpatrick
}
\note{
Ancestry analysis using a global genomic map. %%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
plotAMids
}
\examples{
\dontrun{

library(AncestryMapper)

HGDP_References <- system.file('inst', 'extdata',
                               'HGDP_References.txt',
                                package='AncestryMapper')

HGDP_500SNPs <- system.file('inst', 'extdata',
                            'HGDP_500SNPs.ped',
                             package='AncestryMapper')

HGDP_Phenotypes <- system.file('inst', 'extdata',
                               'HGDP_Phenotypes.txt',
                                package='AncestryMapper')

genetic.distance <- calculateAMids(pedtxtFile=HGDP_500SNPs,
                                   fileReferences=HGDP_References)
plotAMids(AMids=genetic.distance, phenoFile=HGDP_Phenotypes)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{genetic distance}
\keyword{ancestry mapping}
\keyword{population genetics}
