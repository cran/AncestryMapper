\name{plotAMids}
\alias{plotAMids}
\title{Visualise genetic distances}
\description{
plotAMids is used to visualise the relationship amongst individuals and the 51 HGDP references.
}
\usage{
plotAMids(AMids,
          phenoFile = "",
          columnPlot = "I",
          quantilePlot = yes,
          colorPlot = "BlBrewer",
          sepLinesPop = yes,
          sepLineRef = yes,
          plotIndNames = 'no',
          legColor = yes,
          legRef = yes,
          legPheno = yes,
          legAxisPop = yes)
}
\arguments{
  \item{AMids}{
Dataframe of genetic distances calculated by calculateAMids.
}
  \item{phenoFile}{
Optional file with phenotype information for each individual. Columns should be named 'UNIQID' 'Pheno_1' 'Colors_1' and 'Order'. If provided each individual is assigned a color code to provide information on the phenotype. If "Colors_" are not provided the function will randomly assign colors to each phenotype (see HGDP_Phenotypes). Three different phenotypes are possible, in order of size: e.g., Pheno_1 should be population, Pheno_2 Continental Region, Pheno_3 the dataset. Order is used to order the individuals in the plot.
}
  \item{columnPlot}{
Takes values 'I' or 'C'. 'I' is the default option. 
'I' plots the normalised euclidean distances whereas 'C' plots the crude distances.
}
  \item{quantilePlot}{
Takes values 'yes' or 'no'. 'yes' is the default option. 
If columnPlot is 'C', 'yes' will plot the quantiles, 'no' will plot the raw values
}
  \item{colorPlot}{
Colors for the AMids. Possible choices are 'RedBl', 'RedBlGr' and 'BLBrewer'. 
The user can also provide a vector of clors.
}
  \item{sepLinesPop}{
Takes values 'yes' or 'no'. The default is 'yes'. 
If 'yes', a line demarcating populations is plotted.
}
  \item{sepLineRef}{
Takes values 'yes' or 'no'. The default is 'yes'. 
If 'yes', a line demarcating continental regions for the HGDP populations is plotted.
}
  \item{plotIndNames}{
Takes values 'yes' or 'no'. The default is 'no'. 
If 'yes', the individual ids are plotted on the left axis.
}
  \item{legColor}{
Takes values 'yes' or 'no'. The default is 'yes'. 
If 'yes', the legend colour will be plotted in the top left.
}
  \item{legRef}{
Takes values 'yes' or 'no'. The default is 'yes'. If 'yes', indications of continental region will be plotted.
}
  \item{legPheno}{
%%     ~~Describe \code{legPheno} here~~
}
  \item{legAxisPop}{
Takes values 'yes' or 'no'. The default is 'yes'. 
If 'yes' and phenotype information is provided, this will be plotted on the right axis of the plot.
}
}
\details{
plotAMids takes as input the dataframe of genetic distances returned by calculateAMids. For assignment of phenotype the user should provide a file with colors and phenotype.
}
\value{
Return a plot of the relationship of a sample population to the HGDP references to the R plotting device.
}
\references{
Magalh\~aes et al. (2012) PLOS One (accepted)
}
\author{
Tiago R. Magalh\~aes, Darren J. Fitzpatrick
}
\note{
Ancestry analysis using a global genomic map.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
calculateAMids %% ~~objects to See Also as \code{\link{help}}, ~~~
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
                               'HGDP_References.txt',
                               package='AncestryMapper')

genetic.distance <- calculateAMids(pedtxtFile=HGDP_500SNPs,
                                   fileReferences=HGDP_References)
plotAMids(AMids=genetic.distance, phenoFile=HGDP_Phenotypes)}
}
\keyword{ genetic distance }
\keyword{ ancestry mapping }% __ONLY ONE__ keyword per line
\keyword{ population genetics }
