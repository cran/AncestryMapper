####### calculateAMids #######

calculateAMids <- function(pedtxtFile, fileReferences){

### Calculate sumPed from scratch, ped, map

FUN.sumPedG <- function(row){

	row <- as.numeric(row)
### Treat 0s and NAs

row.na <- which(is.na(row))
row.0 <- which(row==0)
row.change <- c(row.na, row.0)

### When row is NA or O - random substitution of 1 or 2

	row[row.change] <- sample(c(1,2), length(row.change), replace=T) 
	line.i <- matrix(unlist(row), length(row)/2, 2, byrow=T)
	line.i <- line.i[,1] + line.i[,2]
	return(line.i)
	}

### Import ped, faster with readLines

	ped <- readLines(pedtxtFile)

### strplit, each vector by ' '

	ped <- lapply(ped, function(row) unlist(strsplit(row, split=' ')))

### Individual names

	indNames <- unlist(lapply(ped, function(row) row[2])) ### row 2 because ped file is Id
### Check duplicates; only looks at the individual ids not fam ids

	if(sum(duplicated(indNames))!=0) stop("There are duplicate individuals in the input file.")

### Extract genotypes, remove pedigree informations

	ped <- lapply(ped, function(row) row[7:length(row)])

### Check only 1,2 and NAs

FUN.checkVectors <- function(row){

	pedUnique <- unique(row)
	if(!identical(setdiff(pedUnique, c(0,1,2,NA)), character(0)))
	stop("Input file can have 0, 1, 2 and NA values only.")
	return(NULL)
	}

### Check each vector has only 0, 1, 2 and NAs

	checkRows <- lapply(ped, FUN.checkVectors)

### Sum each genotype per snp; add snps; convert 2 cols into one

	sumPed <- lapply(ped, FUN.sumPedG)
	print('Summation of genotypes complete.')
	sumPed <- as.data.frame(do.call('rbind', sumPed))
	row.names(sumPed) <- indNames

### function to calculate distances; euclidean divided by number of snps
### maximum of distance 2, minimum 0

FUN.dist2Ref <- function(row){

	row <- as.numeric(row)
	distOut <- sqrt(as.vector(rowMeans(referencesMatrix^2) + mean(row^2) - 2 * referencesMatrix %*% row/ncol(referencesMatrix)))
	distOut <- round(distOut, 4)
	return(distOut)
	}

### Calculate indices: Assign 100 to most similar reference, 0 to least similar

FUN.indexes <- function(row.in){

	range.in <- diff(range(row.in))
	index <- round((row.in-min(row.in))/range.in*100,3)
	index <- abs(index-100)
	return(index)
	}

### Load the references file

	references <- read.table(fileReferences, header=T, as.is=T)

### References matrix: 51 references
### Error if references not in ped file format

	referencesMatrix <- as.matrix(sumPed[references[,'Reference'],])

	if(nrow(referencesMatrix)!=51)

		warning(paste('Number of references is usually 51.'), immediate.=T)

### Distances of every individual to references
### CORE of the script

	count <- 0
	dist2Ref <- apply(sumPed, 1, FUN.dist2Ref)
	print('Distances of Samples to References Computed.')
	dist2Ref <- as.data.frame(t(dist2Ref))

### Names; same as references

	names(dist2Ref) <- paste('C_', references$Pop)

### Indices; min:0 max:100

	Indexes <- data.frame(t(apply(dist2Ref, 1, FUN.indexes)))
	names(Indexes) <- gsub('C_', 'I_', names(Indexes))
	distAll <- data.frame(dist2Ref, Indexes)
	distAll <- data.frame(Id=row.names(distAll), distAll) 

### Return Distances

	return(distAll)
}

####### plotAMids #######

plotAMids <- function(AMids,
			  phenoFile='',
			  columnPlot='I',
			  quantilePlot='yes',
			  colorPlot='BlBrewer',
			  sepLinesPop='yes',
			  sepLineRef='yes',
			  plotIndNames='no',
			  legColor='yes',
			  legRef='yes',
			  legPheno='yes',
			  legAxisPop='yes'){

aid <- AMids

### Merge data from phenotype file if present

if(phenoFile!=''){
	Pheno <- read.table(phenoFile, header=T, as.is=T, comment.char='')
	aid <- merge(aid, Pheno, by.x='Id', by.y='UNIQID', all.x=T)
	aid <- aid[order(aid$Order),]
	}

aid <- aid[nrow(aid):1,]

### Colors taken from RColorBrewer and contrast
RedBl <- c('#2400D9','#191DF7','#2957FF','#3D87FF','#57B0FF',
	   '#75D3FF','#99EBFF','#BDF9FF','#EBFFFF','#FFFFEB',
	   '#FFF2BD','#FFD699','#FFAC75','#FF7857','#FF3D3D',
	   '#F72836','#D91630','#A60021')
RedBl <- rev(RedBl)
RedBlGr <- c(RedBl[1:9], gray.colors(100-length(RedBl)), RedBl[10:18])
BlBrewer <- c('#FFFFD9','#EDF8B1','#C7E9B4','#7FCDBB','#41B6C4',
	      '#1D91C0','#225EA8','#253494','#081D58','black')

### colImage depending on choice from user

colImage <- switch(colorPlot, BlBrewer=BlBrewer, RedBl=RedBl, BlBrewer=BlBrewer,RedBlGr)
plotI <- aid[grep(paste('^',columnPlot,'_', sep=''), names(aid), value=T)]
dimnames(plotI)[[1]] <- aid$Id

if(identical(quantilePlot,'yes')){
	breaks.in <- quantile(unlist(plotI), probs=seq(0,1, length.out=(length(colImage)+1)))
	breaks.in[1] <- breaks.in[1]-1
	breaks.in[length(breaks.in)] <- breaks.in[length(breaks.in)]+1
	}

if(identical(columnPlot,'I')){
	breaks.in <- c(seq(0,95,length.out=length(colImage)-1),99.9,101)
	breaks.in[1] <- -1
	}

### Produce layout: user input on what to plot

matPlot <- rep(0,18) ## 3 rows of 6 columns: 18 layouts
matPlot[7] <- 1
PhenoCols <- grep('Pheno_', names(aid), value=T)

if(identical(legPheno,'yes'))
	matPlot[8:10] <- switch(length(PhenoCols)+1, rep(0,3), c(2,0,0), c(2,3,0), 2:4)
	matPlot[1] <- switch(legColor, yes=max(matPlot)+1, 0)
	matPlot[13] <- switch(legRef, yes=max(matPlot)+1, 0)
	matLayout <- matrix(matPlot, ncol=6, byrow=T)
	widths.mat <- c(7,rep(0,4),0.5)
if(matLayout[2,2]!=0) widths.mat[2] <- 0.2
if(matLayout[2,3]!=0) widths.mat[3] <- 0.2
if(matLayout[2,4]!=0) widths.mat[4] <- 0.2
if(identical(legAxisPop,'yes')) widths.mat[5] <- 0.6
heights.mat <- c(0,6,0)
if(matLayout[1,1]!=0) heights.mat[1] <- 0.5
if(matLayout[3,1]!=0) heights.mat[3] <- 0.5

layout(matLayout,widths=widths.mat, heights=heights.mat)
mar1 <- c(5,3,3,0)
if(identical(plotIndNames,'yes')) mar1[2] <- 8 ## space to plot Individuals
par(mar=mar1)
if(quantilePlot=='yes'|columnPlot=='I')
	image(x=1:(ncol(plotI)+1), y=1:(nrow(plotI)+1), t(plotI), col=colImage, axes=F, ann=F, breaks=breaks.in)
if(columnPlot=='C'&quantilePlot=='no')
	image(x=1:(ncol(plotI)+1), y=1:(nrow(plotI)+1), t(plotI), col=colImage, axes=F, ann=F)
axis(1, at=seq(ncol(plotI))+0.5, labels=gsub('C_|I_','',names(plotI)), cex.axis=T, las=2, cex.axis=0.9)
if(identical(plotIndNames,'yes'))
	axis(2,at=seq(nrow(plotI))+0.5, labels=dimnames(plotI)[[1]],las=2,cex.axis=0.8)

### Separation between populations

colSepPop <- 'pink'
if(identical(grep('Pheno_1',names(aid),value=T),'Pheno_1')){
	sep.pop <- which(!duplicated(aid$Pheno_1))
if(identical(sepLinesPop,'yes')) abline(h=sep.pop, col=colSepPop, lwd=0.5)
	}

### Mark reference lines

sepLinesHGDP <- c(6, 10, 18, 25, 44, 49)
if(identical(sepLineRef,'yes')) abline(v=sepLinesHGDP+1, col='black')

### Image for phenotypes

mar2 <- mar1
mar2[c(2,4)] <- 0
par(mar=mar2)
cexBarPlot <- 0.7

### Make Colors if not existing; for each Pheno_
### Set1 of package RColorBrewer
colPheno <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00',
	      '#FFFF33','#A65628','#F781BF','#999999') 
### Colors of phenotype
PhenoCol <- grep('Pheno_',names(aid),value=T)
for(i in PhenoCol){
	col.i <- gsub('Pheno_','Colors_',i)
	if(!identical(grep(col.i,names(aid),value=T),col.i)){
		aid[,col.i] <- as.factor(aid[,i])
		levels(aid[,col.i]) <- rep(colPheno,length(levels(aid[,col.i])))[1:length(levels(aid[,col.i]))]
		aid[,col.i] <- as.character(aid[,col.i])
		}
	}

if(identical(grep('Pheno_3',names(aid),value=T),'Pheno_3')){
	barplot(rep(1, nrow(aid)), col=aid[,'Colors_3'], beside=F, border=NA, space=0, horiz=T, axes=F, yaxs='i')
if(identical(sepLinesPop,'yes')) abline(h=sep.pop, col=colSepPop, lwd=0.5)
	axis(1,0.5,'Pheno 3',cex.axis=cexBarPlot,las=2,xpd=NA)
}

if(identical(grep('Pheno_2',names(aid),value=T),'Pheno_2')){
	barplot(rep(1, nrow(aid)), col=aid[,'Colors_2'], beside=F, border=NA, space=0, horiz=T, axes=F, yaxs='i')
if(identical(sepLinesPop,'yes')) abline(h=sep.pop, col=colSepPop, lwd=0.5)
	axis(1,0.5,'Pheno 2',cex.axis=cexBarPlot,las=2,,xpd=NA)
	}

if(identical(grep('Pheno_1',names(aid),value=T),'Pheno_1')){
	barplot(rep(1, nrow(aid)), col=aid[,'Colors_1'], beside=F, border=NA, space=0, horiz=T, axes=F, yaxs='i')
if(identical(sepLinesPop,'yes')) abline(h=sep.pop, col=colSepPop, lwd=0.5)
	axis(1,0.5,'Pheno 1',cex.axis=cexBarPlot,las=2,xpd=NA)
	}

if(identical(grep('Pheno_1',names(aid),value=T),'Pheno_1')){
	phenoLeg <- data.frame(Pheno_1=aid$Pheno_1, Index=seq(nrow(aid)))
	phenoLeg <- phenoLeg[!duplicated(phenoLeg$Pheno_1) | !duplicated(phenoLeg$Pheno_1, fromLast=T) ,]
	phenoLeg$time <- 1:2
	phenoLeg <- reshape(phenoLeg, idvar='Pheno_1', v.names='Index', timevar='time', direction='wide')
	phenoLeg$loc <- (phenoLeg$Index.1+phenoLeg$Index.2)/2
	axis(4,at=phenoLeg$loc, labels=phenoLeg$Pheno_1, xpd=NA, las=2, cex.axis=0.9) ## 0.7 / 1.2
	}

### Plot top legend: color gradient
mar3 <- mar1
mar3[c(1,3)] <- 0
par(mar=mar3)
plot(1,1, t='n', axes=F, ann=F,xlim=c(0,length(colImage)*3),xaxs='i')
x1 <- seq(length(colImage))
rect(x1,0.5,x1+1,0.8,col=colImage,xpd=NA)
text(x1[1],0.3,'low',xpd=NA)
text(x1[length(x1)],0.3,'high',xpd=NA)

### Many lines of text; place text at the bottom and sides to indicate continental block

popTextDown <- c(sepLinesHGDP,51)

### separation of big regions in HGDP: Afr, MiddleEast, Eur, CSA, EA, America, Oceania

popText <- c('Africa', 'ME', 'Europe', 'C S Asia', 'East Asia', 'Amer', 'Oce')
col.text <- 'black'
par2 <- mar1
par2[c(1,3)] <- 0
plot(1,1,xlim=c(0,51),ylim=c(1,100), xaxs='i',axes=F,ann=F,t='n')
sepPopLine <- 0.5
lwd.i <- 2
y.seg.pos <- 100
y.text.pos <- 70
cex.i <- 1
x1 <- c(sepPopLine,popTextDown[-length(popTextDown)]+sepPopLine)
x2 <- c(popTextDown-sepPopLine)
segments(x1, y.seg.pos, x2, y.seg.pos, col=col.text, lwd=lwd.i)
textAll <- ((x1+x2)/2)
text(textAll, y.text.pos, popText,cex=cex.i, col=col.text)

}

