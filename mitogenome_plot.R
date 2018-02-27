#########################################################
#### Plotting circular genomes in R using {circlize} ####
#########################################################

# Get the tools
#install.packages("circlize")
#install.packages("seqinr")

### Example of getting set up to run function:


#seq = read.fasta(file.choose()) # Read in sequence data
#seq = seq[[1]][1:16544] # and pull out just the sequence stored in a character vector

#plus = read.csv(file.choose()) # read in plus coords and colors
#minus = read.csv(file.choose()) # read in plus coords and colors

#coverage = read.csv(file.choose()) # read in coverage per base
#coverage = coverage[,1] # make sure its a vector


mitoplotR <- function (seq, plus, minus, coverage){

library(circlize)
circos.clear() 

# For GC content across mitogenome we need to first define some funcitons These are fromhttps://www.cs.us.es/~fran/students/julian/introduction/introduction.html
# Thanks, Julian and Francisco.
  
                          getMultinomialModel <- function(sequence) {
                            frequencies <- getFrequency(sequence)
                            sequence.length <- length(sequence)
                            model <- frequencies / sequence.length
                            names(model) <- names(frequencies)
                            return(model)
                          }
                          
                          getFrequency <- function(sequence, alphabet = c("a","c","g","t"))
                          {
                            freq <- table(sequence)[alphabet]
                            names(freq) <- alphabet
                            # evaluates the value for every symbol of the alphabet
                            for(bp in alphabet){
                              if(is.na(freq[bp])) freq[bp] <- 0
                            }
                            return(freq)
                          }
                          
                          getSequenceLocalData <- function(sequence,window.length,offset)
                          {
                            # initializes the window, sequence and output parameters
                            window.lowest <- 1
                            window.highest <- window.length
                            sequence.length <- length(sequence)
                            result <- double()
                            
                            # while loop to slide the window along the sequence
                            while(window.highest <= sequence.length)
                            {
                              # computes the local base composition, and then adds the GC content
                              # and the position to result
                              local.data <- getMultinomialModel(sequence[window.lowest:window.highest])
                              local.data["GC"] <- local.data["c"] + local.data["g"]
                              local.data["pos"] <- window.lowest
                              result <- rbind(result,local.data)
                              
                              # slides the window updating its parameters by adding offset
                              window.lowest <- window.lowest + offset
                              window.highest <- window.highest + offset
                            }
                            
                            return(result)
                          }

#### Plus strand
  # Change classes
    plus$Color= as.character(plus$Color)
    plus$Start=as.numeric(plus$Start)
    plus$Stop=as.numeric(plus$Stop)
  # hack things around for whatever order you want - here oriented to COI.
    plusOrd = as.character(plus$Name)
    plusFac = plus$Name
    plusFac = factor(plusFac, levels = c(plusOrd))
  # Make all of this into a matrix for the plotting
    plusY = cbind(plus$Start, plus$Stop)
    row.names(plusY) = plus$Name
    plusZ = as.matrix(plusY)

# plot plus strand
  circos.par(cell.padding=c(0,0,0,0),start.degree = 90, gap.degree=0, track.height = 0.1)
  circos.initialize(factors = plusFac, xlim = plusZ )
  circos.track(ylim=c(0,1),track.height = 0.1, bg.border=c(plus$Color), bg.col=c(plus$Color))
  circos.info()
  circos.clear() # and clear it

#### Minus strand
    minus$Color= as.character(minus$Color)
    minus$Start=as.numeric(minus$Start)
    minus$Stop=as.numeric(minus$Stop)
    minusOrd = as.character(minus$Name)
    minusFac = minus$Name
    minusFac = factor(minusFac, levels = c(minusOrd))
  
    minusY = cbind(minus$Start, minus$Stop)
    row.names(minusY) = minus$Name
    minusZ = as.matrix(minusY)

  par(new = TRUE) # as the circlize doc says this is magic
  circos.par(cell.padding=c(0,0,0,0),gap.degree=0,start.degree = 90,"canvas.xlim" = c(-1.1, 1.1), "canvas.ylim" = c(-1.1, 1.1))
  circos.initialize(factors = minusFac, xlim = minusZ)
  circos.track(ylim = c(0, 1), track.height=0.1, bg.border=c(minus$Color), bg.col=c(minus$Color))
  circos.clear() #clear again


#### Coverage

  coverage = log(coverage)-5
  
# need some other things for plotting region and circos lines:
    bases = c(1:length(coverage))
    coverageFac = sample(letters[1], length(coverage), replace = TRUE)

  par(new = TRUE) 
  circos.par(cell.padding=c(0,0.2,0,0.2), gap.degree=0,start.degree = 90, "canvas.xlim" = c(-1.20, 1.20), "canvas.ylim" = c(-1.20, 1.20))
  circos.initialize(factors = coverageFac, xlim = c(1, length(coverage)))
  circos.trackPlotRegion(factors = coverageFac, ylim = c(0, 7), track.height = 0.1, bg.border="white",bg.col="white") 
  circos.lines(bases, coverage, sector.index = "a", type="h", col="grey") # you may get 'notes' here but everything looks goo don't worry
  circos.clear() # clear again


#### GC Content ###
  gcCont = getSequenceLocalData(seq,10, 100)
# again for circos.lines() we need:
  seqFac = sample(letters[1], nrow(gcCont), replace = TRUE)
  position = c(1:nrow(gcCont))

# Plot 
par(new = TRUE)
circos.par(cell.padding=c(0,0.2,0,0.2), gap.degree=0,start.degree = 90, "canvas.xlim" = c(-1.4, 1.4), "canvas.ylim" = c(-1.4, 1.4))
circos.initialize(factors = seqFac, xlim = c(1, nrow(gcCont)))
circos.trackPlotRegion(factors = seqFac, ylim = c(0, 1), track.height = 0.1, bg.border="white",bg.col="white") 
circos.lines(position, gcCont[,5],  sector.index = "a", type="l", area=T, border="white", col="grey")

} # close function


mitoplotR(seq,plus,minus,coverage)
