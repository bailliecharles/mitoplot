#########################################################
#### Plotting circular genomes in R using {circlize} ####
#########################################################

# Get the tools
install.packages("circlize")
library(circlize)
install.packages("seqinr")
library(seqinr)

##### Mitogenome feature tracks #######
### Plus ###
# data simply needs to be in a data frame - I addedd a column for colors with the HEX codes (colorbrewer Dark2), but you could cbind a vector
# after reading in the data. Cols are Color, feature start position (bp), feature stop position (bp), name of feature,
# like so:
# Name	Start	Stop	Color
# cox1	1	    1533	#7570b3
# Gap1	1534	1581	white
# cox2	1582	2265	#7570b3
# ...   ...   ...   ...   and so on.
df = read.csv(file.choose(), header=T)

# check and sort out classes
sapply(df,class)
df$Color= as.character(df$Color)
df$Start=as.numeric(df$Start)
df$Stop=as.numeric(df$Stop)
# hack things around for whatever order you want - here oriented to COI.
ord = as.character(df$Name)
myFac = df$Name
myFac = factor(myFac, levels = c(ord))
myFac
# Make all of this into a matrix for the plotting
y = cbind(df$Start,df$Stop)
row.names(y) = df$Name
z = as.matrix(y)

# plot
circos.par(cell.padding=c(0,0,0,0),start.degree = 90, gap.degree=0, track.height = 0.1)
circos.initialize(factors = myFac, xlim = z )
circos.track(ylim=c(0,1),track.height = 0.1, bg.border=c(df$Color), bg.col=c(df$Color))
circos.info()
# make sure you clear the track before moving on!
circos.clear()

### Minus ###
# do the same for the minus strand
df2 = read.csv(file.choose(), header=T)

sapply(df2,class)
df2$Color= as.character(df2$Color)
df2$Start=as.numeric(df2$Start)
df2$Stop=as.numeric(df2$Stop)
ord2 = as.character(df2$Name)
myFac2 = df2$Name
myFac2 = factor(myFac2, levels = c(ord2))


y2 = cbind(df2$Start,df2$Stop)
row.names(y2) = df2$Name
z2 = as.matrix(y2)
nrow(z2)
nlevels(df2$Name)

par(new = TRUE) # <- as the circlize doc says this is magic

circos.par(cell.padding=c(0,0,0,0),gap.degree=0,start.degree = 90,"canvas.xlim" = c(-1.1, 1.1), "canvas.ylim" = c(-1.1, 1.1))
circos.initialize(factors = myFac2, xlim = z2)
circos.track(ylim = c(0, 1), track.height=0.1, bg.border=c(df2$Color), bg.col=c(df2$Color))
# clear again
circos.clear()


### Coverage ###
# I had to log transform and standardise because coverage was through the roof
# Data
basecov = read.table(file.choose(), header=F)
base = log(basecov[,3])
base2 = base-5

# need some other things for plotting region and circos lines:
positions = c(1:nrow(basecov))
factors = sample(letters[1], nrow(basecov), replace = TRUE)

### Plot
par(new = TRUE) 
circos.par(cell.padding=c(0,0.2,0,0.2), gap.degree=0,start.degree = 90, "canvas.xlim" = c(-1.20, 1.20), "canvas.ylim" = c(-1.20, 1.20))
circos.initialize(factors = factors, xlim = c(1, nrow(basecov)))
circos.trackPlotRegion(factors = factors, ylim = c(0, 7), track.height = 0.1, bg.border="white",bg.col="white") 
circos.lines(positions, base2, sector.index = "a", type="h", col="grey") # you may get 'notes' here but everything looks goo don't worry
# clear again
circos.clear()
          

### GC Content ###
# Read in sequence data
cer = read.fasta(file.choose())
seq = cer[[1]]

### I then borrowed the gc sliding window funcs from https://www.cs.us.es/~fran/students/julian/introduction/introduction.html
# where there are some other cool R things. Thanks, Julian and Francisco.

getMultinomialModel <- function(sequence) {
  frequencies <- getFrequency(sequence)
  sequence.length <- length(sequence)
  model <- frequencies / sequence.length
  names(model) <- names(frequencies)
  return(model)
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

gc_cnt = getSequenceLocalData(seq,10, 100)

# again for circos.lines() we need:
factors2 = sample(letters[1], nrow(gc_cnt), replace = TRUE)
positions2 = c(1:nrow(gc_cnt))

# Plot 
par(new = TRUE)
circos.par(cell.padding=c(0,0.2,0,0.2), gap.degree=0,start.degree = 90, "canvas.xlim" = c(-1.4, 1.4), "canvas.ylim" = c(-1.4, 1.4))
circos.initialize(factors = factors2, xlim = c(1, nrow(gc_cnt)))
circos.trackPlotRegion(factors = factors2, ylim = c(0, 1), track.height = 0.1, bg.border="white",bg.col="white") 
circos.lines(positions2, gc_cnt[,5],  sector.index = "a", type="l", area=T, border="white", col="grey")


