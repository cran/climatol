pkgname <- "climatol"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('climatol')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("dahstat")
### * dahstat

flush(stderr()); flush(stdout())

### Name: dahstat
### Title: Statistical summaries of the homogenized data
### Aliases: dahstat
### Keywords: datagen

### ** Examples

  #This function only works with files. In order to run the example,
  #you should run the examples of the homogen function first.
  #Once the results from homogen have been generated, you can run:
  ## Not run: dahstat("Tmin", 1956, 2005, 1971, 2000)
  #See the resulting monthly mean data for the period 1971-2000
  #in the generated file 'Tmin_1971-2000.med'.



cleanEx()
nameEx("dd2m")
### * dd2m

flush(stderr()); flush(stdout())

### Name: dd2m
### Title: Compute monthly data from daily series
### Aliases: dd2m
### Keywords: datagen

### ** Examples

  #This function only works with files. In order to run the example,
  #you should run the examples of the homogen function first.
  #Once the results from homogen have been generated, you can run:
  ## Not run: dd2m("Tmax", 2007, 2010, ini='2007-01-01')
  #See the resulting monthly mean data in the generated file
  #'Tmin-m_1971-2000.dah'.
  #To take advantage of all the missing daily data estimations, run:
  ## Not run: dd2m("Tmax", 2007, 2010, ini='2007-01-01',nmin=0)
  #(The last quarter of 2010 is still NA (not available), because
  #there were no original data it that period).



cleanEx()
nameEx("diagwl")
### * diagwl

flush(stderr()); flush(stdout())

### Name: diagwl
### Title: Walter & Lieth climatic diagram
### Aliases: diagwl
### Keywords: hplot

### ** Examples

  data(datcli)
  diagwl(datcli,est="Example station",alt=100,per="1961-90",mlab="en")



cleanEx()
nameEx("homogen")
### * homogen

flush(stderr()); flush(stdout())

### Name: homogen
### Title: Automatic homogenization of climatological series
### Aliases: homogen
### Keywords: datagen ts manip

### ** Examples

  #As this function only works with files, you must uncompress the
  #example data in your working directory first. (The example files
  #can be obtained from http://webs.ono.com/climatol/climatol-dat.zip).
  #Afterwards, we can run the example:
  homogen("Tmin", 1956, 2005)
  #See the resulting four output files "Tmin_1956-2005.dah",
  #"Tmin_1956-2005.esh", "Tmin_1956-2005.txt" and "Tmin_1956-2005.pdf".
  #
  #Another example with daily data, but only for filling missing data:
  homogen('Tmax', 2007, 2010, nm=0, tVt=0, snhtt=0, ini='2007-01-01')
  #The results will have been saved to four Tmax-d_2007-2010.* files.



cleanEx()
nameEx("leerdm")
### * leerdm

flush(stderr()); flush(stdout())

### Name: leerdm
### Title: Read the climatological series to homogenize
### Aliases: leerdm
### Keywords: internal IO

### ** Examples

  #As this function only works with files, you must uncompress the
  #example data in your working directory first. (The example files
  #can be obtained from http://webs.ono.com/climatol/climatol-dat.zip).
  #After uncompressing them, we can run:
  leerdm("Tmin", 1956, 2005)
  #data are now in the memory space, in objects dat, est.c, ne and nd.



cleanEx()
nameEx("outrename")
### * outrename

flush(stderr()); flush(stdout())

### Name: outrename
### Title: Rename homogen's output files
### Aliases: outrename
### Keywords: utilities

### ** Examples

  #After having run the first example of the homogen function, you can try:
  ## Not run: outrename("Tmin", 1956, 2005, "old")
  #The previous four output files will have '-old' appended to their base name.



cleanEx()
nameEx("rosavent")
### * rosavent

flush(stderr()); flush(stdout())

### Name: rosavent
### Title: Wind-rose plot
### Aliases: rosavent
### Keywords: hplot

### ** Examples

  data("windfr")
  rosavent(windfr, 4, 4, ang=-3*pi/16, main="Annual windrose")



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
