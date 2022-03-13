
data.j <- read.csv("dataset.csv",sep=",",header=T)
dim(data.j)
names(data.j)

library(maps)
library(maptools)
library(spatstat)

## Using the maps package in R which also contains a shapefile for the US
us.map <- map("usa",plot=F)
## These are the vectors with the longitudes and latitudes for the points in the boundary
us.bdries.x <- us.map$x
range(us.bdries.x,na.rm=T)
#[1] -124.68134  -67.00742
us.bdries.y <- us.map$y
range(us.bdries.y,na.rm=T)
#[1] 25.12993 49.38323

## The polygon for the US is made of different disconnected components, each with its own
## name. Here we are figuring out the number of disconnected components.
US.names <- us.map$names
length(US.names)
#[1] 10

## Names of the different components of the US polygon
US.names


## This is to determine where the lon and lat coordinates for the main part of the US are
length(which(is.na(us.bdries.x)))
# [1] 9
length(which(is.na(us.bdries.y)))
#[1] 9
which(is.na(us.bdries.x))
#[1] 6887 6924 6955 6972 6983 7152 7170 7188 7208
which(is.na(us.bdries.y))
#[1] 6887 6924 6955 6972 6983 7152 7170 7188 7208

## Plotting just the main part of the US without any island
plot(us.bdries.x[1:6887],us.bdries.y[1:6887],type="l",col="black",xlab="Longitude",ylab="Latitude")

## Creating the matrix with the longitude and latitude of the boundaries
bdries.us <- matrix(cbind(us.bdries.x[1:6886],us.bdries.y[1:6886]),nrow=6886,ncol=2)
## Creating the window of the spatial point pattern.
win.usa <- owin(poly=bdries.us)


names(data.j)
lon.shootings <- data.j$lon
lat.shootings <- data.j$lat
## Creating the matrix with the actual coordinates of the "events" (the shootings)
X <- matrix(cbind(lon.shootings,lat.shootings),nrow=length(lon.shootings),ncol=2)

## Creating the ppp oject: some shootings are outside the boundaries (8). You
## might look at them one by one.
ppp.j <- as.ppp(X,win.usa)

