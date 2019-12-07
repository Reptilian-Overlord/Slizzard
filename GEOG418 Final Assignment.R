##Sidney Forstrom V00852415
##GEOG418 Final Assignment 
##SubSample:150

##Libraries


install.packages("sf")
install.packages("plyr")
install.packages("dplyr")
install.packages("spdep")
install.packages("GISTools")
install.packages("raster")
install.packages("maptools")
install.packages("rgdal")
install.packages("spatstat")
install.packages("sp")
install.packages("spatstat")
install.packages("tmap")
install.packages("gstat")
install.packages("spgwr")
install.packages("rgeos")


library(sf)
library(plyr)
library(dplyr)
library(spdep)
library(GISTools)
library(raster)
library(maptools)
library(rgdal)
library(spatstat)
library(sp)
library(spatstat)
library(tmap)
library(gstat)
library(spgwr)
library(rgeos)

##Probably going to need to add more packages to get descriptive stats for NND


#Set working directory
dir <- "C:/Program Files/R/GEOG418/Final_Proj/master"
setwd(dir)

#Reading in particulate matter dataset
pm25 <- read.csv("PM25.csv") #Read in PM2.5 data
#Select only columns 1 and 2
pm25 <- pm25[,1:2]
#Change the column names 
colnames(pm25) <- c("POSTALCODE", "PM25")
pm25 <- na.omit(pm25)

#Reading in postal code shapefile
postalcodes <- shapefile("./BC_Postal_Codes.shp") #Read in related postal code data

#Reading in dissemination tract and income data
income <- read.csv("Income.csv") #Read in census income data  
colnames(income) <- c("DAUID", "Income") #Select only ID and Income columns
census.tracts <- shapefile("./BC_DA.shp") #Read in dissemination tract shapefile
income.tracts <- merge(census.tracts,income, by = "DAUID") #Merge income and dissemination data
nrow(income.tracts) #Determine the number of columns in the dataframe
income.tracts <- income.tracts[!is.na(income.tracts$Income),]


#Create choropleth map of income
med.income <- income.tracts$Income
shades <- auto.shading(med.income, n=6, cols = brewer.pal(6, 'Oranges'))
choropleth(income.tracts, med.income, shades) #map the data with associated colours
choro.legend(3864000, 1965000, shades) #add a legend (you might need to change the location)

#Select postal codes that fall within dissemination tracts)
postalcodes <- intersect(postalcodes,income.tracts)
plot(postalcodes) #See what the data looks like spatially
head(postalcodes) #See what the data looks like in tabular form

##income.tracts <- spTransform(income.tracts, CRS("+init=epsg:4326"))

#Join PM2.5 data with postal code data
pm25.spatial <- merge(postalcodes,pm25,by = "POSTALCODE")

#Aggregate the PM2.5 values in each DA in order to have a single value per DA. Here we aggregate based on the mean.
pm25.aggregate <- aggregate((as.numeric(pm25.spatial$PM25)/10)~pm25.spatial$DAUID,FUN=max)

#Re-join aggregated data to the income.tracts layer.
colnames(pm25.aggregate) <- c("DAUID", "PM25AGG") #Select only ID and Income columns
income.pm25 <- merge(income.tracts,pm25.aggregate, by = "DAUID") #Merge income and dissemination data

#Re-join aggregated data to the pm25.spatial points layer.
pm25.points.aggregate <- merge(pm25.spatial, pm25.aggregate, by = "DAUID")

#Create a subsample of the datapoints provided in the PM2.5 dataset using the sample n provided on CourseSpaces

sampleSize=150
spSample <- pm25.points.aggregate[sample(1:length(pm25.points.aggregate),sampleSize),]

#Create a grid called grd to use in your interpolation
# Create an empty grid where n is the total number of cells
grd <- as.data.frame(spsample(spSample, "regular", n=5000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object
proj4string(grd) <- proj4string(spSample)

crs(grd)



################################################################################################################



##Moran's I



#####Weighted Matrix
crd.nb <- poly2nb(income.tracts)
crd.net <- nb2lines(crd.nb,coords=coordinates(income.tracts))##Converts neighbour matrix into plot


tm_shape(income.tracts) + tm_borders(col='lightgrey') + 
  tm_shape(crd.net) + tm_lines(col='red')

##This graphic shows us the neighbours of the polygons. We are using queen weight in red line (queen = default)

crd.nb2 <- poly2nb(income.tracts, queen = FALSE)#we are using the rook weight for the yellow line one
crd.net2 <- nb2lines(crd.nb2,coords=coordinates(income.tracts))


tm_shape(income.tracts) + tm_borders(col='lightgrey') + 
  tm_shape(crd.net) + tm_lines(col='blue', lwd = 2) + #blue lines are corner connections with polygons
  tm_shape(crd.net2) + tm_lines(col='yellow', lwd = 2)

########################

crd.lw <- nb2listw(crd.nb, zero.policy = TRUE, style = "W")

print.listw(crd.lw, zero.policy = TRUE) #These tell us our list matrix

#####Moran's I
mi <- moran.test(income.tracts$Income, crd.lw, zero.policy = TRUE) ##This is our whole moran's i test.

mi ##shows our results


moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values)) ##shows us what the full range of mi test can be for this data
  
}

moran.range(crd.lw)


mI <- mi$estimate[[1]] ##moran's i value

eI <- mi$estimate[[2]] ##expected i value

var <- mi$estimate[[3]]

z <- (mI - eI) / (sqrt(var))


########################  

lisa.test <- localmoran(income.tracts$Income, crd.lw) ##Very similar equation to moran's i

income.tracts$Ii <- lisa.test[,1] ##Lisa i value
income.tracts$E.Ii<- lisa.test[,2] ##Expected lisa i value
income.tracts$Var.Ii<- lisa.test[,3] ##Lisa variance value
income.tracts$Z.Ii<- lisa.test[,4] ##lisa provides a z value, we don't have to calculate for each polygon
income.tracts$P<- lisa.test[,5]  ##Will give us a bunch of descreptive statistics for each polygon



map_LISA <- tm_shape(income.tracts) + ##make a map out of our crd data, moran's i, and lisa
  tm_polygons(col = "Ii", 
              title = "Local Moran's I", 
              style = "fisher", 
              palette = "viridis", n = 6) 


map_LISA


####Kriging

f.0 <- as.formula(PM25AGG ~ 1) 


var.smpl <- variogram(f.0, spSample, cloud = FALSE) #, cutoff=1000000, width=89900)
dat.fit  <- fit.variogram(var.smpl, fit.ranges = TRUE, fit.sills = TRUE,
                          vgm(model="Exp"))
plot(var.smpl, dat.fit)

# Define the trend model
f.0 <- as.formula(PM25AGG~ 1) 

# Perform the krige interpolation (note the use of the variogram model
# created in the earlier step)
dat.krg <- krige( f.0, spSample, grd, dat.fit) ##This is what interpolates

# Convert kriged surface to a raster object for clipping
r <- raster(dat.krg)
r.m <- mask(r, income.tracts)

# Plot the map
tm_shape(r.m) + 
  tm_raster(n=10, palette="-RdBu",  
            title="Predicted Ozone \n(in ppm)") +
  tm_shape(spSample) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)

r   <- raster(dat.krg, layer="var1.var")
r.m <- mask(r, income.tracts)

tm_shape(r.m) + 
  tm_raster(n=7, palette ="-RdBu",
            title="Variance map \n(in squared ppm)") +tm_shape(spSample) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)

r   <- sqrt(raster(dat.krg, layer="var1.var")) * 1.96
r.m <- mask(r, income.tracts)

tm_shape(r.m) + 
  tm_raster(n=7, palette ="-RdBu",
            title="95% CI map \n(in ppm)") +tm_shape(spSample) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)





################################################################################################################






#These steps will help you combine the outputs from your spatial interpolation with your income data.

#If you have too many cells, you can reduce the number by aggregating values
#step.1 <- aggregate(r, fact= 2, fun=mean)
#plot(step.1)

#Convert the raster dataset to points

#step.2 <-  rasterToPoints(step.1,fun=NULL, spatial=FALSE, crs=spSample)
#step.2 <- as.data.frame(step.2) #convert the point dataset to a spatial dataframe
#Coords <- step.2[,c("x", "y")]  #assign coordinates to a new object
#crs <- crs(census.tracts) #utilize an existing projection
#step.3 <- SpatialPointsDataFrame(coords = Coords, data = step.2, proj4string = crs) #create a spatial points dataframe
#step.4 <- aggregate(x=step.3,by=income.tracts, FUN=mean) #aggregate points into census tracts
#step.5 <- intersect(step.4,income.tracts)  #get the intersection of step.4 with the income.tracts dataset (this will take a while) 

step.5a <- extract(r, income.tracts, fun = mean, sp = TRUE)

#plot(step.5a)

#You are now ready to perform a regression




################################################################################################################



######Linear Regression##########
#Let's say your dataset with both PM2.5 and Income are stored in a dataset called pm.income.poly.

pm.income.poly <- step.5a
#View(step.5@data)

#Plot income and PM2.5 from the pm.income.poly dataset you created
plot(pm.income.poly$Income~pm.income.poly$layer)

#Notice that there are a lot of 0's in this dataset. If you decide to remove them, use the following line:
pm.income.poly <- pm.income.poly[!is.na(pm.income.poly$layer),] ##Remove NA values

#Now plot the data again
plot(pm.income.poly$Income~pm.income.poly$layer)

#Perform a linear regression on the two variables. You should decide which one is dependent.
lm.model <- lm(pm.income.poly$Income~pm.income.poly$layer)
#Add the regression model to the plot you created
abline(lm.model)
#Get the summary of the results
summary(lm.model)

#You want to determine if the model residuals are spatially clustered. 
#First obtain the residuals from the model
model.resids <- as.data.frame(residuals.lm(lm.model))
#Then add the residuals to your spatialpolygon dataframe
pm.income.poly$residuals <- residuals.lm(lm.model)
#Observe the result to make sure it looks correct
head(pm.income.poly)

#Now, create choropleth map of residuals
resids <- pm.income.poly$residuals
shades <- auto.shading(resids, n=6, cols = brewer.pal(6, 'Greens'))
choropleth(income.tracts, resids, shades) #map the data with associated colours
choro.legend("left", "bottom", shades) #add a legend (you might need to change the location)



################################################################################################################



##MoransI2



#####Weighted Matrix
crd.nb.reg <- poly2nb(pm.income.poly)
crd.net.reg <- nb2lines(crd.nb.reg,coords=coordinates(pm.income.poly))##Converts neighbour matrix into plot


tm_shape(pm.income.poly) + tm_borders(col='lightgrey') + 
  tm_shape(crd.net.reg) + tm_lines(col='red')

##This graphic shows us the neighbours of the polygons. We are using queen weight in red line (queen = default)

crd.nb2.reg <- poly2nb(pm.income.poly, queen = FALSE)#we are using the rook weight for the yellow line one
crd.net2.reg <- nb2lines(crd.nb2.reg,coords=coordinates(pm.income.poly))

plot(pm.income.poly)
plot(step.5a)


tm_shape(pm.income.poly) + tm_borders(col='lightgrey') + 
  tm_shape(crd.net.reg) + tm_lines(col='blue', lwd = 2) + #blue lines are corner connections with polygons
  tm_shape(crd.net2.reg) + tm_lines(col='yellow', lwd = 2)

########################

crd.lw.reg <- nb2listw(crd.nb.reg, zero.policy = TRUE, style = "W")

print.listw(crd.lw.reg, zero.policy = TRUE) #These tell us our list matrix

#####Moran's I
mi.reg <- moran.test(pm.income.poly$Income, crd.lw.reg, zero.policy = TRUE) ##This is our whole moran's i test.

mi.reg ##shows our results


moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values)) ##shows us what the full range of mi test can be for this data
  
}

moran.range(crd.lw.reg)


mI.reg <- mi$estimate[[1]] ##moran's i value

eI.reg <- mi$estimate[[2]] ##expected i value

var.reg <- mi$estimate[[3]]

z.reg <- (mI - eI) / (sqrt(var))


########################  

lisa.test <- localmoran(pm.income.poly$Income, crd.lw.reg) ##Very similar equation to moran's i

pm.income.poly$Ii <- lisa.test[,1] ##Lisa i value
pm.income.poly$E.Ii<- lisa.test[,2] ##Expected lisa i value
pm.income.poly$Var.Ii<- lisa.test[,3] ##Lisa variance value
pm.income.poly$Z.Ii<- lisa.test[,4] ##lisa provides a z value, we don't have to calculate for each polygon
pm.income.poly$P<- lisa.test[,5]  ##Will give us a bunch of descreptive statistics for each polygon



map_LISA <- tm_shape(pm.income.poly) + ##make a map out of our crd data, moran's i, and lisa
  tm_polygons(col = "Ii", 
              title = "Local Moran's I", 
              style = "fisher", 
              palette = "viridis", n = 6) 


map_LISA



################################################################################################################


####Geographically Weighted Regression

#Let's say you are continuing with your data from the regression analysis. 

#The first thing you need to do is to add the polygon coordinates to the spatialpolygondataframe.
#You can obtain the coordinates using the "coordinates" function from the sp library
pm.income.poly.coords <- sp::coordinates(pm.income.poly)
#Observe the result
head(pm.income.poly.coords)
#Now add the coordinates back to the spatialpolygondataframe
pm.income.poly$X <- pm.income.poly.coords[,1]
pm.income.poly$Y <- pm.income.poly.coords[,2]
head(pm.income.poly)

##pm.income.poly <- pm.income.poly[!is.na(pm.income.poly$PM25),] ##remove NA values
View(pm.income.poly@data) ##View data of this spatial data frame


###Determine the bandwidth for GWR: this will take a while
GWRbandwidth <- gwr.sel(pm.income.poly$Income~pm.income.poly$layer, 
                        data=pm.income.poly, coords=cbind(pm.income.poly$X,pm.income.poly$Y),adapt=T) 

###Perform GWR on the two variables with the bandwidth determined above
###This will take a looooooong while
gwr.model = gwr(pm.income.poly$Income~pm.income.poly$layer, ###anywhere there is PM25, change to layer
                data=pm.income.poly, coords=cbind(pm.income.poly$X,pm.income.poly$Y), 
                adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 

#Print the results of the model
gwr.model

#Look at the results in detail
results<-as.data.frame(gwr.model$SDF)
head(results)

#Now for the magic. Let's add our local r-square values to the map
pm.income.poly$localr <- results$localR2

#Create choropleth map of r-square values
local.r.square <- pm.income.poly$localr
shades <- auto.shading(local.r.square, n=6, cols = brewer.pal(6, 'Oranges'))
choropleth(income.tracts, local.r.square, shades) #map the data with associated colours
choro.legend("bottomleft", "bottomleft", shades) #add a legend (you might need to change the location)

#Time for more magic. Let's map the coefficients
pm.income.poly$coeff <- results$pm.income.poly.PM25

#Create choropleth map of the coefficients
local.coefficient <- pm.income.poly$coeff
shades <- auto.shading(local.coefficient, n=6, cols = brewer.pal(6, 'Oranges'))
choropleth(income.tracts, local.coefficient, shades) #map the data with associated colours
choro.legend("bottomleft", "bottomleft", shades) #add a legend (you might need to change the location)



################################################################################################################


##Nearest Neighbour Distance
###NEAREST NEIGHBOUR
nearestNeighbour <- nndist(results)

##Convert the nearestNeighbor object into a dataframe.
nearestNeighbour=as.data.frame(as.numeric(nearestNeighbour))
##Change the column name to "Distance"
colnames(nearestNeighbour) = "Distance"

##gArea(results)


##Calculate the nearest neighbor statistic to test for a random spatial distribution.
#mean nearest neighbour

N <- 150

nnd = sum(nearestNeighbour) / N

#mean nearest neighbour for random spatial distribution


studyArea <- gArea(income.tracts)

pointDensity <- N/studyArea

r.nnd = 1 / (2*sqrt(pointDensity))

d.nnd = 1.07453 / sqrt(pointDensity)

R = nnd / r.nnd

SE.NND <- 0.26136 / sqrt(N*pointDensity)

znnd = (nnd - r.nnd) / SE.NND
