### Precipitation Interpolation

##Read the  data and preliminar analysis

#Just read data

d <- read.csv("C:/Users/Lumila/Desktop/Hijmans workshop/precipitation.csv")
head(d)

#Compute annual precipitation plot

d$prec <- rowSums(d[, c(6:17)])
plot(sort(d$prec), ylab='Annual precipitation (mm)', las=1, xlab='Stations')

# Agregar mapa de California con los counties

library(sp)
dsp <- SpatialPoints(d[,4:3], proj4string=CRS("+proj=longlat +datum=NAD83"))
dsp <- SpatialPointsDataFrame(dsp, d)
CA <- readRDS("C:/Users/Lumila/Desktop/Hijmans workshop/counties.rds")

# define groups for mapping

cuts <- c(0,200,300,500,1000,3000)

# set up a palette of interpolated colors

blues <- colorRampPalette(c('yellow', 'orange', 'blue', 'dark blue'))
pols <- list("sp.polygons", CA, fill = "lightgray")
spplot(dsp, 'prec', cuts=cuts, col.regions=blues(5), sp.layout=pols, pch=20, cex=2)

#Transformar lat/long en coordenadas flat usando coordinate reference system for California (“Teale Albers”)

TA <- CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0")
library(rgdal)
dta <- spTransform(dsp, TA)
cata <- spTransform(CA, TA)

##Interpolar los valores de precipitacion, usando la media de todas las observaciones

RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}
null <- RMSE(mean(dsp$prec), dsp$prec)
null

#Uso de poligonos de proximidad para interpolar variables categoricas (tambien se podria usar nearest neighbour)

library(dismo)
v <- voronoi(dta)
plot(v)

#Como se ve raro, ahora limitarlo a California

ca <- aggregate(cata)
vca <- intersect(v, ca)
spplot(vca, 'prec', col.regions=rev(get_col_regions()))

#Ahora se puede rasterizar

r <- raster(cata, res=10000)
vr <- rasterize(vca, r, 'prec')
plot(vr)

##Neighbour Interpolation (Multiple neighbours)

#Matriz de distancia de todos los puntos entre si
# control points
cp <- rasterToPoints(vr)
# distance matrix
d <- pointDistance(cp[, 1:2], dta, lonlat=FALSE)
nrow(dta)

nrow(cp)

# not symmetric!
dim(d)

d[1:5,1:5]

#Encontrar 5 vecinos mas proximos a cada punto

nn <- 5
ngb <- t(apply(d, 1, function(x) order(x)[1:nn]))

#Checkear si tiene sentido

plot(cata)
points(cp[1, 1:2, drop=FALSE], col='blue', pch='x', cex=2)
points(dta[ngb[1,], ], col='red', pch=20)
points(cp[nrow(cp), 1:2, drop=FALSE], col='blue', pch='x', cex=2)
points(dta[ngb[nrow(cp),], ], col='red', pch=20)

#Mapa con pares 

pairs <- cbind(rep(1:nrow(ngb), nn), as.vector(ngb))
values <- dta$prec[pairs[,2]]
pn <- tapply(values, pairs[,1], mean)

#Asignar a raster

nnr <- r
nnr[!is.na(vr)] <- pn
plot(nnr)

