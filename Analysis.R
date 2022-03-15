######load libraries########
library(raster)
library(fields)
library(truncnorm)
library(gridExtra)
library(fBasics)
library(gbm)
library(ggplot2) 
library(grDevices)
library(rgeos)
library(extraDistr)

######define functions########
# Link function g(.)
link <- function(nu){ptruncnorm(nu, a=0, b=Inf, mean = 0, sd = 1)}
# Inverse of that link function g^-1()
ilink <- function(nu){qtruncnorm(nu, a=0, b=Inf, mean = 0, sd = 1)}
# First-order neighborhood matrix from a RasterLayer object (from Hefley et al. 2017)
neighborhood <- function(raster){
  nn <- matrix(,length(raster[]),4)
  for(i in 1:dim(nn)[1]){
    loc <- adjacent(raster,i)[,2]
    ln <- loc[which((loc+1)==i)]
    rn <- loc[which((loc-1)==i)]
    bn <- loc[which((loc-dim(raster)[2])==i)]
    tn <- loc[which((loc+dim(raster)[2])==i)]
    nn[i,1] <- if(length(ln)>0){ln}else{0}
    nn[i,2] <- if(length(rn)>0){rn}else{0}
    nn[i,3] <- if(length(bn)>0){bn}else{0}
    nn[i,4] <- if(length(tn)>0){tn}else{0}
  }
  nn
}
# Propagator matrix for plain diffusion PDE (from Hefley et al. 2017)
propagator.plain <- function(NN,mu,lambda,dx,dy,dt){
  H <- matrix(0,dim(NN)[1],dim(NN)[1])
  for(i in 1:dim(H)[1]){
    if(length(which(NN[i,]>0))==4){
      H[i,i] <- 1-2*mu[i]*(dt/dx^2 + dt/dy^2) + dt*lambda[i]
      H[i,NN[i,1]] <- dt/dx^2*mu[i]
      H[i,NN[i,2]] <- dt/dx^2*mu[i]
      H[i,NN[i,3]] <- dt/dy^2*mu[i]
      H[i,NN[i,4]] <- dt/dy^2*mu[i]}
  }
  H
}

# RasterStack of c(s,t) (from Hefley et al. 2017)
calc.c <- function(H,c0,t.steps,t.keep){
  c.all <- c0
  c.all[] <- H%*%c0[]
  c.all <- stack(mget(rep("c.all",t.steps)))
  for(t in 2:t.steps){
    c.all[[t]][] <- H%*%c.all[[t-1]][]
  }
  c.all[[t.keep]]
}

#time since introduction (Zero truncate beta-binomial)
t.sample <- function(t) {
  success <- FALSE
  while (!success) {
    t <-rbbinom(1, 10, alpha = 1, beta = 0.5)
    success <- t > 0
  }
  return(t)
}

#load in landscape and cell reference rasters
development <- raster("develop.tif")
forest <- raster("forest.tif")
river <- raster("river.tif")
cell <- raster("cell.tif")

#load in positives
all.dat <- read.csv("AllDeerExamined.csv")
all.dat.sub.18 <- all.dat[(all.dat$Year==2018 &  all.dat$Y_COORD > 214002.6 & all.dat$Y_COORD < 391030 & all.dat$X_COORD > 483000 & all.dat$X_COORD < 660027.4),]
all.dat.sub.19 <- all.dat[(all.dat$Year==2019 & all.dat$Y_COORD > 214002.6 & all.dat$Y_COORD < 391030 & all.dat$X_COORD > 483000 & all.dat$X_COORD < 660027.4),]

#extract xy coords for 2017 positives
pos.clust.17 <- all.dat[(all.dat$Year==2017 & all.dat$CWDInd==1 &  all.dat$Y_COORD > 214002.6 & all.dat$Y_COORD < 391030 & all.dat$X_COORD > 483000 & all.dat$X_COORD < 660027.4),]
pos.clustxy <- cbind(pos.clust.17$X_COORD, pos.clust.17$Y_COORD)

#extract xy coords for 2015,2016, 2017 positives
pos.clust.bern <- all.dat[(all.dat$CWDInd==1 &  all.dat$Y_COORD > 214002.6 & all.dat$Y_COORD < 391030 & all.dat$X_COORD > 483000 & all.dat$X_COORD < 660027.4),]
pos.clust.bern <- pos.clust.bern[pos.clust.bern$Year %in% c(2015,2016,2017), ]
pos.clust.bern.xy <- cbind(pos.clust.bern$X_COORD, pos.clust.bern$Y_COORD)



#create shapefile of all positives
pos.all.df <- all.dat[(all.dat$CWDInd==1),]
coordinates(pos.all.df)=~X_COORD+Y_COORD
proj4string(pos.all.df)<- CRS("+proj=omerc +lat_0=45.30916666666666 +lonc=-86 +alpha=337.25556 +k=0.9996 +x_0=2546731.496 +y_0=-4354009.816 +no_uoff +ellps=GRS80 +units=m +no_defs")
pos.all.df.proj <- spTransform(pos.all.df, CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))

pos.clust.pred = subset(pos.all.df.proj, Year %in% c(2015,2016,2017))
pos.clust.fore = subset(pos.all.df.proj, Year %in% c(2018,2019))


#crop cell layer to study area extent
s <- extent(483000,660027.4, 214002.6,391030)
cell <- crop(cell,s)
cell <- setValues(cell,1:ncell(cell))
nsite<-ncell(cell)

#crop landscape layers
development <- crop(development,s)
forest <- crop(forest,s)
river <- crop(river,s)
river[is.na(river)] <- 0

##all rasters to vectors
river.v <- as.vector(river)
river.v[is.na(river.v)] <- 0
forest.v <- as.vector(forest)
development.v <- as.vector(development)
intercept.v <- rep(1,12100)
w <- cbind(river.v, forest.v, development.v)
w <- as.matrix(w)
z <- w
z <- as.matrix(z)

#stack spatial covariates for analysis
spatial.covariates <- stack(cell,river,forest,development)
spatial.covariates$intercept <- 1
spatial.covariates <- subset(spatial.covariates, c(1, 5, 2, 3, 4))

#number of Monte Carlo samples to draw
num.samples <- 1e6
#progress bar
pb <- txtProgressBar(min=1,max=num.samples,style=3)

#matrices for storing data
#2017
Y.pred <- matrix(, ncell(cell), num.samples)
#2018
Y.1 <- matrix(, ncell(cell), num.samples)
#2019
Y.2 <- matrix(, ncell(cell), num.samples)

#run simulations
for(i in 1:num.samples){
  #random vector for diffusion rate according to covariates
  alpha <- c(rnorm(1 ,1.721e+01, 0.098), rnorm(1, 6.254e-01, 0.233),rnorm(1, -1.921e-01, 0.053), rnorm(1, -1.297e-01, 0.12))
  #random Vector for growth rate according to covariates
  gamma <- c(rnorm(1, 6.914e-02,0.04), rnorm(1, 2.189e-01, 0.30), rnorm(1, 2.087e-01, 0.05), rnorm(1, 7.007e-02, 0.08))
  #random vector for demographic effects
  beta <- c(rnorm(1, 8.1e-1, 0.06), rnorm(1, 4.2e-1, 0.03))
  #pick random positive as point of introduction
  d <- SpatialPoints(matrix(pos.clustxy[sample(nrow(pos.clustxy), 1), ], nrow=1, ncol=2))
  #pick random time since 2017 (+2 is to forecast out to 2019 positives)
  t <- t.sample()+2
  #growth rate
  lambda <- raster(, nrows = nsite^0.5, ncols = nsite^0.5, xmn = s[1], xmx = s[2], ymn = s[3],
                   ymx = s[4], crs = CRS)
  lambda[] <- model.matrix(~w) %*% gamma
  #diffusion rate
  mu <- raster(, nrows = nsite^0.5, ncols = nsite^0.5, xmn = s[1], xmx = s[2], ymn = s[3],
               ymx = s[4], crs = CRS)
  mu[] <- exp(model.matrix(~z) %*% alpha)
  # Scaling factor
  us.fact <- 5
  # Diffusion rate for homogenized pde
  mu.bar <- aggregate(mu, fact = us.fact, na.rm = TRUE, fun = function(x, na.rm) {
    1/mean(1/x)})
  # Growth rate for homogenized pde
  lambda.bar <- mu.bar * aggregate(lambda/mu, fact = us.fact, na.rm = TRUE, FUN = mean)
  # First-order neighborhood matrix
  NN <- neighborhood(mu.bar)
  # Propagator matrix
  H <- propagator.plain(NN = NN, mu = mu.bar[], lambda = lambda.bar[], dx=8046.7,dy = 8046.7,dt=1/4)
  #create raster with cells equal to homogenized pde
  crand <- raster(, nrows = nsite^0.5/us.fact, ncols = nsite^0.5/us.fact, xmn = s[1], xmx = s[2], ymn = s[3],  ymx = s[4], crs = CRS)
  crand[] <- 1:ncell(crand)
  #vector location of random introduction
  point <- extract(crand, d)
  #small but nonzero disease intensity number across grid
  theta.all <- rep(1.0e-6, length(crand))
  #initial intensity
  p.star <- ilink(runif(1,1/45,1/15))
  theta.point <- p.star/(exp((1*beta[1]) + (4 * beta[2])))
  theta.all[point] <- theta.point
  ##define initial conditions of homogenized pde
  c0 <- raster(, nrows = nsite^0.5/us.fact, ncols = nsite^0.5/us.fact, xmn = s[1], xmx = s[2], ymn = s[3], ymx = s[4], crs = CRS)
  c0[] <- extract(mu, SpatialPoints(mu.bar))
  c0[] <- values(c0) * theta.all
  #homoegenized pde landscape effects up to 2019
  c.all <- calc.c(H,c0,t*4,1:(t*4))
  u.all <- disaggregate(c.all, us.fact)/mu
  #probability of infection up to 2019
  p <- link((exp((1*beta[1]) + (4 * beta[2]))) * vec(u.all[]))
  p.all <- u.all
  p.all[] <- p
  #extract 2019 probability of infection values
  y.2 <- p.all[[t*4]]
  #extract 2018 probability of infection values
  y.1 <- p.all[[(t-1)*4]]
  #extract 2017 probability of infection values
  y.pred <- p.all[[(t-2)*4]]
  #store data from each iteration
  Y.2[, i] <- values(y.2) 
  Y.1[, i] <- values(y.1) 
  Y.pred[, i] <- values(y.pred) 
  setTxtProgressBar(pb,value=i)
}
close(pb)

#Expected values (Derived quantities) for 2017, 2018, 2019
E.y.17 <- rowMeans(Y.pred)
E.y.18 <- rowMeans(Y.1)
E.y.19 <- rowMeans(Y.2)

#Plot Results
#Define coordinate reference system
coord.ll <- "+proj=omerc +lat_0=45.30916666666666 +lonc=-86 +alpha=337.25556 +k=0.9996 +x_0=2546731.496 +y_0=-4354009.816 +no_uoff +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#2017 Raster 
E.y.rast.17<- raster(, nrows = nsite^0.5, ncols = nsite^0.5, xmn = s[1], xmx = s[2], ymn = s[3], ymx = s[4], crs = coord.ll)
E.y.rast.17[] <- E.y.17

#2018 Raster 
E.y.rast.18<- raster(, nrows = nsite^0.5, ncols = nsite^0.5, xmn = s[1], xmx = s[2], ymn = s[3], ymx = s[4], crs = coord.ll)
E.y.rast.18[] <- E.y.18

#2019 Raster 
E.y.rast.19<- raster(, nrows = nsite^0.5, ncols = nsite^0.5, xmn = s[1], xmx = s[2], ymn = s[3], ymx = s[4], crs = coord.ll)
E.y.rast.19[] <- E.y.19

#Bernoulli Deviance Calculation

#Bayesian Approach
#2018 accuracy
detection.sim.18 <- extract(E.y.rast.18, all.dat.sub.18[,c(10,11)])
bern.out.18 <- sum(-2*dbinom(all.dat.sub.18[,8], 1, detection.sim.18, log=TRUE))
#2019 accuracy
detection.sim.19 <- extract(E.y.rast.19, all.dat.sub.19[,c(10,11)])
bern.out.19 <- sum(-2*dbinom(all.dat.sub.19[,8], 1, detection.sim.19, log=TRUE))
#overall 
bern.tot.real <- bern.out.18 + bern.out.19

#Rule-Based Approach
#Create buffer surrounding known positives
pos.clust.in <- SpatialPointsDataFrame(coords = pos.clust.bern.xy, data = pos.clust.bern,
                               proj4string = CRS(coord.ll))
pos.clust.buff= subset(pos.clust.in, Year %in% c(2015,2016,2017))
circle <- buffer(pos.clust.buff, width=16090, dissolve=TRUE)
v <- extract(cell, circle)
v.vec <- unlist(v)
bern.vec <- 1:12100

#2018 accuracy
bern.rule.18 <- replace(bern.vec, v.vec, max(E.y.18))
bern.rule.18 <- replace(bern.rule.18, bern.rule.18>max(E.y.18), 1e-6)
bern.rule.rast.18 <- setValues(cell,bern.rule.18)
detection.rule.18 <- extract(bern.rule.rast.18, all.dat.sub.18[,c(10,11)])
bern.out.rule.18 <- sum(-2*dbinom(all.dat.sub.18[,8], 1, detection.rule.18, log=TRUE))

#2019 accuracy
bern.rule.19 <- replace(bern.vec, v.vec, max(E.y.19))
bern.rule.19 <- replace(bern.rule.19, bern.rule.19>max(E.y.19), 1e-6)
bern.rule.rast.19 <- setValues(cell,bern.rule.19)
detection.rule.19 <- extract(bern.rule.rast.19, all.dat.sub.19[,c(10,11)])
bern.out.rule.19 <- sum(-2*dbinom(all.dat.sub.19[,8], 1, detection.rule.19, log=TRUE))

#total for simulation
bern.out.rule <- bern.out.rule.18 + bern.out.rule.19


###################Figures#############################

##Import and project layers
coords <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
bounding.box <- shapefile("MI_BoundingBox.shp")
bounding.box.proj <- spTransform(bounding.box, CRS=coords)
mi.diss <- shapefile("MI_diss.shp") 
midiss.proj <- spTransform(mi.diss, CRS=coords)
counties <- shapefile("Counties_v17a.shp")
counties.proj <- spTransform(counties, CRS=coords)
river<- shapefile("LargeRivers.shp")
river.proj <- spTransform(river, CRS=coords)
circle.proj <- spTransform(circle, crs("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))

####################################
#prepare spatial objects for plots
mi <- list("sp.polygons", midiss.proj , col='black', lwd=3)
co <- list("sp.polygons", counties.proj, col='gray', lwd=1)
pos.pred <- list("sp.points", pos.clust.pred, col='black', pch=19, lwd=2)
pos.fore <- list("sp.points", pos.clust.fore, col='red4', pch=3, lwd=4)
bound <- list("sp.polygons", bounding.box.proj, col="black", lwd=3, first=FALSE)
boundg <- list("sp.polygons", bounding.box.proj, col="gray50", lwd=3, first=FALSE)
North <- list("SpatialPolygonsRescale", layout.north.arrow(type=1), 
              offset = c(-89,42.65), scale = .5)
riv <- list("sp.polygons", river.proj, col='gray40', lwd=.5)
tenmile <- list("sp.polygons", circle.proj, col='black', lwd=3, lty=1, first=FALSE)

#100 km scale
scale = list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(-89.5,42.25), scale = 1.21, fill=c("transparent","black"))
text1 = list("sp.text", c(-89.5,42.15), "0")
text2 = list("sp.text", c(-88.29,42.15), "100 km")

#figure 2A
par(font.axis = 2, font=2, cex=1.5, cex.axis=1.5)
spplot(midiss.proj, fill=NA, colorkey=FALSE, xlim=c(-90.41829 -0.4, -82.41348 +0.4), ylim=c(41.69613 -0.4, 48.26269 +0.4), sp.layout=list(co,pos.fore,bound,pos.pred, North, scale, text1, text2), scales = list(draw = T, cex=2), key=list(title = "CWD Detections", cex.title = 2, x = .8, y = .95, corner=c(1.1,1),points=list(col=c("black", "red4"), pch=c(19,3), cex=2, lwd=c(2,4)), text=list(c("2015-2017", "2018-2019"), cex=2)))

#figure 2B
result.cv.proj.17 <- projectRaster(E.y.rast.17, crs="+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")
breaks=seq(0,max(E.y.19), by=0.0001)

#Color palette
colfunc <- colorRampPalette(c("white","brown1","brown2","brown4")) 
color.grad <- c(colfunc(length(breaks)))

#plot
spplot(result.cv.proj.17, colorkey=F,  
       ylim=c(42.45952 -0.4, 44.07051 +0.4), xlim=c(-86.21033 -0.4, -84.00013 +0.4), 
       sp.layout=list(bound, riv, co), at = breaks,
       scales = list(draw = T, cex=2), col.regions = color.grad,
       plot.margin=unit(c(1,15,1.5,1.2),"cm"))
       
#Figure 2C
result.cv.proj.18 <- projectRaster(E.y.rast.18, crs="+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")

#plot
spplot(result.cv.proj.18, colorkey=F,  
       ylim=c(42.45952 -0.4, 44.07051 +0.4), xlim=c(-86.21033 -0.4, -84.00013 +0.4), 
       sp.layout=list(bound, riv, co), at = breaks,
       scales = list(draw = T, cex=2), col.regions = color.grad,
       plot.margin=unit(c(1,15,1.5,1.2),"cm"))

#figure 2D
result.cv.proj.19 <- projectRaster(E.y.rast.19, crs="+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")

spplot(result.cv.proj.19, colorkey=F,  
       ylim=c(42.45952 -0.4, 44.07051 +0.4), xlim=c(-86.21033 -0.4, -84.00013 +0.4), 
       sp.layout=list(bound, riv, co), at = breaks,
       scales = list(draw = T, cex=2), col.regions = color.grad,
       plot.margin=unit(c(1,15,1.5,1.2),"cm"))

###Figure 3
quant <- quantile(E.y.17, probs=c(0.5, 0.75, 0.95))

#50% percentile
E.y.rast.2<- raster(, nrows = nsite^0.5, ncols = nsite^0.5, xmn = s[1], xmx = s[2], ymn = s[3], ymx = s[4], crs = coord.ll)
E.y.rast.2[] <- E.y.17
E.y.rast.2[E.y.rast.2 < quant[1]]=-Inf
r.2 <- E.y.rast.2 > -Inf
r.proj.2 <- projectRaster(r.2, crs="+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")
r.proj.2[r.proj.2 == 0] <- NA

# convert to polygon
mean.2 <- rasterToPolygons(r.proj.2, dissolve=TRUE)
m.2 <- gUnaryUnion(mean.2)
mean.ci.2 <- list("sp.polygons", m.2, col='black', lwd=3, lty=1, first=FALSE)

#75 percentile
E.y.rast.3 <- raster(, nrows = nsite^0.5, ncols = nsite^0.5, xmn = s[1], xmx = s[2], ymn = s[3], ymx = s[4], crs = coord.ll)
E.y.rast.3[] <- E.y.17
E.y.rast.3[E.y.rast.3 < quant[2]]=-Inf
r.3 <- E.y.rast.3 > -Inf
r.proj.3 <- projectRaster(r.3, crs="+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")
r.proj.3[r.proj.3 == 0] <- NA

#convert to polygon
mean.3 <- rasterToPolygons(r.proj.3, dissolve=TRUE)
m.3 <- gUnaryUnion(mean.3)
mean.ci.3 <- list("sp.polygons", m.3, col='black', lwd=3, lty=2, first=FALSE)

#95 percentile
E.y.rast.4 <- raster(, nrows = nsite^0.5, ncols = nsite^0.5, xmn = s[1], xmx = s[2], ymn = s[3], ymx = s[4], crs = coord.ll)
E.y.rast.4[] <- E.y.17
E.y.rast.4[E.y.rast.4 < quant[3]]=-Inf
r.4 <- E.y.rast.4 > -Inf
r.proj.4 <- projectRaster(r.4, crs="+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")
r.proj.4[r.proj.4 == 0] <- NA

#convert to polygon
mean.4 <- rasterToPolygons(r.proj.4, dissolve=TRUE)
m.4 <- gUnaryUnion(mean.4)
mean.ci.4 <- list("sp.polygons", m.4, col='black', lwd=3, lty=3, first=FALSE)


#blank raster for plot
plot.rast <- raster(, nrows = nsite^0.5, ncols = nsite^0.5, xmn = s[1], xmx = s[2], ymn = s[3], ymx = s[4], crs="+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")
plot.rast[] <- 0

#Figure 3A
spplot(plot.rast, fill=NA, colorkey=FALSE, xlim=c(-86.21033 -0.4, -84.00013 +0.4), ylim=c(42.45952 -0.4, 44.07051 +0.4), 
       sp.layout=list(boundg, riv, co, pos.fore, pos.pred, mean.ci.2, mean.ci.3, mean.ci.4 ), scales = list(draw = T, cex=2), key=list(title = "", cex.title = 2, x = .99, y = 1.055, corner=c(1.1,1), lines=list(col=c("black", "black", "black", "black", "red4"), pch=c(NA, NA, NA, 19, 3), type=c("l","l", "l","p", "p"), lty=c(1,2,3,1,1), lwd=c(4,4,4,4,4)), text=list(c("50th","75th", "95th","CWD+ 2015-17", "CWD+ 2018-19"), cex=c(2,2,2,2,2))))
       grid.text(label="A", x=0.28, y=0.88, gp=gpar(fontsize=30, col="black", fontface="bold")) 
       
#figure 3B
spplot(plot.rast, fill=NA, colorkey=FALSE, xlim=c(-86.21033 -0.4, -84.00013 +0.4), ylim=c(42.45952 -0.4, 44.07051 +0.4), 
       sp.layout=list(bound, riv, co, pos.fore, pos.pred, tenmile), scales = list(draw = T, cex=2),
       key=list(title = "", cex.title = 2, x = .95, y = .995, corner=c(1.1,1), lines=list(col=c("black"), lty=c(1), lwd=4), text=list(c("16.09 km"), cex=2)))
       grid.text(label="B", x=0.28, y=0.88, gp=gpar(fontsize=30, col="black", fontface="bold")) 

##Reference##
##Hefley, T. J., Hooten, M. B., Russell, R. E., Walsh, D. P., and Powell, J. A. (2017). When mechanism matters: Bayesian forecasting using models of ecological diffusion. Ecology Letters, 20:640–650.##
