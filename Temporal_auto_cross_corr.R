#***************************#
# Intro to Sp_T Variography
#***************************#

library(spacetime)
data(air) # air quality pm10 measurements voer Germany
ls()

library("dplyr")
library("tidyr")


#===========#
# Variography
#===========#

#============================================#
# Temporal autocorrelation & cross correlation
#============================================#

str(air)
# num [1:70, 1:4383] NA NA NA NA NA NA NA NA NA NA ...
#- attr(*, "dimnames")=List of 2
#..$ : chr [1:70] "DESH001" "DENI063" "DEUB038" "DEBE056" ...

if(!exists("rural")) {
  rural <- STFDF(stations, dates, data.frame(PM10 = as.vector(air)))
}

rr <- rural[, "2005::2010"]
str(rr) # Formal class 'STFDF' [package "spacetime"] with 4 slots

rts <- as(rr, "xts")
str(rts)

unsel <- which(apply(rts, 2, function(x) all(is.na(x)))) # by col
r5to10 <- rr[-unsel, ]
summary(r5to10)
# Length  Class   Mode 
# 96778  STFDF     S4 

str(r5to10)


## arbitraily selet 4 stations
rn <- row.names(r5to10@sp)[4:7]


#----------------------------#
#  Univariate autocorrelation
#-----------------------------#

par(mfrow = c(2, 2))
par(mai = rep(0.5, 4))
for(i in rn) {
  # one station i per plot, univariate
  acf(na.omit(r5to10[i, ]), main = i)
}
par(mfrow = c(1,1))


#-------------------------#
# Auto & cross correlation 
#-------------------------#

## arbitraily selet 4 stations
rn <- row.names(r5to10@sp)[4:7]

acf(na.omit(as(r5to10[rn, ], "xts")))
# auto-corrleation for lag-0 is always 1
# cross-correlation for lag-0 is NOT alway 1
# cross-correlation for lag-h can be asymmetric
  # corr(Z(sA, t), Z(sB, t+h)) = corr(Z(sB, t), Z(sA, t+h))
  # != corr(Z(sA, t), Z(sB, t-h))

# the plot does not show strong asymmetrix, 
# cross-correlation yet quite strong, similar to auto-correlation


# draw a 7*7 plot
acf(na.omit(as(r5to10[4:10, ], "xts")))
  # notice DESH & DESN04 with nearly no cross-correlation
  # have to do spatial distance betw two stations. 


#------------------------------#
# Examine dist btw two stations
#------------------------------#

library(sp)

str(r5to10[4:10, ]) # Formal class 'STFDF' [package "spacetime"] with 4 slots

print(spDists(r5to10[4:10, ]@sp), digits = 3)



#================================#
# Spatial correlation, variograms
#================================#

# sample 100 time instances randomly
dim(r5to10) 
#   space      time variables 
#     53      1826         1 
rs <- sample(dim(r5to10)[2], 100)
# select these instances as a SpatialPointsDataFrame
# and add a time index to them
head(rs) # [1]  751 1350 1155  338  689  963


# then bind them together in a single SpatialPointsDataFrame
# which has time index ti:

str(rs) # int [1:100] 751 1350 1155 338 689 963 1780 1068 841 1065


a <- r5to10[, 1]# formal class SpatialPointsDataFrame 5 slots
b <- r5to10[, 2]
str(a)
str(b)

c <- r5to10[1, ]
str(c)
# An ‘xts’ object on 2005-01-01/2009-12-31 containing:
# Data: num [1:1826, 1:2] 16.7 13.9 16.6 16 21.2 ...
#- attr(*, "dimnames")=List of 2
#..$ : NULL
#..$ : chr [1:2] "PM10" "timeIndex"


lst <- lapply(rs, function(i) {x <- r5to10[, i]; x$ti <- i; rownames(x@coords) = NULL; x})
# large list 100 elements


pts <- do.call("rbind", lst) # formal class SpatialPointsDataFrame 5 slots
str(pts)
tail(pts@data$ti)
# Formal class 'SpatialPointsDataFrame' [package "sp"] with 5 slots
#..@ data       :'data.frame':	5300 obs. of  2 variables:
#  .. ..$ PM10: num [1:5300] 12.21 12.25 NA 9.92 9.26 ...
#.. ..$ ti  : int [1:5300] 751 751 751 751 751 751 751 751 751 751 ...
#..@ coords.nrs : num(0) 
#..@ coords     : num [1:5300, 1:2] 9.59 9.69 9.79 13.65 13.23 ...
#.. ..- attr(*, "dimnames")=List of 2
#.. .. ..$ : NULL
#.. .. ..$ : chr [1:2] "coords.x1" "coords.x2"
#..@ bbox       : num [1:2, 1:2] 6.28 47.81 14.79 54.92
#.. ..- attr(*, "dimnames")=List of 2
#.. .. ..$ : chr [1:2] "coords.x1" "coords.x2"
#.. .. ..$ : chr [1:2] "min" "max"
#..@ proj4string:Formal class 'CRS' [package "sp"] with 1 slot
#.. .. ..@ projargs: chr "+proj=longlat +datum=WGS84 +no_defs"

# 53 * 100 rows of coords PM10 ti

#=============================#
# compute the pooled variogram
#=============================#

library(gstat)

pts$PM10

V <- variogram(PM10 ~ ti,  pts[!is.na(pts$PM10), ])
plot(V)
str(V)


fit.variogram(v, vgm)

vgm # fit a variogram model
vmod <- fit.variogram(V, vgm(1, "Exp", 200, 1))
str(vmod)
vmod 
#   model    psill    range
#1   Exp 43.75494 69.24329

# a rather poor fit as only 53 stations although
# the time resolutionis rich 1862 days

dim(r5to10)
#    space      time variables 
#      53      1826         1 


vv <- variogram(PM10~1, r5to10, width = 20, cutoff = 200, tlags = 0:5)
plot(vv)
