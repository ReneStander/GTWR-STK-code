# Example 1
# This script contains the code to simulate the data for Example 1
# And to apply the GTWR-STK method

# Packages

library(spatstat)
library(dplyr)
library(sf)
library(sp)
library(spacetime)
library(stpp)
library(ggplot2)
library(forecast)
library(reshape2)
library(spdep)
library(gstat)
library(units)
library(readxl)

# Import other R scripts
source("../Functions/timeseries_sim.R")
source("../Functions/import_dpt.R", chdir = T)

# Import shape file
dc <- read_sf("../Shapefiles/DC_SA_2011.shp")

# Simulate the first time step
domain <- dc # 52 areas
domain$id <- 1:nrow(domain)

domain <- domain |> st_transform(32735)

(n1 <- poly2nb(st_geometry(domain), queen = F))
list1 <- nb2listw(n1)

ctd <- st_centroid(domain)
ctd <- st_set_crs(ctd, st_crs(domain))

set.seed(17)
copy <- domain

simulation <- domain |> 
  st_geometry() |> 
  st_sample(1000) |> 
  st_as_sf()

copy$points <- domain |> st_geometry() |>
  st_contains(simulation) |>
  as.matrix() |> 
  rowSums()

tr <- sample(1:52, size = 1)
tr2 <- sample(n1[[tr]],size = 2)
(regions4trend <- c(tr,tr2))

reg1 <- setdiff(1:nrow(copy), regions4trend)
tr <- sample(reg1, size = 1)
tr2 <- sample(n1[[tr]],size = 1)
(regions4trend2 <- c(tr,tr2))

reg <- setdiff(reg1, regions4trend2)
(regions4ma <- sample(reg, size = ceiling(length(reg)/2), replace = F))

# regions 4 ar
regg <- setdiff(reg, regions4ma)
(regions4ar <- sample(regg, size = length(regg), replace = F))

t <- 5
ts <- matrix(NA, nrow = nrow(copy), ncol =  t)

ts[regions4trend,] <- lin_trend(copy[regions4trend,]$points, 
                                times = t, 
                                var_at = 9, 
                                beta = 10)

ts[regions4trend2,] <- lin_trend(copy[regions4trend2,]$points, 
                                 times = t, 
                                 var_at = 9, 
                                 beta = 15)

(ts[regions4ma,] <- ma1(copy[regions4ma,]$points, 
                        times = t, 
                        var_at = 9,
                        theta = 0.8))
ts[regions4ar,] <- ar1(copy[regions4ar,]$points, 
                       times = t, 
                       var_at = 9,
                       theta0 = 2,
                       phi = 0.8)

shape2 <- data.frame(id = 1:nrow(copy),
                     time = "Time 1",
                     points_ba = copy$points)

dates <- paste0("Time ",(1:(t+1)))

for(j in 2:length(dates)){
  d <- data.frame(id = 1:nrow(copy),time = dates[j])
  d$points_ba <- ts[,(j-1)]
  shape2 <- rbind(shape2, d)
}

write.csv(shape2, file = "data_ex1.csv")

dat <- shape2
gtwr_data <- dat |> dplyr::select(id, time, points_ba)
gtwr_data$lag <- as.numeric(unname(sapply(dat$time, function(x) strsplit(x, " ")[[1]][2])))

gtwr_data <- gtwr_data |> filter(lag != max(gtwr_data$lag))

betas_gtwr <- tibble(id = numeric(0), time = numeric(0), beta = numeric(0))
X <- rep(1, nrow(gtwr_data))

for(i in 1:nrow(gtwr_data)){
  dist <- st_distance(ctd[gtwr_data$id[i],"geometry"],
                      ctd[gtwr_data$id,"geometry"]) |> drop_units()
  dist <- dist/1000
  dist2 <- as.numeric(dist^2)
  lagg <- (gtwr_data$lag[i] - gtwr_data$lag)^2
  w <- diag(exp(-(dist2)/400)*exp(-(lagg)/2))
  beta <- solve(t(X)%*%w%*%X)%*%t(X)%*%w%*%gtwr_data$points_ba
  inp <- c(gtwr_data$id[i],
           gtwr_data$time[i],
           beta)
  betas_gtwr <- rbind(betas_gtwr, inp)
  names(betas_gtwr) <- c("id", "time", "beta")
  
}

betas_gtwr$beta <- as.numeric(betas_gtwr$beta)
betas_gtwr$residuals <- gtwr_data$points_ba - betas_gtwr$beta

ids <- unique(betas_gtwr$id)
new <- data.frame(id = ids)
new$time <- paste0("Time ",t + 1)
new$beta <- NA
new$residuals <- NA

for(i in ids){
  d <- betas_gtwr |> filter(id == i)
  ts <- ts(d$beta)
  HW1 <- ets(y=ts)
  HW1.pred <- predict(HW1, 2, prediction.interval = TRUE, level=0.95)
  new$beta[new$id == i] <- HW1.pred$mean[1]
}

betas_gtwr <- rbind(betas_gtwr, new)
write.csv(betas_gtwr, file = "gtwr_ex2.csv")

# Kriging
dom.sp <- domain |> as("Spatial")
union <- domain |> st_union() |> as("Spatial")

x_range <- range(union@bbox[1,])
y_range <- range(union@bbox[2,])

resolution <- 50000
grid <- expand.grid(x = seq(from = x_range[1], to = x_range[2], by = resolution),
                    y = seq(from = y_range[1], to = y_range[2], by = resolution))
coordinates(grid) <- ~x + y
proj4string(grid) <- CRS("+init=epsg:32735") 

mask <- grid[dom.sp]

dat2 <- betas_gtwr

dm <- domain |> st_geometry() |> st_as_sf()
dm$id <- as.character(1:nrow(dm))
dat2 <- left_join(dat2, dm)

resid <- dat2 |> filter(time != paste0("Time ",t + 1))

pnts1 <- spsample(dom.sp, n = 150, "random") |> st_as_sf()
pnts2 <- spsample(dom.sp, n = 150, "nonaligned") |> st_as_sf()
pnts <- rbind(pnts1, pnts2) |> st_geometry()
ctd.g <- ctd |> st_geometry()
pnts1 <- c(pnts, ctd.g)

points.in <- st_contains(st_geometry(domain), pnts1)
samp <- st_as_sf(pnts1)

for(j in 1:t){
  eval(parse(text = paste0("samp$sim_p",j," <- NA")))
  
  for(i in 1:length(points.in)){
    d <- points.in[[i]]
    dd <- dat2 |> filter(time == paste0("Time ",j))
    eval(parse(text = paste0("samp$sim_p",j,"[d] <- dd$residuals[i]")))
  }
}

dates <- paste0("Time ",(1:t))

st.long <- samp |> select(x, sim_p1)
st.long$time <- dates[1]
st.long$id <- 1:nrow(samp)
names(st.long)[1] <- "sim"

for(i in 2:t){
  d <- eval(parse(text = paste0("samp |> select(x, sim_p", i, ")")))
  d$time <- dates[i]
  d$id <- 1:nrow(d)
  names(d)[1] <- "sim"
  
  st.long <- rbind(st.long,d)
}

st.long$id <- as.factor(st.long$id)

st.long <- st.long |>
  mutate(tt = recode(time, 
                     "Time 1" = "2010-01-01", 
                     "Time 2" = "2010-01-02", 
                     "Time 3" = "2010-01-03", 
                     "Time 4" = "2010-01-04", 
                     "Time 5" = "2010-01-05"))

st.long$TIME <- as.Date(st.long$tt, "%Y-%m-%d")
st.sp <- st.long |> as("Spatial")

SP <- SpatialPoints(st.sp@coords,CRS("+init=epsg:32735")) 
DF <- data.frame(sim=st.sp$sim) 
TM <- st.sp$TIME 

timeDF <- STIDF(SP,TM,data=DF) 

var <- variogramST(sim~1,data=timeDF,tunit = "days") 

productSum <- vgmST("productSum", 
                    space = vgm(300,"Sph", 20000, 20),
                    time = vgm(300,"Sph", 20000, 0.5),
                    k = 0.001)

productSum_Vgm <- fit.StVariogram(var, productSum,
                                  method = "L-BFGS-B")

dates_1 <- paste0("2010-01-0",(1:(t+1)))

tm.grid <- as.Date(dates_1[1:(t+1)],"%Y-%m-%d")
grid.ST <- STF(mask,tm.grid) 

pred_SM <- krigeST(sim~1,
                   data=timeDF, 
                   modelList=productSum_Vgm, 
                   newdata=grid.ST)

uns <- t(unstack(pred_SM))
mask_SM <- mask
mask_SM$pred <- uns[,(t+1)]

pred.final <- mask_SM |> st_as_sf()
pred.in <- st_contains(domain$geometry,pred.final)

domain$res <- NA

for(i in 1:length(pred.in)){
  d <- pred.in[[i]]
  domain$res[i] <- mean(pred.final$pred[d])
}

krig <- betas_gtwr |> filter(time == "Time 6")
krig$residuals <- domain$res
krig$residuals[is.na(krig$residuals)] <- 0


resid <- resid |> dplyr::select(-x)

final <- rbind(resid, krig)

pred <- final |> filter(time == paste0("Time ",t + 1))
pred$yhat <- pred$beta + pred$residuals

domain$id <- as.character(domain$id)
shape2$id <- as.character(shape2$id)

pred <- left_join(pred,domain) |> st_as_sf()

orig <- shape2 |> filter(time == paste0("Time ",t + 1))
orig <- left_join(orig, domain) |> st_as_sf()

















