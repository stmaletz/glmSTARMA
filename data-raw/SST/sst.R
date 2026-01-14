# Install STRbook if necessary
# https://github.com/andrewzm/STRbook
# install_github("andrewzm/STRbook")
library("STRbook")
library("Matrix")

data(SST_df)
SST <- subset(SST_df, lon >= 160 & lon <= 240)
SST <- subset(SST, lat >= -29 & lat <= 29)

times <- seq(from = as.Date("1970-01-01"), to = as.Date("2002-12-01"), by = "m")
times <- format(times, "%b %Y")

ts <- matrix(NA, nrow = nrow(unique(cbind(SST$lon, SST$lat))), ncol = length(times))
for(i in seq_along(times)){
  temp <- subset(SST, date == times[i])
  temp <- temp[order(temp$lon, temp$lat),]
  ts[, i] <- temp$sst
}
colnames(ts) <- times


locations <- subset(SST, select = c("lat", "lon"))
locations <- unique(locations)
locations <- locations[order(locations$lon, locations$lat),]
rownames(locations) <- NULL


# Neighbors up to order 4 based on  queen contiguity
W <- as.matrix(dist(locations, method = "maximum"))
W1 <- 1 * (W == 2)
W1 <- W1 / colSums(W1)
W2 <- 1 * (W == 4)
W2 <- W2 / colSums(W2)
W3 <- 1 * (W == 6)
W3 <- W3 / colSums(W3)
W4 <- 1 * (W == 8)
W4 <- W4 / colSums(W4)

# use sparse matrices to save memory
W1s <- as(W1, "dgCMatrix")
W2s <- as(W2, "dgCMatrix")
W3s <- as(W3, "dgCMatrix")
W4s <- as(W4, "dgCMatrix")

W_queen <- list(as(diag(1230), "dgCMatrix"), W1s, W2s, W3s, W4s)
#W_queen <- list(diag(1230), W1, W2, W3, W4)


## Directed neighbors
W_n <- W_s <- W_e <- W_w <- matrix(0, nrow = 1230, ncol = 1230)
for(i in seq(nrow(locations))){
  
  north <- which(locations$lat == locations$lat[i] + 2 &
                   locations$lon == locations$lon[i])
  if(length(north) > 0)
    W_n[i, north] <- 1
  
  south <- which(locations$lat == locations$lat[i] - 2 &
                   locations$lon == locations$lon[i])
  if(length(south) > 0)
    W_s[i, south] <- 1
  
  east  <- which(locations$lat == locations$lat[i] &
                   locations$lon == locations$lon[i] + 2)
  if(length(east) > 0)
    W_e[i, east] <- 1
  
  west  <- which(locations$lat == locations$lat[i] &
                   locations$lon == locations$lon[i] - 2)
  if(length(west) > 0)
    W_w[i, west] <- 1
}

W_directed <- list(
  as(diag(1230), "dgCMatrix"),             
  north = as(W_n, "dgCMatrix"),
  east  = as(W_e, "dgCMatrix"),
  south = as(W_s, "dgCMatrix"),
  west  = as(W_w, "dgCMatrix")
)
# W_directed <- list(diag(1230), north = W_n, east = W_e, south = W_s, west = W_w)

SST <- ts
save(SST, locations, W_directed, file = "../sst.rda", compress = "xz")

