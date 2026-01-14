# Chickenpox cases

chickenpox <- read.csv("hungary_chickenpox.csv")
mat <- as.matrix(chickenpox[, -1])  # Date entfernen
colnames(mat) <- colnames(chickenpox)[-1]
rownames(mat) <- chickenpox$Date
mat <- t(mat)
mat <- mat[order(rownames(mat)),]
chickenpox <- mat

# Spatial Weight Matrix

chickenpox_adjacency <- read.csv("hungary_county_edges.csv")
W <- matrix(0, nrow = 20, ncol = 20)
rn <- row.names(mat)
row.names(W) <- rn
colnames(W) <- rn
for(i in seq(20)){
  cha <- subset(chickenpox_adjacency, name_1 == rn[i])
  W[i, rn %in% cha$name_2] <- 1
}
diag(W) <- 0
W <- W / rowSums(W)
W_hungary <- list(diag(20), W)

# Population in Hungary
population_hungary <- read.csv2("stadat-nep0034-22.1.2.1-hu.csv", skip = 1)
population_hungary <- tail(population_hungary, 30) #Data is seperated in Male, Female, Total. We want total
population_hungary <- subset(population_hungary, grepl("capital|county", Level.of.territorial.units))

# Replace region names with row-names of 'chickenpox'
population_hungary$Name.of.territorial.units <- c("BUDAPEST", "PEST", "FEJER", "KOMAROM",
                                 "VESZPREM", "GYOR", "VAS", "ZALA", "BARANYA",
                                 "SOMOGY", "TOLNA", "BORSOD", "HEVES", "NOGRAD",
                                 "HAJDU", "JASZ", "SZABOLCS", "BACS", "BEKES", "CSONGRAD")


population_hungary <- population_hungary[order(population_hungary$Name.of.territorial.units),]
population_hungary <- subset(population_hungary, select = c("Name.of.territorial.units", paste0("X", 2005:2015)))

y <- do.call("cbind", population_hungary[-1])
y <- matrix(as.numeric(gsub(" ", "", y)), nrow = 20)
row.names(y) <- population_hungary$Name.of.territorial.units

date <- colnames(chickenpox)
date <- as.Date(date, format = "%d/%m/%Y")
table(format(date, "%Y"))


population_hungary <- apply(y, 1, function(x){
  unlist(mapply(function(a, b, d) return(a + seq(from = 0, length.out = d) * b / d), a = x[-length(x)],
                b = diff(x), d = c(52, 52, 53, 52, 52, 52, 52, 53, 52, 52)))
})

population_hungary <- t(population_hungary) / 100000
colnames(population_hungary) <- colnames(chickenpox)

## Save data in /data directory
save(chickenpox, W_hungary, population_hungary, file = "../chickenpox.rda", compress = "xz")

