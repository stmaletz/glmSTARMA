library("data.table")
library("stringdist")
library("spdep")
library("dplyr")
library("sf")

# shape files are retrieved from
# https://github.com/ostojanovic/BSTIM
# probably preprocessed
# orginal source is 
# Bundeamt für Kartographie und Geodäsie (BKG) (2018)
# Don't know which version

# Read in shape file and merge regions that were merged in 2021
shape <- read_sf("germany_county_shapes.json")
shape$RKI_NameDE <- gsub("\u0096", " ", shape$RKI_NameDE)
shape$RKI_NameDE <- gsub("-", " ", shape$RKI_NameDE)
# merge eisenach with wartburgkreis
any(shape$RKI_NameDE == "SK Eisenach")
any(shape$RKI_NameDE == "LK Wartburgkreis")
wartburg_new <- shape %>% 
  filter(RKI_NameDE %in% c("SK Eisenach", "LK Wartburgkreis")) %>% 
  summarise(
    RKI_NameDE = "LK Wartburgkreis",
    RKI_ID = "16063",
    .groups = "drop"
  )
shape <- shape %>% 
  filter(!RKI_NameDE %in% c("SK Eisenach", "LK Wartburgkreis"))
shape <- bind_rows(shape, wartburg_new)


# Import rota counts
# These were downloaded from the RKI (survstat.rki.de) as weekly reports in csv format
# for the years 2001 to 2024 
files <- list.files("raw", recursive = TRUE, full.names = TRUE, pattern = ".csv")
years <- lapply(files, function(file){
  read.csv2(file, fileEncoding = "UTF-16LE", skip = 1, sep = "\t", row.names = 1)
})
names(years) <- sub(".*/(\\d{4})/.*", "\\1", files)
for(i in seq_along(years)){
  years[[i]]$jahr_woche <- paste0(names(years)[i], "-", row.names(years[[i]]))
}

years <- rbindlist(years)
jahr_woche <- years$jahr_woche
years$jahr_woche <- NULL
years$Unbekannt <- NULL
years <- as.matrix(years)
years[is.na(years)] <- 0
years <- t(years)
colnames(years) <- jahr_woche

regionen <- row.names(years)
regionen <- gsub("\\.", " ", regionen)
regionen <- gsub("a d ", "a.d. ", regionen)
regionen <- gsub(" a ", " a. ", regionen)
regionen <- gsub(" i ", " i. ", regionen)
regionen <- gsub("i d ", "i.d. ", regionen)
regionen <- trimws(regionen)
matches <- amatch(shape$RKI_NameDE, regionen, maxDist = 10)  
years <- years[matches,]
row.names(years) <- shape$RKI_NameDE

# Generate feature indicating former GDR regions
# Berlin is special case, since it was owned by both GDR and BRD (value 0.5)
## Overview about counties
# ID;Bundesland
# 01;Schleswig-Holstein
# 02;Hamburg
# 03;Niedersachsen
# 04;Bremen
# 05;Nordrhein-Westfalen
# 06;Hessen
# 07;Rheinland-Pfalz
# 08;Baden-Württemberg
# 09;Bayern
# 10;Saarland
# 11;Berlin -> owned by GDR and BRD
# 12;Brandenburg -> GDR
# 13;Mecklenburg-Vorpommern -> GDR
# 14;Sachsen -> GDR
# 15;Sachsen-Anhalt -> GDR
# 16;Thüringen -> GDR

shape$former_gdr <- 0.5 * (substr(shape$RKI_ID, 1, 2) == "11")
shape$former_gdr[substr(shape$RKI_ID, 1, 2) %in% c("12", "13", "14", "15", "16")] <- 1
gdr_feature <- shape$former_gdr


# Get information about population per county over time
# Note: population data in above repository is outdated (only until 2018)
# Downloaded full data from Statistisches Bundesamt (2025-12-08)
# https://www-genesis.destatis.de/datenbank/online/statistic/12411/table/12411-0015/table-toolbar

# We need to split data from Berlin in 12 parts (subregions), we do this by the proportion of population
# obtained from the 2022 census
# https://ergebnisse.zensus2022.de/datenbank/online/statistic/1000A/table/1000A-0000


census <- read.csv2("census_berlin.csv", skip = 4)
census$Personen <- NULL
census$e <- NULL
census$Anzahl <- NULL
names(census)[3] <- "Einwohner"
census <- subset(census, substr(DG, 1, 2) == "11")
census$prop <- census$Einwohner / sum(census$Einwohner)
census$DG <- gsub("0000000000", "0", census$DG)


# Read rest of population data
population <- read.csv2("12411-0015_de.csv", skip = 5)
names(population)[1:2] <- c("ID", "Name")
population <- subset(population, select = !grepl("X\\.", names(population)))
population[-(1:2)] <- lapply(population[-(1:2)], as.numeric)
names(population)[3:27] <- gsub("31.12.", "", names(population)[3:27])


fuse_regions <- function(old, new, when, pop = population){
  sub_data <- subset(pop, ID %in% old)
  sums <- sapply(sub_data[-(1:2)], sum)
  years <- 2000:2024
  pop[pop$ID == new, which(years < when) + 2] <- sums[years < when]
  pop <- pop[!(pop$ID %in% old),]
  return(pop)
}


# goettingen osterode
population <- fuse_regions(c("03152", "03156"), "03159", 2016)
# hannover
population <- fuse_regions(c("03201", "03253"), "03241", 2001)
# aachen
population <- fuse_regions(c("05354", "05313"), "05334", 2009)
# mecklenburg-vorpommern
population <- fuse_regions(c("13053", "13051"), "13072", 2011)
population <- fuse_regions(c("13005", "13061", "13057"), "13073", 2011)
population <- fuse_regions(c("13006", "13058"), "13074", 2011)
population <- fuse_regions(c("13054", "13060"), "13076", 2011)
population <- fuse_regions(c("13002", "13056", "13055", "13052"), "13071", 2011)
population <- fuse_regions(c("13001", "13059", "13062"), "13075", 2011)
# sachsen
population <- fuse_regions(c("14161"), "14511", 2008)
population <- fuse_regions(c("14191", "14171", "14188", "14181"), "14521", 2008)
population <- fuse_regions(c("14375", "14177", "14182"), "14522", 2008)
population <- fuse_regions(c("14178", "14166"), "14523", 2008)
population <- fuse_regions(c("14173", "14167", "14193"), "14524", 2008)
population <- fuse_regions(c("14262"), "14612", 2008)
population <- fuse_regions(c("14272", "14292", "14264"), "14625", 2008)
population <- fuse_regions(c("14263", "14286", "14284"), "14626", 2008)
population <- fuse_regions(c("14280", "14285"), "14627", 2008)
population <- fuse_regions(c("14287", "14290"), "14628", 2008)
population <- fuse_regions(c("14365"), "14713", 2008)
population <- fuse_regions(c("14379", "14383"), "14729", 2008)
population <- fuse_regions(c("14374", "14389"), "14730", 2008)
# sachsen-anhalt
population <- fuse_regions(c("15202"), "15002", 2007)
population <- fuse_regions(c("15303"), "15003", 2007)
population <- fuse_regions(c("15370"), "15081", 2007)
population <- fuse_regions(c("15355", "15362"), "15083", 2007)
population <- fuse_regions(c("15268", "15256"), "15084", 2007)
population <- fuse_regions(c("15260", "15266"), "15087", 2007)
population <- fuse_regions(c("15265", "15261"), "15088", 2007)
population <- fuse_regions(c("15363"), "15090", 2007)

population[population$ID == "15001", 3:9] <- population[population$ID == "15101", 3:9] + 0.2 * population[population$ID == "15151", 3:9]
population[population$ID == "15082", 3:9] <- population[population$ID == "15154", 3:9] + population[population$ID == "15159", 3:9] + 0.45 * population[population$ID == "15151", 3:9]
population[population$ID == "15085", 3:9] <- population[population$ID == "15357", 3:9] + population[population$ID == "15364", 3:9] + population[population$ID == "15369", 3:9] + 0.05 * population[population$ID == "15352", 3:9]
population[population$ID == "15086", 3:9] <- population[population$ID == "15358", 3:9] + 0.05 * population[population$ID == "15151", 3:9]
population[population$ID == "15089", 3:9] <- population[population$ID == "15153", 3:9] + population[population$ID == "15367", 3:9] + 0.95 * population[population$ID == "15352", 3:9]
population[population$ID == "15091", 3:9] <- population[population$ID == "15171", 3:9] + 0.3 * population[population$ID == "15151", 3:9]
population <- subset(population, !(ID %in% c("15101", "15151", "15153", "15154", "15159", "15171", "15352", "15357", "15358", "15367", "15369", "15364")))
# thueringen
population[population$ID == "16063", 3:23] <- population[population$ID == "16056", 3:23] + population[population$ID == "16063", 3:23]
population <- subset(population, !(ID %in% c("16056")))
population <- population[1:400,]
#berlin
berlin <- as.data.frame(matrix(NA, nrow = 12, ncol = 27))
names(berlin) <- names(population)
berlin$ID <- census$DG
berlin$Name <- census$Deutschland
for(i in 3:27){
  berlin[[i]] <- census$prop * population[[i]][population$ID == "11000"]
}
population <- rbind(population, berlin)
population <- subset(population, ID != "11000")
population$Name <- shape$RKI_NameDE[match(population$ID, shape$RKI_ID)]
population <- population[match(shape$RKI_NameDE, population$Name),]


pop <- as.matrix(population[, 3:27])

# interpolate population
population_germany <- apply(pop, 1, function(x){
  unlist(mapply(function(a, b, d){
    return(seq(from = a, to = b, length.out = d))}, a = x[-length(x)],
                b = x[-1], d = c(52, 52, 52, 53, 52, 52, 52, 52, 53, 52, 52, 52, 52, 52, 53, 52, 52, 52, 52, 53, 52, 52, 52, 52)))
})

population_germany <- t(population_germany / 100000)



nb <- poly2nb(shape, row.names = shape$RKI_NameDE)
W1 <- nb2mat(nb, style="W")
# row.names(W1) <- colnames(W1) <- shape$RKI

W2 <- 1 * (W1 > 0)
W2 <- W2 %*% W2
diag(W2) <- 0
W2[W1 > 0] <- 0
W2 <- W2 / colSums(W2)


## name the data set properly
rota <- years
gdr_feature
population_germany
W_germany <- list(diag(411), W1, W2)


save(rota, gdr_feature, population_germany, W_germany, file = "../../data/rota.rda", compress = "xz")
