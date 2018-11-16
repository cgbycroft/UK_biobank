#################
# Convert to mapping coordinates
#################

# Clare actually ran this whole thing on her laptop, as it requires the correct rgeos libraries etc. This is just for records.

# THESE MAPS WERE DOWNLOADED FROM: http://www.gadm.org/download

setwd('/Users/clare/Documents/ClareDPhil/DPhil/UKBiobank_V2/testMirror/QC-Scripts/R')


GBR = readRDS('~/Downloads/GBR_adm0.rds')
IRL = readRDS('~/Downloads/IRL_adm0.rds')
IMN = readRDS('~/Downloads/IMN_adm0.rds')

GBR = readRDS('~/Downloads/GBR_adm0.rds')
IRL = readRDS('~/Downloads/IRL_adm0.rds')
IMN = readRDS('~/Downloads/IMN_adm0.rds')

# Merge into one big country!
GBR2 <- spChFIDs(GBR, paste("GBR", row.names(GBR), sep="."))
IRL2 <- spChFIDs(IRL, paste("IRL", row.names(IRL), sep="."))
IMN2 <- spChFIDs(IMN, paste("IMN", row.names(IMN), sep="."))
all = rbind(GBR,IRL2,IMN2)

# Simplify the polygons
allSimple <- gSimplify(all,tol=0.01,topologyPreserve=T)
allSimpledf  <- SpatialPolygonsDataFrame(allSimple,data=all@data,match.ID=T)

# Convert to british national grid
bng = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs' 

allSimpledf_bng = spTransform(allSimpledf,CRSobj = CRS(bng))

british_isles_sp_map_lonLat = allSimpledf
british_isles_sp_map_bng = allSimpledf_bng

save(british_isles_sp_map_lonLat,british_isles_sp_map_bng,file="GBR_IRL_IMN_sp_maps.RData")

