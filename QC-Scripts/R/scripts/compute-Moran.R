#############
# Function to compute Moran's i statistic on som emapped continous variable
#############

#install.packages("~/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/spdep_0.6-13.tar.gz",repos=NULL,type="source")
# Installed some dependencies from source files.

library(spdep)

getSpatialGrid <- function(myMap,cellSize,extent=NULL){
    
    print("Making spatial grid object...")
    
    if(is.null(extent)) bb <- bbox(myMap) else bb=extent

    cs <- c(1,1)*cellSize # cell size (in units of proj4string(espMapTotalProj). i.e meters)
                                        #cs <- c(1,1)*100000 # cell size (in units of proj4string(espMapTotalProj). i.e meters)

                                        # 1 ft = 3.28084 m;  1degree X-direction = 93606.42 meters
    cc <- bb[, 1] + (cs/2)  # cell offset
    cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
    grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
    grd

    sp_grd2 <- SpatialGridDataFrame(grd,
                                    data=data.frame(id=1:prod(cd)),
                                    proj4string=CRS(proj4string(myMap)))
    # how many grid points?
    dim(sp_grd2@data)

    sp_grd2Points = sp_grd2; gridded(sp_grd2Points) = FALSE  # convert grid to points (centres of grids)

                                        # in the borders
#    Notinthesea = which(!is.na(over(sp_grd2,myMap)[,1]))

#    sp_grd3 <- sp_grd2[Notinthesea,]
    sp_grd3 <- sp_grd2
    return(sp_grd3)
}
    

gridMeans <- function(coords,values,mySpatialGrid){

    myValuesPoints = SpatialPointsDataFrame(coords, as.data.frame(values), coords.nrs = numeric(0), proj4string = CRS(proj4string(mySpatialGrid)),bbox = NULL)

    #print("Aggregating values within each spatial grids...")

    agg = aggregate(myValuesPoints, mySpatialGrid, mean)

    return(agg)
}


getWeights <- function(myMap,agg,NotintheseaOrNA,method="nearestNeighbours",nNeighbours=NULL,sdGaussian=NULL){
    
                                        # get distance between centres of the grid. These determine the weights of the moran statistic

    if(!method %in% c("distance")){
                                        # K nearest neighbours binary.
        aggPoints = agg[NotintheseaOrNA,]; gridded(aggPoints) = FALSE
                                        #    print(agg)
        print(dim(aggPoints@coords))
        
        if( method == "nearestNeighbours" ){
            
            print(paste0("Computing weights for ",length(NotintheseaOrNA),"x",length(NotintheseaOrNA)," matrix, using ",nNeighbours," neighbours in statistic..."))
            print(date())

            w = knn2nb(knearneigh(aggPoints@coords,k=nNeighbours))
            weightsList = nb2listw(w)
        }

        if( method == "Gaussian" ){
            
            print(paste0("Computing weights for ",length(NotintheseaOrNA),"x",length(NotintheseaOrNA)," matrix, using ",sdGaussian," as sd in weights..."))
            print(date())
            
            d <- as.matrix(dist(coordinates(agg)[NotintheseaOrNA,]))
            w = dnorm(d/10000,sd=sdGaussian) # a gaussian decay function in units of 10km

            # truncate anything further than 1000Km away 
            tol = dnorm(100,sd=sdGaussian)
            w[w<tol] = 0
            diag(w) <- 0
            
            weightsList = mat2listw(w)
        }
        if( method == "Gaussian2" ){
            
            print(paste0("Computing weights for ",length(NotintheseaOrNA),"x",length(NotintheseaOrNA)," matrix, using ",sdGaussian," as sd in weights..."))
            print(date())

            # Set SD to be proportional to the number of points per grid.
            d <- as.matrix(dist(coordinates(agg)[NotintheseaOrNA,]))
            w = dnorm(d/10000,sd=max(nPointsPerGrid,na.rm=TRUE)*sdGaussian/nPointsPerGrid) # a gaussian decay function in units of 10km

            # truncate anything further than 1000Km away 
            tol = dnorm(100,sd=sdGaussian)
            w[w<tol] = 0
            diag(w) <- 0
            
            weightsList = mat2listw(w)
        }

        
    } else {

        print(paste0("Computing weights for ",length(NotintheseaOrNA),"x",length(NotintheseaOrNA)," matrix, using basic distances..."))
        print(date())

        w <- 1/as.matrix(dist(coordinates(agg)[NotintheseaOrNA,])) # 1/distance
        diag(w) <- 0
        weightsList = mat2listw(w) # This is super slow! but hopefully not tooooo slow.
        
    }

    return(weightsList)
}



myMoran <- function( myMap=uk0,values,coords,mySpatialGrid=NULL,agg=NULL,cellSize=10000,colors=magma,baseMap,nNeighbours=NULL,weightsList=NULL,makeMap=TRUE){
    
# NOTE: cellSize is in the units of the map (e.g meters)
# coords = cbind(x,y)
# values = PCs[["PC1"]]

     # remove any NAs    
    notNA = (!is.na(rowSums(coords[,1:2]))) & (!is.na(values))

    print(paste0("Computing statistics using ",sum(notNA)," samples with non-missing data."))
    
    coords = coords[notNA,]
    values = values[notNA]
    
    if(is.null(mySpatialGrid)) mySpatialGrid = getSpatialGrid(myMap,cellSize,extent=rbind(range(coords[,1]),range(coords[,2])))

                                        # Get gridded means (this takes a while, but makes the next bit faster)
                                        #    cellSize = 10000
    if(is.null(agg)) {
        print(paste0("Aggregating over cells of size ",mySpatialGrid@grid@cellsize[1],"..."))
        print(date())
        agg = gridMeans(coords,values,mySpatialGrid)
    }

                                        # only apply the test to non-NA areas.
    Notinthesea = which(!is.na(over(agg,myMap)[,1]))
    NotintheseaOrNA = intersect(Notinthesea,which(!is.na(agg@data[,1])))

    
    if(is.null(weightsList)) weightsList = getWeights(myMap,agg,NotintheseaOrNA,method="nearestNeighbours",nNeighbours,sdGaussian) else print("Using supplied weightsList")
            
    myValues = agg@data$values[NotintheseaOrNA]
    
    print("Computing moran statistics...")
    print(date())

    moranResults = moran.test(myValues,weightsList)
    print("Done computing moran statistics...")
    print(date())

    if(makeMap){
                                        # Plot the aggregated values
        print("Making map...")

        fixedLims = range(myValues)
        colorCuts =  seq(fixedLims[1],fixedLims[2],length.out=200)
        colorVector = colour.scale(colorCuts,colourSet=colors,nBreaks=200,fixedLims=fixedLims)
        
        gridMap = spplot(agg[NotintheseaOrNA,],zcol="values",col.regions=colorVector,cuts=200,at=colorCuts)
        map = update(baseMap,col.regions="transparent",main=paste("PC",i)) + gridMap
        
    } else {
        map=NULL
    }
    
    return(list("map"=map,"moranResults"=moranResults))
    
}



    
