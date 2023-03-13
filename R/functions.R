## create modified extract function

extract.mod <- function(
    raster.obj, # The raster 
    spdf.obj    # The spatial points data frame
){
  
  ### Attempt blind extraction
  ext <- raster::extract(raster.obj, spdf.obj)
  
  ### Add to spdf 
  spdf.obj@data$extract.col <- ext
  
  ### Subset spdf to include only the NA points
  spdf.sub <- subset(spdf.obj, is.na(extract.col))
  
  ### Turn raster into spdf
  raster.spdf <- as(raster.obj, "SpatialPointsDataFrame")
  
  ### Get matrix of distances
  dists <- gDistance(spdf.sub, raster.spdf, byid = T)
  
  ### Find minimums within columns (the spdf sites are in the columns)
  dists.min <- apply(dists, 2, min)
  
  ### Loop to find the corresponding rows for the minima
  n.elements <- length(dists.min)
  min.rows <- vector("numeric", n.elements)
  
  for(i in 1:n.elements){
    
    min.rows[i] <- which(dists[,i] == dists.min[i])
    
  }
  
  ### Extract raster values at these rows
  na.nearest.values <- raster.spdf@data[min.rows, 1]
  
  ### Join back to spdf
  spdf.obj@data$extract.col[is.na(spdf.obj@data$extract.col)] <- na.nearest.values
  
  ### Return spdf
  return(spdf.obj)
  
}

