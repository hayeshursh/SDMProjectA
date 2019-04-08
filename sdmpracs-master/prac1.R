####################################################################################
###### This code was created by Jen Cruz as part of the species distribution ######
#####   modeling class, co-taught with Benjamin Zuckerberg!                  #######
##### Practical 1: preparing and visualizing species and predictor data      ######
###################################################################################


########## clean workspace and load required packages ####################
###########################################################################
#####clean workspace to improve efficiency: ###
rm(list = ls() ) 
#set working dir() 
gc() #releases memory

# Note that you only need to install each package once ###
install.packages( "dplyr" ) 
install.packages( "tidyr" )
install.packages( "rgdal" )
install.packages( "raster" )
install.packages( "rasterVis" )
install.packages( "prism" )
install.packages( "psych" ) 

####### load relevant packages ###
#library( tidyr ) #combines dataframes
library( dplyr ) #manipulates dataframes
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( tidyr ) #spread and other dataframe functions
library( ggplot2 ) #fancy plots
#library( lubridate ) #easy date adjustments and calculations.
library( rgdal ) #imports projects transforms spatial data. Requires sp.
#library( rgeos ) #spatial manipulations e.g. buffer, interpolate.
#library( maptools ) #reads and manipulates ESRI shapefiles.
library( raster ) #manipulates rasters
library( rasterVis ) #visualises rasters
library( prism ) #extracts weather data from PRIS

# those additional packages, not being used in this exercise may still be useful when you are manipulating your own data. #
########## end of package loading ###########

########################################################################################
#######################    import relevant data    #####################################
#########################################################################################
# Set your working directory. In our case we use the current directory, but you should #
# adjust this to wherever you have your code and data:
workdir <- getwd() 

############### Importing and checking species data #################################
# The species data comes from the latest Wisconsin Bird Atlas (2016-2017): #
# http://www.uwgb.edu/birds/wbba/index.htm #
# See prac1notes.docx for details on both species and methodology. 

# Import species and site lcations from a .csv file that includes #
# presence/absence observations and easting and northings for site locations #
spp_df <- read.csv( file = paste( workdir, "/data/sppdata.csv", sep = "" ), #specify file name
                    header = TRUE, #keeps column labels
                    strip.white =TRUE ) #remove white spaces 
# Check that the data were imported the way you wanted it, i.e. no empty rows or columns #
glimpse( spp_df )
# View first and last rows of the dataframe:
head( spp_df ); tail( spp_df )
# During our in-class practicals we are going to be working with the Hermit thrush data.# 
# Data for the Eastern Meadowlark is available for those that need data for the #
# take home exercises. Note that you will have to modify your predictors to fit #
# the ecology of Eastern Meadowlarks in the take home exercises (our use your own data) #

# Here is where you are going to have to clean up your species data if necessary ###

# View how many presence and absences were observed 
table( spp_df$heth )
# What does this tell you about your sample size?

# View spatial distribution of your observations: 
ggplot( data = spp_df, aes( x = x, y  = y, group = factor(heth), color = factor(heth) ) ) + 
        geom_point()  

# What can you see regarding the observed spatial distribution for hermit thrush?
# What might it tell you about the ecology? 
# What other checks do you need to make for your own data? Is the sampling biased? 

# We need to match our predictor data to our site locations. We thus need to #
# extract locations from our imported dataframe. #
# Site locations were recorded using NAD_1983_HARN_Wisconsin_TM. We used this to define #
# out coordinate system first: #
siteproj <- sp::CRS(  "+init=epsg:3071" )
# An easy way to define coordinate systems is by using their epsg code. #
# You can find a searchable list of codes at http://www.spatialreference.org #

# Define site locations by extracting eastings and northings from the dataframe, #
# converting them to spatial points with the correct coordinate system #
sites <- SpatialPoints( spp_df[, c( 'x', 'y' ) ], proj4string = siteproj )
# View attributes and first rows:
class( sites ); head( sites )
# Check accuracy of your site locations. Here we are just plotting all points to see
# that they fall inside Wisconsin. 
# Alternatively, you could check if your sites are within your study area polygon. 
plot( sites )
# Are they within your study region? Mistakes are common when transcribing x, y locations

###################################################################################
######### Import predictor data ###################################################
##################################################################################
# Several datasets are available at the national level and need to be clipped to #
# the relevant study region. In particular, raster files can take a lot of space. #
# Our study area is Wisconsin so lets get a shapefile of the State first. #

# As you probably know, there are multiple ways of achieving the same outcome in R. #
# Sometimes it's a matter of what you know, but sometimes different approaches have #
# different advantages. Importing a shape file with rgdal also imports its projection #
# You can also import shapefiles using maptools. However, that won't import projections. #
# Import Wisconsin shapefile including its projection using rgdal #
WI <- rgdal::readOGR( paste( workdir, "/data/WI_state_outline.shp", sep="" ) )
# Note that we specify which package the function comes from. This is often good practice # 
# when using functions with common names like select(), which are used by multiple packages. #

# View file details:
summary( WI )
# Quick plot:
plot( WI )
# Check if projection matches site location data:
proj4string( WI )
# Convert projection so that it does match:
WI <- spTransform( WI, proj4string( sites ) )
# Note, a common mistake when working with spatial data is forgetting to match their projections.#
# This can introduce significant errors in your analysis. It is particularly common when you #
# are sourcing your data from multiple locations. #

# Our dataset contains presence absence records but in the coming weeks we will compare #
# analyses using presence/absence against those using presence/background data only. #
# We thus generate random background random points but keep them separate from our original #
# dataset #
# We will generate the same number of background points as presence observations
# Number of presence observations for hermit thrush:
M <- table( spp_df$heth )[2]
# Generate background points within Wisconsin using a random sampling scheme
bkgrd.pnts <- sp::spsample( WI, n = M, type = 'random' )
# View
bkgrd.pnts 
# Plot results
points( bkgrd.pnts )

# Was this a reasonable number of background points? What sampling scheme would you use for your
# own data? #
# More ideas on how to choose bakcground points can be found in: 
# Jarnevich et al. (2017) Minimizing effects of methodological decisions on interpretation and #
# prediction in species distribution studies: An example with background selection. #
# Ecological Modelling, 363, 48-56.#

######### Manipulating landcover data ##############
## We use landcover data from the National Geospatial Data Asset (NGDA) Land Use Land Cover #
# https://www.mrlc.gov/nlcd01_data.php #
# We downloaded data directly from the website and extracted them into our data folder.
# We are interested in landcover surrounding our sites only. So we import the complete file #
# first and then extract landcover values around our sites (using a buffer). We do this by #
# adapting code from http://mbjoseph.github.io/2014/11/08/nlcd.html #

# STOP ##
# SKIP THIS SECTION OF BRINGING IN NATIONAL LAND COVER DATA (proceed to the import of land cover data for WI)
# This is only for importing land cover for the whole United States

# Define file location and name:
filename <- paste( workdir, #this is our working directory
                   "/data/nlcd_2001_landcover_2011ed/", #this is the folder path for our file
                  "nlcd_2001_landcover_2011_edition_2014_10_10.img", #this is the actual file name
                  sep="" ) #this pastes those 3 components together with no spaces
# Import NLCD (landcover) data file using the raster package:
NLCD <- raster::raster( filename )
# Plot:
plot( NLCD )

# Check if projection matches survey location data:
proj4string( NLCD )

# Convert shapefile of WI to landcover projection so that it does match:
WI_proj <- spTransform( WI, proj4string( NLCD ) )

## Crop raster file to Wisconsin only:
NLCD_WI <- raster::crop( NLCD, WI_proj ) 

#Plot land cover for Wisconsin
plot(NLCD_WI)

#Write out
writeRaster( NLCD_WI, filename = paste( workdir, "/data/landcover", sep = "" ), 
             format = "HFA" )

# IMPORT LAND COVER FOR WISCONSIN ONLY

# Read in land cover data just for Wisconsin
# Define file location and name:
filename2 <- paste( workdir, #this is our working directory
                   "/data/", #this is the folder path for our file
                   "landcover.img", #this is the actual file name
                   sep="" ) #this pastes those 3 components together with no spaces
# Import NLCD (landcover) data file using the raster package:
NLCD_WI <- raster::raster( filename2 )

# Convert locations to landcover projection so that it does match:
sites_transf <- spTransform( sites, proj4string( NLCD_WI ) )
# Note, we create a second object where site locations are converted to the landcover projection. #
# Can you think of why we did it this way, rather than converting landcover data to the site projection? #
# Do the same for your backgrond points
bkgrd_transf <- spTransform( bkgrd.pnts, proj4string( NLCD_WI ) )
# Define buffer size in meters
buf <- 200

# Extract landcover data for each site, given buffer size #
### NOTE THIS FUNCTION WILL TAKE AT LEAST 15 MINUTES TO RUN ###
Landcover <- raster::extract( NLCD_WI, sites_transf, buffer = buf )
# We have to repeat it also for the random points:
LC_bkgrd <-  raster::extract( NLCD_WI, bkgrd_transf, buffer = buf )

# Create a function that summarizes the proportion of each cover type for each data point:
# the function requires the landcover raster as input and the site ids
summarize_landcover <- function( lcraster, siteid ){
  # Summarize landcover as a proportions at each point
  summ <- lapply( lcraster, function( x ){
          prop.table( table( x ) )
        } )

  # Generate landcover number to name conversions:
  num.codes <- unique( unlist( lcraster ) )
  cover.names <- NLCD_WI@data@attributes[[1]]$Land.Cover.Class[ num.codes + 1 ]
  levels( cover.names )[1] <- NA # first level is ""
  conversions <- data.frame( num.codes, cover.names )
  conversions <- na.omit( conversions )
  conversions <- conversions[ order( conversions$num.codes ), ]
  # Why do we need to do this section? #

  # Summarize results in a data frame:
  lcdf <- data.frame( id = rep( siteid, lapply( summ, length ) ),
                   cover = names( unlist( summ ) ),
                   percent = unlist( summ ) )
  # Create cover name column:
  lcdf$cover2 <- lcdf$cover
  levels( lcdf$cover2 ) <- conversions$cover.names #ensure to match factor levels
  return ( lcdf )
} #end of function

# Run function for our site locations:
NLCD_df <- summarize_landcover( Landcover, spp_df$id )
# Run function for background data points:
bkgrd_df <- summarize_landcover( lcraster = LC_bkgrd, siteid = 1:M )
##delete any temp files created by raster:
raster::removeTmpFiles()

# View first rows:
head( NLCD_df )
head( bkgrd_df )
#Check cover levels found in our sites:
levels( NLCD_df$cover2 ); levels( NLCD_df$cover )
levels( bkgrd_df$cover2 ); levels( bkgrd_df$cover )
# The number of landcover categories required is species dependent. Note that #
# the more categories you use, the more data you require. Think whether you can #
# combine some prior to your analysis. For our Hermit thrush, which is a forest-nesting #
# species, we combine some categories as follows: #
# "1", Developed: Developed, Medium Intensity; Developed, High Intensity
# "2", Open: Developed, Open Space; Barren Land; Hay/Pasture
# "3", Herb: Herbaceuous; Emergent Herbaceuous Wetlands

# Sometimes modifying dataframes can be tricky and will require trial and error. So you don't #
# have to run your entire code again, create new, modified dataframes. However, #
# be careful that you are not filling your entire workspace with intermitent, large dataframes. #

# Create new dataframe from NLCD cover dataframe using dplyr and spread your cover categories
cover_df <- NLCD_df %>% dplyr::select( -cover2 ) %>% group_by( id ) %>% #remove cover 2 and group rows by id
            tidyr::spread( cover, value = percent ) #spread cover to wide format
# If you haven't encountered piping ( %>% ), dplyr and tidyr before, I recommend you look them up as a #
# great way to manipulate, summarize and combine multiple tasks #

# Rename cover categories     
covernames <- c( "id", "Water", "DevOpen", "DevLow", "DevMed", "DevHigh", "Barren",
                           "Deciduous", "Evergreen", "Mixed", "Shrub", "Herbaceuous", 
                          "Hay", "Crops", "WoodyWetland", "HerbaceuousWetland"  )
colnames( cover_df ) <- covernames
cover_df <- cover_df %>% #add new columns for combined landcover categories:
            dplyr::mutate( Developed = sum( DevLow, DevMed, DevHigh, na.rm = TRUE ), 
                           Open = sum( DevOpen, Barren, Hay, na.rm = TRUE ),
                           Herb = sum( Herbaceuous, HerbaceuousWetland , na.rm = TRUE ), 
                           Forest = sum( Deciduous, Evergreen, Mixed, na.rm = TRUE ) ) %>% 
            dplyr::select( -Water, -DevLow, -DevMed, -DevHigh, -DevOpen, -Barren, -Hay, -Herbaceuous, 
                           -HerbaceuousWetland ) #remove categories not being used
# Note that we are adding a general forest variable that combines all subcategories, but leaving those in.
# Note that we are removing Water as we don't think water drives the distribution of our species. 

# View
head( cover_df )
# Replace NA values with 0:
cover_df[ is.na(cover_df) ] <- 0
# Check that there are no empty categories
colSums( cover_df )
glimpse( cover_df )

# Repeat process for your background points but keeping the same categories that you chose for
# your data
bkgrdcover_df <- bkgrd_df %>% dplyr::select( -cover2 ) %>% group_by( id ) %>% #remove cover 2 and group rows by id
  tidyr::spread( cover, value = percent ) #spread cover to wide format
colnames( bkgrdcover_df ) <- covernames
bkgrdcover_df <- bkgrdcover_df %>% #add new columns for combined landcover categories:
  dplyr::mutate( Developed = sum( DevLow, DevMed, DevHigh, na.rm = TRUE ), 
                 Open = sum( DevOpen, Barren, Hay, na.rm = TRUE ),
                 Herb = sum( Herbaceuous, HerbaceuousWetland , na.rm = TRUE ), 
                 Forest = sum( Deciduous, Evergreen, Mixed, na.rm = TRUE ) ) %>% 
  dplyr::select( -Water, -DevLow, -DevMed, -DevHigh, -DevOpen, -Barren, -Hay, -Herbaceuous, 
                 -HerbaceuousWetland ) #remove categories not being used
# View
head( bkgrdcover_df )
# Replace NA values with 0:
bkgrdcover_df[ is.na( bkgrdcover_df ) ] <- 0
# Check that there are no empty categories
colSums( bkgrdcover_df )

#### end of landcover data import #####

################## Manipulating climate data from PRISM #################
# We downladed temperature and rainfall data from PRISM website: #
#  http://prism.oregonstate.edu/
# using the PRISMextract.R code. Data were saved in the data folder. #
# Here we use the downloaded data to match it to our site locations. #
# Minimum temperature normals:
# Define file names
mintfiles <- paste( workdir, "/data/PRISM_Jun-Jul_minT_30yrnorm", sep="" )
# Set path of prism data you want to manipulate:
options( prism.path = mintfiles )

# Stack all minimum temperature files extracted (for Jun-Jul):
MinTAll <- prism_stack( ls_prism_data()[ , 1 ] ) ##raster file
# Convert to same coordinate system as our species data:
MinTWI <- raster::projectRaster( MinTAll, crs = proj4string( sites ) )
## Crop raster file to Wisconsin only:
MinTWI <- raster::crop( MinTWI, WI ) 
# Get seasonal mean of minimum temperature normals for the breeding season (Jun-Jul):
MinT <- stackApply( MinTWI, 1:1, fun = mean ) 
# 1:1 just assigns same index to both files so we get the combined seasonal mean
# Label resulting column:
names( MinT ) <- c( "minT" )
# Save the cropped raster so that it can be used in later practicals
writeRaster( MinT, filename = paste( workdir, "/data/MinT.tif", sep = "" ), 
             format = "GTiff" )

# Plot seasonal (Jun-Jul) mean for minimum temperature 30yr normals:
levelplot( MinT )  + latticeExtra::layer( sp.polygons( WI, col = 'white', 
                                                           lwd = 2 ) )
# Remove intermitent, larger rasters
rm( "MinTAll", "MinTWI" )

# Now that we have seasonal data for Wisconsin we need to match them to 
# our site locations #
# Extract value for each site location:
minTvals <- raster::extract( MinT, sites )
# Why did we not use the buffer around each site location this time?
#View 
head( minTvals ); length( minTvals )

# Extract value for each background location:
minTvals.bkgrd <- raster::extract( MinT, bkgrd.pnts )

# Rainfall normals:
# Define file names
rainfiles <- paste( workdir, "/data/PRISM_Jun-Jul_rain_30yrnorm", sep="" )
# Set path of prism data you want to manipulate:
options( prism.path = rainfiles )

# Stack all files extracted (for Jun-Jul):
RainAll <- prism_stack( ls_prism_data()[ , 1 ] ) ##raster file
# Convert to same coordinate system as our species data:
RainWI <- raster::projectRaster( RainAll, crs = proj4string( sites ) )
## Crop raster file to Wisconsin only:
RainWI <- raster::crop( RainWI, WI ) 
# Get seasonal mean of rainfall normals for the breeding season (Jun-Jul):
Rain <- stackApply( RainWI, 1:1, fun = mean ) 
# 1:1 just assigns same index to both files so we get the combined seasonal mean
# Label resulting column:
names( Rain ) <- c( "rain" )
# Save the cropped raster so that it can be used in later practicals
writeRaster( Rain, filename = paste( workdir, "/data/Rain.tif", sep = "" ), 
             format = "GTiff" )

# Plot seasonal (Jun-Jul) mean for minimum temperature 30yr normals:
levelplot( Rain )  + latticeExtra::layer( sp.polygons( WI, col = 'white', 
                                                       lwd = 2 ) )
# Remove intermitent, larger rasters
rm( "RainAll", "RainWI" )

# Now that we have seasonal data for Wisconsin we need to match them to 
# our site locations #
# Extract value for each site location:
rainvals <- raster::extract( Rain, sites )
# View 
head( rainvals ); length( rainvals )
# Extract value for each background location:
rainvals.bkgrd <- raster::extract( Rain, bkgrd.pnts )
# View 
head( rainvals.bkgrd ); length( rainvals.bkgrd )

## Delete any temp files created by raster:
raster::removeTmpFiles()
##### end of weather data import #######

############### combine species and predictor data  ############################
# Join species and landcover dataframes
alldata <- left_join( spp_df, cover_df, by = 'id' ) %>%
            arrange( id ) #make sure is sorter by site id
# Add weather data:
alldata$MinT <- minTvals
alldata$Rain <- rainvals
# View
head( alldata )

##### save as .csv file so that it can be used by another script:
write.csv( x = alldata, file = paste( workdir, "/data/alldata.csv", sep="" ), 
           row.names = FALSE )

# Repeat the process for your background points
bkgrddata <- as.data.frame( bkgrd.pnts@coords[,c("x", "y" ) ] )
bkgrddata$id <- 1:M
bkgrddata <- left_join( bkgrddata, bkgrdcover_df, by = "id" )
# Add weather data:
bkgrddata$MinT <- minTvals.bkgrd
bkgrddata$Rain <- rainvals.bkgrd
# View
head( bkgrddata )
##### save as .csv file so that it can be used by another script:
write.csv( x = bkgrddata, file = paste( workdir, "/data/bkgrddata.csv", sep="" ), 
           row.names = FALSE )

########## end of data combine ##############
############# Checking predictors ##########################
# We start by evaluating the correlation amongst our predictors and also against
# our site locations
library( psych ) # we only load this package now because we are just using it briefly
alldata %>% ungroup() %>% #ensure data is not grouped
  dplyr::select( -id ) %>% #remove columns that we don't want
  pairs.panels() #create correlation plot
detach( "package:psych", unload=TRUE ) # we remove the package from the workspace

# What can you see from this analysis? Are there any variables that are highly correlated?

# We then evaluate the spread and distribution of predictors:
prednames <- c( "Forest", "Deciduous", "Evergreen", "Mixed", "Shrub", "Crops", "WoodyWetland", "Developed",
                "Open", "Herb", "MinT", "Rain" )
par( mfrow = c( 4, 3 ) )
for( n in 1:length( prednames ) ){
  hist( alldata[ , prednames[ n ] ], main = prednames[ n ], xlab = "" )
  print( c( prednames[ n ], min( alldata[ , prednames[ n ] ], na.rm = TRUE ), 
   max( alldata[ , prednames[ n ] ], na.rm = TRUE ) ) )
}

# Any concerns regarding the spread or distribution of the predictors?

# We now want to visualize some of our ecological predictors spatially
ggplot( data = alldata ) + geom_point( aes( x = x, y  = y, color = Forest ) )
ggplot( data = alldata ) + geom_point( aes( x = x, y  = y, color = Deciduous ) )
ggplot( data = alldata ) + geom_point( aes( x = x, y  = y, color = MinT ) )
ggplot( data = alldata ) + geom_point( aes( x = x, y  = y, color = Rain ) )

# How else do you plot spatial data?

#### end of predictor assessment #####

######################## END OF SCRIPT ################################################
########################################################################################
