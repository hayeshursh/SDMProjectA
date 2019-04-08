####################################################################################
###### This code was created by Jen Cruz as part of the species distribution ######
#####   modeling class, co-taught with Ben Zuckerberg.                       #######
##### Practical 2: Fitting species distribution models for presence/absence   ######
####      and presence/background datasets using the sdm package             ######
###################################################################################

########## clean workspace and load required packages ####################
###########################################################################
#####clean workspace to improve efficiency: ###
rm(list = ls() ) 
#set working dir() 
gc() #releases memory

##### install relevant packages ####
# For our analysis we will use the sdm package
install.packages( "sdm" ) 

####### load relevant packages ###
#library( tidyr ) #combines dataframes
library( dplyr ) #manipulates dataframes
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( tidyr ) #spread and other dataframe functions
#library( ggplot2 ) #fancy plots
#library( lubridate ) #easy date adjustments and calculations.
#library( rgdal ) #imports projects transforms spatial data. Requires sp.
# library( rgeos ) #spatial manipulations e.g. buffer, interpolate.
#library( maptools ) #reads and manipulates ESRI shapefiles.
#library( raster ) #manipulates rasters
#library( rasterVis ) #visualises rasters
library( sdm ) #Species distribution analysis package.
# sdm relies on other packages for analysis. Install those associated packages:
installAll( )

########## end of package loading ###########

###########################################################################################
########    loading data and defining required objects           ##########################
###########################################################################################
# Set working directory
workdir <- getwd()
# Import .csv file containing combined species and predictor data:
alldata <- read.csv( file = paste( workdir, "/data/alldata.csv", sep="" ) )
# View imported data
head( alldata )
glimpse( alldata )
# Import .csv file containing background points
bkgrddata <- read.csv( file = paste( workdir, "/data/bkgrddata.csv", sep="" ) )
# View
head( bkgrddata )

# We want to see any potential benefits of using presence/absence data over just #
# presence records. To do this we will analyze presence/absence data (alldata) and #
# presence/background data. #
# For this second analysis we need to combine our presence records with our #
# background points. #
# For homework with your own data, just choose the analysis that matches your data type.#

# We start by adding a column to the background points that highlight the species was #
# not observed there. #
# Add response column to background points:
bkgrddata$heth <- rep( 0, dim( bkgrddata )[1] )
# We then join presence observations with background points. # 
# First we find common columns for both dataframes:
common_cols <- intersect( colnames( bkgrddata ), colnames( alldata ) )
# Then join them based on those common columns (this removes the column for the second bird species)
PBdata <- rbind(
          alldata[ which( alldata$heth == 1 ), common_cols ], #Select presence only and common cols
          bkgrddata[ , common_cols ] #select common columns
          )
# View
head( PBdata ); tail( PBdata )

### Setting general vectors: ####
# Total number of sites sampled in presence/absence data:
M <- max( alldata$id )
# in presence background data
N <- dim( PBdata )[1]

# We also want to split our data into our training and testing sets: # 
# Define number of data points to use for testing:
TP <- 1000
# Select which rows to use for testing for our presence/absence data:
t.parows <- sample( x = 1:M, size = TP, replace = FALSE )
# Select which rows to use for testing for our presence/background data:
t.pbrows <- sample( x = 1:N, size = TP, replace = FALSE )

###########################################################################################
######## Preparing required response and covariate data for analyses  #####################
###########################################################################################
# Predictor data is commonly standardized prior to analysis so that effect sizes are directly #
# comparable. We want to also keep the original values to aid interpretation of results. #

#For the presence/absence dataframe:
# Create a dataframe to contain standardize covariates:
pa.data <- alldata %>% dplyr::select( -eame ) #remove data for Eastern Meadowlark
# View
head( pa.data )
# Define predictor columns that require standardizing:
covcols <- which( !names( pa.data ) %in% c("id", "x", "y", "heth" ) )
# Define predictor names
covnames <- colnames( pa.data )[covcols]
# Check that you selected the right columns
head( pa.data[ , covcols ] )
# Standardize each column
pa.data[ , covcols ] <- apply( pa.data[ , covcols ], MARGIN = 2, scale )
# Note that scale standardizes continuous predictors. How do we standardize categorical ones?
# or binomial ones?

# Plot predictor distributions:
par( mfrow = c( 4, 3 ) )
for (n in covcols ){
  hist( pa.data[ ,n ], breaks = 8, main = colnames( pa.data )[n] )
}
# Check correlations amongst covariates
cor( pa.data[ , covcols ] )
# Create training dataset
train.padata <- pa.data[ -t.parows, ]
# View
head( train.padata ); dim( train.padata )
# Create testing dataset:
test.padata <- pa.data[ t.parows, ]
# View
head( test.padata ); dim( test.padata )

##### save our dataframes as .csv files so that they can be used by other scripts:
write.csv( x = train.padata, file = paste( workdir, "/data/train_padata.csv", sep="" ), 
           row.names = FALSE )
write.csv( x = test.padata, file = paste( workdir, "/data/test_padata.csv", sep="" ), 
           row.names = FALSE )

#For presence/background dataframe:
# Select predictor columns
pb.data <- PBdata 
# View
head( pb.data ); dim( pb.data )
# Define predictor columns that require standardizing:
bk.covcols <- which( !names( pb.data ) %in% c("id", "x", "y", "heth" ) )
# Define predictor names
bk.covnames <- colnames( pb.data )[ bk.covcols ]
# Check that you selected the right columns
head( pb.data[ , bk.covcols ] )
# Standardize each column
pb.data[ , bk.covcols ] <- apply( pb.data[ , bk.covcols ], MARGIN = 2, scale )
# Plot predictor distributions:
par( mfrow = c( 4, 3 ) )
for (n in bk.covcols ){
  hist( pb.data[ ,n ], breaks = 8, main = colnames( pb.data )[n] )
}
# Check correlations amongst covariates
cor( pb.data[ , bk.covcols ] )
# Create training dataset
train.pbdata <- pb.data[ -t.pbrows, ]
# View
head( train.pbdata ); dim( train.pbdata )
# Create testing dataset:
test.pbdata <- pb.data[ t.pbrows, ]
# View
head( test.pbdata ); dim( test.pbdata )

##### save our dataframes as .csv files so that they can be used by other scripts:
write.csv( x = train.pbdata, file = paste( workdir, "/data/train_pbdata.csv", sep="" ), 
           row.names = FALSE )
write.csv( x = test.pbdata, file = paste( workdir, "/data/test_pbdata.csv", sep="" ), 
           row.names = FALSE )

####### end #############
###################################################################################################
########                     Analyzing our data                   #################################
###################################################################################################

########## Presence/absence data analysis ####################
# We analyze our data using the sdm package. For details on what this package can do read:
# Naimi, B. & Araújo, M.B. (2016) sdm: a reproducible and extensible R platform for species #
# distribution modelling. Ecography, 39, 368-375. #

# To analyze data using sdm we first need to get data ready for analysis using sdmData(). # 
# By using the formula option we can specify which predictors we #
# want to include, and can incorporate multiple species as responses (on the left-hand side of the equation) #
# We can specify site locations, with coords(x + y), and the dataframe that contains our data. #
# We could actually provide rasters as our predictors, and implement several approaches for #
# estimating a testing dataset (including Cross-Validation). #
# In our example, we do those outside of the package, can you think of reasons why?

# define our data, response and predictors:
hethd <- sdmData( formula = heth ~ Deciduous + Evergreen + Mixed + Forest + Shrub + Crops + 
                    WoodyWetland + Developed + Open + Herb + MinT + Rain + coords(x + y), 
                  train = train.padata, test = test.padata )
# View data details:
hethd

# Now that we've inputed the data we are ready to run the analysis. The sdm package allows the #
# concurrent use of multiple methods. Make sure that you select methods that are suitable for #
# your data type. In this example with presence/absence data we choose three common methods: #
# glm: generalized linear model, rf: random forest, and brt: boosted regression trees. 
m1 <- sdm( heth ~., data = hethd, methods = c( "glm", "rf", "brt" ) )
# Since we already defined our predictors in our data object, we can use . to tell the function #
# to include them all. #
# Note that by doing this we are including covariates that are highly correlated #
# The authors claim the sdm package can use a variance inflation factor (VIF) to measures how much a #
# predictor is explained by alternative predictors. If it doesn't provide additional information it is removed #
# The authors claim this can be done in sdm using a stepwise procedure #

# What are alternative approaches to deal with collinearity in our predictors? #

# View results
m1
getModelInfo( m1 )
# What does these results tell us? #

# What about our predictors? Which are important in modeling the Hermit thrush distribution? #
# The sdm package can estimates which variables were important for each method using our training data:
vi1.1 <- getVarImp( m1, id = 1, wtest = "training" )
vi1.2 <- getVarImp( m1, id = 2, wtest = "training" )
vi1.3 <- getVarImp( m1, id = 3, wtest = "training" )
# Plot results
par( mfrow = c( 3, 1 ) )
plot( vi1.1, "auc", main = "glm" )
plot( vi1.2, "auc", main = "rf" )
plot( vi1.3, "auc", main = "brt" )
# According to these results, which variables are influencing model results?

# Model evaluation:
# Can our model discriminate whether a species is present at a site or not?
# We assess this by estimating, the plot Receiver Operating Characteristic (ROC) Curve (AUC):
roc( m1)
# What is AUC a measure of? What do the x and y axes represent? 
# You may want to refresh your memory by re-reading:
# Jiménez-Valverde, A. (2012) Insights into the area under the receiver operating #
# characteristic curve (AUC) as a discrimination measure in species distribution modelling. #
# Global Ecology and Biogeography, 21, 498-507. #

# How do the AUC values compare to other measures of misclassification? sdm package also allows us to 
# evaluate the rate of misclassification using other techniques. 
# Here we calculate: #
# the true skill statistic (TSS), the Receiver Operating Characteristic Curve (AUC), #
# the sensitivity (proportion of presences correctly predicted as such), #
# and specificity (proportion of absences correctly predicted as such) #
getEvaluation( m1, stat = c( "TSS", "AUC",  "Sensitivity", "Specificity" ), opt = "max(se+sp)" )
# Note that model performance was evaluated on our independent, testing data #
# For threshold-based procedures we can select the criteria for optimising the threshold #
# here we chose the option to maximise the true negative and positive rates #
# Based on these results, which model would you choose as better? #

###### end ##########
########## Presence/background data analysis ####################

# We start our presence/background analysis using the sdm package again. #
# sdm has internal methods for developing background data but since we already have some, #
# we provide background records directly. We also provide our test data but note that the package can also derive that #
# using a variety of methods including Cross-Validation:
hethbd <- sdmData( formula = heth ~ Deciduous + Evergreen + Mixed + Forest + Shrub + Crops + 
                     WoodyWetland + Developed + Open + Herb + MinT + Rain + coords(x+y), 
                    train = train.pbdata[ which( train.pbdata$heth == 1), ], 
                    bg = train.pbdata[ which( train.pbdata$heth == 0), ],
                   test = test.pbdata )
hethbd
# Run model
bm1 <- sdm( heth ~., data = hethbd, methods = c( "glm", "rf", "brt" ) )
# View results:
bm1
getModelInfo( bm1 )

# Determining significant predictors: 
vib1.1 <- getVarImp( bm1, id = 1, wtest = "training" )
vib1.2 <- getVarImp( bm1, id = 2, wtest = "training" )
vib1.3 <- getVarImp( bm1, id = 3, wtest = "training" )
# Plot results
par( mfrow = c( 3, 1 ) )
plot( vib1.1, "auc", main = "glm" )
plot( vib1.2, "auc", main = "rf" )
plot( vib1.3, "auc", main = "brt" )
# According to these results, which variables are influencing model results?
# Did both types of analyses (presence/absence vs presence/background select similar predictors? #

# Model evaluation:
# Evaluating missclassification rates of the model:
# Plot Receivor Operating Characteristic (ROC) Curve (AUC):
roc( bm1 )
# Is this ROC appropriate for presence/background data? 
# What should the x axes represent instead? 

getEvaluation( bm1, stat = c( "TSS", "Sensitivity", "Specificity" ), opt = "max(se+sp)" )
# Based on these results, which method would you use?
# How do those results compare to the presence/absence analysis:
getEvaluation( m1, stat = c( "TSS", "Sensitivity", "Specificity" ), opt = "max(se+sp)" )

# Which one do you trust more? Are our measures of missclassification enough to guide our decision? #
# What else can help us decide? #
# What other model evaluation procedures are we missing? What other things could we have tested for? #

# We end by saving our workspace so that we can continue working with our results in future prac classes:
save.image( "D:/RProjects/SDM/prac2results.RData" )

################## end of analyses ##################################

######################## END OF SCRIPT ################################################
########################################################################################