####################################################################################
###### This code was created by Jen Cruz as part of the species distribution ######
#####   modeling class, co-taught with Ben Zuckerberg.                       #######
##### Practical 4: Fitting species distribution models for presence/absence   ######
######  data with a glm and an occupancy model in JAGs (Bayesian)            ######
###################################################################################

########## clean workspace and load required packages ####################
###########################################################################
#####clean workspace to improve efficiency: ###
rm(list = ls() ) 
#set working dir() 
gc() #releases memory

# We are going to be learning how to use JAGs to run a glm and an occupancy #
# model under a Bayesian framework. For this, you need to first downloand and #
# install a copy of JAGs on your computer from here: #
# https://sourceforge.net/projects/mcmc-jags/
# You'll also need jagsUI which allows R to interact with JAGs #
install.packages( "jagsUI" ) 

# To get details on how occupancy models work in a Bayesian framework, you can read: #
# Royle, J.A. & Kéry, M. (2007) A bayesian state-space formulation of dynamic occupancy #
# models. Ecology, 88, 1813-1823. #

# For variations on the basic framework including multi-season, multi-state and dynamic #
# occupancy models I recommend you get Andrew Royle's or Marc Kery's books. They #
# provide all you need to teach yourself how to run these models. #
# Briefly, if you are wondering why bother with Bayesian models: they allow you #
# to easily include hierarchical structures (e.g. an observation submodel, random effects), #
# missing data (those times that you couldn't get to your sites) and low sample sizes (those #
# working with endangered species may relate). Moving from binary to poisson to normal structures #
# is also really easy so you can go from estimating occupancy to abundance with very small #
# changes to your code. #

####### load relevant packages ###
#library( tidyr ) #combines dataframes
library( dplyr ) #manipulates dataframes
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( tidyr ) #spread and other dataframe functions
library( ggplot2 ) #fancy plots
#library( lubridate ) #easy date adjustments and calculations.
#library( rgdal ) #imports projects transforms spatial data. Requires sp.
# library( rgeos ) #spatial manipulations e.g. buffer, interpolate.
#library( maptools ) #reads and manipulates ESRI shapefiles.
#library( raster ) #manipulates rasters
#library( rasterVis ) #visualises rasters
#library( sdm ) #Species distribution analysis package.
library( jagsUI ) #To run JAGs
########## end of package loading ###########

############## add functions ##########################
expit <- function( p ){
  return( exp( p ) / ( exp( p ) + 1 ) )
}
logit <- function( x ){
    log( x / ( 1 - x ) ) }

############
########## multiple plots function###############
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#############

###########################################################################################
########    loading data and defining required objects           ##########################
###########################################################################################
# Set working directory
workdir <- getwd()

# We will be using our presence/absence training and testing datasets that we created #
# in prac 2 #
#Import from csv files:
train.padata <- read.csv( file = paste( workdir, "/data/train_padata.csv", sep="" ) )
# View
head( train.padata )
test.padata <- read.csv( file = paste( workdir, "/data/test_padata.csv", sep="" ) )
# View
head( test.padata )

# We also import our complete dataset with original predictor values (not standardized) #
# which we will need later to help us interpret our results. #
alldata <- read.csv( file = paste( workdir, "/data/alldata.csv", sep="" ) )
# View
head( alldata )

# The presence/absence observations for hermit thrushes that we have been using so far #
# were actually derived by pooling the results from three separate repeat surveys at #
# each site. These replicate data actually allows us to estimate how good we were at #
# finding our species (probability of detection). A failure to observe a species  #
# does not necessarily mean that the species was not present. Imperfect detection #
# can bias our estimates of species distributions. Here, we will compare how our #
# distribution changes when we estimate it using a model that ignores imperfect #
# detection (glm) and one that accounts for it (occupancy model). #

# We import our detailed data file containing observations for each of three #
# replicate surveys as well as the date of survey (reported as julian date) #
# and hour of survey. These last two predictors can affect our detection. #
detdata <- read.csv( file = paste( workdir, "/data/detdata.csv", sep="" ) )
# View
head( detdata )

# We had split our original species dataset in two (training and testing sets) #
# we now have to match our detection records to each dataset:
train.padata <- inner_join( train.padata, detdata, by = "id" )
head( train.padata )
test.padata <- inner_join( test.padata, detdata, by = "id" )
head( test.padata )

#### We use our training data for analyses #

# We define which predictors we want in our ecological model:
econames <- c( "Shrub", "Crops", "WoodyWetland", 
               "Developed", "Open", "Herb", "Forest", "MinT", "Rain" )  
# We then extract those columns from our dataframe:
ecocovs <- as.matrix( train.padata[ , econames ]  )
# View
head( ecocovs )
# We extract predictors for the observation model as a matrices #
# First julian dates for each survey:
julmat <- as.matrix( train.padata[ , c( "jul1", "jul2", "jul3" ) ] )
# Scale each column:
julmat <- apply( julmat, MARGIN = 2, scale )
# View
head( julmat )
# Then for hour of survey:
hrmat <- as.matrix( train.padata[ , c( "hr1", "hr2", "hr3" ) ] )
# Scale each column:
hrmat <- apply( hrmat, MARGIN = 2, scale )
# View
head( hrmat )

# For our glm we need a vector containing our observations across the three surveys #
# i.e. the summary observations that we have been using in previous classes #
ys <- train.padata$heth 

# For occupancy we need a matrix containing our observations:
y_obs <- train.padata %>% dplyr::select( yobs1:yobs3 )
#convert to matrix:
y_obs <- as.matrix( y_obs )
#remove column names from matrix
colnames( y_obs ) <- NULL
# View
head( y_obs )
# if we want to see differences in observations amongst surveys:
apply( y_obs, MARGIN = 2, table )

# We also need some summary values
# Number of sites
M <- dim( train.padata )[1]
# Number of ecological predictors:
xno <- dim( ecocovs )[2]
# Number of replicate surveys
J <- 3

# JAGs relies on Multi-Chain MonteCarlo (MCMC) algorithms to update model parameter #
# estimates #
# So we also need to define how long we will run our models, how many chains, and how #
# many records we will keep. #

ni <- 20000; nt <- 5; nb <- 5000; nc <- 3 #iterations, thinning, burnin, chains


####### end #############
###################################################################################################
########                     Analyzing our data                   #################################
###################################################################################################

###################     Generalized linear model              ###################
############## Specfy model in bugs language:  #####################
sink( "glm1.txt" )
cat( "
     model{
     
      #priors
      #define intercepts as mean prob:
      int.psi <- logit( mean.psi )
      mean.psi ~ dbeta( 2, 2 ) #mean occupancy probability
      # coefficients for predictors
      for( n in 1:xno ){ 
        beta[ n ] ~ dnorm( 0, 0.05 ) #using strong priors
      }
     
      #ecological model:
      for ( i in 1:M ){  #loop over sites
        ys[ i ] ~ dbern( psi[ i ] ) #latent, true occupancy 
        #model for occupancy probability:
        logit( psi[ i ] ) <- int.psi + 
                          #matrix multiplication of coefficients with predictors
                          inprod( beta[1:xno ], ecocovs[ i, 1:xno ] ) 
      } # close i loop
     
     } #model close
     
     ", fill = TRUE )

sink()
################ end of model specification  #####################################
modelname <- "glm1.txt"

# Provide initial values
inits <- function(){ list( beta = rnorm( n = xno) ) }
#parameters monitored
params <- c( "mean.psi" #mean prob of occupancy
             , "beta" #coefficients in ecological model
             , "psi" #estimated probability of occupancy at each site
)

# Define data to input into model:
str( win.data <- list( ys = ys, M = M, J = J , 
                       ecocovs = ecocovs, xno = xno
) )

#call JAGS to run the model:
GLM1 <- jags( data = win.data, inits = inits, 
              params, model.file = modelname, 
              n.chains = nc, n.thin = nt, n.iter = ni, 
              n.burnin = nb, parallel = TRUE ) 

#### end of glm ####

##################################################################################################
###################            Occupancy model              ###################
############## Specfy model in bugs language:  #####################
sink( "om1.txt" )
cat( "
     model{
     
      #priors
      #for ecological model:    
      #define intercepts as mean prob:
      int.psi <- logit( mean.psi )
      mean.psi ~ dbeta( 2, 2 ) #mean occupancy probability
      # coefficients for predictors
      for( n in 1:xno ){ 
        beta[ n ] ~ dnorm( 0, 0.05 ) #using strong priors
      }

      #priors for detection models:
      #for occupancy
      int.p <- logit( mean.p ) #intercept
      mean.p ~ dbeta( 2, 2 ) #mean detection probability
      #coefficients for detection predictors
      for( n in 1:3 ){ 
        alpha[ n ] ~ dnorm( 0, 0.05 ) #using strong priors
      }
     
      #ecological model:
      for ( i in 1:M ){  #loop over sites
        z[ i ] ~ dbern( psi[ i ] ) #latent, true occupancy 
        #model for occupancy probability:
        logit( psi[ i ] ) <- int.psi + 
                          #matrix multiplication of coefficients with predictors
                          inprod( beta[1:xno ], ecocovs[ i, 1:xno ] ) 
      } # close i loop

      #observation model:
      for ( i in 1:M ) { #loop over sites 
        for( j in 1:J ) { #loop over replicate surveys
          #relationship with date of survey modelled as quadratic relationship
          #and with hour as linear
          logit( p[ i, j ] ) <- int.p + alpha[ 1 ] * julmat[ i, j ] +
                          alpha[ 2 ] * pow( julmat[ i, j ], 2 ) + 
                          alpha[ 3 ] * hrmat[ i, j ]
          #observations depend on occupancy and detection
          y_obs[ i, j ] ~ dbern( z[ i ] * p[ i, j ] ) 
        } #close j loop
      } #close i loop

   } #model close
     
     ", fill = TRUE )

sink()
################ end of model specification  #####################################
modelname <- "om1.txt"

#create initial values
inits <- function(){ list( z = ys ) }

#parameters monitored
params <- c( "mean.psi" #mean prob of occupancy
             ,"mean.p" #mean prob of detection
             , "beta" #coefficients in ecological model
             , "alpha" #coefficients in detection model
             , "z" # estimated occupancy at each site
             , "psi" #estimated probability of occupancy at each site
            )

# Define data to input into model:
str( win.data <- list( y_obs = y_obs, M = M, J = J , 
                       ecocovs = ecocovs, xno = xno,
                       julmat = julmat, hrmat = hrmat
                      ) )

#call JAGS to run the model:
moc1 <- jags( win.data, inits = inits, params, modelname, 
            n.chains = nc, n.thin = nt, n.iter = ni, 
            n.burnin = nb, parallel = TRUE ) 

######### end of occupancy model ######

# Save your results:
save.image( paste( workdir, "/prac4results.RData", sep="" ) )
#####################################################################################
########################### PLOTTING RESULTS #########################################
### General summaries ####
summary( GLM1 )

# Evaluating convergence:
par( mfrow = c( 2, 2 ) )
traceplot( GLM1, parameters = c( 'beta') )
traceplot( moc1, parameters = c( 'beta') )

# Parameter estimates:
par( mfrow = c( 2, 2 ) )
whiskerplot( GLM1, parameters = c( "beta" ) )
whiskerplot( moc1, parameters = c( "beta" ) )
whiskerplot( GLM1, parameters = c( "mean.psi" ) )
whiskerplot( moc1, parameters = c( "mean.psi" ) )
whiskerplot( moc1, parameters = c( "alpha" ) )
whiskerplot( moc1, parameters = c( "mean.p" ) )

### end #####
# rename your model a generic name so you don't have to change the code between models
mr <- moc1

###################################################################################
################# plot predicted relationships ###########
################################################################################
############ Ecological partial relationships #########
# Define parameter rangers for important predictors:
econames
# Forest
forest <- seq( min( alldata[ , 'Forest'], na.rm = TRUE ), 
               max( alldata[ , 'Forest'], na.rm = TRUE ), 
               length.out = 100 )
sclfor <- scale( forest )
# Minimum temperature
mint <- seq( min( alldata[ , 'MinT' ], na.rm = TRUE ), 
                max( alldata[ , 'MinT' ], na.rm = TRUE ), 
                length.out = 100 )
sclmint <- scale( mint )
# Rain
rain <- seq( min( alldata[ , 'Rain'], na.rm = TRUE ), max( alldata[ , 'Rain'], na.rm = TRUE ), 
                 length.out = 100 )
sclrain <- scale( rain )
# Intercept
int <- rep( 1, 100 )

# Estimate mean relationship with predictor while keeping all others constant at their mean
#multiply coefficient matrix (for all iterations) against transposed predictor matrix using matrix algebra
psi.for <- cbind( logit( mr$sims.list$mean.psi ), mr$sims.list$beta[ , 7 ] ) %*% 
            t( cbind( int, sclfor ) ) 
#calculate its mean and get inverse logit
psi.form <- expit( apply( psi.for, MARGIN = 2, FUN = mean ) )
#calculate its 95% CIs and get its inverse logit:
psi.forCI <- expit( apply( psi.for, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975) ) )
#get a matrix (nrows= no of iteractions from mcmc, ncols = number of values to predict over )
psi.mint <- cbind( logit( mr$sims.list$mean.psi ), mr$sims.list$beta[ , 8 ] ) %*% 
            t( cbind( int, sclmint ) )
#calculate its mean and get inverse logit
psi.mintm <- expit( apply( psi.mint, MARGIN = 2, FUN = mean ) )
#calculate its 95% CIs and get its inverse logit:
psi.mintCI <- expit( apply( psi.mint, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975) ) )

psi.rain <- cbind( logit( mr$sims.list$mean.psi ), mr$sims.list$beta[ , 9 ] ) %*% 
            t( cbind( int, sclrain ) )
#calculate its mean and get inverse logit
psi.rainm <- expit( apply( psi.rain, MARGIN = 2, FUN = mean ) )
#calculate its 95% CIs and get its inverse logit:
psi.rainCI <- expit( apply( psi.rain, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975) ) )

#combine predicted estimates into a dataframe
ecorships <- data.frame( sclfor, forest, sclmint, mint, sclrain, rain, 
                          psi.form, t( psi.forCI ) , 
                          psi.mintm, t( psi.mintCI ), 
                          psi.rainm, t( psi.rainCI )
)
#view
head( ecorships )

#plot using ggplot
ecop <- ggplot( data = ecorships ) + theme_classic() +
        theme( legend.position = "none", 
         text = element_text( size = 18 ), 
         axis.line = element_line( size = 1.3 ) ) + #, 
        ylab( "Probability of ocupancy" )

#par( mfrow = c( 3, 1 ), ask = F , mar = c(3,4,2,2) )
forp <- ecop +  xlab( "Forest (%)" ) +
  geom_line( aes( x = forest, y = psi.form ), size = 1.2 ) +
  geom_ribbon( alpha = 0.4, aes( x = forest, ymin = X2.5. , ymax = X97.5. ) ) 

mintp <- ecop + xlab( "Minimum temperature (breeding season)" ) +
  geom_line( aes( x = mint, y = psi.mintm ), size = 1.2 ) +
  geom_ribbon( alpha = 0.4, aes( x = mint, ymin = X2.5..1, ymax = X97.5..1 ) )

rainp <- ecop + xlab( "Total rain (breeding season)" ) +
  geom_line( aes( x = rain, y = psi.rainm ), size = 1.2 ) +
  geom_ribbon( alpha = 0.4, aes( x = rain, ymin = X2.5..2, ymax = X97.5..2 ) ) 

multiplot( forp, mintp, rainp, cols = 2 )

############ end of partial ecological relationships ####

########## partial detection relationships ###########
# Julian day of survey
jul <- seq( min( detdata[ , 'jul1' ], na.rm = TRUE ), max( detdata[ , 'jul1' ], na.rm = TRUE ), 
             length.out = 100 )
scljul <- scale( jul )
# Hour of survey
hr <- seq( min( detdata[ , 'hr1' ], na.rm = TRUE ), 
           max( detdata[ , 'hr1' ], na.rm = TRUE ), 
           length.out = 100 )
sclhr <- scale( hr )

# Estimate mean relationship with predictor while keeping all others constant at their mean
#multiply coefficient matrix (for all iterations) against transposed predictor matrix using matrix algebra
#get a matrix (nrows= no of iteractions from mcmc, ncols = number of values to predict over )
p.jul <- cbind( logit( mr$sims.list$mean.p ), mr$sims.list$alpha[ , 1 ], 
                mr$sims.list$alpha[ , 2 ] ) %*% t( cbind( int, scljul, scljul^2 ) )
#calculate its mean and get inverse logit
p.julm <- expit( apply( p.jul, MARGIN = 2, FUN = mean ) )
#calculate its 95% CIs and get its inverse logit:
p.julCI <- expit( apply( p.jul, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975) ) )

# Now for hour
p.hr <- cbind( logit( mr$sims.list$mean.p ), mr$sims.list$alpha[ , 3 ] ) %*% 
        t( cbind( int, sclhr ) ) 
#calculate its mean and get inverse logit
p.hrm <- expit( apply( p.hr, MARGIN = 2, FUN = mean ) )
#calculate its 95% CIs and get its inverse logit:
p.hrCI <- expit( apply( p.hr, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975) ) )

#combine estimates into a dataframe
detrships <- data.frame( sclhr, hr, scljul, jul, 
                         p.hrm, t( p.hrCI ) , 
                         p.julm, t( p.julCI )
)
#view
head( detrships )

#plot using ggplot
detp <- ggplot( data = detrships ) + theme_classic() +
  theme( legend.position = "none", 
         text = element_text( size = 18 ), 
         axis.line = element_line( size = 1.3 ) ) + #, 
  ylab( "Probability of detection" )

#par( mfrow = c( 3, 1 ), ask = F , mar = c(3,4,2,2) )
hrp <- detp +  xlab( "Hour of survey" ) +
  geom_line( aes( x = hr, y = p.hrm ), size = 1.2 ) +
  geom_ribbon( alpha = 0.4, aes( x = hr, ymin = X2.5. , ymax = X97.5. ) ) 

julp <- detp + xlab( "Julian day of survey" ) +
  geom_line( aes( x = jul, y = p.julm ), size = 1.2 ) +
  geom_ribbon( alpha = 0.4, aes( x = jul, ymin = X2.5..1, ymax = X97.5..1 ) )

multiplot( hrp, julp, cols = 2 )

############ end of partial detection relationships ####

################################################################################
################## Create prediction maps across sites #################

# We want to compare the estimated probability of occupancy (psi) and estimated #
#occupancy (z) against our summed observations of presence at each site. # 

###### First for our training set (within sample validation) #####
# We estimate psi and z in JAGs and keep them as parameters to trace so that we can #
# easily plot them spatially #

# To plot them in ggplot we combine results into a dataframe. Here we focus on mean values#
# but from the partial plots we also know that the 95% CIs are also easily extracted #
# and plotted 
resdf <- data.frame( mr$mean$psi, train.padata$heth, train.padata$x, 
                     train.padata$y )
colnames( resdf ) <- c( "psi", "naive.z", "x", "y" )
# We can then create plots that display our results spatially #
# for the estimated probabilities of presence at our training sites: 
psip <- ggplot( data = resdf ) + geom_point( aes( x = x, y  = y, color = psi ), 
                                                  size = 2, alpha = 0.6  ) + 
  scale_color_gradient( low = "skyblue1", high = "orangered4", name = "Occupancy probability" )
# You can find more color options here: http://sape.inf.usi.ch/quick-reference/ggplot2/colour

# For the actual observations of presence at those sites:
naivezp <- ggplot( data = resdf ) + geom_point( aes( x = x, y  = y, color = as.factor( naive.z ) ), 
                                                  size = 2, alpha = 0.6 ) + 
  scale_color_manual( values =c( "skyblue1", "orangered4" ), name = "Naive occupancy" )

multiplot( psip, naivezp, cols = 2 )

##### FOR OUR OCCUPANCY MODEL ONLY:
resdf$z <- round( mr$mean$z )
head( resdf )
# We also have estimates of occupancy at each site:
zp <- ggplot( data = resdf ) + geom_point( aes( x = x, y  = y, color = as.factor( z ) ), 
                                           size = 2, alpha = 0.6 ) + 
  scale_color_manual( values = c( "skyblue1", "orangered4" ), name = "Estimated occupancy" )

# We can  plot the differences between the estimated and naive occupancies to get #
# a clearer picture of where our model improved inference:
difp <- ggplot( data = resdf ) + geom_point( aes( x = x, y  = y, 
        color = as.factor( z - naive.z ) ), size = 2, alpha = 0.6 ) + 
  scale_color_manual( values = c( "grey80", "red3" ), name = "Differences" )

multiplot( psip, zp, naivezp, difp, cols = 2 )

##### end of within sample plots #####

#### For our testing dataset (out-of-sample prediction) ######

# We now want to predict outside of our sites. Here we use the testing sites #
# but we can also use grid cells for our study area #
# Psi and z are site specific, so to predict to new sites we need to do it outside of #
# JAGs. we demonstrate how to estimate psi. You can also estimate z but that #
# requires further assumptions (threshold values at which you classify psi as z=1 ) #

# To predict psi, we use a similar approach to our partial predicted relationships #
# except that we use all ecological predictors in the model, rather than focusing on a single one #
head( test.padata )

# store transposed intercept and parameter values:
parvals <- t( cbind( rep( 1,  dim( test.padata )[1] ), test.padata[, econames ] ) )
# Multiply by mean coefficient estimates (but could use sims.list instead if you want #
# CIs same as we did with the partial prediction plots) #
test.psi <- expit( c( logit( mr$mean$mean.p ), mr$mean$beta ) %*% 
              parvals )

# Create dataframe for ggplot
ggdata <- test.padata
ggdata$psi <- t( test.psi )
head( ggdata )
# Plot
ggplot( data =  ggdata ) + geom_point( aes( x = x, y  = y, 
              color = psi ), size = 2 ) + 
            scale_color_gradient( low = "skyblue1", high = "orangered4", 
              name = "Occupancy probability" )

###############     end of map plots ##################################

######################## END OF SCRIPT ################################################
########################################################################################