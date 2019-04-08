####################################################################################
###### This code was created by Jen Cruz as part of the species distribution ######
#####   modeling class, co-taught with Ben Zuckerberg.                       #######
##### Practical 3: Fitting species distribution models for presence/absence   ######
######            data using the glm() and the dismo package                  ######
###################################################################################

########## clean workspace and load required packages ####################
###########################################################################
#####clean workspace to improve efficiency: ###
rm(list = ls() ) 
#set working dir() 
gc() #releases memory

####### load relevant packages ###
library( ggplot2 ) #fancy plots
library( dismo ) #Species distribution analysis package.
library( randomForest ) #Random forest analysis.
########## end of package loading ###########

###########################################################################################
########    loading data and defining required objects           ##########################
###########################################################################################
# Rather than reloading our datasets, this time we load the workspace we saved in prac2 #
# what are the advantages/disadvantages of doing that? #

load( "D:/RProjects/SDM/prac2results.RData" )

###################    functions #############################################
# R allows us to create functions for processes that we want to repeat on multiple objects #

# We will use a logistic model for one of our analyses below. Here we create a function #
# that allows us to go from the logit to the real scale. How else could we do this in R? #
expit <- function( p ){
  return( exp( p ) / ( exp( p ) + 1 ) )
}

# For your own analysis, can you think of any functions that can help you streamline your #
# work? 

############
###########################################################################################
########    data analysis           ##########################
###########################################################################################

# In prac2 we explored the sdm package for analysing the distributions of Hermit thrush in Wisconsin. #
# Although it claimed to be very powerful, we found that some of the options are yet #
# to be implemented. Here we go back to more established packages that allow sdm analyses. #
# We focus on our presence-absence data only, but you can adjust to your data needs accordingly #
# for the homework exercise. Pick only one analysis for that. #

##### Start of GLM anlysis #########
# Generalized linear models are easily fit using the glm() function that comes loaded in R. #
# We want to compare results from glm() with those obtained last week in sdm() using the glm method. #
# Last week we were missing predictor coefficients. These are easily obtained in glm(). #
# Last week we incorporated all predictors because the package claimed to be able to deal with collinearity. #
# Here, we deal with collinearity by comparing two models (so highly correlated predictors are not included #
# in the same model). This is a more standard approach. #

# Our first glm excludes Forest and includes other forest-related categories:
glm1 <- glm( heth ~ Deciduous + Evergreen + Mixed  + Shrub + Crops + 
               WoodyWetland + Developed + Open + Herb + MinT + Rain, family = binomial,
             data = train.padata )
# View output
summary( glm1 )
# Now we can actually see coefficient estimates for our predictors #
# Which predictors are important? How do we assess this?
# The anova function can provide some additional clues:
anova( glm1 )
# What does it tell us?

# Our second glm replaces detailed forest categories with an overall forest value. #
glm2 <- glm( heth ~ Forest + Shrub + Crops + 
               WoodyWetland + Developed + Open + Herb + MinT + Rain, family = binomial,
             data = train.padata )
# View output
summary( glm2 )
anova( glm2 )
# How do the two models compare? #
# Which variables appear important? How do the coefficients vary between models?

# Which model should we use? ###

# We can compare our models using AIC. When is AIC unsuitable for model comparisons? #
# If you want to learn more about model evaluation and selection techniques read: #
# Hooten, M.B. & Hobbs, N.T. (2015) A guide to Bayesian model selection for ecologists. #
# Ecological Monographs, 85, 3-28. #

# Which model would you select based on AIC? Why?

# Model evaluation:
# We can use dismo to assess misclassification rates, against our test data, similarly to last week #
# For model 1:
evglm1 <- dismo::evaluate( p = test.padata[ which( test.padata$heth == 1 ), covnames ], #predictors for presence records
                 a = test.padata[ which( test.padata$heth == 0 ), covnames ], #predictors for absence records
                 model = glm1 )
# For model 2:
evglm2 <- dismo::evaluate( p = test.padata[ which( test.padata$heth == 1 ), covnames ], #predictors for presence records
                           a = test.padata[ which( test.padata$heth == 0 ), covnames ], #predictors for absence records
                           model = glm2 )

# Compare model results:
par( mfrow = c( 3, 2 ) )
plot( evglm1, "ROC" )
plot( evglm2, "ROC" )
boxplot( evglm1 )
boxplot( evglm2 )
density( evglm1 )
density( evglm2 )

# Can these results help you discriminate which model to choose? #
# What else could you use to help your decision? #

# What about our biological understanding? Do the relationships between our response and #
# important predictors make sense? #

# We can plot response curves to get a better idea of the relationships with predictors. #
par( mfrow = c( 3, 3 ) )
termplot( glm1 ) #remeber to hit return after running this line
termplot( glm2 )#remeber to hit return after running this line
# What do the y axes represent?#

# Based on the pieces of information collected above, let's choose a model as a top choice #

# We probably want a better idea of how our predictors relate to the probability of presence. #
# We can create our own response curves. We only do this for predictors of interest #
# Here we demonstrate using Forest and MinT using glm2. #
# For this we need the range of predictor values we want to predict our responses over. #
# Remember our dataframe containing our unscaled predictors? This was loaded when we re-loaded the workspace:
head( alldata )

# We use it to get the range of values for our predictors of interest #
# We work with 100 values:
predrows <- 100
forvals <- seq( min( alldata[ ,"Forest" ] ), max( alldata[ ,"Forest" ] ), 
                       length.out = predrows  )
mintvals <- seq( min( alldata[ ,"MinT" ] ), max( alldata[ ,"MinT" ] ), 
                length.out = predrows  )

# we want to use the appropriate coefficients to create our partial relationships:
coefficients( glm2 )
int.glm2 <- coefficients( glm2 )[1]
forcoef <- coefficients( glm2 )[2]
mintcoef <- coefficients( glm2 )[9]

#We then calculate partial responses while assuming all other predictors are kept at their #
# mean values i.e zero:
forresp <- expit( int.glm2 + forcoef * scale( forvals ) )
mintresp <-expit( int.glm2 + mintcoef * scale( mintvals ) )
# Plot partial response curves against actual predictor values (not scaled) so they are easy to #
# interpret: #
par( mfrow = c( 2,1 ) )
plot( forvals, forresp, xlab = "Forest", ylab = "Prob of presence" )
plot( mintvals, mintresp, xlab = "Minimum, temperature", ylab = "Prob of presence" )

# What do these response plot tell us in how hermit thrush responds to those predictors?

#### end of glm analysis #####

#### Start of random forest analysis #####
# Random forest is a machine learning option, which is gaining recent popularity in the sdm field. #
# Here, we replicate our randomforest analysis using the randomForest and dismo packages. #
# Machine learning approaches use regularization and are not as sensitive to collinearity so #
# we are able to include all predictors and let the model decide which are important. #

# We start by defining our model including all predictors:
rfmod <- heth ~ Forest + Deciduous + Evergreen + Mixed  + Shrub + Crops + 
          WoodyWetland + Developed + Open + Herb + MinT + Rain
# Then we run random forest analysis:
rfm1 <- randomForest( rfmod, data = train.padata )
# Variable importance can be assessed through:
importance( rfm1 )
# Which variables are being selected as important? Do they match the biology of the species?

# We can more easily plot the response curves for important predictors:
par( mfrow = c( 2, 1 ) )
randomForest::partialPlot( rfm1, pred.data = train.padata, x.var = Forest )
randomForest::partialPlot( rfm1, pred.data = train.padata, x.var = MinT )
# Are these predictors the same as those chosen by the analysis last week? #
# How do these relationships compare to those chosen by the top GLM model above?

# We can evaluate misclassication rates for randomForest also using dismo:
ev.rfm1 <- dismo::evaluate( p = test.padata[ which( test.padata$heth == 1 ), covnames ], #predictors for presence records
                            a = test.padata[ which( test.padata$heth == 0 ), covnames ], #predictors for absence records
                            model = rfm1 )
# Plot model evaluation results:
par( mfrow = c( 3, 1 ) )
plot( ev.rfm1, "ROC" )
boxplot( ev.rfm1 )
density( ev.rfm1 )

# Again what are some of the issues in relying on AUC as the only approach to model evaluation? #
# What else should we be trying? #

################################################################################################
############### Predicting our species distribution base on model results #####################
###############################################################################################

# We can evaluate the output of our model spatially in two  ways. First, by assessing how our #
# model estimates relate to our site observations (using our training data). Second, by assessing # 
# how our model predics to new locations. We could predict the probability of presence to the #
# entire state of Wisconsin by creating rasters of our predictors at our site resolution and extent #
# so 200 m grid cells to the state of Wisconsin. We would then standardize those (scale) in the same #
# way that we did with our site predictors (in prac1) so that we can apply our model to these new input predictors. #
# However, remember that extracting landcover data for our <10,000 sites took >25min to run.#
# NLCD comes at a 30 x 30 m resolution, even if you crop it to Wisconsin we are still talking #
# about several million cells. Remember this when attempting to do this for your own work. #
# You may find faster options outside of R. Have any of you used those? #

# So instead, we evaluate how our model predicts to new locations using our testing data. #
# for which we already extracted predictor values (in prac1). # These locations also have #
# actual observations of Hermit thrush, so we can actually assess how well our model does #
# at predicting those. #

# We use our top glm model as an example and take advantage of the predict() function #
# We need to include the same predictors as we used in the model:
glm2names <- covnames[ 4:length( covnames ) ]
# We can then use the predict.glm to obtain our estimated values for each site used in the model #
# and append to our training dataset
train.padata$psi.glm2 <- predict.glm( glm2, newdata = train.padata[ ,glm2names ], type = "response" ) 
# View
head( train.padata )
# We do the same to test our model against the new sites (from test dataset):
test.padata$psi.glm2 <- predict.glm( glm2, newdata = test.padata[ ,glm2names ], type = "response" ) 

# Note that dismo also has a predict function that you can adapt if you do random forest #
# instead. The sdm package also has a predict function if you want to go that route also. #

# We can then create plots that display our results spatially #
# for the estimated probabilities of presence at our training sites: 
p1 <- ggplot( data = train.padata ) + geom_point( aes( x = x, y  = y, color = psi.glm2 ), 
        size = 2, alpha = 0.6  ) + 
        scale_color_gradient( low = "skyblue1", high = "orangered4", name = "Presence probability" )
# You can find more color options here: http://sape.inf.usi.ch/quick-reference/ggplot2/colour

# For the actual observations of presence at those sites:
p2 <- ggplot( data = train.padata ) + geom_point( aes( x = x, y  = y, color = as.factor( heth ) ), 
        size = 2, alpha = 0.6 ) + 
        scale_color_manual( values = c( "skyblue1", "orangered4" ), name = "Observed presence" )
# For the predicted probabilities of presence at your new testing sites:
p3 <- ggplot( data = test.padata ) + geom_point( aes( x = x, y  = y, color = psi.glm2 ), 
      size = 2, alpha = 0.6 ) + 
      scale_color_gradient( low = "skyblue1", high = "orangered4", name = "Presence probability" )
# For your observations at the testing sites:
p4 <- ggplot( data = test.padata ) + geom_point( aes( x = x, y  = y, color = as.factor( heth ) ), 
      size = 2, alpha = 0.6 ) + 
      scale_color_manual( values = c( "skyblue1", "orangered4" ), name = "Observed presence" )

# An easy way to plot multiple ggplots that come from different datasets we can use the #
# multiplot function below
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

# Then run it to plot your graphs in the same screen: 
multiplot( cols = 2, p1, p3, p2, p4 )

# What do these results tell you? How well did the model estimate the data we used to run it? #
# How well did it predict in new areas? # 
# Is the distribution what you expect for Hermit thrush in Wisconsin? #

######################## END OF SCRIPT ################################################
########################################################################################