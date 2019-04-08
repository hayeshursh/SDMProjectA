##############################################################################
###   Code created by Jen Cruz to download PRISM weather data        ########
###               http://prism.oregonstate.edu/                      ########
#### We download 30yr normals (1981-2010) at 800 m resolution for    #########
#### the summer breeding season (Jun-Jul)                            ########
##############################################################################

########### clean workspace and load required packages ####################
###########################################################################

#####clean workspace to improve efficiency: ###
rm(list = ls() ) 
gc() #releases memory
# Note that you only need to install each package once ###
install.packages( "prism" ) 
####### load relevant packages ###
library( prism ) #to download PRISM data

######### end ######################################################

######################################################################
####               Download PRISM data                      ##########
##### 800 m grid resolution for the lower 48 states of USA    ########
##########
# Set path to store data. In our case we use the data subfolder #
# withing our project directory, but you should adjust this to #
# wherever you want to store your data #
datadir <- paste( getwd(), #points to working project directory
                  "/data/", #points to data subfolder 
                  sep="" ) #joins path with no spaces

######### Downloading minimum temperature normals:
# Define file names
mintfiles <- paste( datadir, "PRISM_Jun-Jul_minT_30yrnorm", sep="" )
# Set as prism path:
options( prism.path = mintfiles )
# Download data
get_prism_normals( type = 'tmin', resolution = "800m", mon = 6:7 )

######### Downloading rainfall normals:
# Define file names
rainfiles <- paste( datadir, "PRISM_Jun-Jul_rain_30yrnorm", sep="" )
# Set as prism path:
options( prism.path = rainfiles )
# Download data
get_prism_normals( type = 'ppt', resolution = "800m", mon = 6:7 )

#look at downloaded prism files product name:
head( ls_prism_data( name = TRUE ) )
#look at downloaded prism file absolute path:
head( ls_prism_data( absPath = TRUE ),30 )
####end data download #####

######################## END OF SCRIPT ################################################
########################################################################################
