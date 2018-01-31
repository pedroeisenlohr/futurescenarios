#####################################################
#### NICHE MODELS FORECASTED FOR FUTURE SCENARIOS ###
#####################################################


### By: Pedro V. Eisenlohr - UNEMAT, Alta Floresta, MT, Brazil ###
### Contact: pedro.eisenlohr@unemat.br ###########################


### Acknowledgments ###
## Guarino Colli (Universidade de Brasília) 
## Diogo Bezerra Rocha (Instituto de Pesquisas Jardim Botânico do Rio de Janeiro)
## My students of LabEc (Ecology Lab, Alta Floresta-MT, Brazil)



###########################
## SET WORKING DIRECTORY ##
###########################

# Each user should adjust this!
setwd("C:/Models")
getwd()


### Whenever it's necessary:
# Parallel processing
cl <- makeCluster(detectCores()) # number of cores in computer
registerDoParallel(cl)
getDoParWorkers()

#Enhance memory limit:
memory.limit(1000000) #or another suitable memory size (in KB)




###############################
## INSTALL AND LOAD THE REQUIRED PACKAGES ##
###############################

install.packages("biomod2", dep=T)
install.packages("car", dep=T)
install.packages("colorRamps", dep=T)
install.packages("corrplot", dep=T)
install.packages("dismo", dep=T)
install.packages("doParallel", dep=T)
install.packages("dplyr", dep=T)
install.packages("maps", dep=T)
install.packages("maptools", dep=T)
install.packages("plotKML", dep=T)
install.packages("PresenceAbsence", dep=T)
install.packages("raster", dep=T)
install.packages("rgdal", dep=T)
install.packages("sdm", dep=T)
install.packages("SDMTools", dep=T)
install.packages("sqldf", dep=T)
install.packages("testthat", dep=T)
install.packages("usdm", dep=T)
install.packages("FactoMineR", dep=T)

library(biomod2)
library(car)
library(colorRamps)
library(corrplot)
library(dismo)
library(doParallel)
library(dplyr)
library(maps)
library(maptools)
library(plotKML)
library(PresenceAbsence)
library(raster)
library(rgdal)
library(sdm)
library(SDMTools)
library(sqldf)
library(testthat)
library(usdm)
library(FactoMineR)
library(foreach)


####################################
###### IMPORTING DATA ##############
####################################

# Current climate
bioclim <- raster::getData("worldclim", var = "bio", res = 5)
bio <-stack(bioclim)

# Biotic data
spp<-read.table(file.choose(),row.names=1,header=T,sep=",")
dim(spp)
edit(spp)

# Remove duplicates
dups2 <- duplicated(spp[, 1:2])
sum(dups2) # how many duplicates?
spp <- spp[!dups2, 1:2] # remove duplicates
dim(spp)

# Selecting spatially unique records #
mask <- bio[[1]]
cell <- cellFromXY(mask, spp[,1:2]) # get the cell number for each point
cell
dup <- duplicated(cbind(spp[,1:2],cell))
spp <- spp[!dup, ]# select the records that are not duplicated
dim(spp)


#######################################
## FIRST STEP OF COLLINEARITY CHECK ##
#######################################

presvals <- extract(bio, spp)
dim(presvals)
edit(presvals)

# PCA:
pca <- PCA(presvals,graph=FALSE)
plot(pca, choix="var")

# Variance Inflation Factors (VIF) of environmental variables
v1 <- vifcor(presvals, th=0.8) #correlation
v1
### Check if no variable presents VIF>10


#######################################
## SECOND STEP OF COLLINEARITY CHECK ##
#######################################

## Crop WorldClim layers
#  *********************

# Visualize boundaries
data(wrld_simpl)
plot(wrld_simpl, xlim = c(-100, -30), ylim = c(-60, 25), axes = TRUE, col = "light grey")

# Set boundaries
ext <- extent(-100, -30, -60, 25) #Define here the extention of your projection area.

# Crop layers (present)
bio.crop1 <- crop(bio, ext)
bio.crop1
plot(bio.crop1)

# Create a new object excluding the collinear variables:
bio.crop2 <- exclude(bio.crop1,v1)
bio.crop2
bio.crop2 <-stack(bio.crop2)

# Alternatively, exclude variables of your choice:
#bio.crop2 <- dropLayer(bio.crop1,c(1,2,3,4,5,6,8,9,10,12,13,14,17,18,19)) ###Include only variables that should be discarded.
#bio.crop2 #Check with caution!

# Select 10000 random points from mask
mask <- bio.crop2$bio7 ###Any variable present in bio.crop2
rnd.points <- randomPoints(mask, 10000)
plot(!is.na(mask), legend = F)
points(rnd.points, cex = 0.5)

# Principal Components Analysis (PCA) of environmental variables
env.data <- extract(bio.crop2, rnd.points)
pca.env.data <- princomp(env.data, cor = T)
biplot(pca.env.data, pc.biplot = T)

# Variance Inflation Factors (VIF) of environmental variables
v1 <- vifcor(bio.crop2, th=0.8) #correlation
v1
### Check if no variable presents VIF>10

# Subset environmental stack
env.selected <- exclude(bio.crop2, v1) #exclude collinear variables identified with vifcor 
env.selected
plot(env.selected)

# Name variables (only the selected variables!)
names(env.selected) <- c("bio7", "bio11", "bio15", "bio16")
names(env.selected)



################################################
## GENERATING OTHER REQUIRED OBJECTS FOR SDM ###
################################################

# Convert dataset to SpatialPointsDataFrame (only presences)
myRespXY <- spp[,c("long","lat")]
# Creating occurrence data object
occurrence.resp <-  rep(1, length(myRespXY$long))



############################################
## FIT SPECIES DISTRIBUTION MODELS - SDMS ##
############################################

########## Prepare data ##########

### For CTA, GBM and RF algorithms:
### for example, if you have n=100 ocurrence records
sppBiomodData.PA.equal <- BIOMOD_FormatingData(
	resp.var = occurrence.resp,
	expl.var = env.selected,
	resp.xy = myRespXY,
	resp.name = "Occurrence",
	PA.nb.rep = 10,
	PA.nb.absences = 100,
	PA.strategy = "sre", #Read manual with caution
	PA.sre.quant = 0.025)
sppBiomodData.PA.equal

# For the other 7 algorithms:
#Don't modify anything here, except, if you wish, PA.strategy and/or PA.sre.quant:
sppBiomodData.PA.10000 <- BIOMOD_FormatingData(
	resp.var = occurrence.resp,
	expl.var = env.selected,
	resp.xy = myRespXY,
	resp.name = "Occurrence",
	PA.nb.rep = 10,
	PA.nb.absences = 10000,
	PA.strategy = "sre",
	PA.sre.quant = 0.025)
sppBiomodData.PA.10000



# MaxEnt .jar
  jar <- paste0(system.file(package = "dismo"), "/java/maxent.jar")
  if (file.exists(jar) != T) {
    url = "http://biodiversityinformatics.amnh.org/open_source/maxent/maxent.php?op=download"
    download.file(url, dest = "maxent.zip", mode = "wb")
    unzip("maxent.zip", files = "maxent.jar", exdir = system.file("java", package = "dismo"))
    unlink("maxent.zip")
    warning("Maxent into directory")
  } 
system.file("java", package = "dismo")

# Here, define the directory where MaxEnt has been installed, exactly as defined above
myBiomodOption <- BIOMOD_ModelingOptions(MAXENT.Phillips = list(path_to_maxent.jar="C:/Documents/R/win-library/3.3/dismo/java"))




#################
### Modeling ####
#################

sppModelOut.PA.equal <- BIOMOD_Modeling(sppBiomodData.PA.equal, 
	models = c("GBM", "CTA", "RF"), 
	models.options = myBiomodOption, 
	NbRunEval = 10,
	DataSplit = 70, #percentage of data used to train
	Prevalence = NULL, 
	VarImport = 0,
	models.eval.meth = c("TSS","ROC","SR","POD","ACCURACY","BIAS"),
	SaveObj = TRUE,
	rescal.all.models = FALSE,
	do.full.models = FALSE,
	modeling.id = "spp_current")
sppModelOut.PA.equal


sppModelOut.PA.10000 <- BIOMOD_Modeling( 
	sppBiomodData.PA.10000, 
	models = c("GLM","GAM","ANN","SRE","FDA","MARS","MAXENT.Phillips"), 
	models.options = myBiomodOption, 
	NbRunEval = 10,
	DataSplit = 70, #percentage of data used to train
	Prevalence = NULL, 
	VarImport = 0,
	models.eval.meth = c("TSS","ROC","SR","POD","ACCURACY","BIAS"),
	SaveObj = TRUE,
	rescal.all.models = FALSE,
	do.full.models = FALSE,
	modeling.id = "spp_current")
sppModelOut.PA.10000



###################################
## EVALUATE MODELS USING BIOMOD2 ##
###################################

# Get evaluations
sppModelEval.PA.equal <- get_evaluations(sppModelOut.PA.equal)
sppModelEval.PA.equal
write.table(sppModelEval.PA.equal, "EvaluationsAll_1.csv")
sppModelEval.PA.10000 <- get_evaluations(sppModelOut.PA.10000)
sppModelEval.PA.10000
write.table(sppModelEval.PA.10000, "EvaluationsAll_2.csv")


# Get summaries (mean and std.dev.) of model evaluation - 1
sdm.models1 <- c("GBM","CTA","RF") #3 models
sdm.models1
eval.methods1 <- c("TSS","ROC","SR","POD","ACCURACY","BIAS") #6 evaluation methods
eval.methods1

means.i <- numeric(0)
means.j <- numeric(6)
for (i in 1:3){
	for (j in 1:6){
	means.j[j] <- mean(sppModelEval.PA.equal[paste(eval.methods1[j]),"Testing.data",paste(sdm.models1[i]),,])
	}
	means.i <- c(means.i, means.j)
}

summary.eval.equal <- data.frame(rep(sdm.models1,each=6), rep(eval.methods1,3), means.i)
names(summary.eval.equal) <- c("Model", "Method", "Mean")
summary.eval.equal
write.table(summary.eval.equal,"Models1_Evaluation_Mean.csv")

sd.i <- numeric(0)
sd.j <- numeric(6)
for (i in 1:3){
	for (j in 1:6){
	sd.j[j] <- sd(sppModelEval.PA.equal[paste(eval.methods1[j]),"Testing.data",paste(sdm.models1[i]),,])
	}
	sd.i <- c(sd.i, sd.j)
}

summary.eval.equal <- data.frame(rep(sdm.models1,each=6), rep(eval.methods1,3), sd.i)
names(summary.eval.equal) <- c("Model", "Method", "SD")
summary.eval.equal
write.table(summary.eval.equal,"Models1_Evaluation_SD.csv")


# Get summaries (mean and std.dev.) of model evaluation - 2
sdm.models2 <- c("GLM","GAM","ANN","SRE","MARS","MAXENT.Phillips","FDA") #7 models
sdm.models2
eval.methods2 <- c("TSS","ROC","SR","POD","ACCURACY","BIAS") #6 evaluation methods
eval.methods2

means.i <- numeric(0)
means.j <- numeric(6)
for (i in 1:7){
	for (j in 1:6){
	means.j[j] <- mean(sppModelEval.PA.10000[paste(eval.methods2[j]),"Testing.data",paste(sdm.models2[i]),,], na.rm=T)
	}
	means.i <- c(means.i, means.j)
}

summary.eval.10000 <- data.frame(rep(sdm.models2,each=6), rep(eval.methods2,7), means.i)
names(summary.eval.10000) <- c("Model", "Method", "Mean")
summary.eval.10000
write.table(summary.eval.10000,"Models2_Evaluation_Mean.csv")

sd.i <- numeric(0)
sd.j <- numeric(6)
for (i in 1:7){
	for (j in 1:6){
	sd.j[j] <- sd(sppModelEval.PA.10000[paste(eval.methods2[j]),"Testing.data",paste(sdm.models2[i]),,])
	}
	sd.i <- c(sd.i, sd.j)
}

summary.eval.10000 <- data.frame(rep(sdm.models2,each=6), rep(eval.methods2,7), sd.i)
names(summary.eval.10000) <- c("Model", "Method", "SD")
summary.eval.10000
write.table(summary.eval.10000,"Models2_Evaluation_SD.csv")




###############################################################
############ Which algorithms have been retained? #############
###############################################################
### Answering such a question is very important from now on ###
###############################################################
### When all algorithms are called, please choose only ######
############## those you have selected ########################
###############################################################



###############################
## PRODUCE MODEL PROJECTIONS ##
###############################

spp.projections_1 <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env.selected,
	proj.name = "Cur1_current",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections_2 <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env.selected,
	proj.name = "Cur2_current",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")


# Stack projections
### Define the directory where proj_Cur1_Occurrence.grd is:
projections_1 <-stack("C:/Models/Current/proj_Cur1_current/proj_Cur1_current_Occurrence.grd")
names(projections_1)

### Define the directory where proj_Cur2_Occurrence.grd is:
projections_2 <-stack("C:/Models/Current/proj_Cur2_current/proj_Cur2_current_Occurrence.grd")
names(projections_2)


### Mean models for each algorithm:
projections.RF.all <- subset(projections_1, grep("RF", names(projections_1)))
projections.RF.mean <- mean(projections.RF.all)/10

projections.GBM.all <-subset(projections_1, grep("GBM", names(projections_1)))
projections.GBM.mean <- mean(projections.GBM.all)/10

projections.CTA.all <-subset(projections_1,grep("CTA", names(projections_1)))
projections.CTA.mean <- mean(projections.CTA.all)/10

projections.GLM.all <-subset(projections_2,grep("GLM", names(projections_2)))
projections.GLM.mean <- mean(projections.GLM.all)/10

projections.GAM.all <-subset(projections_2,grep("GAM", names(projections_2)))
projections.GAM.mean <- mean(projections.GAM.all)/10

projections.ANN.all <- subset(projections_2,grep("ANN", names(projections_2)))
projections.ANN.mean <- mean(projections.ANN.all)/10

projections.SRE.all <- subset(projections_2,grep("SRE", names(projections_2)))
projections.SRE.mean <- mean(projections.SRE.all)/10

projections.MARS.all <- subset(projections_2,grep("MARS", names(projections_2)))
projections.MARS.mean <- mean(projections.MARS.all)/10

projections.FDA.all <- subset(projections_2,grep("FDA", names(projections_2)))
projections.FDA.mean <- mean(projections.FDA.all)/10

projections.MAXENT.all <- subset(projections_2,grep("MAXENT.Phillips", names(projections_2)))
projections.MAXENT.mean <- mean(projections.MAXENT.all)/10


##############################################################
######## AVERAGE ENSEMBLE AMONG SELECTED ALGORITHMS ##########
##############################################################

projections.all.mean <- mean(projections.RF.mean + projections.GBM.mean +
	projections.CTA.mean + projections.GLM.mean +
	projections.MARS.mean + projections.MAXENT.mean)
windows(w=6, h=6)
plot(projections.all.mean, col = matlab.like(100), main = "Ensemble - Current Climate", las = 1)
plot(wrld_simpl, add = TRUE, col="transparent", border="white", lwd = 0.5)
writeRaster(projections.all.mean, filename="Ensemble - Current Climate.tif", format="GTiff")
writeRaster(projections.all.mean, filename="Ensemble - Current Climate.asc", format="ascii")
plot(projections.all.mean)


#Converting to binary maps:
projections.binary <- BinaryTransformation(projections.all.mean, th) #Replace 'th' by the threshold value you wish to apply.
class(projections.binary)
summary(values(projections.binary))
writeRaster(projections.binary, filename="Current Climate_BINARY.tif", format="GTiff")
writeRaster(projections.binary, filename="Current Climate_BINARY.asc", formato="ascii")
plot(projections.binary)




## Get WorldClim layers for future (with rcp 45 and 85 to compare) - GCMs: AC, BC, CC, CN, GS, HD, HG, HE, IN, IP, MI, MR, MC, MP, MG, NO
#  ************************************************************************************************

        
##########################
###### 2050 ##############
##########################

bio50.45_AC <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 45, model = "AC", year = 50)
bio50.85_AC <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 85, model = "AC", year = 50)

bio50.45_BC <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 45, model = "BC", year = 50)
bio50.85_BC <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 85, model = "BC", year = 50)

bio50.45_CC <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 45, model = "CC", year = 50)
bio50.85_CC <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 85, model = "CC", year = 50)

bio50.45_CN <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 45, model = "CN", year = 50)
bio50.85_CN <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 85, model = "CN", year = 50)

bio50.45_GS <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 45, model = "GS", year = 50)
bio50.85_GS <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 85, model = "GS", year = 50)

bio50.45_HD <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 45, model = "HD", year = 50)
bio50.85_HD <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 85, model = "HD", year = 50)

bio50.45_HG <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 45, model = "HG", year = 50)
bio50.85_HG <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 85, model = "HG", year = 50)

bio50.45_HE <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 45, model = "HE", year = 50)
bio50.85_HE <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 85, model = "HE", year = 50)

bio50.45_IN <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 45, model = "IN", year = 50)
bio50.85_IN <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 85, model = "IN", year = 50)

bio50.45_IP <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 45, model = "IP", year = 50)
bio50.85_IP <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 85, model = "IP", year = 50)

bio50.45_MI <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 45, model = "MI", year = 50)
bio50.85_MI <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 85, model = "MI", year = 50)

bio50.45_MR <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 45, model = "MR", year = 50)
bio50.85_MR <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 85, model = "MR", year = 50)

bio50.45_MC <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 45, model = "MC", year = 50)
bio50.85_MC <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 85, model = "MC", year = 50)

bio50.45_MP <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 45, model = "CC", year = 50)
bio50.85_MP <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 85, model = "CC", year = 50)

bio50.45_MG <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 45, model = "CC", year = 50)
bio50.85_MG <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 85, model = "CC", year = 50)

bio50.45_NO <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 45, model = "CC", year = 50)
bio50.85_NO <- raster::getData("CMIP5", var = "bio", res = 5, rcp = 85, model = "CC", year = 50)


# Crop WorldClim layers
bio50.45.crop_AC <- crop(bio50.45_AC, ext)
bio50.85.crop_AC <- crop(bio50.85_AC, ext)

bio50.45.crop_BC <- crop(bio50.45_BC, ext)
bio50.85.crop_BC <- crop(bio50.85_BC, ext)

bio50.45.crop_CC <- crop(bio50.45_CC, ext)
bio50.85.crop_CC <- crop(bio50.85_CC, ext)

bio50.45.crop_CN <- crop(bio50.45_CN, ext)
bio50.85.crop_CN <- crop(bio50.85_CN, ext)

bio50.45.crop_GS <- crop(bio50.45_GS, ext)
bio50.85.crop_GS <- crop(bio50.85_GS, ext)

bio50.45.crop_HD <- crop(bio50.45_HD, ext)
bio50.85.crop_HD <- crop(bio50.85_HD, ext)

bio50.45.crop_HG <- crop(bio50.45_HG, ext)
bio50.85.crop_HG <- crop(bio50.85_HG, ext)

bio50.45.crop_HE <- crop(bio50.45_HE, ext)
bio50.85.crop_HE <- crop(bio50.85_HE, ext)

bio50.45.crop_IN <- crop(bio50.45_IN, ext)
bio50.85.crop_IN <- crop(bio50.85_IN, ext)

bio50.45.crop_IP <- crop(bio50.45_IP, ext)
bio50.85.crop_IP <- crop(bio50.85_IP, ext)

bio50.45.crop_MI <- crop(bio50.45_MI, ext)
bio50.85.crop_MI <- crop(bio50.85_MI, ext)

bio50.45.crop_MR <- crop(bio50.45_MR, ext)
bio50.85.crop_MR <- crop(bio50.85_MR, ext)

bio50.45.crop_MC <- crop(bio50.45_MC, ext)
bio50.85.crop_MC <- crop(bio50.85_MC, ext)

bio50.45.crop_MP <- crop(bio50.45_MP, ext)
bio50.85.crop_MP <- crop(bio50.85_MP, ext)

bio50.45.crop_MG <- crop(bio50.45_MG, ext)
bio50.85.crop_MG <- crop(bio50.85_MG, ext)

bio50.45.crop_NO <- crop(bio50.45_NO, ext)
bio50.85.crop_NO <- crop(bio50.85_NO, ext)


# Stack environmental layers
environment50.45_AC <- stack(bio50.45.crop_AC)
environment50.85_AC <- stack(bio50.85.crop_AC)

environment50.45_BC <- stack(bio50.45.crop_BC)
environment50.85_BC <- stack(bio50.85.crop_BC)
  
environment50.45_CC <- stack(bio50.45.crop_CC)
environment50.85_CC <- stack(bio50.85.crop_CC)

environment50.45_CN <- stack(bio50.45.crop_CN)
environment50.85_CN <- stack(bio50.85.crop_CN)

environment50.45_GS <- stack(bio50.45.crop_GS)
environment50.85_GS <- stack(bio50.85.crop_GS)

environment50.45_HD <- stack(bio50.45.crop_HD)
environment50.85_HD <- stack(bio50.85.crop_HD)

environment50.45_HG <- stack(bio50.45.crop_HG)
environment50.85_HG <- stack(bio50.85.crop_HG)

environment50.45_HE <- stack(bio50.45.crop_HE)
environment50.85_HE <- stack(bio50.85.crop_HE)

environment50.45_IN <- stack(bio50.45.crop_IN)
environment50.85_IN <- stack(bio50.85.crop_IN)

environment50.45_IP <- stack(bio50.45.crop_IP)
environment50.85_IP <- stack(bio50.85.crop_IP)

environment50.45_MI <- stack(bio50.45.crop_MI)
environment50.85_MI <- stack(bio50.85.crop_MI)

environment50.45_MR <- stack(bio50.45.crop_MR)
environment50.85_MR <- stack(bio50.85.crop_MR)

environment50.45_MC <- stack(bio50.45.crop_MC)
environment50.85_MC <- stack(bio50.85.crop_MC)

environment50.45_MP <- stack(bio50.45.crop_MP)
environment50.85_MP <- stack(bio50.85.crop_MP)

environment50.45_MG <- stack(bio50.45.crop_MG)
environment50.85_MG <- stack(bio50.85.crop_MG)

environment50.45_NO <- stack(bio50.45.crop_NO)
environment50.85_NO <- stack(bio50.85.crop_NO)





# Name variables 

#### WARNING!!! #########
#### WARNING!!! #########
### THE SEQUENCE OF THE NAMES BELOW MUST BE EXACTLY THE SAME AS DEFINED ABOVE ####
### SOMETIMES THE SEQUENCE IS "bio1", "bio10" etc. IF THIS IS THE CASE, PLEASE MAKE ADJUSTMENTS BELOW ###########################
#################################################################################################################################


names(environment50.45_AC) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.85_AC) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.45_BC) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.85_BC) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.45_CC) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.85_CC) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.45_CN) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.85_CN) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.45_GS) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.85_GS) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.45_HD) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.85_HD) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.45_HG) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.85_HG) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.45_HE) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.85_HE) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.45_IN) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.85_IN) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.45_IP) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.85_IP) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.45_MI) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.85_MI) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.45_MR) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.85_MR) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.45_MC) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.85_MC) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.45_MP) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.85_MP) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")
names(environment50.45_MG) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.85_MG) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")
names(environment50.45_NO) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

names(environment50.85_NO) <- c("bio1","bio2","bio3","bio4","bio5",
					"bio6","bio7","bio8","bio9","bio10",
					"bio11","bio12","bio13","bio14","bio15",
					"bio16","bio17","bio18","bio19")

# Subset environmental stack for future scenarios (Here, you should included only your selected variables)
env50.45.selected_AC <- subset(environment50.45_AC, c("bio7", "bio11", "bio15", "bio16"))
env50.85.selected_AC <- subset(environment50.85_AC, c("bio7", "bio11", "bio15", "bio16"))

env50.45.selected_BC <- subset(environment50.45_BC, c("bio7", "bio11", "bio15", "bio16"))
env50.85.selected_BC <- subset(environment50.85_BC, c("bio7", "bio11", "bio15", "bio16"))

env50.45.selected_CC <- subset(environment50.45_CC, c("bio7", "bio11", "bio15", "bio16"))
env50.85.selected_CC <- subset(environment50.85_CC, c("bio7", "bio11", "bio15", "bio16"))

env50.45.selected_CN <- subset(environment50.45_CN, c("bio7", "bio11", "bio15", "bio16"))
env50.85.selected_CN <- subset(environment50.85_CN, c("bio7", "bio11", "bio15", "bio16"))

env50.45.selected_GS <- subset(environment50.45_GS, c("bio7", "bio11", "bio15", "bio16"))
env50.85.selected_GS <- subset(environment50.85_GS, c("bio7", "bio11", "bio15", "bio16"))

env50.45.selected_HD <- subset(environment50.45_HD, c("bio7", "bio11", "bio15", "bio16"))
env50.85.selected_HD <- subset(environment50.85_HD, c("bio7", "bio11", "bio15", "bio16"))

env50.45.selected_HG <- subset(environment50.45_HG, c("bio7", "bio11", "bio15", "bio16"))
env50.85.selected_HG <- subset(environment50.85_HG, c("bio7", "bio11", "bio15", "bio16"))

env50.45.selected_HE <- subset(environment50.45_HE, c("bio7", "bio11", "bio15", "bio16"))
env50.85.selected_HE <- subset(environment50.85_HE, c("bio7", "bio11", "bio15", "bio16"))

env50.45.selected_IN <- subset(environment50.45_IN, c("bio7", "bio11", "bio15", "bio16"))
env50.85.selected_IN <- subset(environment50.85_IN, c("bio7", "bio11", "bio15", "bio16"))

env50.45.selected_IP <- subset(environment50.45_IP, c("bio7", "bio11", "bio15", "bio16"))
env50.85.selected_IP <- subset(environment50.85_IP, c("bio7", "bio11", "bio15", "bio16"))

env50.45.selected_MI <- subset(environment50.45_MI, c("bio7", "bio11", "bio15", "bio16"))
env50.85.selected_MI <- subset(environment50.85_MI, c("bio7", "bio11", "bio15", "bio16"))

env50.45.selected_MR <- subset(environment50.45_MR, c("bio7", "bio11", "bio15", "bio16"))
env50.85.selected_MR <- subset(environment50.85_MR, c("bio7", "bio11", "bio15", "bio16"))

env50.45.selected_MC <- subset(environment50.45_MC, c("bio7", "bio11", "bio15", "bio16"))
env50.85.selected_MC <- subset(environment50.85_MC, c("bio7", "bio11", "bio15", "bio16"))

env50.45.selected_MP <- subset(environment50.45_MP, c("bio7", "bio11", "bio15", "bio16"))
env50.85.selected_MP <- subset(environment50.85_MP, c("bio7", "bio11", "bio15", "bio16"))

env50.45.selected_MG <- subset(environment50.45_MG, c("bio7", "bio11", "bio15", "bio16"))
env50.85.selected_MG <- subset(environment50.85_MG, c("bio7", "bio11", "bio15", "bio16"))

env50.45.selected_NO <- subset(environment50.45_NO, c("bio7", "bio11", "bio15", "bio16"))
env50.85.selected_NO <- subset(environment50.85_NO, c("bio7", "bio11", "bio15", "bio16"))


############################################
## Model projection for the future - 2050 ##
############################################

### Listing all objects produced for future scenarios:
scenario.list<-list(env50.45.selected_AC, env50.85.selected_AC,
				env50.45.selected_BC, env50.85.selected_BC,
				env50.45.selected_CC, env50.85.selected_CC,
				env50.45.selected_CN, env50.85.selected_CN,
				env50.45.selected_GS, env50.85.selected_GS,
				env50.45.selected_HD, env50.85.selected_HD,
				env50.45.selected_HG, env50.85.selected_HG,
				env50.45.selected_HE, env50.85.selected_HE,
				env50.45.selected_IN, env50.85.selected_IN,
				env50.45.selected_IP, env50.85.selected_IP,
				env50.45.selected_MI, env50.85.selected_MI,
				env50.45.selected_MR, env50.85.selected_MR,
				env50.45.selected_MC, env50.85.selected_MC,
				env50.45.selected_MP, env50.85.selected_MP,
				env50.45.selected_MG, env50.85.selected_MG,
				env50.45.selected_NO, env50.85.selected_NO)		

names(scenario.list)<-c("env50.45_AC", "env50.85_AC", 
				"env50.45_BC", "env50.85_BC",
				"env50.45_CC", "env50.85_CC",
				"env50.45_CN", "env50.85_CN",
				"env50.45_GS", "env50.85_GS",
				"env50.45_HD", "env50.85_HD",
				"env50.45_HG", "env50.85_HG",
				"env50.45_HE", "env50.85_HE",
				"env50.45_IN", "env50.85_IN",
				"env50.45_IP", "env50.85_IP",
				"env50.45_MI", "env50.85_MI",
				"env50.45_MR", "env50.85_MR",
				"env50.45_MC", "env50.85_MC",
				"env50.45_MP", "env50.85_MP",
				"env50.45_MG", "env50.85_MG",
				"env50.45_NO", "env50.85_NO")
				


######################################
## 2050 - optimistic scenario (rcp 45)
######################################

spp.projections.2050.rcp45_1_AC <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.45.selected_AC,
	proj.name = "2050.rcp45_1_AC",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp45_2_AC <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.45.selected_AC,
	proj.name = "2050.rcp45_2_AC",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp45_1_BC <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.45.selected_BC,
	proj.name = "2050.rcp45_1_BC",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp45_2_BC <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.45.selected_BC,
	proj.name = "2050.rcp45_2_BC",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp45_1_CC <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.45.selected_CC,
	proj.name = "2050.rcp45_1_CC",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp45_2_CC <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.45.selected_CC,
	proj.name = "2050.rcp45_2_CC",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp45_1_CN <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.45.selected_CN,
	proj.name = "2050.rcp45_1_CN",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp45_2_CN <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.45.selected_CN,
	proj.name = "2050.rcp45_2_CN",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp45_1_GS <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.45.selected_GS,
	proj.name = "2050.rcp45_1_GS",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp45_2_GS <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.45.selected_GS,
	proj.name = "2050.rcp45_2_GS",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp45_1_HD <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.45.selected_HD,
	proj.name = "2050.rcp45_1_HD",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp45_2_HD <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.45.selected_HD,
	proj.name = "2050.rcp45_2_HD",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp45_1_HG <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.45.selected_HG,
	proj.name = "2050.rcp45_1_HG",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp45_2_HG <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.45.selected_HG,
	proj.name = "2050.rcp45_2_HG",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp45_1_HE <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.45.selected_HE,
	proj.name = "2050.rcp45_1_HE",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp45_2_HE <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.45.selected_HE,
	proj.name = "2050.rcp45_2_HE",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp45_1_IN <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.45.selected_IN,
	proj.name = "2050.rcp45_1_IN",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp45_2_IN <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.45.selected_IN,
	proj.name = "2050.rcp45_2_IN",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp45_1_IP <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.45.selected_IP,
	proj.name = "2050.rcp45_1_IP",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp45_2_IP <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.45.selected_IP,
	proj.name = "2050.rcp45_2_IP",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp45_1_MI <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.45.selected_MI,
	proj.name = "2050.rcp45_1_MI",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp45_2_MI <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.45.selected_MI,
	proj.name = "2050.rcp45_2_MI",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp45_1_MR <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.45.selected_MR,
	proj.name = "2050.rcp45_1_MR",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp45_2_MR <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.45.selected_MR,
	proj.name = "2050.rcp45_2_MR",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp45_1_MC <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.45.selected_MC,
	proj.name = "2050.rcp45_1_MC",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp45_2_MC <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.45.selected_MC,
	proj.name = "2050.rcp45_2_MC",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp45_1_MP <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.45.selected_MP,
	proj.name = "2050.rcp45_1_MP",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp45_2_MP <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.45.selected_MP,
	proj.name = "2050.rcp45_2_MP",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp45_1_MG <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.45.selected_MG,
	proj.name = "2050.rcp45_1_MG",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp45_2_MG <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.45.selected_MG,
	proj.name = "2050.rcp45_2_MG",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp45_1_NO <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.45.selected_NO,
	proj.name = "2050.rcp45_1_NO",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp45_2_NO <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.45.selected_NO,
	proj.name = "2050.rcp45_2_NO",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")


# Stack projections
projections.2050.rcp45_1_AC <- stack("C:/Models/Occurrence/proj_2050.rcp45_1_AC/proj_2050.rcp45_1_AC_Occurrence.grd")
names(projections.2050.rcp45_1_AC)
projections.2050.rcp45_2_AC <- stack("C:/Models/Occurrence/proj_2050.rcp45_2_AC/proj_2050.rcp45_2_AC_Occurrence.grd")
names(projections.2050.rcp45_2_AC)

projections.2050.rcp45_1_BC <- stack("C:/Models/Occurrence/proj_2050.rcp45_1_BC/proj_2050.rcp45_1_BC_Occurrence.grd")
names(projections.2050.rcp45_1_BC)
projections.2050.rcp45_2_BC <- stack("C:/Models/Occurrence/proj_2050.rcp45_2_BC/proj_2050.rcp45_2_BC_Occurrence.grd")
names(projections.2050.rcp45_2_BC)

projections.2050.rcp45_1_CC <- stack("C:/Models/Occurrence/proj_2050.rcp45_1_CC/proj_2050.rcp45_1_CC_Occurrence.grd")
names(projections.2050.rcp45_1_CC)
projections.2050.rcp45_2_CC <- stack("C:/Models/Occurrence/proj_2050.rcp45_2_CC/proj_2050.rcp45_2_CC_Occurrence.grd")
names(projections.2050.rcp45_2_CC)

projections.2050.rcp45_1_CN <- stack("C:/Models/Occurrence/proj_2050.rcp45_1_CN/proj_2050.rcp45_1_CN_Occurrence.grd")
names(projections.2050.rcp45_1_CN)
projections.2050.rcp45_2_CN <- stack("C:/Models/Occurrence/proj_2050.rcp45_2_CN/proj_2050.rcp45_2_CN_Occurrence.grd")
names(projections.2050.rcp45_2_CN)

projections.2050.rcp45_1_GS <- stack("C:/Models/Occurrence/proj_2050.rcp45_1_GS/proj_2050.rcp45_1_GS_Occurrence.grd")
names(projections.2050.rcp45_1_GS)
projections.2050.rcp45_2_GS <- stack("C:/Models/Occurrence/proj_2050.rcp45_2_GS/proj_2050.rcp45_2_GS_Occurrence.grd")
names(projections.2050.rcp45_2_GS)

projections.2050.rcp45_1_HD <- stack("C:/Models/Occurrence/proj_2050.rcp45_1_HD/proj_2050.rcp45_1_HD_Occurrence.grd")
names(projections.2050.rcp45_1_HD)
projections.2050.rcp45_2_HD <- stack("C:/Models/Occurrence/proj_2050.rcp45_2_HD/proj_2050.rcp45_2_HD_Occurrence.grd")
names(projections.2050.rcp45_2_HD)

projections.2050.rcp45_1_HG <- stack("C:/Models/Occurrence/proj_2050.rcp45_1_HG/proj_2050.rcp45_1_HG_Occurrence.grd")
names(projections.2050.rcp45_1_HG)
projections.2050.rcp45_2_HG <- stack("C:/Models/Occurrence/proj_2050.rcp45_2_HG/proj_2050.rcp45_2_HG_Occurrence.grd")
names(projections.2050.rcp45_2_HG)

projections.2050.rcp45_1_HE <- stack("C:/Models/Occurrence/proj_2050.rcp45_1_HE/proj_2050.rcp45_1_HE_Occurrence.grd")
names(projections.2050.rcp45_1_HE)
projections.2050.rcp45_2_HE <- stack("C:/Models/Occurrence/proj_2050.rcp45_2_HE/proj_2050.rcp45_2_HE_Occurrence.grd")
names(projections.2050.rcp45_2_HE)

projections.2050.rcp45_1_IN <- stack("C:/Models/Occurrence/proj_2050.rcp45_1_IN/proj_2050.rcp45_1_IN_Occurrence.grd")
names(projections.2050.rcp45_1_IN)
projections.2050.rcp45_2_IN <- stack("C:/Models/Occurrence/proj_2050.rcp45_2_IN/proj_2050.rcp45_2_IN_Occurrence.grd")
names(projections.2050.rcp45_2_IN)

projections.2050.rcp45_1_IP <- stack("C:/Models/Occurrence/proj_2050.rcp45_1_IP/proj_2050.rcp45_1_IP_Occurrence.grd")
names(projections.2050.rcp45_1_IP)
projections.2050.rcp45_2_IP <- stack("C:/Models/Occurrence/proj_2050.rcp45_2_IP/proj_2050.rcp45_2_IP_Occurrence.grd")
names(projections.2050.rcp45_2_IP)

projections.2050.rcp45_1_MI <- stack("C:/Models/Occurrence/proj_2050.rcp45_1_MI/proj_2050.rcp45_1_MI_Occurrence.grd")
names(projections.2050.rcp45_1_MI)
projections.2050.rcp45_2_MI <- stack("C:/Models/Occurrence/proj_2050.rcp45_2_MI/proj_2050.rcp45_2_MI_Occurrence.grd")
names(projections.2050.rcp45_2_MI)

projections.2050.rcp45_1_MR <- stack("C:/Models/Occurrence/proj_2050.rcp45_1_MR/proj_2050.rcp45_1_MR_Occurrence.grd")
names(projections.2050.rcp45_1_MR)
projections.2050.rcp45_2_MR <- stack("C:/Models/Occurrence/proj_2050.rcp45_2_MR/proj_2050.rcp45_2_MR_Occurrence.grd")
names(projections.2050.rcp45_2_MR)

projections.2050.rcp45_1_MC <- stack("C:/Models/Occurrence/proj_2050.rcp45_1_MC/proj_2050.rcp45_1_MC_Occurrence.grd")
names(projections.2050.rcp45_1_MC)
projections.2050.rcp45_2_MC <- stack("C:/Models/Occurrence/proj_2050.rcp45_2_MC/proj_2050.rcp45_2_MC_Occurrence.grd")
names(projections.2050.rcp45_2_MC)

projections.2050.rcp45_1_MP <- stack("C:/Models/Occurrence/proj_2050.rcp45_1_MP/proj_2050.rcp45_1_MP_Occurrence.grd")
names(projections.2050.rcp45_1_MP)
projections.2050.rcp45_2_MP <- stack("C:/Models/Occurrence/proj_2050.rcp45_2_MP/proj_2050.rcp45_2_MP_Occurrence.grd")
names(projections.2050.rcp45_2_MP)

projections.2050.rcp45_1_MG <- stack("C:/Models/Occurrence/proj_2050.rcp45_1_MG/proj_2050.rcp45_1_MG_Occurrence.grd")
names(projections.2050.rcp45_1_MG)
projections.2050.rcp45_2_MG <- stack("C:/Models/Occurrence/proj_2050.rcp45_2_MG/proj_2050.rcp45_2_MG_Occurrence.grd")
names(projections.2050.rcp45_2_MG)

projections.2050.rcp45_1_NO <- stack("C:/Models/Occurrence/proj_2050.rcp45_1_NO/proj_2050.rcp45_1_NO_Occurrence.grd")
names(projections.2050.rcp45_1_NO)
projections.2050.rcp45_2_NO <- stack("C:/Models/Occurrence/proj_2050.rcp45_2_NO/proj_2050.rcp45_2_NO_Occurrence.grd")
names(projections.2050.rcp45_2_NO)





### Average consensus models for each algorithm:
#AC
projections.RF.all.2050.rcp45_AC <- subset(projections.2050.rcp45_1_AC, grep("RF", names(projections.2050.rcp45_1_AC)))
projections.RF.mean.2050.rcp45_AC <- mean(projections.RF.all.2050.rcp45_AC)/10

projections.GBM.all.2050.rcp45_AC <-subset(projections.2050.rcp45_1_AC, grep("GBM", names(projections.2050.rcp45_1_AC)))
projections.GBM.mean.2050.rcp45_AC <- mean(projections.GBM.all.2050.rcp45_AC)/10

projections.CTA.all.2050.rcp45_AC <-subset(projections.2050.rcp45_1_AC,grep("CTA", names(projections.2050.rcp45_1_AC)))
projections.CTA.mean.2050.rcp45_AC <- mean(projections.CTA.all.2050.rcp45_AC)/10

projections.GLM.all.2050.rcp45_AC <-subset(projections.2050.rcp45_2_AC,grep("GLM", names(projections.2050.rcp45_2_AC)))
projections.GLM.mean.2050.rcp45_AC <- mean(projections.GLM.all.2050.rcp45_AC)/10

projections.GAM.all.2050.rcp45_AC <-subset(projections.2050.rcp45_2_AC,grep("GAM", names(projections.2050.rcp45_2_AC)))
projections.GAM.mean.2050.rcp45_AC <- mean(projections.GAM.all.2050.rcp45_AC)/10

projections.ANN.all.2050.rcp45_AC <- subset(projections.2050.rcp45_2_AC,grep("ANN", names(projections.2050.rcp45_2_AC)))
projections.ANN.mean.2050.rcp45_AC <- mean(projections.ANN.all.2050.rcp45_AC)/10

projections.SRE.all.2050.rcp45_AC <- subset(projections.2050.rcp45_2_AC,grep("SRE", names(projections.2050.rcp45_2_AC)))
projections.SRE.mean.2050.rcp45_AC <- mean(projections.SRE.all.2050.rcp45_AC)/10

projections.MARS.all.2050.rcp45_AC <- subset(projections.2050.rcp45_2_AC,grep("MARS", names(projections.2050.rcp45_2_AC)))
projections.MARS.mean.2050.rcp45_AC <- mean(projections.MARS.all.2050.rcp45_AC)/10

projections.FDA.all.2050.rcp45_AC <- subset(projections.2050.rcp45_2_AC,grep("FDA", names(projections.2050.rcp45_2_AC)))
projections.FDA.mean.2050.rcp45_AC <- mean(projections.FDA.all.2050.rcp45_AC)/10

projections.MAXENT.all.2050.rcp45_AC <- subset(projections.2050.rcp45_2_AC,grep("MAXENT.Phillips", names(projections.2050.rcp45_2_AC)))
projections.MAXENT.mean.2050.rcp45_AC <- mean(projections.MAXENT.all.2050.rcp45_AC)/10


#BC
projections.RF.all.2050.rcp45_BC <- subset(projections.2050.rcp45_1_BC, grep("RF", names(projections.2050.rcp45_1_BC)))
projections.RF.mean.2050.rcp45_BC <- mean(projections.RF.all.2050.rcp45_BC)/10

projections.GBM.all.2050.rcp45_BC <-subset(projections.2050.rcp45_1_BC, grep("GBM", names(projections.2050.rcp45_1_BC)))
projections.GBM.mean.2050.rcp45_BC <- mean(projections.GBM.all.2050.rcp45_BC)/10

projections.CTA.all.2050.rcp45_BC <-subset(projections.2050.rcp45_1_BC,grep("CTA", names(projections.2050.rcp45_1_BC)))
projections.CTA.mean.2050.rcp45_BC <- mean(projections.CTA.all.2050.rcp45_BC)/10

projections.GLM.all.2050.rcp45_BC <-subset(projections.2050.rcp45_2_BC,grep("GLM", names(projections.2050.rcp45_2_BC)))
projections.GLM.mean.2050.rcp45_BC <- mean(projections.GLM.all.2050.rcp45_BC)/10

projections.GAM.all.2050.rcp45_BC <-subset(projections.2050.rcp45_2_BC,grep("GAM", names(projections.2050.rcp45_2_BC)))
projections.GAM.mean.2050.rcp45_BC <- mean(projections.GAM.all.2050.rcp45_BC)/10

projections.ANN.all.2050.rcp45_BC <- subset(projections.2050.rcp45_2_BC,grep("ANN", names(projections.2050.rcp45_2_BC)))
projections.ANN.mean.2050.rcp45_BC <- mean(projections.ANN.all.2050.rcp45_BC)/10

projections.SRE.all.2050.rcp45_BC <- subset(projections.2050.rcp45_2_BC,grep("SRE", names(projections.2050.rcp45_2_BC)))
projections.SRE.mean.2050.rcp45_BC <- mean(projections.SRE.all.2050.rcp45_BC)/10

projections.MARS.all.2050.rcp45_BC <- subset(projections.2050.rcp45_2_BC,grep("MARS", names(projections.2050.rcp45_2_BC)))
projections.MARS.mean.2050.rcp45_BC <- mean(projections.MARS.all.2050.rcp45_BC)/10

projections.FDA.all.2050.rcp45_BC <- subset(projections.2050.rcp45_2_BC,grep("FDA", names(projections.2050.rcp45_2_BC)))
projections.FDA.mean.2050.rcp45_BC <- mean(projections.FDA.all.2050.rcp45_BC)/10

projections.MAXENT.all.2050.rcp45_BC <- subset(projections.2050.rcp45_2_BC,grep("MAXENT.Phillips", names(projections.2050.rcp45_2_BC)))
projections.MAXENT.mean.2050.rcp45_BC <- mean(projections.MAXENT.all.2050.rcp45_BC)/10


#CC
projections.RF.all.2050.rcp45_CC <- subset(projections.2050.rcp45_1_CC, grep("RF", names(projections.2050.rcp45_1_CC)))
projections.RF.mean.2050.rcp45_CC <- mean(projections.RF.all.2050.rcp45_CC)/10

projections.GBM.all.2050.rcp45_CC <-subset(projections.2050.rcp45_1_CC, grep("GBM", names(projections.2050.rcp45_1_CC)))
projections.GBM.mean.2050.rcp45_CC <- mean(projections.GBM.all.2050.rcp45_CC)/10

projections.CTA.all.2050.rcp45_CC <-subset(projections.2050.rcp45_1_CC,grep("CTA", names(projections.2050.rcp45_1_CC)))
projections.CTA.mean.2050.rcp45_CC <- mean(projections.CTA.all.2050.rcp45_CC)/10

projections.GLM.all.2050.rcp45_CC <-subset(projections.2050.rcp45_2_CC,grep("GLM", names(projections.2050.rcp45_2_CC)))
projections.GLM.mean.2050.rcp45_CC <- mean(projections.GLM.all.2050.rcp45_CC)/10

projections.GAM.all.2050.rcp45_CC <-subset(projections.2050.rcp45_2_CC,grep("GAM", names(projections.2050.rcp45_2_CC)))
projections.GAM.mean.2050.rcp45_CC <- mean(projections.GAM.all.2050.rcp45_CC)/10

projections.ANN.all.2050.rcp45_CC <- subset(projections.2050.rcp45_2_CC,grep("ANN", names(projections.2050.rcp45_2_CC)))
projections.ANN.mean.2050.rcp45_CC <- mean(projections.ANN.all.2050.rcp45_CC)/10

projections.SRE.all.2050.rcp45_CC <- subset(projections.2050.rcp45_2_CC,grep("SRE", names(projections.2050.rcp45_2_CC)))
projections.SRE.mean.2050.rcp45_CC <- mean(projections.SRE.all.2050.rcp45_CC)/10

projections.MARS.all.2050.rcp45_CC <- subset(projections.2050.rcp45_2_CC,grep("MARS", names(projections.2050.rcp45_2_CC)))
projections.MARS.mean.2050.rcp45_CC <- mean(projections.MARS.all.2050.rcp45_CC)/10

projections.FDA.all.2050.rcp45_CC <- subset(projections.2050.rcp45_2_CC,grep("FDA", names(projections.2050.rcp45_2_CC)))
projections.FDA.mean.2050.rcp45_CC <- mean(projections.FDA.all.2050.rcp45_CC)/10

projections.MAXENT.all.2050.rcp45_CC <- subset(projections.2050.rcp45_2_CC,grep("MAXENT.Phillips", names(projections.2050.rcp45_2_CC)))
projections.MAXENT.mean.2050.rcp45_CC <- mean(projections.MAXENT.all.2050.rcp45_CC)/10


#CN
projections.RF.all.2050.rcp45_CN <- subset(projections.2050.rcp45_1_CN, grep("RF", names(projections.2050.rcp45_1_CN)))
projections.RF.mean.2050.rcp45_CN <- mean(projections.RF.all.2050.rcp45_CN)/10

projections.GBM.all.2050.rcp45_CN <-subset(projections.2050.rcp45_1_CN, grep("GBM", names(projections.2050.rcp45_1_CN)))
projections.GBM.mean.2050.rcp45_CN <- mean(projections.GBM.all.2050.rcp45_CN)/10

projections.CTA.all.2050.rcp45_CN <-subset(projections.2050.rcp45_1_CN,grep("CTA", names(projections.2050.rcp45_1_CN)))
projections.CTA.mean.2050.rcp45_CN <- mean(projections.CTA.all.2050.rcp45_CN)/10

projections.GLM.all.2050.rcp45_CN <-subset(projections.2050.rcp45_2_CN,grep("GLM", names(projections.2050.rcp45_2_CN)))
projections.GLM.mean.2050.rcp45_CN <- mean(projections.GLM.all.2050.rcp45_CN)/10

projections.GAM.all.2050.rcp45_CN <-subset(projections.2050.rcp45_2_CN,grep("GAM", names(projections.2050.rcp45_2_CN)))
projections.GAM.mean.2050.rcp45_CN <- mean(projections.GAM.all.2050.rcp45_CN)/10

projections.ANN.all.2050.rcp45_CN <- subset(projections.2050.rcp45_2_CN,grep("ANN", names(projections.2050.rcp45_2_CN)))
projections.ANN.mean.2050.rcp45_CN <- mean(projections.ANN.all.2050.rcp45_CN)/10

projections.SRE.all.2050.rcp45_CN <- subset(projections.2050.rcp45_2_CN,grep("SRE", names(projections.2050.rcp45_2_CN)))
projections.SRE.mean.2050.rcp45_CN <- mean(projections.SRE.all.2050.rcp45_CN)/10

projections.MARS.all.2050.rcp45_CN <- subset(projections.2050.rcp45_2_CN,grep("MARS", names(projections.2050.rcp45_2_CN)))
projections.MARS.mean.2050.rcp45_CN <- mean(projections.MARS.all.2050.rcp45_CN)/10

projections.FDA.all.2050.rcp45_CN <- subset(projections.2050.rcp45_2_CN,grep("FDA", names(projections.2050.rcp45_2_CN)))
projections.FDA.mean.2050.rcp45_CN <- mean(projections.FDA.all.2050.rcp45_CN)/10

projections.MAXENT.all.2050.rcp45_CN <- subset(projections.2050.rcp45_2_CN,grep("MAXENT.Phillips", names(projections.2050.rcp45_2_CN)))
projections.MAXENT.mean.2050.rcp45_CN <- mean(projections.MAXENT.all.2050.rcp45_CN)/10


#GS
projections.RF.all.2050.rcp45_GS <- subset(projections.2050.rcp45_1_GS, grep("RF", names(projections.2050.rcp45_1_GS)))
projections.RF.mean.2050.rcp45_GS <- mean(projections.RF.all.2050.rcp45_GS)/10

projections.GBM.all.2050.rcp45_GS <-subset(projections.2050.rcp45_1_GS, grep("GBM", names(projections.2050.rcp45_1_GS)))
projections.GBM.mean.2050.rcp45_GS <- mean(projections.GBM.all.2050.rcp45_GS)/10

projections.CTA.all.2050.rcp45_GS <-subset(projections.2050.rcp45_1_GS,grep("CTA", names(projections.2050.rcp45_1_GS)))
projections.CTA.mean.2050.rcp45_GS <- mean(projections.CTA.all.2050.rcp45_GS)/10

projections.GLM.all.2050.rcp45_GS <-subset(projections.2050.rcp45_2_GS,grep("GLM", names(projections.2050.rcp45_2_GS)))
projections.GLM.mean.2050.rcp45_GS <- mean(projections.GLM.all.2050.rcp45_GS)/10

projections.GAM.all.2050.rcp45_GS <-subset(projections.2050.rcp45_2_GS,grep("GAM", names(projections.2050.rcp45_2_GS)))
projections.GAM.mean.2050.rcp45_GS <- mean(projections.GAM.all.2050.rcp45_GS)/10

projections.ANN.all.2050.rcp45_GS <- subset(projections.2050.rcp45_2_GS,grep("ANN", names(projections.2050.rcp45_2_GS)))
projections.ANN.mean.2050.rcp45_GS <- mean(projections.ANN.all.2050.rcp45_GS)/10

projections.SRE.all.2050.rcp45_GS <- subset(projections.2050.rcp45_2_GS,grep("SRE", names(projections.2050.rcp45_2_GS)))
projections.SRE.mean.2050.rcp45_GS <- mean(projections.SRE.all.2050.rcp45_GS)/10

projections.MARS.all.2050.rcp45_GS <- subset(projections.2050.rcp45_2_GS,grep("MARS", names(projections.2050.rcp45_2_GS)))
projections.MARS.mean.2050.rcp45_GS <- mean(projections.MARS.all.2050.rcp45_GS)/10

projections.FDA.all.2050.rcp45_GS <- subset(projections.2050.rcp45_2_GS,grep("FDA", names(projections.2050.rcp45_2_GS)))
projections.FDA.mean.2050.rcp45_GS <- mean(projections.FDA.all.2050.rcp45_GS)/10

projections.MAXENT.all.2050.rcp45_GS <- subset(projections.2050.rcp45_2_GS,grep("MAXENT.Phillips", names(projections.2050.rcp45_2_GS)))
projections.MAXENT.mean.2050.rcp45_GS <- mean(projections.MAXENT.all.2050.rcp45_GS)/10



#HD
projections.RF.all.2050.rcp45_HD <- subset(projections.2050.rcp45_1_HD, grep("RF", names(projections.2050.rcp45_1_HD)))
projections.RF.mean.2050.rcp45_HD <- mean(projections.RF.all.2050.rcp45_HD)/10

projections.GBM.all.2050.rcp45_HD <-subset(projections.2050.rcp45_1_HD, grep("GBM", names(projections.2050.rcp45_1_HD)))
projections.GBM.mean.2050.rcp45_HD <- mean(projections.GBM.all.2050.rcp45_HD)/10

projections.CTA.all.2050.rcp45_HD <-subset(projections.2050.rcp45_1_HD,grep("CTA", names(projections.2050.rcp45_1_HD)))
projections.CTA.mean.2050.rcp45_HD <- mean(projections.CTA.all.2050.rcp45_HD)/10

projections.GLM.all.2050.rcp45_HD <-subset(projections.2050.rcp45_2_HD,grep("GLM", names(projections.2050.rcp45_2_HD)))
projections.GLM.mean.2050.rcp45_HD <- mean(projections.GLM.all.2050.rcp45_HD)/10

projections.GAM.all.2050.rcp45_HD <-subset(projections.2050.rcp45_2_HD,grep("GAM", names(projections.2050.rcp45_2_HD)))
projections.GAM.mean.2050.rcp45_HD <- mean(projections.GAM.all.2050.rcp45_HD)/10

projections.ANN.all.2050.rcp45_HD <- subset(projections.2050.rcp45_2_HD,grep("ANN", names(projections.2050.rcp45_2_HD)))
projections.ANN.mean.2050.rcp45_HD <- mean(projections.ANN.all.2050.rcp45_HD)/10

projections.SRE.all.2050.rcp45_HD <- subset(projections.2050.rcp45_2_HD,grep("SRE", names(projections.2050.rcp45_2_HD)))
projections.SRE.mean.2050.rcp45_HD <- mean(projections.SRE.all.2050.rcp45_HD)/10

projections.MARS.all.2050.rcp45_HD <- subset(projections.2050.rcp45_2_HD,grep("MARS", names(projections.2050.rcp45_2_HD)))
projections.MARS.mean.2050.rcp45_HD <- mean(projections.MARS.all.2050.rcp45_HD)/10

projections.FDA.all.2050.rcp45_HD <- subset(projections.2050.rcp45_2_HD,grep("FDA", names(projections.2050.rcp45_2_HD)))
projections.FDA.mean.2050.rcp45_HD <- mean(projections.FDA.all.2050.rcp45_HD)/10

projections.MAXENT.all.2050.rcp45_HD <- subset(projections.2050.rcp45_2_HD,grep("MAXENT.Phillips", names(projections.2050.rcp45_2_HD)))
projections.MAXENT.mean.2050.rcp45_HD <- mean(projections.MAXENT.all.2050.rcp45_HD)/10



#HG
projections.RF.all.2050.rcp45_HG <- subset(projections.2050.rcp45_1_HG, grep("RF", names(projections.2050.rcp45_1_HG)))
projections.RF.mean.2050.rcp45_HG <- mean(projections.RF.all.2050.rcp45_HG)/10

projections.GBM.all.2050.rcp45_HG <-subset(projections.2050.rcp45_1_HG, grep("GBM", names(projections.2050.rcp45_1_HG)))
projections.GBM.mean.2050.rcp45_HG <- mean(projections.GBM.all.2050.rcp45_HG)/10

projections.CTA.all.2050.rcp45_HG <-subset(projections.2050.rcp45_1_HG,grep("CTA", names(projections.2050.rcp45_1_HG)))
projections.CTA.mean.2050.rcp45_HG <- mean(projections.CTA.all.2050.rcp45_HG)/10

projections.GLM.all.2050.rcp45_HG <-subset(projections.2050.rcp45_2_HG,grep("GLM", names(projections.2050.rcp45_2_HG)))
projections.GLM.mean.2050.rcp45_HG <- mean(projections.GLM.all.2050.rcp45_HG)/10

projections.GAM.all.2050.rcp45_HG <-subset(projections.2050.rcp45_2_HG,grep("GAM", names(projections.2050.rcp45_2_HG)))
projections.GAM.mean.2050.rcp45_HG <- mean(projections.GAM.all.2050.rcp45_HG)/10

projections.ANN.all.2050.rcp45_HG <- subset(projections.2050.rcp45_2_HG,grep("ANN", names(projections.2050.rcp45_2_HG)))
projections.ANN.mean.2050.rcp45_HG <- mean(projections.ANN.all.2050.rcp45_HG)/10

projections.SRE.all.2050.rcp45_HG <- subset(projections.2050.rcp45_2_HG,grep("SRE", names(projections.2050.rcp45_2_HG)))
projections.SRE.mean.2050.rcp45_HG <- mean(projections.SRE.all.2050.rcp45_HG)/10

projections.MARS.all.2050.rcp45_HG <- subset(projections.2050.rcp45_2_HG,grep("MARS", names(projections.2050.rcp45_2_HG)))
projections.MARS.mean.2050.rcp45_HG <- mean(projections.MARS.all.2050.rcp45_HG)/10

projections.FDA.all.2050.rcp45_HG <- subset(projections.2050.rcp45_2_HG,grep("FDA", names(projections.2050.rcp45_2_HG)))
projections.FDA.mean.2050.rcp45_HG <- mean(projections.FDA.all.2050.rcp45_HG)/10

projections.MAXENT.all.2050.rcp45_HG <- subset(projections.2050.rcp45_2_HG,grep("MAXENT.Phillips", names(projections.2050.rcp45_2_HG)))
projections.MAXENT.mean.2050.rcp45_HG <- mean(projections.MAXENT.all.2050.rcp45_HG)/10



#HE
projections.RF.all.2050.rcp45_HE <- subset(projections.2050.rcp45_1_HE, grep("RF", names(projections.2050.rcp45_1_HE)))
projections.RF.mean.2050.rcp45_HE <- mean(projections.RF.all.2050.rcp45_HE)/10

projections.GBM.all.2050.rcp45_HE <-subset(projections.2050.rcp45_1_HE, grep("GBM", names(projections.2050.rcp45_1_HE)))
projections.GBM.mean.2050.rcp45_HE <- mean(projections.GBM.all.2050.rcp45_HE)/10

projections.CTA.all.2050.rcp45_HE <-subset(projections.2050.rcp45_1_HE,grep("CTA", names(projections.2050.rcp45_1_HE)))
projections.CTA.mean.2050.rcp45_HE <- mean(projections.CTA.all.2050.rcp45_HE)/10

projections.GLM.all.2050.rcp45_HE <-subset(projections.2050.rcp45_2_HE,grep("GLM", names(projections.2050.rcp45_2_HE)))
projections.GLM.mean.2050.rcp45_HE <- mean(projections.GLM.all.2050.rcp45_HE)/10

projections.GAM.all.2050.rcp45_HE <-subset(projections.2050.rcp45_2_HE,grep("GAM", names(projections.2050.rcp45_2_HE)))
projections.GAM.mean.2050.rcp45_HE <- mean(projections.GAM.all.2050.rcp45_HE)/10

projections.ANN.all.2050.rcp45_HE <- subset(projections.2050.rcp45_2_HE,grep("ANN", names(projections.2050.rcp45_2_HE)))
projections.ANN.mean.2050.rcp45_HE <- mean(projections.ANN.all.2050.rcp45_HE)/10

projections.SRE.all.2050.rcp45_HE <- subset(projections.2050.rcp45_2_HE,grep("SRE", names(projections.2050.rcp45_2_HE)))
projections.SRE.mean.2050.rcp45_HE <- mean(projections.SRE.all.2050.rcp45_HE)/10

projections.MARS.all.2050.rcp45_HE <- subset(projections.2050.rcp45_2_HE,grep("MARS", names(projections.2050.rcp45_2_HE)))
projections.MARS.mean.2050.rcp45_HE <- mean(projections.MARS.all.2050.rcp45_HE)/10

projections.FDA.all.2050.rcp45_HE <- subset(projections.2050.rcp45_2_HE,grep("FDA", names(projections.2050.rcp45_2_HE)))
projections.FDA.mean.2050.rcp45_HE <- mean(projections.FDA.all.2050.rcp45_HE)/10

projections.MAXENT.all.2050.rcp45_HE <- subset(projections.2050.rcp45_2_HE,grep("MAXENT.Phillips", names(projections.2050.rcp45_2_HE)))
projections.MAXENT.mean.2050.rcp45_HE <- mean(projections.MAXENT.all.2050.rcp45_HE)/10



#IN
projections.RF.all.2050.rcp45_IN <- subset(projections.2050.rcp45_1_IN, grep("RF", names(projections.2050.rcp45_1_IN)))
projections.RF.mean.2050.rcp45_IN <- mean(projections.RF.all.2050.rcp45_IN)/10

projections.GBM.all.2050.rcp45_IN <-subset(projections.2050.rcp45_1_IN, grep("GBM", names(projections.2050.rcp45_1_IN)))
projections.GBM.mean.2050.rcp45_IN <- mean(projections.GBM.all.2050.rcp45_IN)/10

projections.CTA.all.2050.rcp45_IN <-subset(projections.2050.rcp45_1_IN,grep("CTA", names(projections.2050.rcp45_1_IN)))
projections.CTA.mean.2050.rcp45_IN <- mean(projections.CTA.all.2050.rcp45_IN)/10

projections.GLM.all.2050.rcp45_IN <-subset(projections.2050.rcp45_2_IN,grep("GLM", names(projections.2050.rcp45_2_IN)))
projections.GLM.mean.2050.rcp45_IN <- mean(projections.GLM.all.2050.rcp45_IN)/10

projections.GAM.all.2050.rcp45_IN <-subset(projections.2050.rcp45_2_IN,grep("GAM", names(projections.2050.rcp45_2_IN)))
projections.GAM.mean.2050.rcp45_IN <- mean(projections.GAM.all.2050.rcp45_IN)/10

projections.ANN.all.2050.rcp45_IN <- subset(projections.2050.rcp45_2_IN,grep("ANN", names(projections.2050.rcp45_2_IN)))
projections.ANN.mean.2050.rcp45_IN <- mean(projections.ANN.all.2050.rcp45_IN)/10

projections.SRE.all.2050.rcp45_IN <- subset(projections.2050.rcp45_2_IN,grep("SRE", names(projections.2050.rcp45_2_IN)))
projections.SRE.mean.2050.rcp45_IN <- mean(projections.SRE.all.2050.rcp45_IN)/10

projections.MARS.all.2050.rcp45_IN <- subset(projections.2050.rcp45_2_IN,grep("MARS", names(projections.2050.rcp45_2_IN)))
projections.MARS.mean.2050.rcp45_IN <- mean(projections.MARS.all.2050.rcp45_IN)/10

projections.FDA.all.2050.rcp45_IN <- subset(projections.2050.rcp45_2_IN,grep("FDA", names(projections.2050.rcp45_2_IN)))
projections.FDA.mean.2050.rcp45_IN <- mean(projections.FDA.all.2050.rcp45_IN)/10

projections.MAXENT.all.2050.rcp45_IN <- subset(projections.2050.rcp45_2_IN,grep("MAXENT.Phillips", names(projections.2050.rcp45_2_IN)))
projections.MAXENT.mean.2050.rcp45_IN <- mean(projections.MAXENT.all.2050.rcp45_IN)/10


#IP
projections.RF.all.2050.rcp45_IP <- subset(projections.2050.rcp45_1_IP, grep("RF", names(projections.2050.rcp45_1_IP)))
projections.RF.mean.2050.rcp45_IP <- mean(projections.RF.all.2050.rcp45_IP)/10

projections.GBM.all.2050.rcp45_IP <-subset(projections.2050.rcp45_1_IP, grep("GBM", names(projections.2050.rcp45_1_IP)))
projections.GBM.mean.2050.rcp45_IP <- mean(projections.GBM.all.2050.rcp45_IP)/10

projections.CTA.all.2050.rcp45_IP <-subset(projections.2050.rcp45_1_IP,grep("CTA", names(projections.2050.rcp45_1_IP)))
projections.CTA.mean.2050.rcp45_IP <- mean(projections.CTA.all.2050.rcp45_IP)/10

projections.GLM.all.2050.rcp45_IP <-subset(projections.2050.rcp45_2_IP,grep("GLM", names(projections.2050.rcp45_2_IP)))
projections.GLM.mean.2050.rcp45_IP <- mean(projections.GLM.all.2050.rcp45_IP)/10

projections.GAM.all.2050.rcp45_IP <-subset(projections.2050.rcp45_2_IP,grep("GAM", names(projections.2050.rcp45_2_IP)))
projections.GAM.mean.2050.rcp45_IP <- mean(projections.GAM.all.2050.rcp45_IP)/10

projections.ANN.all.2050.rcp45_IP <- subset(projections.2050.rcp45_2_IP,grep("ANN", names(projections.2050.rcp45_2_IP)))
projections.ANN.mean.2050.rcp45_IP <- mean(projections.ANN.all.2050.rcp45_IP)/10

projections.SRE.all.2050.rcp45_IP <- subset(projections.2050.rcp45_2_IP,grep("SRE", names(projections.2050.rcp45_2_IP)))
projections.SRE.mean.2050.rcp45_IP <- mean(projections.SRE.all.2050.rcp45_IP)/10

projections.MARS.all.2050.rcp45_IP <- subset(projections.2050.rcp45_2_IP,grep("MARS", names(projections.2050.rcp45_2_IP)))
projections.MARS.mean.2050.rcp45_IP <- mean(projections.MARS.all.2050.rcp45_IP)/10

projections.FDA.all.2050.rcp45_IP <- subset(projections.2050.rcp45_2_IP,grep("FDA", names(projections.2050.rcp45_2_IP)))
projections.FDA.mean.2050.rcp45_IP <- mean(projections.FDA.all.2050.rcp45_IP)/10

projections.MAXENT.all.2050.rcp45_IP <- subset(projections.2050.rcp45_2_IP,grep("MAXENT.Phillips", names(projections.2050.rcp45_2_IP)))
projections.MAXENT.mean.2050.rcp45_IP <- mean(projections.MAXENT.all.2050.rcp45_IP)/10



#MI
projections.RF.all.2050.rcp45_MI <- subset(projections.2050.rcp45_1_MI, grep("RF", names(projections.2050.rcp45_1_MI)))
projections.RF.mean.2050.rcp45_MI <- mean(projections.RF.all.2050.rcp45_MI)/10

projections.GBM.all.2050.rcp45_MI <-subset(projections.2050.rcp45_1_MI, grep("GBM", names(projections.2050.rcp45_1_MI)))
projections.GBM.mean.2050.rcp45_MI <- mean(projections.GBM.all.2050.rcp45_MI)/10

projections.CTA.all.2050.rcp45_MI <-subset(projections.2050.rcp45_1_MI,grep("CTA", names(projections.2050.rcp45_1_MI)))
projections.CTA.mean.2050.rcp45_MI <- mean(projections.CTA.all.2050.rcp45_MI)/10

projections.GLM.all.2050.rcp45_MI <-subset(projections.2050.rcp45_2_MI,grep("GLM", names(projections.2050.rcp45_2_MI)))
projections.GLM.mean.2050.rcp45_MI <- mean(projections.GLM.all.2050.rcp45_MI)/10

projections.GAM.all.2050.rcp45_MI <-subset(projections.2050.rcp45_2_MI,grep("GAM", names(projections.2050.rcp45_2_MI)))
projections.GAM.mean.2050.rcp45_MI <- mean(projections.GAM.all.2050.rcp45_MI)/10

projections.ANN.all.2050.rcp45_MI <- subset(projections.2050.rcp45_2_MI,grep("ANN", names(projections.2050.rcp45_2_MI)))
projections.ANN.mean.2050.rcp45_MI <- mean(projections.ANN.all.2050.rcp45_MI)/10

projections.SRE.all.2050.rcp45_MI <- subset(projections.2050.rcp45_2_MI,grep("SRE", names(projections.2050.rcp45_2_MI)))
projections.SRE.mean.2050.rcp45_MI <- mean(projections.SRE.all.2050.rcp45_MI)/10

projections.MARS.all.2050.rcp45_MI <- subset(projections.2050.rcp45_2_MI,grep("MARS", names(projections.2050.rcp45_2_MI)))
projections.MARS.mean.2050.rcp45_MI <- mean(projections.MARS.all.2050.rcp45_MI)/10

projections.FDA.all.2050.rcp45_MI <- subset(projections.2050.rcp45_2_MI,grep("FDA", names(projections.2050.rcp45_2_MI)))
projections.FDA.mean.2050.rcp45_MI <- mean(projections.FDA.all.2050.rcp45_MI)/10

projections.MAXENT.all.2050.rcp45_MI <- subset(projections.2050.rcp45_2_MI,grep("MAXENT.Phillips", names(projections.2050.rcp45_2_MI)))
projections.MAXENT.mean.2050.rcp45_MI <- mean(projections.MAXENT.all.2050.rcp45_MI)/10



#MR
projections.RF.all.2050.rcp45_MR <- subset(projections.2050.rcp45_1_MR, grep("RF", names(projections.2050.rcp45_1_MR)))
projections.RF.mean.2050.rcp45_MR <- mean(projections.RF.all.2050.rcp45_MR)/10

projections.GBM.all.2050.rcp45_MR <-subset(projections.2050.rcp45_1_MR, grep("GBM", names(projections.2050.rcp45_1_MR)))
projections.GBM.mean.2050.rcp45_MR <- mean(projections.GBM.all.2050.rcp45_MR)/10

projections.CTA.all.2050.rcp45_MR <-subset(projections.2050.rcp45_1_MR,grep("CTA", names(projections.2050.rcp45_1_MR)))
projections.CTA.mean.2050.rcp45_MR <- mean(projections.CTA.all.2050.rcp45_MR)/10

projections.GLM.all.2050.rcp45_MR <-subset(projections.2050.rcp45_2_MR,grep("GLM", names(projections.2050.rcp45_2_MR)))
projections.GLM.mean.2050.rcp45_MR <- mean(projections.GLM.all.2050.rcp45_MR)/10

projections.GAM.all.2050.rcp45_MR <-subset(projections.2050.rcp45_2_MR,grep("GAM", names(projections.2050.rcp45_2_MR)))
projections.GAM.mean.2050.rcp45_MR <- mean(projections.GAM.all.2050.rcp45_MR)/10

projections.ANN.all.2050.rcp45_MR <- subset(projections.2050.rcp45_2_MR,grep("ANN", names(projections.2050.rcp45_2_MR)))
projections.ANN.mean.2050.rcp45_MR <- mean(projections.ANN.all.2050.rcp45_MR)/10

projections.SRE.all.2050.rcp45_MR <- subset(projections.2050.rcp45_2_MR,grep("SRE", names(projections.2050.rcp45_2_MR)))
projections.SRE.mean.2050.rcp45_MR <- mean(projections.SRE.all.2050.rcp45_MR)/10

projections.MARS.all.2050.rcp45_MR <- subset(projections.2050.rcp45_2_MR,grep("MARS", names(projections.2050.rcp45_2_MR)))
projections.MARS.mean.2050.rcp45_MR <- mean(projections.MARS.all.2050.rcp45_MR)/10

projections.FDA.all.2050.rcp45_MR <- subset(projections.2050.rcp45_2_MR,grep("FDA", names(projections.2050.rcp45_2_MR)))
projections.FDA.mean.2050.rcp45_MR <- mean(projections.FDA.all.2050.rcp45_MR)/10

projections.MAXENT.all.2050.rcp45_MR <- subset(projections.2050.rcp45_2_MR,grep("MAXENT.Phillips", names(projections.2050.rcp45_2_MR)))
projections.MAXENT.mean.2050.rcp45_MR <- mean(projections.MAXENT.all.2050.rcp45_MR)/10



#MC
projections.RF.all.2050.rcp45_MC <- subset(projections.2050.rcp45_1_MC, grep("RF", names(projections.2050.rcp45_1_MC)))
projections.RF.mean.2050.rcp45_MC <- mean(projections.RF.all.2050.rcp45_MC)/10

projections.GBM.all.2050.rcp45_MC <-subset(projections.2050.rcp45_1_MC, grep("GBM", names(projections.2050.rcp45_1_MC)))
projections.GBM.mean.2050.rcp45_MC <- mean(projections.GBM.all.2050.rcp45_MC)/10

projections.CTA.all.2050.rcp45_MC <-subset(projections.2050.rcp45_1_MC,grep("CTA", names(projections.2050.rcp45_1_MC)))
projections.CTA.mean.2050.rcp45_MC <- mean(projections.CTA.all.2050.rcp45_MC)/10

projections.GLM.all.2050.rcp45_MC <-subset(projections.2050.rcp45_2_MC,grep("GLM", names(projections.2050.rcp45_2_MC)))
projections.GLM.mean.2050.rcp45_MC <- mean(projections.GLM.all.2050.rcp45_MC)/10

projections.GAM.all.2050.rcp45_MC <-subset(projections.2050.rcp45_2_MC,grep("GAM", names(projections.2050.rcp45_2_MC)))
projections.GAM.mean.2050.rcp45_MC <- mean(projections.GAM.all.2050.rcp45_MC)/10

projections.ANN.all.2050.rcp45_MC <- subset(projections.2050.rcp45_2_MC,grep("ANN", names(projections.2050.rcp45_2_MC)))
projections.ANN.mean.2050.rcp45_MC <- mean(projections.ANN.all.2050.rcp45_MC)/10

projections.SRE.all.2050.rcp45_MC <- subset(projections.2050.rcp45_2_MC,grep("SRE", names(projections.2050.rcp45_2_MC)))
projections.SRE.mean.2050.rcp45_MC <- mean(projections.SRE.all.2050.rcp45_MC)/10

projections.MARS.all.2050.rcp45_MC <- subset(projections.2050.rcp45_2_MC,grep("MARS", names(projections.2050.rcp45_2_MC)))
projections.MARS.mean.2050.rcp45_MC <- mean(projections.MARS.all.2050.rcp45_MC)/10

projections.FDA.all.2050.rcp45_MC <- subset(projections.2050.rcp45_2_MC,grep("FDA", names(projections.2050.rcp45_2_MC)))
projections.FDA.mean.2050.rcp45_MC <- mean(projections.FDA.all.2050.rcp45_MC)/10

projections.MAXENT.all.2050.rcp45_MC <- subset(projections.2050.rcp45_2_MC,grep("MAXENT.Phillips", names(projections.2050.rcp45_2_MC)))
projections.MAXENT.mean.2050.rcp45_MC <- mean(projections.MAXENT.all.2050.rcp45_MC)/10



#MP
projections.RF.all.2050.rcp45_MP <- subset(projections.2050.rcp45_1_MP, grep("RF", names(projections.2050.rcp45_1_MP)))
projections.RF.mean.2050.rcp45_MP <- mean(projections.RF.all.2050.rcp45_MP)/10

projections.GBM.all.2050.rcp45_MP <-subset(projections.2050.rcp45_1_MP, grep("GBM", names(projections.2050.rcp45_1_MP)))
projections.GBM.mean.2050.rcp45_MP <- mean(projections.GBM.all.2050.rcp45_MP)/10

projections.CTA.all.2050.rcp45_MP <-subset(projections.2050.rcp45_1_MP,grep("CTA", names(projections.2050.rcp45_1_MP)))
projections.CTA.mean.2050.rcp45_MP <- mean(projections.CTA.all.2050.rcp45_MP)/10

projections.GLM.all.2050.rcp45_MP <-subset(projections.2050.rcp45_2_MP,grep("GLM", names(projections.2050.rcp45_2_MP)))
projections.GLM.mean.2050.rcp45_MP <- mean(projections.GLM.all.2050.rcp45_MP)/10

projections.GAM.all.2050.rcp45_MP <-subset(projections.2050.rcp45_2_MP,grep("GAM", names(projections.2050.rcp45_2_MP)))
projections.GAM.mean.2050.rcp45_MP <- mean(projections.GAM.all.2050.rcp45_MP)/10

projections.ANN.all.2050.rcp45_MP <- subset(projections.2050.rcp45_2_MP,grep("ANN", names(projections.2050.rcp45_2_MP)))
projections.ANN.mean.2050.rcp45_MP <- mean(projections.ANN.all.2050.rcp45_MP)/10

projections.SRE.all.2050.rcp45_MP <- subset(projections.2050.rcp45_2_MP,grep("SRE", names(projections.2050.rcp45_2_MP)))
projections.SRE.mean.2050.rcp45_MP <- mean(projections.SRE.all.2050.rcp45_MP)/10

projections.MARS.all.2050.rcp45_MP <- subset(projections.2050.rcp45_2_MP,grep("MARS", names(projections.2050.rcp45_2_MP)))
projections.MARS.mean.2050.rcp45_MP <- mean(projections.MARS.all.2050.rcp45_MP)/10

projections.FDA.all.2050.rcp45_MP <- subset(projections.2050.rcp45_2_MP,grep("FDA", names(projections.2050.rcp45_2_MP)))
projections.FDA.mean.2050.rcp45_MP <- mean(projections.FDA.all.2050.rcp45_MP)/10

projections.MAXENT.all.2050.rcp45_MP <- subset(projections.2050.rcp45_2_MP,grep("MAXENT.Phillips", names(projections.2050.rcp45_2_MP)))
projections.MAXENT.mean.2050.rcp45_MP <- mean(projections.MAXENT.all.2050.rcp45_MP)/10



#MG
projections.RF.all.2050.rcp45_MG <- subset(projections.2050.rcp45_1_MG, grep("RF", names(projections.2050.rcp45_1_MG)))
projections.RF.mean.2050.rcp45_MG <- mean(projections.RF.all.2050.rcp45_MG)/10

projections.GBM.all.2050.rcp45_MG <-subset(projections.2050.rcp45_1_MG, grep("GBM", names(projections.2050.rcp45_1_MG)))
projections.GBM.mean.2050.rcp45_MG <- mean(projections.GBM.all.2050.rcp45_MG)/10

projections.CTA.all.2050.rcp45_MG <-subset(projections.2050.rcp45_1_MG,grep("CTA", names(projections.2050.rcp45_1_MG)))
projections.CTA.mean.2050.rcp45_MG <- mean(projections.CTA.all.2050.rcp45_MG)/10

projections.GLM.all.2050.rcp45_MG <-subset(projections.2050.rcp45_2_MG,grep("GLM", names(projections.2050.rcp45_2_MG)))
projections.GLM.mean.2050.rcp45_MG <- mean(projections.GLM.all.2050.rcp45_MG)/10

projections.GAM.all.2050.rcp45_MG <-subset(projections.2050.rcp45_2_MG,grep("GAM", names(projections.2050.rcp45_2_MG)))
projections.GAM.mean.2050.rcp45_MG <- mean(projections.GAM.all.2050.rcp45_MG)/10

projections.ANN.all.2050.rcp45_MG <- subset(projections.2050.rcp45_2_MG,grep("ANN", names(projections.2050.rcp45_2_MG)))
projections.ANN.mean.2050.rcp45_MG <- mean(projections.ANN.all.2050.rcp45_MG)/10

projections.SRE.all.2050.rcp45_MG <- subset(projections.2050.rcp45_2_MG,grep("SRE", names(projections.2050.rcp45_2_MG)))
projections.SRE.mean.2050.rcp45_MG <- mean(projections.SRE.all.2050.rcp45_MG)/10

projections.MARS.all.2050.rcp45_MG <- subset(projections.2050.rcp45_2_MG,grep("MARS", names(projections.2050.rcp45_2_MG)))
projections.MARS.mean.2050.rcp45_MG <- mean(projections.MARS.all.2050.rcp45_MG)/10

projections.FDA.all.2050.rcp45_MG <- subset(projections.2050.rcp45_2_MG,grep("FDA", names(projections.2050.rcp45_2_MG)))
projections.FDA.mean.2050.rcp45_MG <- mean(projections.FDA.all.2050.rcp45_MG)/10

projections.MAXENT.all.2050.rcp45_MG <- subset(projections.2050.rcp45_2_MG,grep("MAXENT.Phillips", names(projections.2050.rcp45_2_MG)))
projections.MAXENT.mean.2050.rcp45_MG <- mean(projections.MAXENT.all.2050.rcp45_MG)/10



#NO
projections.RF.all.2050.rcp45_NO <- subset(projections.2050.rcp45_1_NO, grep("RF", names(projections.2050.rcp45_1_NO)))
projections.RF.mean.2050.rcp45_NO <- mean(projections.RF.all.2050.rcp45_NO)/10

projections.GBM.all.2050.rcp45_NO <-subset(projections.2050.rcp45_1_NO, grep("GBM", names(projections.2050.rcp45_1_NO)))
projections.GBM.mean.2050.rcp45_NO <- mean(projections.GBM.all.2050.rcp45_NO)/10

projections.CTA.all.2050.rcp45_NO <-subset(projections.2050.rcp45_1_NO,grep("CTA", names(projections.2050.rcp45_1_NO)))
projections.CTA.mean.2050.rcp45_NO <- mean(projections.CTA.all.2050.rcp45_NO)/10

projections.GLM.all.2050.rcp45_NO <-subset(projections.2050.rcp45_2_NO,grep("GLM", names(projections.2050.rcp45_2_NO)))
projections.GLM.mean.2050.rcp45_NO <- mean(projections.GLM.all.2050.rcp45_NO)/10

projections.GAM.all.2050.rcp45_NO <-subset(projections.2050.rcp45_2_NO,grep("GAM", names(projections.2050.rcp45_2_NO)))
projections.GAM.mean.2050.rcp45_NO <- mean(projections.GAM.all.2050.rcp45_NO)/10

projections.ANN.all.2050.rcp45_NO <- subset(projections.2050.rcp45_2_NO,grep("ANN", names(projections.2050.rcp45_2_NO)))
projections.ANN.mean.2050.rcp45_NO <- mean(projections.ANN.all.2050.rcp45_NO)/10

projections.SRE.all.2050.rcp45_NO <- subset(projections.2050.rcp45_2_NO,grep("SRE", names(projections.2050.rcp45_2_NO)))
projections.SRE.mean.2050.rcp45_NO <- mean(projections.SRE.all.2050.rcp45_NO)/10

projections.MARS.all.2050.rcp45_NO <- subset(projections.2050.rcp45_2_NO,grep("MARS", names(projections.2050.rcp45_2_NO)))
projections.MARS.mean.2050.rcp45_NO <- mean(projections.MARS.all.2050.rcp45_NO)/10

projections.FDA.all.2050.rcp45_NO <- subset(projections.2050.rcp45_2_NO,grep("FDA", names(projections.2050.rcp45_2_NO)))
projections.FDA.mean.2050.rcp45_NO <- mean(projections.FDA.all.2050.rcp45_NO)/10

projections.MAXENT.all.2050.rcp45_NO <- subset(projections.2050.rcp45_2_NO,grep("MAXENT.Phillips", names(projections.2050.rcp45_2_NO)))
projections.MAXENT.mean.2050.rcp45_NO <- mean(projections.MAXENT.all.2050.rcp45_NO)/10




##############################################################
####### AVERAGE ENSEMBLE AMONG SELECTED ALGORITHMS ###########
##############################################################
#GCMs: AC, BC, CC, CN, GF***, GS, HD, HG, HE, IN, IP, MI, MR, MC, MP, MG, NO


#AC
projections.all.mean.2050.rcp45_AC <- mean(projections.RF.mean.2050.rcp45_AC + projections.GBM.mean.2050.rcp45_AC +
	projections.CTA.mean.2050.rcp45_AC + projections.GLM.mean.2050.rcp45_AC + projections.GAM.mean.2050.rcp45_AC +
	projections.ANN.mean.2050.rcp45_AC + projections.SRE.mean.2050.rcp45_AC +
	projections.MARS.mean.2050.rcp45_AC + projections.FDA.mean.2050.rcp45_AC + projections.MAXENT.mean.2050.rcp45_AC)
writeRaster(projections.all.mean.2050.rcp45_AC, filename="Future Climate - 2050_rcp4.5_AC.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp45_AC, filename="Future Climate - 2050_rcp4.5_AC.asc", format="ascii")

#BC
projections.all.mean.2050.rcp45_BC <- mean(projections.RF.mean.2050.rcp45_BC + projections.GBM.mean.2050.rcp45_BC +
	projections.CTA.mean.2050.rcp45_BC + projections.GLM.mean.2050.rcp45_BC + projections.GAM.mean.2050.rcp45_BC +
	projections.ANN.mean.2050.rcp45_BC + projections.SRE.mean.2050.rcp45_BC +
	projections.MARS.mean.2050.rcp45_BC + projections.FDA.mean.2050.rcp45_BC + projections.MAXENT.mean.2050.rcp45_BC)
writeRaster(projections.all.mean.2050.rcp45_BC, filename="Future Climate - 2050_rcp4.5_BC.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp45_BC, filename="Future Climate - 2050_rcp4.5_BC.asc", format="ascii")

#CC
projections.all.mean.2050.rcp45_CC <- mean(projections.RF.mean.2050.rcp45_CC + projections.GBM.mean.2050.rcp45_CC +
	projections.CTA.mean.2050.rcp45_CC + projections.GLM.mean.2050.rcp45_CC + projections.GAM.mean.2050.rcp45_CC +
	projections.ANN.mean.2050.rcp45_CC + projections.SRE.mean.2050.rcp45_CC +
	projections.MARS.mean.2050.rcp45_CC + projections.FDA.mean.2050.rcp45_CC + projections.MAXENT.mean.2050.rcp45_CC)
writeRaster(projections.all.mean.2050.rcp45_CC, filename="Future Climate - 2050_rcp4.5_CC.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp45_CC, filename="Future Climate - 2050_rcp4.5_CC.asc", format="ascii")

#CN
projections.all.mean.2050.rcp45_CN <- mean(projections.RF.mean.2050.rcp45_CN + projections.GBM.mean.2050.rcp45_CN +
	projections.CTA.mean.2050.rcp45_CN + projections.GLM.mean.2050.rcp45_CN + projections.GAM.mean.2050.rcp45_CN +
	projections.ANN.mean.2050.rcp45_CN + projections.SRE.mean.2050.rcp45_CN +
	projections.MARS.mean.2050.rcp45_CN + projections.FDA.mean.2050.rcp45_CN + projections.MAXENT.mean.2050.rcp45_CN)
writeRaster(projections.all.mean.2050.rcp45_CN, filename="Future Climate - 2050_rcp4.5_CN.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp45_CN, filename="Future Climate - 2050_rcp4.5_CN.asc", format="ascii")

#GS
projections.all.mean.2050.rcp45_GS <- mean(projections.RF.mean.2050.rcp45_GS + projections.GBM.mean.2050.rcp45_GS +
	projections.CTA.mean.2050.rcp45_GS + projections.GLM.mean.2050.rcp45_GS + projections.GAM.mean.2050.rcp45_GS +
	projections.ANN.mean.2050.rcp45_GS + projections.SRE.mean.2050.rcp45_GS +
	projections.MARS.mean.2050.rcp45_GS + projections.FDA.mean.2050.rcp45_GS + projections.MAXENT.mean.2050.rcp45_GS)
writeRaster(projections.all.mean.2050.rcp45_GS, filename="Future Climate - 2050_rcp4.5_GS.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp45_GS, filename="Future Climate - 2050_rcp4.5_GS.asc", format="ascii")

#HD
projections.all.mean.2050.rcp45_HD <- mean(projections.RF.mean.2050.rcp45_HD + projections.GBM.mean.2050.rcp45_HD +
	projections.CTA.mean.2050.rcp45_HD + projections.GLM.mean.2050.rcp45_HD + projections.GAM.mean.2050.rcp45_HD +
	projections.ANN.mean.2050.rcp45_HD + projections.SRE.mean.2050.rcp45_HD +
	projections.MARS.mean.2050.rcp45_HD + projections.FDA.mean.2050.rcp45_HD + projections.MAXENT.mean.2050.rcp45_HD)
writeRaster(projections.all.mean.2050.rcp45_HD, filename="Future Climate - 2050_rcp4.5_HD.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp45_HD, filename="Future Climate - 2050_rcp4.5_HD.asc", format="ascii")

#HG
projections.all.mean.2050.rcp45_HG <- mean(projections.RF.mean.2050.rcp45_HG + projections.GBM.mean.2050.rcp45_HG +
	projections.CTA.mean.2050.rcp45_HG + projections.GLM.mean.2050.rcp45_HG + projections.GAM.mean.2050.rcp45_HG +
	projections.ANN.mean.2050.rcp45_HG + projections.SRE.mean.2050.rcp45_HG +
	projections.MARS.mean.2050.rcp45_HG + projections.FDA.mean.2050.rcp45_HG + projections.MAXENT.mean.2050.rcp45_HG)
writeRaster(projections.all.mean.2050.rcp45_HG, filename="Future Climate - 2050_rcp4.5_HG.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp45_HG, filename="Future Climate - 2050_rcp4.5_HG.asc", format="ascii")

#HE
projections.all.mean.2050.rcp45_HE <- mean(projections.RF.mean.2050.rcp45_HE + projections.GBM.mean.2050.rcp45_HE +
	projections.CTA.mean.2050.rcp45_HE + projections.GLM.mean.2050.rcp45_HE + projections.GAM.mean.2050.rcp45_HE +
	projections.ANN.mean.2050.rcp45_HE + projections.SRE.mean.2050.rcp45_HE +
	projections.MARS.mean.2050.rcp45_HE + projections.FDA.mean.2050.rcp45_HE + projections.MAXENT.mean.2050.rcp45_HE)
writeRaster(projections.all.mean.2050.rcp45_HE, filename="Future Climate - 2050_rcp4.5_HE.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp45_HE, filename="Future Climate - 2050_rcp4.5_HE.asc", format="ascii")

#IN
projections.all.mean.2050.rcp45_IN <- mean(projections.RF.mean.2050.rcp45_IN + projections.GBM.mean.2050.rcp45_IN +
	projections.CTA.mean.2050.rcp45_IN + projections.GLM.mean.2050.rcp45_IN + projections.GAM.mean.2050.rcp45_IN +
	projections.ANN.mean.2050.rcp45_IN + projections.SRE.mean.2050.rcp45_IN +
	projections.MARS.mean.2050.rcp45_IN + projections.FDA.mean.2050.rcp45_IN + projections.MAXENT.mean.2050.rcp45_IN)
writeRaster(projections.all.mean.2050.rcp45_IN, filename="Future Climate - 2050_rcp4.5_IN.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp45_IN, filename="Future Climate - 2050_rcp4.5_IN.asc", format="ascii")

#IP
projections.all.mean.2050.rcp45_IP <- mean(projections.RF.mean.2050.rcp45_IP + projections.GBM.mean.2050.rcp45_IP +
	projections.CTA.mean.2050.rcp45_IP + projections.GLM.mean.2050.rcp45_IP + projections.GAM.mean.2050.rcp45_IP +
	projections.ANN.mean.2050.rcp45_IP + projections.SRE.mean.2050.rcp45_IP +
	projections.MARS.mean.2050.rcp45_IP + projections.FDA.mean.2050.rcp45_IP + projections.MAXENT.mean.2050.rcp45_IP)
writeRaster(projections.all.mean.2050.rcp45_IP, filename="Future Climate - 2050_rcp4.5_IP.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp45_IP, filename="Future Climate - 2050_rcp4.5_IP.asc", format="ascii")

#MI
projections.all.mean.2050.rcp45_MI <- mean(projections.RF.mean.2050.rcp45_MI + projections.GBM.mean.2050.rcp45_MI +
	projections.CTA.mean.2050.rcp45_MI + projections.GLM.mean.2050.rcp45_MI + projections.GAM.mean.2050.rcp45_MI +
	projections.ANN.mean.2050.rcp45_MI + projections.SRE.mean.2050.rcp45_MI +
	projections.MARS.mean.2050.rcp45_MI + projections.FDA.mean.2050.rcp45_MI + projections.MAXENT.mean.2050.rcp45_MI)
writeRaster(projections.all.mean.2050.rcp45_MI, filename="Future Climate - 2050_rcp4.5_MI.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp45_MI, filename="Future Climate - 2050_rcp4.5_MI.asc", format="ascii")

#MR
projections.all.mean.2050.rcp45_MR <- mean(projections.RF.mean.2050.rcp45_MR + projections.GBM.mean.2050.rcp45_MR +
	projections.CTA.mean.2050.rcp45_MR + projections.GLM.mean.2050.rcp45_MR + projections.GAM.mean.2050.rcp45_MR +
	projections.ANN.mean.2050.rcp45_MR + projections.SRE.mean.2050.rcp45_MR +
	projections.MARS.mean.2050.rcp45_MR + projections.FDA.mean.2050.rcp45_MR + projections.MAXENT.mean.2050.rcp45_MR)
writeRaster(projections.all.mean.2050.rcp45_MR, filename="Future Climate - 2050_rcp4.5_MR.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp45_MR, filename="Future Climate - 2050_rcp4.5_MR.asc", format="ascii")

#MC
projections.all.mean.2050.rcp45_MC <- mean(projections.RF.mean.2050.rcp45_MC + projections.GBM.mean.2050.rcp45_MC +
	projections.CTA.mean.2050.rcp45_MC + projections.GLM.mean.2050.rcp45_MC + projections.GAM.mean.2050.rcp45_MC +
	projections.ANN.mean.2050.rcp45_MC + projections.SRE.mean.2050.rcp45_MC +
	projections.MARS.mean.2050.rcp45_MC + projections.FDA.mean.2050.rcp45_MC + projections.MAXENT.mean.2050.rcp45_MC)
writeRaster(projections.all.mean.2050.rcp45_MC, filename="Future Climate - 2050_rcp4.5_MC.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp45_MC, filename="Future Climate - 2050_rcp4.5_MC.asc", format="ascii")

#MP
projections.all.mean.2050.rcp45_MP <- mean(projections.RF.mean.2050.rcp45_MP + projections.GBM.mean.2050.rcp45_MP +
	projections.CTA.mean.2050.rcp45_MP + projections.GLM.mean.2050.rcp45_MP + projections.GAM.mean.2050.rcp45_MP +
	projections.ANN.mean.2050.rcp45_MP + projections.SRE.mean.2050.rcp45_MP +
	projections.MARS.mean.2050.rcp45_MP + projections.FDA.mean.2050.rcp45_MP + projections.MAXENT.mean.2050.rcp45_MP)
writeRaster(projections.all.mean.2050.rcp45_MP, filename="Future Climate - 2050_rcp4.5_MP.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp45_MP, filename="Future Climate - 2050_rcp4.5_MP.asc", format="ascii")

#MG
projections.all.mean.2050.rcp45_MG <- mean(projections.RF.mean.2050.rcp45_MG + projections.GBM.mean.2050.rcp45_MG +
	projections.CTA.mean.2050.rcp45_MG + projections.GLM.mean.2050.rcp45_MG + projections.GAM.mean.2050.rcp45_MG +
	projections.ANN.mean.2050.rcp45_MG + projections.SRE.mean.2050.rcp45_MG +
	projections.MARS.mean.2050.rcp45_MG + projections.FDA.mean.2050.rcp45_MG + projections.MAXENT.mean.2050.rcp45_MG)
writeRaster(projections.all.mean.2050.rcp45_MG, filename="Future Climate - 2050_rcp4.5_MG.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp45_MG, filename="Future Climate - 2050_rcp4.5_MG.asc", format="ascii")

#NO
projections.all.mean.2050.rcp45_NO <- mean(projections.RF.mean.2050.rcp45_NO + projections.GBM.mean.2050.rcp45_NO +
	projections.CTA.mean.2050.rcp45_NO + projections.GLM.mean.2050.rcp45_NO + projections.GAM.mean.2050.rcp45_NO +
	projections.ANN.mean.2050.rcp45_NO + projections.SRE.mean.2050.rcp45_NO +
	projections.MARS.mean.2050.rcp45_NO + projections.FDA.mean.2050.rcp45_NO + projections.MAXENT.mean.2050.rcp45_NO)
writeRaster(projections.all.mean.2050.rcp45_NO, filename="Future Climate - 2050_rcp4.5_NO.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp45_NO, filename="Future Climate - 2050_rcp4.5_NO.asc", format="ascii")


#Ensemble 2050 - rcp4.5
#GCMs: AC, BC, CC, CN, GS, HD, HG, HE, IN, IP, MI, MR, MC, MP, MG, NO
ensemble2050rcp4.5 <- mean(projections.all.mean.2050.rcp45_AC, projections.all.mean.2050.rcp45_BC,
					projections.all.mean.2050.rcp45_CC, projections.all.mean.2050.rcp45_CN,
					projections.all.mean.2050.rcp45_GS, projections.all.mean.2050.rcp45_HD,
					projections.all.mean.2050.rcp45_HG, projections.all.mean.2050.rcp45_HE,
					projections.all.mean.2050.rcp45_IN, projections.all.mean.2050.rcp45_IP,
					projections.all.mean.2050.rcp45_MI, projections.all.mean.2050.rcp45_MR,
					projections.all.mean.2050.rcp45_MC, projections.all.mean.2050.rcp45_MP,
					projections.all.mean.2050.rcp45_MG, projections.all.mean.2050.rcp45_NO)
writeRaster(ensemble2050rcp4.5, filename="Ensemble - Future Climate - 2050_rcp4.5.tif", format="GTiff")
writeRaster(ensemble2050rcp4.5, filename="Ensemble - Future Climate - 2050_rcp4.5.asc", format="ascii")
plot(ensemble2050rcp4.5)

#Converting to binary maps:
projections.binary_2050rcp4.5 <- BinaryTransformation(ensemble2050rcp4.5, th) #Replace 'th' by the threshold value you wish to apply
class(projections.binary_2050rcp4.5)
summary(values(projections.binary_2050rcp4.5))
writeRaster(projections.binary_2050rcp4.5, filename="Ensemble - Future Climate - 2050_rcp4.5_BINARY.tif", format="GTiff")
writeRaster(projections.binary_2050rcp4.5, filename="Ensemble - Future Climate - 2050_rcp4.5_BINARY.asc", formato="ascii")
plot(projections.binary_2050rcp4.5)





#####################################
## Model projection for the future ##
#####################################

## 2050 - pessimistic scenario (rcp 85)

spp.projections.2050.rcp85_1_AC <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_AC,
	proj.name = "2050.rcp85_1_AC",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp85_2_AC <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_AC,
	proj.name = "2050.rcp85_2_AC",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp85_1_BC <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_BC,
	proj.name = "2050.rcp85_1_BC",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp85_2_BC <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_BC,
	proj.name = "2050.rcp85_2_BC",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp85_1_CC <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_CC,
	proj.name = "2050.rcp85_1_CC",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp85_2_CC <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_CC,
	proj.name = "2050.rcp85_2_CC",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp85_1_CN <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_CN,
	proj.name = "2050.rcp85_1_CN",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp85_2_CN <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_CN,
	proj.name = "2050.rcp85_2_CN",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp85_1_GS <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_GS,
	proj.name = "2050.rcp48_1_GS",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp85_2_GS <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_GS,
	proj.name = "2050.rcp85_2_GS",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp85_1_HD <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_HD,
	proj.name = "2050.rcp85_1_HD",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp85_2_HD <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_HD,
	proj.name = "2050.rcp85_2_HD",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp85_1_HG <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_HG,
	proj.name = "2050.rcp85_1_HG",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp85_2_HG <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_HG,
	proj.name = "2050.rcp85_2_HG",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp85_1_HE <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_HE,
	proj.name = "2050.rcp85_1_HE",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp85_2_HE <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_HE,
	proj.name = "2050.rcp85_2_HE",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp85_1_IN <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_IN,
	proj.name = "2050.rcp85_1_IN",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp85_2_IN <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_IN,
	proj.name = "2050.rcp85_2_IN",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp85_1_IP <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_IP,
	proj.name = "2050.rcp85_1_IP",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp85_2_IP <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_IP,
	proj.name = "2050.rcp85_2_IP",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp85_1_MI <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_MI,
	proj.name = "2050.rcp85_1_MI",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp85_2_MI <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_MI,
	proj.name = "2050.rcp85_2_MI",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp85_1_MR <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_MR,
	proj.name = "2050.rcp85_1_MR",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp85_2_MR <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_MR,
	proj.name = "2050.rcp85_2_MR",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp85_1_MC <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_MC,
	proj.name = "2050.rcp85_1_MC",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp85_2_MC <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_MC,
	proj.name = "2050.rcp85_2_MC",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp85_1_MP <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_MP,
	proj.name = "2050.rcp85_1_MP",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp85_2_MP <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_MP,
	proj.name = "2050.rcp85_2_MP",
	selected.models = "all",
	output.format = ".grd")

spp.projections.2050.rcp85_1_MG <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_MG,
	proj.name = "2050.rcp85_1_MG",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp85_2_MG <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_MG,
	proj.name = "2050.rcp85_2_MG",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections.2050.rcp85_1_NO <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env50.85.selected_NO,
	proj.name = "2050.rcp85_1_NO",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")
spp.projections.2050.rcp85_2_NO <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env50.85.selected_NO,
	proj.name = "2050.rcp85_2_NO",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")




# Stack projections
projections.2050.rcp85_1_AC <- stack("C:/Models/Occurrence/proj_2050.rcp85_1_AC/proj_2050.rcp85_1_AC_Occurrence.grd")
names(projections.2050.rcp85_1_AC)
projections.2050.rcp85_2_AC <- stack("C:/Models/Occurrence/proj_2050.rcp85_2_AC/proj_2050.rcp85_2_AC_Occurrence.grd")
names(projections.2050.rcp85_2_AC)

projections.2050.rcp85_1_BC <- stack("C:/Models/Occurrence/proj_2050.rcp85_1_BC/proj_2050.rcp85_1_BC_Occurrence.grd")
names(projections.2050.rcp85_1_BC)
projections.2050.rcp85_2_BC <- stack("C:/Models/Occurrence/proj_2050.rcp85_2_BC/proj_2050.rcp85_2_BC_Occurrence.grd")
names(projections.2050.rcp85_2_BC)

projections.2050.rcp85_1_CC <- stack("C:/Models/Occurrence/proj_2050.rcp85_1_CC/proj_2050.rcp85_1_CC_Occurrence.grd")
names(projections.2050.rcp85_1_CC)
projections.2050.rcp85_2_CC <- stack("C:/Models/Occurrence/proj_2050.rcp85_2_CC/proj_2050.rcp85_2_CC_Occurrence.grd")
names(projections.2050.rcp85_2_CC)

projections.2050.rcp85_1_CN <- stack("C:/Models/Occurrence/proj_2050.rcp85_1_CN/proj_2050.rcp85_1_CN_Occurrence.grd")
names(projections.2050.rcp85_1_CN)
projections.2050.rcp85_2_CN <- stack("C:/Models/Occurrence/proj_2050.rcp85_2_CN/proj_2050.rcp85_2_CN_Occurrence.grd")
names(projections.2050.rcp85_2_CN)

projections.2050.rcp85_1_GS <- stack("C:/Models/Occurrence/proj_2050.rcp85_1_GS/proj_2050.rcp85_1_GS_Occurrence.grd")
names(projections.2050.rcp85_1_GS)
projections.2050.rcp85_2_GS <- stack("C:/Models/Occurrence/proj_2050.rcp85_2_GS/proj_2050.rcp85_2_GS_Occurrence.grd")
names(projections.2050.rcp85_2_GS)

projections.2050.rcp85_1_HD <- stack("C:/Models/Occurrence/proj_2050.rcp85_1_HD/proj_2050.rcp85_1_HD_Occurrence.grd")
names(projections.2050.rcp85_1_HD)
projections.2050.rcp85_2_HD <- stack("C:/Models/Occurrence/proj_2050.rcp85_2_HD/proj_2050.rcp85_2_HD_Occurrence.grd")
names(projections.2050.rcp85_2_HD)

projections.2050.rcp85_1_HG <- stack("C:/Models/Occurrence/proj_2050.rcp85_1_HG/proj_2050.rcp85_1_HG_Occurrence.grd")
names(projections.2050.rcp85_1_HG)
projections.2050.rcp85_2_HG <- stack("C:/Models/Occurrence/proj_2050.rcp85_2_HG/proj_2050.rcp85_2_HG_Occurrence.grd")
names(projections.2050.rcp85_2_HG)

projections.2050.rcp85_1_HE <- stack("C:/Models/Occurrence/proj_2050.rcp85_1_HE/proj_2050.rcp85_1_HE_Occurrence.grd")
names(projections.2050.rcp85_1_HE)
projections.2050.rcp85_2_HE <- stack("C:/Models/Occurrence/proj_2050.rcp85_2_HE/proj_2050.rcp85_2_HE_Occurrence.grd")
names(projections.2050.rcp85_2_HE)

projections.2050.rcp85_1_IN <- stack("C:/Models/Occurrence/proj_2050.rcp85_1_IN/proj_2050.rcp85_1_IN_Occurrence.grd")
names(projections.2050.rcp85_1_IN)
projections.2050.rcp85_2_IN <- stack("C:/Models/Occurrence/proj_2050.rcp85_2_IN/proj_2050.rcp85_2_IN_Occurrence.grd")
names(projections.2050.rcp85_2_IN)

projections.2050.rcp85_1_IP <- stack("C:/Models/Occurrence/proj_2050.rcp85_1_IP/proj_2050.rcp85_1_IP_Occurrence.grd")
names(projections.2050.rcp85_1_IP)
projections.2050.rcp85_2_IP <- stack("C:/Models/Occurrence/proj_2050.rcp85_2_IP/proj_2050.rcp85_2_IP_Occurrence.grd")
names(projections.2050.rcp85_2_IP)

projections.2050.rcp85_1_MI <- stack("C:/Models/Occurrence/proj_2050.rcp85_1_MI/proj_2050.rcp85_1_MI_Occurrence.grd")
names(projections.2050.rcp85_1_MI)
projections.2050.rcp85_2_MI <- stack("C:/Models/Occurrence/proj_2050.rcp85_2_MI/proj_2050.rcp85_2_MI_Occurrence.grd")
names(projections.2050.rcp85_2_MI)

projections.2050.rcp85_1_MR <- stack("C:/Models/Occurrence/proj_2050.rcp85_1_MR/proj_2050.rcp85_1_MR_Occurrence.grd")
names(projections.2050.rcp85_1_MR)
projections.2050.rcp85_2_MR <- stack("C:/Models/Occurrence/proj_2050.rcp85_2_MR/proj_2050.rcp85_2_MR_Occurrence.grd")
names(projections.2050.rcp85_2_MR)

projections.2050.rcp85_1_MC <- stack("C:/Models/Occurrence/proj_2050.rcp85_1_MC/proj_2050.rcp85_1_MC_Occurrence.grd")
names(projections.2050.rcp85_1_MC)
projections.2050.rcp85_2_MC <- stack("C:/Models/Occurrence/proj_2050.rcp85_2_MC/proj_2050.rcp85_2_MC_Occurrence.grd")
names(projections.2050.rcp85_2_MC)

projections.2050.rcp85_1_MP <- stack("C:/Models/Occurrence/proj_2050.rcp85_1_MP/proj_2050.rcp85_1_MP_Occurrence.grd")
names(projections.2050.rcp85_1_MP)
projections.2050.rcp85_2_MP <- stack("C:/Models/Occurrence/proj_2050.rcp85_2_MP/proj_2050.rcp85_2_MP_Occurrence.grd")
names(projections.2050.rcp85_2_MP)

projections.2050.rcp85_1_MG <- stack("C:/Models/Occurrence/proj_2050.rcp85_1_MG/proj_2050.rcp85_1_MG_Occurrence.grd")
names(projections.2050.rcp85_1_MG)
projections.2050.rcp85_2_MG <- stack("C:/Models/Occurrence/proj_2050.rcp85_2_MG/proj_2050.rcp85_2_MG_Occurrence.grd")
names(projections.2050.rcp85_2_MG)

projections.2050.rcp85_1_NO <- stack("C:/Models/Occurrence/proj_2050.rcp85_1_NO/proj_2050.rcp85_1_NO_Occurrence.grd")
names(projections.2050.rcp85_1_NO)
projections.2050.rcp85_2_NO <- stack("C:/Models/Occurrence/proj_2050.rcp85_2_NO/proj_2050.rcp85_2_NO_Occurrence.grd")
names(projections.2050.rcp85_2_NO)





### Average models for each algorithm:
#AC
projections.RF.all.2050.rcp85_AC <- subset(projections.2050.rcp85_1_AC, grep("RF", names(projections.2050.rcp85_1_AC)))
projections.RF.mean.2050.rcp85_AC <- mean(projections.RF.all.2050.rcp85_AC)/10

projections.GBM.all.2050.rcp85_AC <-subset(projections.2050.rcp85_1_AC, grep("GBM", names(projections.2050.rcp85_1_AC)))
projections.GBM.mean.2050.rcp85_AC <- mean(projections.GBM.all.2050.rcp85_AC)/10

projections.CTA.all.2050.rcp85_AC <-subset(projections.2050.rcp85_1_AC,grep("CTA", names(projections.2050.rcp85_1_AC)))
projections.CTA.mean.2050.rcp85_AC <- mean(projections.CTA.all.2050.rcp85_AC)/10

projections.GLM.all.2050.rcp85_AC <-subset(projections.2050.rcp85_2_AC,grep("GLM", names(projections.2050.rcp85_2_AC)))
projections.GLM.mean.2050.rcp85_AC <- mean(projections.GLM.all.2050.rcp85_AC)/10

projections.GAM.all.2050.rcp85_AC <-subset(projections.2050.rcp85_2_AC,grep("GAM", names(projections.2050.rcp85_2_AC)))
projections.GAM.mean.2050.rcp85_AC <- mean(projections.GAM.all.2050.rcp85_AC)/10

projections.ANN.all.2050.rcp85_AC <- subset(projections.2050.rcp85_2_AC,grep("ANN", names(projections.2050.rcp85_2_AC)))
projections.ANN.mean.2050.rcp85_AC <- mean(projections.ANN.all.2050.rcp85_AC)/10

projections.SRE.all.2050.rcp85_AC <- subset(projections.2050.rcp85_2_AC,grep("SRE", names(projections.2050.rcp85_2_AC)))
projections.SRE.mean.2050.rcp85_AC <- mean(projections.SRE.all.2050.rcp85_AC)/10

projections.MARS.all.2050.rcp85_AC <- subset(projections.2050.rcp85_2_AC,grep("MARS", names(projections.2050.rcp85_2_AC)))
projections.MARS.mean.2050.rcp85_AC <- mean(projections.MARS.all.2050.rcp85_AC)/10

projections.FDA.all.2050.rcp85_AC <- subset(projections.2050.rcp85_2_AC,grep("FDA", names(projections.2050.rcp85_2_AC)))
projections.FDA.mean.2050.rcp85_AC <- mean(projections.FDA.all.2050.rcp85_AC)/10

projections.MAXENT.all.2050.rcp85_AC <- subset(projections.2050.rcp85_2_AC,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_AC)))
projections.MAXENT.mean.2050.rcp85_AC <- mean(projections.MAXENT.all.2050.rcp85_AC)/10


#BC
projections.RF.all.2050.rcp85_BC <- subset(projections.2050.rcp85_1_BC, grep("RF", names(projections.2050.rcp85_1_BC)))
projections.RF.mean.2050.rcp85_BC <- mean(projections.RF.all.2050.rcp85_BC)/10

projections.GBM.all.2050.rcp85_BC <-subset(projections.2050.rcp85_1_BC, grep("GBM", names(projections.2050.rcp85_1_BC)))
projections.GBM.mean.2050.rcp85_BC <- mean(projections.GBM.all.2050.rcp85_BC)/10

projections.CTA.all.2050.rcp85_BC <-subset(projections.2050.rcp85_1_BC,grep("CTA", names(projections.2050.rcp85_1_BC)))
projections.CTA.mean.2050.rcp85_BC <- mean(projections.CTA.all.2050.rcp85_BC)/10

projections.GLM.all.2050.rcp85_BC <-subset(projections.2050.rcp85_2_BC,grep("GLM", names(projections.2050.rcp85_2_BC)))
projections.GLM.mean.2050.rcp85_BC <- mean(projections.GLM.all.2050.rcp85_BC)/10

projections.GAM.all.2050.rcp85_BC <-subset(projections.2050.rcp85_2_BC,grep("GAM", names(projections.2050.rcp85_2_BC)))
projections.GAM.mean.2050.rcp85_BC <- mean(projections.GAM.all.2050.rcp85_BC)/10

projections.ANN.all.2050.rcp85_BC <- subset(projections.2050.rcp85_2_BC,grep("ANN", names(projections.2050.rcp85_2_BC)))
projections.ANN.mean.2050.rcp85_BC <- mean(projections.ANN.all.2050.rcp85_BC)/10

projections.SRE.all.2050.rcp85_BC <- subset(projections.2050.rcp85_2_BC,grep("SRE", names(projections.2050.rcp85_2_BC)))
projections.SRE.mean.2050.rcp85_BC <- mean(projections.SRE.all.2050.rcp85_BC)/10

projections.MARS.all.2050.rcp85_BC <- subset(projections.2050.rcp85_2_BC,grep("MARS", names(projections.2050.rcp85_2_BC)))
projections.MARS.mean.2050.rcp85_BC <- mean(projections.MARS.all.2050.rcp85_BC)/10

projections.FDA.all.2050.rcp85_BC <- subset(projections.2050.rcp85_2_BC,grep("FDA", names(projections.2050.rcp85_2_BC)))
projections.FDA.mean.2050.rcp85_BC <- mean(projections.FDA.all.2050.rcp85_BC)/10

projections.MAXENT.all.2050.rcp85_BC <- subset(projections.2050.rcp85_2_BC,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_BC)))
projections.MAXENT.mean.2050.rcp85_BC <- mean(projections.MAXENT.all.2050.rcp85_BC)/10


#CC
projections.RF.all.2050.rcp85_CC <- subset(projections.2050.rcp85_1_CC, grep("RF", names(projections.2050.rcp85_1_CC)))
projections.RF.mean.2050.rcp85_CC <- mean(projections.RF.all.2050.rcp85_CC)/10

projections.GBM.all.2050.rcp85_CC <-subset(projections.2050.rcp85_1_CC, grep("GBM", names(projections.2050.rcp85_1_CC)))
projections.GBM.mean.2050.rcp85_CC <- mean(projections.GBM.all.2050.rcp85_CC)/10

projections.CTA.all.2050.rcp85_CC <-subset(projections.2050.rcp85_1_CC,grep("CTA", names(projections.2050.rcp85_1_CC)))
projections.CTA.mean.2050.rcp85_CC <- mean(projections.CTA.all.2050.rcp85_CC)/10

projections.GLM.all.2050.rcp85_CC <-subset(projections.2050.rcp85_2_CC,grep("GLM", names(projections.2050.rcp85_2_CC)))
projections.GLM.mean.2050.rcp85_CC <- mean(projections.GLM.all.2050.rcp85_CC)/10

projections.GAM.all.2050.rcp85_CC <-subset(projections.2050.rcp85_2_CC,grep("GAM", names(projections.2050.rcp85_2_CC)))
projections.GAM.mean.2050.rcp85_CC <- mean(projections.GAM.all.2050.rcp85_CC)/10

projections.ANN.all.2050.rcp85_CC <- subset(projections.2050.rcp85_2_CC,grep("ANN", names(projections.2050.rcp85_2_CC)))
projections.ANN.mean.2050.rcp85_CC <- mean(projections.ANN.all.2050.rcp85_CC)/10

projections.SRE.all.2050.rcp85_CC <- subset(projections.2050.rcp85_2_CC,grep("SRE", names(projections.2050.rcp85_2_CC)))
projections.SRE.mean.2050.rcp85_CC <- mean(projections.SRE.all.2050.rcp85_CC)/10

projections.MARS.all.2050.rcp85_CC <- subset(projections.2050.rcp85_2_CC,grep("MARS", names(projections.2050.rcp85_2_CC)))
projections.MARS.mean.2050.rcp85_CC <- mean(projections.MARS.all.2050.rcp85_CC)/10

projections.FDA.all.2050.rcp85_CC <- subset(projections.2050.rcp85_2_CC,grep("FDA", names(projections.2050.rcp85_2_CC)))
projections.FDA.mean.2050.rcp85_CC <- mean(projections.FDA.all.2050.rcp85_CC)/10

projections.MAXENT.all.2050.rcp85_CC <- subset(projections.2050.rcp85_2_CC,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_CC)))
projections.MAXENT.mean.2050.rcp85_CC <- mean(projections.MAXENT.all.2050.rcp85_CC)/10


#CN
projections.RF.all.2050.rcp85_CN <- subset(projections.2050.rcp85_1_CN, grep("RF", names(projections.2050.rcp85_1_CN)))
projections.RF.mean.2050.rcp85_CN <- mean(projections.RF.all.2050.rcp85_CN)/10

projections.GBM.all.2050.rcp85_CN <-subset(projections.2050.rcp85_1_CN, grep("GBM", names(projections.2050.rcp85_1_CN)))
projections.GBM.mean.2050.rcp85_CN <- mean(projections.GBM.all.2050.rcp85_CN)/10

projections.CTA.all.2050.rcp85_CN <-subset(projections.2050.rcp85_1_CN,grep("CTA", names(projections.2050.rcp85_1_CN)))
projections.CTA.mean.2050.rcp85_CN <- mean(projections.CTA.all.2050.rcp85_CN)/10

projections.GLM.all.2050.rcp85_CN <-subset(projections.2050.rcp85_2_CN,grep("GLM", names(projections.2050.rcp85_2_CN)))
projections.GLM.mean.2050.rcp85_CN <- mean(projections.GLM.all.2050.rcp85_CN)/10

projections.GAM.all.2050.rcp85_CN <-subset(projections.2050.rcp85_2_CN,grep("GAM", names(projections.2050.rcp85_2_CN)))
projections.GAM.mean.2050.rcp85_CN <- mean(projections.GAM.all.2050.rcp85_CN)/10

projections.ANN.all.2050.rcp85_CN <- subset(projections.2050.rcp85_2_CN,grep("ANN", names(projections.2050.rcp85_2_CN)))
projections.ANN.mean.2050.rcp85_CN <- mean(projections.ANN.all.2050.rcp85_CN)/10

projections.SRE.all.2050.rcp85_CN <- subset(projections.2050.rcp85_2_CN,grep("SRE", names(projections.2050.rcp85_2_CN)))
projections.SRE.mean.2050.rcp85_CN <- mean(projections.SRE.all.2050.rcp85_CN)/10

projections.MARS.all.2050.rcp85_CN <- subset(projections.2050.rcp85_2_CN,grep("MARS", names(projections.2050.rcp85_2_CN)))
projections.MARS.mean.2050.rcp85_CN <- mean(projections.MARS.all.2050.rcp85_CN)/10

projections.FDA.all.2050.rcp85_CN <- subset(projections.2050.rcp85_2_CN,grep("FDA", names(projections.2050.rcp85_2_CN)))
projections.FDA.mean.2050.rcp85_CN <- mean(projections.FDA.all.2050.rcp85_CN)/10

projections.MAXENT.all.2050.rcp85_CN <- subset(projections.2050.rcp85_2_CN,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_CN)))
projections.MAXENT.mean.2050.rcp85_CN <- mean(projections.MAXENT.all.2050.rcp85_CN)/10


#GS
projections.RF.all.2050.rcp85_GS <- subset(projections.2050.rcp85_1_GS, grep("RF", names(projections.2050.rcp85_1_GS)))
projections.RF.mean.2050.rcp85_GS <- mean(projections.RF.all.2050.rcp85_GS)/10

projections.GBM.all.2050.rcp85_GS <-subset(projections.2050.rcp85_1_GS, grep("GBM", names(projections.2050.rcp85_1_GS)))
projections.GBM.mean.2050.rcp85_GS <- mean(projections.GBM.all.2050.rcp85_GS)/10

projections.CTA.all.2050.rcp85_GS <-subset(projections.2050.rcp85_1_GS,grep("CTA", names(projections.2050.rcp85_1_GS)))
projections.CTA.mean.2050.rcp85_GS <- mean(projections.CTA.all.2050.rcp85_GS)/10

projections.GLM.all.2050.rcp85_GS <-subset(projections.2050.rcp85_2_GS,grep("GLM", names(projections.2050.rcp85_2_GS)))
projections.GLM.mean.2050.rcp85_GS <- mean(projections.GLM.all.2050.rcp85_GS)/10

projections.GAM.all.2050.rcp85_GS <-subset(projections.2050.rcp85_2_GS,grep("GAM", names(projections.2050.rcp85_2_GS)))
projections.GAM.mean.2050.rcp85_GS <- mean(projections.GAM.all.2050.rcp85_GS)/10

projections.ANN.all.2050.rcp85_GS <- subset(projections.2050.rcp85_2_GS,grep("ANN", names(projections.2050.rcp85_2_GS)))
projections.ANN.mean.2050.rcp85_GS <- mean(projections.ANN.all.2050.rcp85_GS)/10

projections.SRE.all.2050.rcp85_GS <- subset(projections.2050.rcp85_2_GS,grep("SRE", names(projections.2050.rcp85_2_GS)))
projections.SRE.mean.2050.rcp85_GS <- mean(projections.SRE.all.2050.rcp85_GS)/10

projections.MARS.all.2050.rcp85_GS <- subset(projections.2050.rcp85_2_GS,grep("MARS", names(projections.2050.rcp85_2_GS)))
projections.MARS.mean.2050.rcp85_GS <- mean(projections.MARS.all.2050.rcp85_GS)/10

projections.FDA.all.2050.rcp85_GS <- subset(projections.2050.rcp85_2_GS,grep("FDA", names(projections.2050.rcp85_2_GS)))
projections.FDA.mean.2050.rcp85_GS <- mean(projections.FDA.all.2050.rcp85_GS)/10

projections.MAXENT.all.2050.rcp85_GS <- subset(projections.2050.rcp85_2_GS,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_GS)))
projections.MAXENT.mean.2050.rcp85_GS <- mean(projections.MAXENT.all.2050.rcp85_GS)/10



#HD
projections.RF.all.2050.rcp85_HD <- subset(projections.2050.rcp85_1_HD, grep("RF", names(projections.2050.rcp85_1_HD)))
projections.RF.mean.2050.rcp85_HD <- mean(projections.RF.all.2050.rcp85_HD)/10

projections.GBM.all.2050.rcp85_HD <-subset(projections.2050.rcp85_1_HD, grep("GBM", names(projections.2050.rcp85_1_HD)))
projections.GBM.mean.2050.rcp85_HD <- mean(projections.GBM.all.2050.rcp85_HD)/10

projections.CTA.all.2050.rcp85_HD <-subset(projections.2050.rcp85_1_HD,grep("CTA", names(projections.2050.rcp85_1_HD)))
projections.CTA.mean.2050.rcp85_HD <- mean(projections.CTA.all.2050.rcp85_HD)/10

projections.GLM.all.2050.rcp85_HD <-subset(projections.2050.rcp85_2_HD,grep("GLM", names(projections.2050.rcp85_2_HD)))
projections.GLM.mean.2050.rcp85_HD <- mean(projections.GLM.all.2050.rcp85_HD)/10

projections.GAM.all.2050.rcp85_HD <-subset(projections.2050.rcp85_2_HD,grep("GAM", names(projections.2050.rcp85_2_HD)))
projections.GAM.mean.2050.rcp85_HD <- mean(projections.GAM.all.2050.rcp85_HD)/10

projections.ANN.all.2050.rcp85_HD <- subset(projections.2050.rcp85_2_HD,grep("ANN", names(projections.2050.rcp85_2_HD)))
projections.ANN.mean.2050.rcp85_HD <- mean(projections.ANN.all.2050.rcp85_HD)/10

projections.SRE.all.2050.rcp85_HD <- subset(projections.2050.rcp85_2_HD,grep("SRE", names(projections.2050.rcp85_2_HD)))
projections.SRE.mean.2050.rcp85_HD <- mean(projections.SRE.all.2050.rcp85_HD)/10

projections.MARS.all.2050.rcp85_HD <- subset(projections.2050.rcp85_2_HD,grep("MARS", names(projections.2050.rcp85_2_HD)))
projections.MARS.mean.2050.rcp85_HD <- mean(projections.MARS.all.2050.rcp85_HD)/10

projections.FDA.all.2050.rcp85_HD <- subset(projections.2050.rcp85_2_HD,grep("FDA", names(projections.2050.rcp85_2_HD)))
projections.FDA.mean.2050.rcp85_HD <- mean(projections.FDA.all.2050.rcp85_HD)/10

projections.MAXENT.all.2050.rcp85_HD <- subset(projections.2050.rcp85_2_HD,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_HD)))
projections.MAXENT.mean.2050.rcp85_HD <- mean(projections.MAXENT.all.2050.rcp85_HD)/10



#HG
projections.RF.all.2050.rcp85_HG <- subset(projections.2050.rcp85_1_HG, grep("RF", names(projections.2050.rcp85_1_HG)))
projections.RF.mean.2050.rcp85_HG <- mean(projections.RF.all.2050.rcp85_HG)/10

projections.GBM.all.2050.rcp85_HG <-subset(projections.2050.rcp85_1_HG, grep("GBM", names(projections.2050.rcp85_1_HG)))
projections.GBM.mean.2050.rcp85_HG <- mean(projections.GBM.all.2050.rcp85_HG)/10

projections.CTA.all.2050.rcp85_HG <-subset(projections.2050.rcp85_1_HG,grep("CTA", names(projections.2050.rcp85_1_HG)))
projections.CTA.mean.2050.rcp85_HG <- mean(projections.CTA.all.2050.rcp85_HG)/10

projections.GLM.all.2050.rcp85_HG <-subset(projections.2050.rcp85_2_HG,grep("GLM", names(projections.2050.rcp85_2_HG)))
projections.GLM.mean.2050.rcp85_HG <- mean(projections.GLM.all.2050.rcp85_HG)/10

projections.GAM.all.2050.rcp85_HG <-subset(projections.2050.rcp85_2_HG,grep("GAM", names(projections.2050.rcp85_2_HG)))
projections.GAM.mean.2050.rcp85_HG <- mean(projections.GAM.all.2050.rcp85_HG)/10

projections.ANN.all.2050.rcp85_HG <- subset(projections.2050.rcp85_2_HG,grep("ANN", names(projections.2050.rcp85_2_HG)))
projections.ANN.mean.2050.rcp85_HG <- mean(projections.ANN.all.2050.rcp85_HG)/10

projections.SRE.all.2050.rcp85_HG <- subset(projections.2050.rcp85_2_HG,grep("SRE", names(projections.2050.rcp85_2_HG)))
projections.SRE.mean.2050.rcp85_HG <- mean(projections.SRE.all.2050.rcp85_HG)/10

projections.MARS.all.2050.rcp85_HG <- subset(projections.2050.rcp85_2_HG,grep("MARS", names(projections.2050.rcp85_2_HG)))
projections.MARS.mean.2050.rcp85_HG <- mean(projections.MARS.all.2050.rcp85_HG)/10

projections.FDA.all.2050.rcp85_HG <- subset(projections.2050.rcp85_2_HG,grep("FDA", names(projections.2050.rcp85_2_HG)))
projections.FDA.mean.2050.rcp85_HG <- mean(projections.FDA.all.2050.rcp85_HG)/10

projections.MAXENT.all.2050.rcp85_HG <- subset(projections.2050.rcp85_2_HG,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_HG)))
projections.MAXENT.mean.2050.rcp85_HG <- mean(projections.MAXENT.all.2050.rcp85_HG)/10



#HE
projections.RF.all.2050.rcp85_HE <- subset(projections.2050.rcp85_1_HE, grep("RF", names(projections.2050.rcp85_1_HE)))
projections.RF.mean.2050.rcp85_HE <- mean(projections.RF.all.2050.rcp85_HE)/10

projections.GBM.all.2050.rcp85_HE <-subset(projections.2050.rcp85_1_HE, grep("GBM", names(projections.2050.rcp85_1_HE)))
projections.GBM.mean.2050.rcp85_HE <- mean(projections.GBM.all.2050.rcp85_HE)/10

projections.CTA.all.2050.rcp85_HE <-subset(projections.2050.rcp85_1_HE,grep("CTA", names(projections.2050.rcp85_1_HE)))
projections.CTA.mean.2050.rcp85_HE <- mean(projections.CTA.all.2050.rcp85_HE)/10

projections.GLM.all.2050.rcp85_HE <-subset(projections.2050.rcp85_2_HE,grep("GLM", names(projections.2050.rcp85_2_HE)))
projections.GLM.mean.2050.rcp85_HE <- mean(projections.GLM.all.2050.rcp85_HE)/10

projections.GAM.all.2050.rcp85_HE <-subset(projections.2050.rcp85_2_HE,grep("GAM", names(projections.2050.rcp85_2_HE)))
projections.GAM.mean.2050.rcp85_HE <- mean(projections.GAM.all.2050.rcp85_HE)/10

projections.ANN.all.2050.rcp85_HE <- subset(projections.2050.rcp85_2_HE,grep("ANN", names(projections.2050.rcp85_2_HE)))
projections.ANN.mean.2050.rcp85_HE <- mean(projections.ANN.all.2050.rcp85_HE)/10

projections.SRE.all.2050.rcp85_HE <- subset(projections.2050.rcp85_2_HE,grep("SRE", names(projections.2050.rcp85_2_HE)))
projections.SRE.mean.2050.rcp85_HE <- mean(projections.SRE.all.2050.rcp85_HE)/10

projections.MARS.all.2050.rcp85_HE <- subset(projections.2050.rcp85_2_HE,grep("MARS", names(projections.2050.rcp85_2_HE)))
projections.MARS.mean.2050.rcp85_HE <- mean(projections.MARS.all.2050.rcp85_HE)/10

projections.FDA.all.2050.rcp85_HE <- subset(projections.2050.rcp85_2_HE,grep("FDA", names(projections.2050.rcp85_2_HE)))
projections.FDA.mean.2050.rcp85_HE <- mean(projections.FDA.all.2050.rcp85_HE)/10

projections.MAXENT.all.2050.rcp85_HE <- subset(projections.2050.rcp85_2_HE,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_HE)))
projections.MAXENT.mean.2050.rcp85_HE <- mean(projections.MAXENT.all.2050.rcp85_HE)/10



#IN
projections.RF.all.2050.rcp85_IN <- subset(projections.2050.rcp85_1_IN, grep("RF", names(projections.2050.rcp85_1_IN)))
projections.RF.mean.2050.rcp85_IN <- mean(projections.RF.all.2050.rcp85_IN)/10

projections.GBM.all.2050.rcp85_IN <-subset(projections.2050.rcp85_1_IN, grep("GBM", names(projections.2050.rcp85_1_IN)))
projections.GBM.mean.2050.rcp85_IN <- mean(projections.GBM.all.2050.rcp85_IN)/10

projections.CTA.all.2050.rcp85_IN <-subset(projections.2050.rcp85_1_IN,grep("CTA", names(projections.2050.rcp85_1_IN)))
projections.CTA.mean.2050.rcp85_IN <- mean(projections.CTA.all.2050.rcp85_IN)/10

projections.GLM.all.2050.rcp85_IN <-subset(projections.2050.rcp85_2_IN,grep("GLM", names(projections.2050.rcp85_2_IN)))
projections.GLM.mean.2050.rcp85_IN <- mean(projections.GLM.all.2050.rcp85_IN)/10

projections.GAM.all.2050.rcp85_IN <-subset(projections.2050.rcp85_2_IN,grep("GAM", names(projections.2050.rcp85_2_IN)))
projections.GAM.mean.2050.rcp85_IN <- mean(projections.GAM.all.2050.rcp85_IN)/10

projections.ANN.all.2050.rcp85_IN <- subset(projections.2050.rcp85_2_IN,grep("ANN", names(projections.2050.rcp85_2_IN)))
projections.ANN.mean.2050.rcp85_IN <- mean(projections.ANN.all.2050.rcp85_IN)/10

projections.SRE.all.2050.rcp85_IN <- subset(projections.2050.rcp85_2_IN,grep("SRE", names(projections.2050.rcp85_2_IN)))
projections.SRE.mean.2050.rcp85_IN <- mean(projections.SRE.all.2050.rcp85_IN)/10

projections.MARS.all.2050.rcp85_IN <- subset(projections.2050.rcp85_2_IN,grep("MARS", names(projections.2050.rcp85_2_IN)))
projections.MARS.mean.2050.rcp85_IN <- mean(projections.MARS.all.2050.rcp85_IN)/10

projections.FDA.all.2050.rcp85_IN <- subset(projections.2050.rcp85_2_IN,grep("FDA", names(projections.2050.rcp85_2_IN)))
projections.FDA.mean.2050.rcp85_IN <- mean(projections.FDA.all.2050.rcp85_IN)/10

projections.MAXENT.all.2050.rcp85_IN <- subset(projections.2050.rcp85_2_IN,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_IN)))
projections.MAXENT.mean.2050.rcp85_IN <- mean(projections.MAXENT.all.2050.rcp85_IN)/10


#IP
projections.RF.all.2050.rcp85_IP <- subset(projections.2050.rcp85_1_IP, grep("RF", names(projections.2050.rcp85_1_IP)))
projections.RF.mean.2050.rcp85_IP <- mean(projections.RF.all.2050.rcp85_IP)/10

projections.GBM.all.2050.rcp85_IP <-subset(projections.2050.rcp85_1_IP, grep("GBM", names(projections.2050.rcp85_1_IP)))
projections.GBM.mean.2050.rcp85_IP <- mean(projections.GBM.all.2050.rcp85_IP)/10

projections.CTA.all.2050.rcp85_IP <-subset(projections.2050.rcp85_1_IP,grep("CTA", names(projections.2050.rcp85_1_IP)))
projections.CTA.mean.2050.rcp85_IP <- mean(projections.CTA.all.2050.rcp85_IP)/10

projections.GLM.all.2050.rcp85_IP <-subset(projections.2050.rcp85_2_IP,grep("GLM", names(projections.2050.rcp85_2_IP)))
projections.GLM.mean.2050.rcp85_IP <- mean(projections.GLM.all.2050.rcp85_IP)/10

projections.GAM.all.2050.rcp85_IP <-subset(projections.2050.rcp85_2_IP,grep("GAM", names(projections.2050.rcp85_2_IP)))
projections.GAM.mean.2050.rcp85_IP <- mean(projections.GAM.all.2050.rcp85_IP)/10

projections.ANN.all.2050.rcp85_IP <- subset(projections.2050.rcp85_2_IP,grep("ANN", names(projections.2050.rcp85_2_IP)))
projections.ANN.mean.2050.rcp85_IP <- mean(projections.ANN.all.2050.rcp85_IP)/10

projections.SRE.all.2050.rcp85_IP <- subset(projections.2050.rcp85_2_IP,grep("SRE", names(projections.2050.rcp85_2_IP)))
projections.SRE.mean.2050.rcp85_IP <- mean(projections.SRE.all.2050.rcp85_IP)/10

projections.MARS.all.2050.rcp85_IP <- subset(projections.2050.rcp85_2_IP,grep("MARS", names(projections.2050.rcp85_2_IP)))
projections.MARS.mean.2050.rcp85_IP <- mean(projections.MARS.all.2050.rcp85_IP)/10

projections.FDA.all.2050.rcp85_IP <- subset(projections.2050.rcp85_2_IP,grep("FDA", names(projections.2050.rcp85_2_IP)))
projections.FDA.mean.2050.rcp85_IP <- mean(projections.FDA.all.2050.rcp85_IP)/10

projections.MAXENT.all.2050.rcp85_IP <- subset(projections.2050.rcp85_2_IP,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_IP)))
projections.MAXENT.mean.2050.rcp85_IP <- mean(projections.MAXENT.all.2050.rcp85_IP)/10



#MI
projections.RF.all.2050.rcp85_MI <- subset(projections.2050.rcp85_1_MI, grep("RF", names(projections.2050.rcp85_1_MI)))
projections.RF.mean.2050.rcp85_MI <- mean(projections.RF.all.2050.rcp85_MI)/10

projections.GBM.all.2050.rcp85_MI <-subset(projections.2050.rcp85_1_MI, grep("GBM", names(projections.2050.rcp85_1_MI)))
projections.GBM.mean.2050.rcp85_MI <- mean(projections.GBM.all.2050.rcp85_MI)/10

projections.CTA.all.2050.rcp85_MI <-subset(projections.2050.rcp85_1_MI,grep("CTA", names(projections.2050.rcp85_1_MI)))
projections.CTA.mean.2050.rcp85_MI <- mean(projections.CTA.all.2050.rcp85_MI)/10

projections.GLM.all.2050.rcp85_MI <-subset(projections.2050.rcp85_2_MI,grep("GLM", names(projections.2050.rcp85_2_MI)))
projections.GLM.mean.2050.rcp85_MI <- mean(projections.GLM.all.2050.rcp85_MI)/10

projections.GAM.all.2050.rcp85_MI <-subset(projections.2050.rcp85_2_MI,grep("GAM", names(projections.2050.rcp85_2_MI)))
projections.GAM.mean.2050.rcp85_MI <- mean(projections.GAM.all.2050.rcp85_MI)/10

projections.ANN.all.2050.rcp85_MI <- subset(projections.2050.rcp85_2_MI,grep("ANN", names(projections.2050.rcp85_2_MI)))
projections.ANN.mean.2050.rcp85_MI <- mean(projections.ANN.all.2050.rcp85_MI)/10

projections.SRE.all.2050.rcp85_MI <- subset(projections.2050.rcp85_2_MI,grep("SRE", names(projections.2050.rcp85_2_MI)))
projections.SRE.mean.2050.rcp85_MI <- mean(projections.SRE.all.2050.rcp85_MI)/10

projections.MARS.all.2050.rcp85_MI <- subset(projections.2050.rcp85_2_MI,grep("MARS", names(projections.2050.rcp85_2_MI)))
projections.MARS.mean.2050.rcp85_MI <- mean(projections.MARS.all.2050.rcp85_MI)/10

projections.FDA.all.2050.rcp85_MI <- subset(projections.2050.rcp85_2_MI,grep("FDA", names(projections.2050.rcp85_2_MI)))
projections.FDA.mean.2050.rcp85_MI <- mean(projections.FDA.all.2050.rcp85_MI)/10

projections.MAXENT.all.2050.rcp85_MI <- subset(projections.2050.rcp85_2_MI,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_MI)))
projections.MAXENT.mean.2050.rcp85_MI <- mean(projections.MAXENT.all.2050.rcp85_MI)/10



#MR
projections.RF.all.2050.rcp85_MR <- subset(projections.2050.rcp85_1_MR, grep("RF", names(projections.2050.rcp85_1_MR)))
projections.RF.mean.2050.rcp85_MR <- mean(projections.RF.all.2050.rcp85_MR)/10

projections.GBM.all.2050.rcp85_MR <-subset(projections.2050.rcp85_1_MR, grep("GBM", names(projections.2050.rcp85_1_MR)))
projections.GBM.mean.2050.rcp85_MR <- mean(projections.GBM.all.2050.rcp85_MR)/10

projections.CTA.all.2050.rcp85_MR <-subset(projections.2050.rcp85_1_MR,grep("CTA", names(projections.2050.rcp85_1_MR)))
projections.CTA.mean.2050.rcp85_MR <- mean(projections.CTA.all.2050.rcp85_MR)/10

projections.GLM.all.2050.rcp85_MR <-subset(projections.2050.rcp85_2_MR,grep("GLM", names(projections.2050.rcp85_2_MR)))
projections.GLM.mean.2050.rcp85_MR <- mean(projections.GLM.all.2050.rcp85_MR)/10

projections.GAM.all.2050.rcp85_MR <-subset(projections.2050.rcp85_2_MR,grep("GAM", names(projections.2050.rcp85_2_MR)))
projections.GAM.mean.2050.rcp85_MR <- mean(projections.GAM.all.2050.rcp85_MR)/10

projections.ANN.all.2050.rcp85_MR <- subset(projections.2050.rcp85_2_MR,grep("ANN", names(projections.2050.rcp85_2_MR)))
projections.ANN.mean.2050.rcp85_MR <- mean(projections.ANN.all.2050.rcp85_MR)/10

projections.SRE.all.2050.rcp85_MR <- subset(projections.2050.rcp85_2_MR,grep("SRE", names(projections.2050.rcp85_2_MR)))
projections.SRE.mean.2050.rcp85_MR <- mean(projections.SRE.all.2050.rcp85_MR)/10

projections.MARS.all.2050.rcp85_MR <- subset(projections.2050.rcp85_2_MR,grep("MARS", names(projections.2050.rcp85_2_MR)))
projections.MARS.mean.2050.rcp85_MR <- mean(projections.MARS.all.2050.rcp85_MR)/10

projections.FDA.all.2050.rcp85_MR <- subset(projections.2050.rcp85_2_MR,grep("FDA", names(projections.2050.rcp85_2_MR)))
projections.FDA.mean.2050.rcp85_MR <- mean(projections.FDA.all.2050.rcp85_MR)/10

projections.MAXENT.all.2050.rcp85_MR <- subset(projections.2050.rcp85_2_MR,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_MR)))
projections.MAXENT.mean.2050.rcp85_MR <- mean(projections.MAXENT.all.2050.rcp85_MR)/10



#MC
projections.RF.all.2050.rcp85_MC <- subset(projections.2050.rcp85_1_MC, grep("RF", names(projections.2050.rcp85_1_MC)))
projections.RF.mean.2050.rcp85_MC <- mean(projections.RF.all.2050.rcp85_MC)/10

projections.GBM.all.2050.rcp85_MC <-subset(projections.2050.rcp85_1_MC, grep("GBM", names(projections.2050.rcp85_1_MC)))
projections.GBM.mean.2050.rcp85_MC <- mean(projections.GBM.all.2050.rcp85_MC)/10

projections.CTA.all.2050.rcp85_MC <-subset(projections.2050.rcp85_1_MC,grep("CTA", names(projections.2050.rcp85_1_MC)))
projections.CTA.mean.2050.rcp85_MC <- mean(projections.CTA.all.2050.rcp85_MC)/10

projections.GLM.all.2050.rcp85_MC <-subset(projections.2050.rcp85_2_MC,grep("GLM", names(projections.2050.rcp85_2_MC)))
projections.GLM.mean.2050.rcp85_MC <- mean(projections.GLM.all.2050.rcp85_MC)/10

projections.GAM.all.2050.rcp85_MC <-subset(projections.2050.rcp85_2_MC,grep("GAM", names(projections.2050.rcp85_2_MC)))
projections.GAM.mean.2050.rcp85_MC <- mean(projections.GAM.all.2050.rcp85_MC)/10

projections.ANN.all.2050.rcp85_MC <- subset(projections.2050.rcp85_2_MC,grep("ANN", names(projections.2050.rcp85_2_MC)))
projections.ANN.mean.2050.rcp85_MC <- mean(projections.ANN.all.2050.rcp85_MC)/10

projections.SRE.all.2050.rcp85_MC <- subset(projections.2050.rcp85_2_MC,grep("SRE", names(projections.2050.rcp85_2_MC)))
projections.SRE.mean.2050.rcp85_MC <- mean(projections.SRE.all.2050.rcp85_MC)/10

projections.MARS.all.2050.rcp85_MC <- subset(projections.2050.rcp85_2_MC,grep("MARS", names(projections.2050.rcp85_2_MC)))
projections.MARS.mean.2050.rcp85_MC <- mean(projections.MARS.all.2050.rcp85_MC)/10

projections.FDA.all.2050.rcp85_MC <- subset(projections.2050.rcp85_2_MC,grep("FDA", names(projections.2050.rcp85_2_MC)))
projections.FDA.mean.2050.rcp85_MC <- mean(projections.FDA.all.2050.rcp85_MC)/10

projections.MAXENT.all.2050.rcp85_MC <- subset(projections.2050.rcp85_2_MC,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_MC)))
projections.MAXENT.mean.2050.rcp85_MC <- mean(projections.MAXENT.all.2050.rcp85_MC)/10



#MP
projections.RF.all.2050.rcp85_MP <- subset(projections.2050.rcp85_1_MP, grep("RF", names(projections.2050.rcp85_1_MP)))
projections.RF.mean.2050.rcp85_MP <- mean(projections.RF.all.2050.rcp85_MP)/10

projections.GBM.all.2050.rcp85_MP <-subset(projections.2050.rcp85_1_MP, grep("GBM", names(projections.2050.rcp85_1_MP)))
projections.GBM.mean.2050.rcp85_MP <- mean(projections.GBM.all.2050.rcp85_MP)/10

projections.CTA.all.2050.rcp85_MP <-subset(projections.2050.rcp85_1_MP,grep("CTA", names(projections.2050.rcp85_1_MP)))
projections.CTA.mean.2050.rcp85_MP <- mean(projections.CTA.all.2050.rcp85_MP)/10

projections.GLM.all.2050.rcp85_MP <-subset(projections.2050.rcp85_2_MP,grep("GLM", names(projections.2050.rcp85_2_MP)))
projections.GLM.mean.2050.rcp85_MP <- mean(projections.GLM.all.2050.rcp85_MP)/10

projections.GAM.all.2050.rcp85_MP <-subset(projections.2050.rcp85_2_MP,grep("GAM", names(projections.2050.rcp85_2_MP)))
projections.GAM.mean.2050.rcp85_MP <- mean(projections.GAM.all.2050.rcp85_MP)/10

projections.ANN.all.2050.rcp85_MP <- subset(projections.2050.rcp85_2_MP,grep("ANN", names(projections.2050.rcp85_2_MP)))
projections.ANN.mean.2050.rcp85_MP <- mean(projections.ANN.all.2050.rcp85_MP)/10

projections.SRE.all.2050.rcp85_MP <- subset(projections.2050.rcp85_2_MP,grep("SRE", names(projections.2050.rcp85_2_MP)))
projections.SRE.mean.2050.rcp85_MP <- mean(projections.SRE.all.2050.rcp85_MP)/10

projections.MARS.all.2050.rcp85_MP <- subset(projections.2050.rcp85_2_MP,grep("MARS", names(projections.2050.rcp85_2_MP)))
projections.MARS.mean.2050.rcp85_MP <- mean(projections.MARS.all.2050.rcp85_MP)/10

projections.FDA.all.2050.rcp85_MP <- subset(projections.2050.rcp85_2_MP,grep("FDA", names(projections.2050.rcp85_2_MP)))
projections.FDA.mean.2050.rcp85_MP <- mean(projections.FDA.all.2050.rcp85_MP)/10

projections.MAXENT.all.2050.rcp85_MP <- subset(projections.2050.rcp85_2_MP,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_MP)))
projections.MAXENT.mean.2050.rcp85_MP <- mean(projections.MAXENT.all.2050.rcp85_MP)/10



#MG
projections.RF.all.2050.rcp85_MG <- subset(projections.2050.rcp85_1_MG, grep("RF", names(projections.2050.rcp85_1_MG)))
projections.RF.mean.2050.rcp85_MG <- mean(projections.RF.all.2050.rcp85_MG)/10

projections.GBM.all.2050.rcp85_MG <-subset(projections.2050.rcp85_1_MG, grep("GBM", names(projections.2050.rcp85_1_MG)))
projections.GBM.mean.2050.rcp85_MG <- mean(projections.GBM.all.2050.rcp85_MG)/10

projections.CTA.all.2050.rcp85_MG <-subset(projections.2050.rcp85_1_MG,grep("CTA", names(projections.2050.rcp85_1_MG)))
projections.CTA.mean.2050.rcp85_MG <- mean(projections.CTA.all.2050.rcp85_MG)/10

projections.GLM.all.2050.rcp85_MG <-subset(projections.2050.rcp85_2_MG,grep("GLM", names(projections.2050.rcp85_2_MG)))
projections.GLM.mean.2050.rcp85_MG <- mean(projections.GLM.all.2050.rcp85_MG)/10

projections.GAM.all.2050.rcp85_MG <-subset(projections.2050.rcp85_2_MG,grep("GAM", names(projections.2050.rcp85_2_MG)))
projections.GAM.mean.2050.rcp85_MG <- mean(projections.GAM.all.2050.rcp85_MG)/10

projections.ANN.all.2050.rcp85_MG <- subset(projections.2050.rcp85_2_MG,grep("ANN", names(projections.2050.rcp85_2_MG)))
projections.ANN.mean.2050.rcp85_MG <- mean(projections.ANN.all.2050.rcp85_MG)/10

projections.SRE.all.2050.rcp85_MG <- subset(projections.2050.rcp85_2_MG,grep("SRE", names(projections.2050.rcp85_2_MG)))
projections.SRE.mean.2050.rcp85_MG <- mean(projections.SRE.all.2050.rcp85_MG)/10

projections.MARS.all.2050.rcp85_MG <- subset(projections.2050.rcp85_2_MG,grep("MARS", names(projections.2050.rcp85_2_MG)))
projections.MARS.mean.2050.rcp85_MG <- mean(projections.MARS.all.2050.rcp85_MG)/10

projections.FDA.all.2050.rcp85_MG <- subset(projections.2050.rcp85_2_MG,grep("FDA", names(projections.2050.rcp85_2_MG)))
projections.FDA.mean.2050.rcp85_MG <- mean(projections.FDA.all.2050.rcp85_MG)/10

projections.MAXENT.all.2050.rcp85_MG <- subset(projections.2050.rcp85_2_MG,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_MG)))
projections.MAXENT.mean.2050.rcp85_MG <- mean(projections.MAXENT.all.2050.rcp85_MG)/10



#NO
projections.RF.all.2050.rcp85_NO <- subset(projections.2050.rcp85_1_NO, grep("RF", names(projections.2050.rcp85_1_NO)))
projections.RF.mean.2050.rcp85_NO <- mean(projections.RF.all.2050.rcp85_NO)/10

projections.GBM.all.2050.rcp85_NO <-subset(projections.2050.rcp85_1_NO, grep("GBM", names(projections.2050.rcp85_1_NO)))
projections.GBM.mean.2050.rcp85_NO <- mean(projections.GBM.all.2050.rcp85_NO)/10

projections.CTA.all.2050.rcp85_NO <-subset(projections.2050.rcp85_1_NO,grep("CTA", names(projections.2050.rcp85_1_NO)))
projections.CTA.mean.2050.rcp85_NO <- mean(projections.CTA.all.2050.rcp85_NO)/10

projections.GLM.all.2050.rcp85_NO <-subset(projections.2050.rcp85_2_NO,grep("GLM", names(projections.2050.rcp85_2_NO)))
projections.GLM.mean.2050.rcp85_NO <- mean(projections.GLM.all.2050.rcp85_NO)/10

projections.GAM.all.2050.rcp85_NO <-subset(projections.2050.rcp85_2_NO,grep("GAM", names(projections.2050.rcp85_2_NO)))
projections.GAM.mean.2050.rcp85_NO <- mean(projections.GAM.all.2050.rcp85_NO)/10

projections.ANN.all.2050.rcp85_NO <- subset(projections.2050.rcp85_2_NO,grep("ANN", names(projections.2050.rcp85_2_NO)))
projections.ANN.mean.2050.rcp85_NO <- mean(projections.ANN.all.2050.rcp85_NO)/10

projections.SRE.all.2050.rcp85_NO <- subset(projections.2050.rcp85_2_NO,grep("SRE", names(projections.2050.rcp85_2_NO)))
projections.SRE.mean.2050.rcp85_NO <- mean(projections.SRE.all.2050.rcp85_NO)/10

projections.MARS.all.2050.rcp85_NO <- subset(projections.2050.rcp85_2_NO,grep("MARS", names(projections.2050.rcp85_2_NO)))
projections.MARS.mean.2050.rcp85_NO <- mean(projections.MARS.all.2050.rcp85_NO)/10

projections.FDA.all.2050.rcp85_NO <- subset(projections.2050.rcp85_2_NO,grep("FDA", names(projections.2050.rcp85_2_NO)))
projections.FDA.mean.2050.rcp85_NO <- mean(projections.FDA.all.2050.rcp85_NO)/10

projections.MAXENT.all.2050.rcp85_NO <- subset(projections.2050.rcp85_2_NO,grep("MAXENT.Phillips", names(projections.2050.rcp85_2_NO)))
projections.MAXENT.mean.2050.rcp85_NO <- mean(projections.MAXENT.all.2050.rcp85_NO)/10





##############################################################
####### AVERAGE ENSEMBLE FOR SELECTED ALGORITHMS #############
##############################################################
#GCMs: AC, BC, CC, CN, GS, HD, HG, HE, IN, IP, MI, MR, MC, MP, MG, NO


#AC
projections.all.mean.2050.rcp85_AC <- mean(projections.RF.mean.2050.rcp85_AC + projections.GBM.mean.2050.rcp85_AC +
	projections.CTA.mean.2050.rcp85_AC + projections.GLM.mean.2050.rcp85_AC + projections.GAM.mean.2050.rcp85_AC +
	projections.ANN.mean.2050.rcp85_AC + projections.SRE.mean.2050.rcp85_AC +
	projections.MARS.mean.2050.rcp85_AC + projections.FDA.mean.2050.rcp85_AC + projections.MAXENT.mean.2050.rcp85_AC)
writeRaster(projections.all.mean.2050.rcp85_AC, filename="Future Climate - 2050_rcp8.5_AC.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp85_AC, filename="Future Climate - 2050_rcp8.5_AC.asc", format="ascii")

#BC
projections.all.mean.2050.rcp85_BC <- mean(projections.RF.mean.2050.rcp85_BC + projections.GBM.mean.2050.rcp85_BC +
	projections.CTA.mean.2050.rcp85_BC + projections.GLM.mean.2050.rcp85_BC + projections.GAM.mean.2050.rcp85_BC +
	projections.ANN.mean.2050.rcp85_BC + projections.SRE.mean.2050.rcp85_BC +
	projections.MARS.mean.2050.rcp85_BC + projections.FDA.mean.2050.rcp85_BC + projections.MAXENT.mean.2050.rcp85_BC)
writeRaster(projections.all.mean.2050.rcp85_BC, filename="Future Climate - 2050_rcp8.5_BC.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp85_BC, filename="Future Climate - 2050_rcp8.5_BC.asc", format="ascii")

#CC
projections.all.mean.2050.rcp85_CC <- mean(projections.RF.mean.2050.rcp85_CC + projections.GBM.mean.2050.rcp85_CC +
	projections.CTA.mean.2050.rcp85_CC + projections.GLM.mean.2050.rcp85_CC + projections.GAM.mean.2050.rcp85_CC +
	projections.ANN.mean.2050.rcp85_CC + projections.SRE.mean.2050.rcp85_CC +
	projections.MARS.mean.2050.rcp85_CC + projections.FDA.mean.2050.rcp85_CC + projections.MAXENT.mean.2050.rcp85_CC)
writeRaster(projections.all.mean.2050.rcp85_CC, filename="Future Climate - 2050_rcp8.5_CC.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp85_CC, filename="Future Climate - 2050_rcp8.5_CC.asc", format="ascii")

#CN
projections.all.mean.2050.rcp85_CN <- mean(projections.RF.mean.2050.rcp85_CN + projections.GBM.mean.2050.rcp85_CN +
	projections.CTA.mean.2050.rcp85_CN + projections.GLM.mean.2050.rcp85_CN + projections.GAM.mean.2050.rcp85_CN +
	projections.ANN.mean.2050.rcp85_CN + projections.SRE.mean.2050.rcp85_CN +
	projections.MARS.mean.2050.rcp85_CN + projections.FDA.mean.2050.rcp85_CN + projections.MAXENT.mean.2050.rcp85_CN)
writeRaster(projections.all.mean.2050.rcp85_CN, filename="Future Climate - 2050_rcp8.5_CN.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp85_CN, filename="Future Climate - 2050_rcp8.5_CN.asc", format="ascii")

#GS
projections.all.mean.2050.rcp85_GS <- mean(projections.RF.mean.2050.rcp85_GS + projections.GBM.mean.2050.rcp85_GS +
	projections.CTA.mean.2050.rcp85_GS + projections.GLM.mean.2050.rcp85_GS + projections.GAM.mean.2050.rcp85_GS +
	projections.ANN.mean.2050.rcp85_GS + projections.SRE.mean.2050.rcp85_GS +
	projections.MARS.mean.2050.rcp85_GS + projections.FDA.mean.2050.rcp85_GS + projections.MAXENT.mean.2050.rcp85_GS)
writeRaster(projections.all.mean.2050.rcp85_GS, filename="Future Climate - 2050_rcp8.5_GS.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp85_GS, filename="Future Climate - 2050_rcp8.5_GS.asc", format="ascii")

#HD
projections.all.mean.2050.rcp85_HD <- mean(projections.RF.mean.2050.rcp85_HD + projections.GBM.mean.2050.rcp85_HD +
	projections.CTA.mean.2050.rcp85_HD + projections.GLM.mean.2050.rcp85_HD + projections.GAM.mean.2050.rcp85_HD +
	projections.ANN.mean.2050.rcp85_HD + projections.SRE.mean.2050.rcp85_HD +
	projections.MARS.mean.2050.rcp85_HD + projections.FDA.mean.2050.rcp85_HD + projections.MAXENT.mean.2050.rcp85_HD)
writeRaster(projections.all.mean.2050.rcp85_HD, filename="Future Climate - 2050_rcp8.5_HD.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp85_HD, filename="Future Climate - 2050_rcp8.5_HD.asc", format="ascii")

#HG
projections.all.mean.2050.rcp85_HG <- mean(projections.RF.mean.2050.rcp85_HG + projections.GBM.mean.2050.rcp85_HG +
	projections.CTA.mean.2050.rcp85_HG + projections.GLM.mean.2050.rcp85_HG + projections.GAM.mean.2050.rcp85_HG +
	projections.ANN.mean.2050.rcp85_HG + projections.SRE.mean.2050.rcp85_HG +
	projections.MARS.mean.2050.rcp85_HG + projections.FDA.mean.2050.rcp85_HG + projections.MAXENT.mean.2050.rcp85_HG)
writeRaster(projections.all.mean.2050.rcp85_HG, filename="Future Climate - 2050_rcp8.5_HG.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp85_HG, filename="Future Climate - 2050_rcp8.5_HG.asc", format="ascii")


#HE
projections.all.mean.2050.rcp85_HE <- mean(projections.RF.mean.2050.rcp85_HE + projections.GBM.mean.2050.rcp85_HE +
	projections.CTA.mean.2050.rcp85_HE + projections.GLM.mean.2050.rcp85_HE + projections.GAM.mean.2050.rcp85_HE +
	projections.ANN.mean.2050.rcp85_HE + projections.SRE.mean.2050.rcp85_HE +
	projections.MARS.mean.2050.rcp85_HE + projections.FDA.mean.2050.rcp85_HE + projections.MAXENT.mean.2050.rcp85_HE)
writeRaster(projections.all.mean.2050.rcp85_HE, filename="Future Climate - 2050_rcp8.5_HE.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp85_HE, filename="Future Climate - 2050_rcp8.5_HE.asc", format="ascii")


#IN
projections.all.mean.2050.rcp85_IN <- mean(projections.RF.mean.2050.rcp85_IN + projections.GBM.mean.2050.rcp85_IN +
	projections.CTA.mean.2050.rcp85_IN + projections.GLM.mean.2050.rcp85_IN + projections.GAM.mean.2050.rcp85_IN +
	projections.ANN.mean.2050.rcp85_IN + projections.SRE.mean.2050.rcp85_IN +
	projections.MARS.mean.2050.rcp85_IN + projections.FDA.mean.2050.rcp85_IN + projections.MAXENT.mean.2050.rcp85_IN)
writeRaster(projections.all.mean.2050.rcp85_IN, filename="Future Climate - 2050_rcp8.5_IN.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp85_IN, filename="Future Climate - 2050_rcp8.5_IN.asc", format="ascii")

#IP
projections.all.mean.2050.rcp85_IP <- mean(projections.RF.mean.2050.rcp85_IP + projections.GBM.mean.2050.rcp85_IP +
	projections.CTA.mean.2050.rcp85_IP + projections.GLM.mean.2050.rcp85_IP + projections.GAM.mean.2050.rcp85_IP +
	projections.ANN.mean.2050.rcp85_IP + projections.SRE.mean.2050.rcp85_IP +
	projections.MARS.mean.2050.rcp85_IP + projections.FDA.mean.2050.rcp85_IP + projections.MAXENT.mean.2050.rcp85_IP)
writeRaster(projections.all.mean.2050.rcp85_IP, filename="Future Climate - 2050_rcp8.5_IP.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp85_IP, filename="Future Climate - 2050_rcp8.5_IP.asc", format="ascii")

#MI
projections.all.mean.2050.rcp85_MI <- mean(projections.RF.mean.2050.rcp85_MI + projections.GBM.mean.2050.rcp85_MI +
	projections.CTA.mean.2050.rcp85_MI + projections.GLM.mean.2050.rcp85_MI + projections.GAM.mean.2050.rcp85_MI +
	projections.ANN.mean.2050.rcp85_MI + projections.SRE.mean.2050.rcp85_MI +
	projections.MARS.mean.2050.rcp85_MI + projections.FDA.mean.2050.rcp85_MI + projections.MAXENT.mean.2050.rcp85_MI)
writeRaster(projections.all.mean.2050.rcp85_MI, filename="Future Climate - 2050_rcp8.5_MI.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp85_MI, filename="Future Climate - 2050_rcp8.5_MI.asc", format="ascii")

#MR
projections.all.mean.2050.rcp85_MR <- mean(projections.RF.mean.2050.rcp85_MR + projections.GBM.mean.2050.rcp85_MR +
	projections.CTA.mean.2050.rcp85_MR + projections.GLM.mean.2050.rcp85_MR + projections.GAM.mean.2050.rcp85_MR +
	projections.ANN.mean.2050.rcp85_MR + projections.SRE.mean.2050.rcp85_MR +
	projections.MARS.mean.2050.rcp85_MR + projections.FDA.mean.2050.rcp85_MR + projections.MAXENT.mean.2050.rcp85_MR)
writeRaster(projections.all.mean.2050.rcp85_MR, filename="Future Climate - 2050_rcp8.5_MR.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp85_MR, filename="Future Climate - 2050_rcp8.5_MR.asc", format="ascii")

#MC
projections.all.mean.2050.rcp85_MC <- mean(projections.RF.mean.2050.rcp85_MC + projections.GBM.mean.2050.rcp85_MC +
	projections.CTA.mean.2050.rcp85_MC + projections.GLM.mean.2050.rcp85_MC + projections.GAM.mean.2050.rcp85_MC +
	projections.ANN.mean.2050.rcp85_MC + projections.SRE.mean.2050.rcp85_MC +
	projections.MARS.mean.2050.rcp85_MC + projections.FDA.mean.2050.rcp85_MC + projections.MAXENT.mean.2050.rcp85_MC)
writeRaster(projections.all.mean.2050.rcp85_MC, filename="Future Climate - 2050_rcp8.5_MC.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp85_MC, filename="Future Climate - 2050_rcp8.5_MC.asc", format="ascii")

#MP
projections.all.mean.2050.rcp85_MP <- mean(projections.RF.mean.2050.rcp85_MP + projections.GBM.mean.2050.rcp85_MP +
	projections.CTA.mean.2050.rcp85_MP + projections.GLM.mean.2050.rcp85_MP + projections.GAM.mean.2050.rcp85_MP +
	projections.ANN.mean.2050.rcp85_MP + projections.SRE.mean.2050.rcp85_MP +
	projections.MARS.mean.2050.rcp85_MP + projections.FDA.mean.2050.rcp85_MP + projections.MAXENT.mean.2050.rcp85_MP)
writeRaster(projections.all.mean.2050.rcp85_MP, filename="Future Climate - 2050_rcp8.5_MP.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp85_MP, filename="Future Climate - 2050_rcp8.5_MP.asc", format="ascii")

#MG
projections.all.mean.2050.rcp85_MG <- mean(projections.RF.mean.2050.rcp85_MG + projections.GBM.mean.2050.rcp85_MG +
	projections.CTA.mean.2050.rcp85_MG + projections.GLM.mean.2050.rcp85_MG + projections.GAM.mean.2050.rcp85_MG +
	projections.ANN.mean.2050.rcp85_MG + projections.SRE.mean.2050.rcp85_MG +
	projections.MARS.mean.2050.rcp85_MG + projections.FDA.mean.2050.rcp85_MG + projections.MAXENT.mean.2050.rcp85_MG)
writeRaster(projections.all.mean.2050.rcp85_MG, filename="Future Climate - 2050_rcp8.5_MG.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp85_MG, filename="Future Climate - 2050_rcp8.5_MG.asc", format="ascii")

#NO
projections.all.mean.2050.rcp85_NO <- mean(projections.RF.mean.2050.rcp85_NO + projections.GBM.mean.2050.rcp85_NO +
	projections.CTA.mean.2050.rcp85_NO + projections.GLM.mean.2050.rcp85_NO + projections.GAM.mean.2050.rcp85_NO +
	projections.ANN.mean.2050.rcp85_NO + projections.SRE.mean.2050.rcp85_NO +
	projections.MARS.mean.2050.rcp85_NO + projections.FDA.mean.2050.rcp85_NO + projections.MAXENT.mean.2050.rcp85_NO)
writeRaster(projections.all.mean.2050.rcp85_NO, filename="Future Climate - 2050_rcp8.5_NO.tif", format="GTiff")
writeRaster(projections.all.mean.2050.rcp85_NO, filename="Future Climate - 2050_rcp8.5_NO.asc", format="ascii")


#Ensemble 2050 - rcp8.5
#GCMs: AC, BC, CC, CN, GS, HD, HG, HE, IN, IP, MI, MR, MC, MP, MG, NO
ensemble2050rcp8.5 <- mean(projections.all.mean.2050.rcp85_AC, projections.all.mean.2050.rcp85_BC,
					projections.all.mean.2050.rcp85_CC, projections.all.mean.2050.rcp85_CN,
					projections.all.mean.2050.rcp85_GS, projections.all.mean.2050.rcp85_HD,
					projections.all.mean.2050.rcp85_HG, projections.all.mean.2050.rcp85_HE,
					projections.all.mean.2050.rcp85_IN, projections.all.mean.2050.rcp85_IP,
					projections.all.mean.2050.rcp85_MI, projections.all.mean.2050.rcp85_MR,
					projections.all.mean.2050.rcp85_MC, projections.all.mean.2050.rcp85_MP,
					projections.all.mean.2050.rcp85_MG, projections.all.mean.2050.rcp85_NO)
writeRaster(ensemble2050rcp8.5, filename="Ensemble - Future Climate - 2050_rcp8.5.tif", format="GTiff")
writeRaster(ensemble2050rcp8.5, filename="Ensemble - Future Climate - 2050_rcp8.5.asc", format="ascii")
plot(ensemble2050rcp8.5)

#Converting to binary maps:
projections.binary_2050rcp8.5 <- BinaryTransformation(ensemble2050rcp8.5, th) #Replace 'th' by the threshold value you wish to apply.
class(projections.binary_2050rcp8.5)
summary(values(projections.binary_2050rcp8.5))
writeRaster(projections.binary_2050rcp8.5, filename="Ensemble - Future Climate - 2050_rcp8.5_BINARY.tif", format="GTiff")
writeRaster(projections.binary_2050rcp8.5, filename="Ensemble - Future Climate - 2050_rcp8.5_BINARY.asc", formato="ascii")
plot(projections.binary_2050rcp8.5)



## End
## ***
