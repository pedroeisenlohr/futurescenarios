#####################################################
#### NICHE MODELS FORECASTED FOR FUTURE SCENARIOS ###
#####################################################

# Released on: 16th July, 2019

# Contact: Pedro V. Eisenlohr (pedro.eisenlohr@unemat.br)

###################### Acknowledgments ##########################
### Dr. Guarino Colli's team of Universidade de Brasília.
### Dr. Diogo Souza Bezerra Rocha (Instituto de Pesquisas Jardim Botânico/RJ).
### My students of Ecology Lab, mainly João Carlos Pires de Oliveira.

#!/usr/bin/env Rscript



#---------------------------------#
# DEFINING THE WORK DIRECTORY ####
#-------------------------------#



#To follow this routine, it is necessary:
# 1) Define a folder named "Habitat Suitability Models" as the working directory.
# 2) Inside Habitat Suitability Models ", include the following folders and subfolders:
# 2.1) Environmental Layers
# 2.1.1) CHELSA.
# 2.1.2) CHELSA_Future.
# 2.2) Shapefiles.
# If you would like the complete contents of these folders, please contact me (pedro.eisenlohr@unemat.br)
# Your spreadsheet needs to contain three columns (sp, lon, lat)
#setwd(choose.dir()) #Defina a pasta Habitat Suitability Models
# if (dir.exists("Habitat Suitability Models") == F) {
#   dir.create("Habitat Suitability Models")
# }
# 
 # setwd(choose.dir()) # escolha a pasta que contem todas as demais pastas dos GCMs e a pasta 'presente'
getwd()
dir() #Dentre as pastas, DEVE haver as pastas com os nomes dos GCMs e uma com nome 'presente'


#----------------------------------------#
# INSTALL AND LIBRARY THE  PACKAGES ####
#--------------------------------------#

#install.packages("biomod2", dep=T)
#install.packages("car", dep=T)
#install.packages("colorRamps", dep=T)
#install.packages("doParallel", dep=T)
#install.packages("dplyr", dep=T)
#install.packages("maps", dep=T)
#install.packages("plotKML", dep=T)
#install.packages("rgdal", dep=T)
#install.packages("sdm", dep=T)
#install.packages("sqldf", dep=T)
#install.packages("testthat", dep=T)
#install.packages("usdm", dep=T)
#install.packages("FactoMineR", dep=T)
#install.packages("foreach", dep=T)
#install.packages("sdmvspecies",dep=T)
#install.packages("filesstrings",dep=T)
# install.packages("filesstrings",dep=T)
#install.packages("githubinstall",dep=T)
# library(githubinstall)
#githubinstall("virtualspecies", force=T) # need the last the virtualspecies package version from Github

library(biomod2)
library(foreach)
library(doParallel)
library(car)
library(maptools)
library(colorRamps)
library(dismo)
library(dplyr)
library(maps)
library(plotKML)
library(rgdal)
library(sdm)
library(sqldf)
library(testthat)
library(usdm)
library(FactoMineR)
library(sdmvspecies)
library(filesstrings)
library(virtualspecies)

if (dir.exists("outputs") == F) {
  dir.create("outputs")
}
##MUDAR DIRETÓRIO DE ARQUIVOS TEMPORÁRIOS DO raster
## define the name of a temp directory where raster tmp files will be stored
raster_tmp_dir <- "raster_tmp"
if (dir.exists("raster_tmp") == F) {
  ## create the directory (somente ao iniciar a modelagem):
  dir.create(raster_tmp_dir,
             showWarnings = F,
             recursive = T)
}

## set raster options
rasterOptions(tmpdir = raster_tmp_dir)
### Sempre que necessário:
# Processamento paralelo #

detectCores()
getDoParWorkers()
cl <- parallel::makeCluster(2)
registerDoParallel(cl)
getDoParWorkers()

# Aumento da alocação de memória
memory.limit(100000000) # ou algum outro valor de memória (em kB)

#Baixar o Maxent (Apenas rode a funcao abaixo se você não tiver baixado o Maxent)
#MaxEnt .jar
jar <- paste0(system.file(package = "dismo"), "/java/maxent.jar")
if (file.exists(jar) != T) {
  url = "http://biodiversityinformatics.amnh.org/open_source/maxent/maxent.php?op=download"
  download.file(url, dest = "maxent.zip", mode = "wb")
  unzip("maxent.zip",
        files = "maxent.jar",
        exdir = system.file("java", package = "dismo"))
  unlink("maxent.zip")
  warning("Maxent foi colocado no diretório")
}
system.file("java", package = "dismo")


#----------------------------------#
# IMPORTANDO E CHECANDO OS DADOS ##
#--------------------------------#

# Importando dados climáticos do presente

bio.crop <-
  list.files(
    "./CHELSA",
    full.names = TRUE,
    pattern = "tif$"
  )
bio.crop
bio.crop <- stack(bio.crop)
bio.crop <-disaggregate(bio.crop, fact=5, fun=mean)
bio.crop <-aggregate(bio.crop, fact=5, fun=mean)
bio.crop<-rescale(stack(bio.crop))
names(bio.crop)
res(bio.crop)
names(bio.crop)<-names(bio.crop)<-paste0("PC",1:6)


# Importando dados bióticos

spp<-read.table(file.choose(), header=T, sep=",")
dim(spp)
head(spp, 10)

table(spp$sp)

especies <- unique(spp$sp)
especies
#------------------------#
#beginning of the modeling ####
#----------------------#

# Creating objects for models calibration
models1<-c("CTA","RF", "GBM")
models2<-c("MAXENT.Phillips", "GAM", "MARS","ANN", "GLM","FDA","MAXENT.Tsuruoka")
n.runs = 3 # number of RUNs
n.algo1 = length(models1)# number of algorithms
n.algo2 = length(models2) #numero de algorithms
n.conj.pa2 = 2 # set of pseudo-absences
env.selected = bio.crop
especie = especies[1] # To model without a loop, remove the '#' of this line and add it to the 'for', 'foreach' and '.packages'
#-------------------------#
#beginning of the loop####
#-----------------------#
# for(especie in especies[1:length(especies)]){
# foreach(especie = especies, # For parallel looping (Multiple Species)
# .packages = c("raster", "biomod2", 'sp', "sdmvspecies", "filesstrings")) %dopar% {
# ini1 = Sys.time()
# criando tabela para uma especie
occs <- spp[spp$sp == especie, c("lon", "lat")]

# nome = strsplit(as.vector(especie), " ")
# especie = paste(nome[[1]][1], nome[[1]][2], sep = ".")

# Selecionado pontos espacialmente únicos #
mask <- env.selected[[1]]
{(cell <-
    cellFromXY(mask, occs[, 1:2])) # get the cell number for each point
  (x<-(cbind(occs[, 1:2], cell)))
  #dup <- duplicated(cbind(occs[, 1:2], cell))
  (dup2 <- duplicated(cbind(cell)))
  xv<-data.frame(x,dup2)
  xv[xv=="TRUE"]<-NA
  (xv<-na.omit(xv))
  xv<-xv[,1:2]
  occs =xv # select the records that are not duplicated
}


#-----------------------------------------------#
# GENERATING OTHER REQUIRED OBJECTS FOR SDM ####
#---------------------------------------------#

# Convert dataset to SpatialPointsDataFrame (only presences)
myRespXY <-
  occs[, c("lon", "lat")] #Caso dê algum erro aqui, veja como você intitulou as colunas da sua matriz.
# Creating occurrence data object
occurrence.resp <-  rep(1, length(myRespXY$lon))


#------------------------------------------#
# FIT SPECIES DISTRIBUTION MODELS - SDMS ####
#----------------------------------------#

try({    
  coord1 = occs
  sp::coordinates(coord1) <- ~ lon + lat
  raster::crs(coord1) <- raster::crs(env.selected)
  
  dist.mean <- mean(sp::spDists(
    x = coord1,
    longlat = T,
    segments = FALSE
  ))
  dist.min = 5
  dist.min <-  min(sp::spDists(x = coord1,
                               longlat = T,
                               segments = F))
  dist.min = 5
  
  write.table(
    c(dist.min, dist.mean),
    paste0('./outputs/', especie,"_", ".csv"),
    row.names = F,
    sep = ","
  )
})
dim(occs)
PA.number <- length(occs[, 1])
PA.number #número de pontos de ocorrência espacialmente únicos

diretorio = paste0("Occurrence.", especie)

##### FORMATING DATA #####

# Preparando para CTA, GBM e RF:
sppBiomodData.PA.equal <- BIOMOD_FormatingData(
  resp.var = occurrence.resp,
  expl.var = env.selected,
  resp.xy = myRespXY,
  resp.name = diretorio,
  PA.nb.rep = n.conj.pa2, #número de datasets de pseudoausências
  PA.nb.absences = PA.number, #= número de pseudoausências = número de pontos espacialmente únicos
  PA.strategy = "disk",
  # PA.sre.quant = 0.10,
  PA.dist.min = dist.min * 1000,
  PA.dist.max = dist.mean * 1000,
  na.rm = TRUE
)

#Preparando para os demais algoritmos:
sppBiomodData.PA.10000 <- BIOMOD_FormatingData(
  resp.var = occurrence.resp,
  expl.var = env.selected,
  resp.xy = myRespXY,
  resp.name = diretorio,
  PA.nb.rep = n.conj.pa2,
  PA.nb.absences = 1000,
  PA.strategy = "disk",
  # PA.sre.quant = 0.10,
  PA.dist.min = dist.min * 1000,
  PA.dist.max = dist.mean * 1000,
  na.rm = TRUE
)
sppBiomodData.PA.10000

myBiomodOption <-
  BIOMOD_ModelingOptions(MAXENT.Phillips = list(path_to_maxent.jar = jar))

# save.image()
#---------------#
# Modeling ####
#-------------#

# Com partição treino x teste:
sppModelOut.PA.equal <- BIOMOD_Modeling(
  sppBiomodData.PA.equal,
  models =models1,
  models.options = NULL,
  NbRunEval = n.runs, #número de repetições para cada algoritmo
  DataSplit = 70,#percentagem de pts para treino.
  Prevalence = 0.5,
  VarImport = 0,#caso queira avaliar a importância das variáveis, mudar para 10 ou 100 permutações
  models.eval.meth = c("TSS", "ROC"),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = FALSE,
  modeling.id = "spp_presente"
)
# import.var.equal<-data.frame(sppModelOut.PA.equal@variables.importances@val)
# names(import.var.equal)<-rep(c('GBM','CTA','RF'),n.runs + n.conj.pa2)
# import.var.equal
# write.table(import.var.equal,
#             paste0("./outputs/", especie, "_", "Var.import.PA.equal.csv"), sep = ',')

sppModelOut.PA.10000 <- BIOMOD_Modeling(
  sppBiomodData.PA.10000,
  models = models2,
  models.options = myBiomodOption,
  NbRunEval = n.runs,  #número de repetições para cada algoritmo
  DataSplit = 70, #percentagem de pts para treino.
  Prevalence = 0.5,
  VarImport = 0, #caso queira avaliar a importância das variáveis, mudar para 10 ou 100 permutações
  models.eval.meth = c("TSS", "ROC"),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = FALSE,
  modeling.id = "spp_presente"
)

# import.var.1000<-data.frame(sppModelOut.PA.10000@variables.importances@val)
# names(import.var.1000)<-rep(c("MAXENT.Phillips", "GLM", "GAM", "ANN", "FDA", "MARS","MAXENT.Tsuruoka"),n.runs + n.conj.pa2)
# import.var.1000
# write.table(import.var.1000,
#             paste0("./outputs/", especie, "_", "Var.import.PA.1000.csv"), sep = ',')

#---------------------------------#
# EVALUATE MODELS USING BIOMOD2 ##
#-------------------------------#

# Sobre as métricas avaliativas,
# ver http://www.cawcr.gov.au/projects/verification/#Methods_for_dichotomous_forecasts


##### Evaluation of Models ####
sppModelEval.PA.equal <-
  get_evaluations(sppModelOut.PA.equal)#GBM, CTA e RF
sppModelEval.PA.equal
write.table(
  sppModelEval.PA.equal,
  paste0("./outputs/", especie, "_", "EvaluationsAll_1.csv")
)


sppModelEval.PA.10000 <-
  get_evaluations(sppModelOut.PA.10000) #Os demais.
sppModelEval.PA.10000
write.table(
  sppModelEval.PA.10000,
  paste0("./outputs/", especie, "_", "EvaluationsAll_2.csv")
)


# Sumarizando as métricas avaliativas
sdm.models1 <-models1
sdm.models1
eval.methods1 <- c("TSS", "ROC") #2 evaluation methods
eval.methods1

##### Eval.1 ####

means.i1 <- numeric(0)
for (i in 1:n.algo1) {
  m1 <-
    sppModelEval.PA.equal[paste(eval.methods1[1]), "Testing.data", paste(sdm.models1[i]), ,]
  means.i1 = c(means.i1, m1) 
}

summary.eval.equal <-
  data.frame(rep(sdm.models1, each =  n.runs*n.conj.pa2),
             rep(1:n.conj.pa2, each = n.runs),
             rep(1:n.runs, n.algo1),
             means.i1)
names(summary.eval.equal) <- c("Model", "PA","Run", "TSS")
summary.eval.equal
write.table(
  summary.eval.equal,
  paste0("./outputs/", especie, "_", "Models1_Evaluation.csv")
)

#----------------------------------------------------------------------------------------#
means.i1 <- numeric(0)
for (i in 1:n.algo1) {
  m1 <-
    sppModelEval.PA.equal[paste(eval.methods1[2]), "Sensitivity", paste(sdm.models1[i]), ,]
  means.i1 = c(means.i1, m1)
}

summary.eval.equal.1 <-
  data.frame(means.i1)
summary.eval.equal.1
(test1<-cbind(summary.eval.equal,summary.eval.equal.1))
names(test1)<-c("Model", "PA","Run","TSS","Se")
test1
#----------------------------------------------------------------------------------------#

means.i1.1 <- numeric(0)
means.j1.1 <- numeric(2)
for (i in 1:n.algo1){
  for (j in 1:2){
    means.j1.1[j] <- mean(sppModelEval.PA.equal[paste(eval.methods1[j]),"Testing.data",paste(sdm.models1[i]),,])
  }
  means.i1.1 <- c(means.i1.1, means.j1.1)
}

summary.eval.equal.mean <- data.frame(rep(sdm.models1,each=j), rep(eval.methods1,i), means.i1.1)
names(summary.eval.equal.mean) <- c("Model", "Method", "Mean")
summary.eval.equal.mean
write.table(summary.eval.equal.mean,
            paste0("./outputs/", especie, "_", "Models1_Evaluation_Mean.csv"))

sd.i1 <- numeric(0)
sd.j1 <- numeric(2)
for (i in 1:n.algo1) {
  for (j in 1:2) {
    sd.j1[j] <-
      sd(sppModelEval.PA.equal[paste(eval.methods1[j]), "Testing.data", paste(sdm.models1[i]), ,])
  }
  sd.i1 <- c(sd.i1, sd.j1)
}

summary.eval.equal.sd <-
  data.frame(rep(sdm.models1, each = 2), rep(eval.methods1, n.algo1), sd.i1)
names(summary.eval.equal.sd) <- c("Model", "Method", "SD")
summary.eval.equal.sd
write.table(
  summary.eval.equal.sd,
  paste0("./outputs/", especie, "_", "Models1_Evaluation_SD.csv")
)


sdm.models2 <-models2 #7 models
sdm.models2
eval.methods2 <- c("TSS", "ROC") #2 evaluation methods
eval.methods2

##### Eval.2 ####

means.i2 <- numeric(0)
for (i2 in 1:n.algo2) {
  m2 <-
    sppModelEval.PA.10000[paste(eval.methods2[1]), "Testing.data", paste(sdm.models2[i2]), ,]
  means.i2 = c(means.i2, m2)
}

summary.eval.10000 <-
  data.frame(rep(sdm.models2, each =  n.runs*n.conj.pa2),
             rep(1:n.conj.pa2, each = n.runs),
             rep(1:n.runs, n.algo2),
             means.i2)
names(summary.eval.10000) <- c("Model", "PA","Run", "TSS")
summary.eval.10000
write.table(
  summary.eval.10000,
  paste0("./outputs/", especie, "_", "Models2_Evaluation.csv")
)

#----------------------------------------------------------------------------------------#
means.i21 <- numeric(0)
for (i21 in 1:n.algo2) {
  m21 <-
    sppModelEval.PA.10000[paste(eval.methods2[2]), "Sensitivity", paste(sdm.models2[i21]), ,]
  means.i21 = c(means.i21, m21)
}

summary.eval.10000.1 <-
  data.frame(means.i21)
summary.eval.10000.1
(test2<-cbind(summary.eval.10000,summary.eval.10000.1))
names(test2)<-c("Model", "PA","Run","TSS","Se")
test2
#----------------------------------------------------------------------------------------#

means.i2.2 <- numeric(0)
means.j2.2 <- numeric(2)
for (i in 1:n.algo2){
  for (j in 1:2){
    means.j2.2[j] <- mean(sppModelEval.PA.10000[paste(eval.methods2[j]),"Testing.data",paste(sdm.models2[i]),,], na.rm = T)
  }
  means.i2.2 <- c(means.i2.2, means.j2.2)
}

summary.eval.10000.mean <- data.frame(rep(sdm.models2,each=j), rep(eval.methods2,i), means.i2.2)
names(summary.eval.10000.mean) <- c("Model", "Method", "Mean")
summary.eval.10000.mean
write.table(summary.eval.10000.mean,
            paste0("./outputs/", especie, "_", "Models2_Evaluation_Mean.csv"))

sd.i2 <- numeric(0)
sd.j2 <- numeric(2)
for (i in 1:n.algo2) {
  for (j in 1:2) {
    sd.j2[j] <-
      sd(sppModelEval.PA.10000[paste(eval.methods2[j]), "Testing.data", paste(sdm.models2[i]), ,])
  }
  sd.i2 <- c(sd.i2, sd.j2)
}

summary.eval.10000.sd <-
  data.frame(rep(sdm.models2, each = 2), rep(eval.methods2, n.algo2), sd.i2)
names(summary.eval.10000.sd) <- c("Model", "Method", "SD")
summary.eval.10000.sd
write.table(
  summary.eval.10000.sd,
  paste0("./outputs/", especie, "_", "Models2_Evaluation_SD.csv")
)



#-----------------------------#
# BUILDING OF PROJECTIONS ####
#---------------------------#

spp.projections_1 <- BIOMOD_Projection(
  modeling.output = sppModelOut.PA.equal,
  new.env = env.selected,
  proj.name = "Cur1_presente",
  selected.models = "all",
  #binary.meth = "ROC",
  output.format = ".grd"
)

spp.projections_2 <- BIOMOD_Projection(
  modeling.output = sppModelOut.PA.10000,
  new.env = env.selected,
  proj.name = "Cur2_presente",
  selected.models = "all",
  #binary.meth = "ROC",
  output.format = ".grd"
)

# save.image()
### Definir diretório onde está o arquivo proj_Cur1_presente_Occurrence.grd
projections_1 <-
  stack(
    paste0(
      "./",
      diretorio,
      "/proj_Cur1_presente/proj_Cur1_presente_Occurrence.",
      especie,
      ".grd"
    )
  )
names(projections_1)
summary.eval.equal_1<-test1
x1<-length(na.omit(summary.eval.equal_1$TSS))
summary.eval.equal_1 <-na.omit(summary.eval.equal_1)
summary.eval.equal_1 = summary.eval.equal_1[order(summary.eval.equal_1$Run),]
summary.eval.equal_1 = summary.eval.equal_1[order(summary.eval.equal_1$PA),]

summary.eval.equal_1$ID = 1:x1

sel = summary.eval.equal_1[summary.eval.equal_1[, "TSS"] > 0.700,]
sel <- na.omit(sel)

projections.1 = (subset(projections_1, sel[, "ID"]))
proj.select1 <- names(projections.1)
### Definir diretório onde está o arquivo proj_Cur2_presente_Occurrence.grd
projections_2 <-
  stack(
    paste0(
      "./",
      diretorio,
      "/proj_Cur2_presente/proj_Cur2_presente_Occurrence.",
      especie,
      ".grd"
    )
  )
names(projections_2)
summary.eval.10000_1<-test2
x2<-length(na.omit(summary.eval.10000_1$TSS))
summary.eval.10000_1 <-na.omit(summary.eval.10000_1)
summary.eval.10000_1 = summary.eval.10000_1[order(summary.eval.10000_1$Run),]
summary.eval.10000_1 = summary.eval.10000_1[order(summary.eval.10000_1$PA),]
summary.eval.10000_1$ID = 1:x2

sel2 = summary.eval.10000_1[summary.eval.10000_1[, "TSS"] > 0.700,]
sel2 <- na.omit(sel2)

projections.2 = (subset(projections_2, sel2[, "ID"]))
proj.select2 <- names(projections.2)
#-----------------------------------------------#
# Mean of the models by algorithm (Present) ####
#---------------------------------------------#
projections.all1 <- stack(projections.1)
 
projections.all2 <- stack(projections.2)


#--------------------------------#
# Ensemble - Current Climate ####
#------------------------------#
all.pres<-stack(projections.1, projections.2)

# Regression

RG<-c("GLM", "GAM", "FDA", "MARS","MAXENT.Tsuruoka")
fam.reg<-stack()
for (l in 1:length(RG)) {
  fam.reg<- stack(fam.reg, subset(all.pres, grep(RG[l], names(all.pres))))
}
fam.reg
fam.reg.m<-mean(fam.reg)
writeRaster(
  fam.reg.m,
  filename = paste0("./outputs/", especie, "_", "Regression - Current Climate.tif"),
  format = "GTiff",
  overwrite = TRUE
)


# Machine Learning

MC<-c("MAXENT.Phillips", "RF", "ANN","GBM", "CTA")
fam.mac<-stack()
for (l in 1:length(MC)) {
  fam.mac<- stack(fam.mac, subset(all.pres, grep(MC[l], names(all.pres))))
}
fam.mac
fam.mac.m<-(mean(fam.mac))
writeRaster(
  fam.mac.m,
  filename = paste0("./outputs/", especie, "_", "Machine - Current Climate.tif"),
  format = "GTiff",
  overwrite = TRUE
)

# All

# try({
projections.all.mean <-
  mean(fam.reg.m,fam.mac.m) / 1000

writeRaster(
  projections.all.mean,
  filename = paste0("./outputs/", especie, "_", "Ensemble - Current Climate.tif"),
  format = "GTiff",
  overwrite = TRUE
)
# })



#--------------------------#
# Scores ROC Threshold ####
#------------------------#

scores_ROC_equal<-subset(sel, select = c(Model, Se))
scores_ROC_equal[scores_ROC_equal=='-Inf']<-NA
scores_ROC_equal[scores_ROC_equal=='Inf']<-NA
scores_ROC_equal<-na.omit(scores_ROC_equal)
write.table(scores_ROC_equal, paste0("./outputs/",especie, "_", "scores_equal_.csv"))


## Evaluation Scores of the  Projections with PA.10000
scores_ROC_10000<-subset(sel2, select = c(Model, Se))
scores_ROC_10000[scores_ROC_10000=='-Inf']<-NA
scores_ROC_10000[scores_ROC_10000=='Inf']<-NA
scores_ROC_10000<-na.omit(scores_ROC_10000)
write.table(scores_ROC_10000, paste0("./outputs/",especie, "_", "scores_10000_.csv"))



#Scores mean
t<-rbind(scores_ROC_equal, scores_ROC_10000)
(score.1<-mean(sel$Se))
(score.2<-mean(sel2$Se))
(score.all<-(mean(cbind(score.1,score.2)/100)))
# write.table(th_mean, paste0("./outputs/",especie, "_", "scores_mean.csv"))
# Regression
fam.reg.d<-NULL
for (l in 1:length(RG)) {
  fam.reg.d<- rbind(fam.reg.d, subset(t, Model== RG[l], select = c(Model, Se)))
}
fam.reg.d.m<-mean(fam.reg.d$Se)

# Machine Learning
fam.mac.d<-NULL
for (l in 1:length(MC)) {
  fam.mac.d<- rbind(fam.mac.d, subset(t, Model== MC[l], select = c(Model, Se)))
}
fam.mac.d.m<-mean(fam.mac.d$Se)

# score mean
(s.m<-mean(fam.reg.d.m,fam.mac.d.m)/100)
#-------------------------------------------------------#
# Binary models by each algorithm (Current Climate) ####
#-----------------------------------------------------#
{th<- function(x,y){
  if("RasterLayer" %in% class(x)){ 
    v<-as.data.frame(x, h=T,xy=F)
    v[v=='0.000']<-NA
    v.l<-na.omit(v)
    (vlen<-length(v.l))
    n<-raster::ncell(x)
    (PR<-vlen/n) # PR
  }else{ 
    cat("x need be raste layer object")
  }
  if("numeric" %in% class(y)){
    (Se<-y) #Sencitivity 0 to 1
    (VDl <- Se-PR)
  }else stop( # VDI
    cat("y need be numeric object"))
  PA <- convertToPA(
    x,
    PA.method = "probability",
    prob.method = "logistic",
    beta = VDl,
    alpha = -0.05,
    plot = T
  )
}
}


#---------------------#          
# Ensenble Binary ####
#-------------------#

Convert.p<-th(projections.all.mean,s.m)
projections.binary.all <- Convert.p$pa.raster
writeRaster(
  projections.binary.all,
  filename = paste0("./outputs/", especie, "_","Ensemble Binary - Current Climate.tif"),
  format = "GTiff",
  overwrite = TRUE
)
# }
# {#---------------------------------#
########## FUTURE #################
#---------------------------------#
#--------------------------------------------------#          
# Loading raster layer with future projections ####
#------------------------------------------------#          
###GCM 1: CanESM2

bio50_CA <-
  list.files(
    "./CanESM2",
    pattern = ".tif$",
    full.names = TRUE
  )
bio50_CA
bio50_CA <- stack(bio50_CA)
bio50_CA <-disaggregate(bio50_CA, fact=5, fun=mean)
bio50_CA <-aggregate(bio50_CA, fact=5, fun=mean)
environment50_CA <- rescale(stack(bio50_CA))
names(environment50_CA)

###GCM 2: CSIRO-Mk3-6-0
bio50_CS <-
  list.files(
    "./CSIRO-Mk3-6-0",
    pattern = ".tif$",
    full.names = TRUE
  )
bio50_CS
bio50_CS <- stack(bio50_CS)
bio50_CS <-disaggregate(bio50_CS, fact=5, fun=mean)
bio50_CS <-aggregate(bio50_CS, fact=5, fun=mean)
environment50_CS <- rescale(stack(bio50_CS))
names(environment50_CS)

###GCM 3: IPSL-CM5A-LR

bio50_IP <-
  list.files(
    "./IPSL-CM5A-LR",
    pattern = ".tif$",
    full.names = TRUE
  )
bio50_IP
bio50_IP <- stack(bio50_IP)
bio50_IP <-disaggregate(bio50_IP, fact=5, fun=mean)
bio50_IP <-aggregate(bio50_IP, fact=5, fun=mean)
environment50_IP <-rescale(stack(bio50_IP))
names(environment50_IP)

# environment50_CA <-
#   subset(environment50_CA, c(names(env.selected)))
# environment50_CS <-
#   subset(environment50_CS, c(names(env.selected)))
# environment50_IP <-
  # subset(environment50_IP, c(names(env.selected)))

#--------------------------------------------#
# Model projection for the future - 2050 ####
#------------------------------------------#

#-----------------------------------------------------#
# Listing all objects produced for future scenarios ##
#---------------------------------------------------#
scenario.list <- list(
  environment50_CA,
  environment50_CS,
  environment50_IP
)

names(scenario.list) <- c(
  "env50_CA",
  "env50_CS",
  "env50_IP"
)

# save.image()
spp.projections.2050_1_CA <- BIOMOD_Projection(
  modeling.output = sppModelOut.PA.equal,
  new.env = environment50_CA,
  proj.name = "2050_1_CA",
  selected.models = proj.select1,
  output.format = ".grd"
)
spp.projections.2050_2_CA <- BIOMOD_Projection(
  modeling.output = sppModelOut.PA.10000,
  new.env = environment50_CA,
  proj.name = "2050_2_CA",
  selected.models = proj.select2,
  output.format = ".grd"
)

spp.projections.2050_1_CS <- BIOMOD_Projection(
  modeling.output = sppModelOut.PA.equal,
  new.env = environment50_CS,
  proj.name = "2050_1_CS",
  selected.models = proj.select1,
  output.format = ".grd"
)
spp.projections.2050_2_CS <- BIOMOD_Projection(
  modeling.output = sppModelOut.PA.10000,
  new.env = environment50_CS,
  proj.name = "2050_2_CS",
  selected.models = proj.select2,
  output.format = ".grd"
)

spp.projections.2050_1_IP <- BIOMOD_Projection(
  modeling.output = sppModelOut.PA.equal,
  new.env = environment50_IP,
  proj.name = "2050_1_IP",
  selected.models = proj.select1,
  output.format = ".grd"
)
spp.projections.2050_2_IP <- BIOMOD_Projection(
  modeling.output = sppModelOut.PA.10000,
  new.env = environment50_IP,
  proj.name = "2050_2_IP",
  selected.models = proj.select2,
  output.format = ".grd"
)


#-----------------------#          
# Stack projections ####
#---------------------#          
projections.2050_1_CA <-
  stack(
    paste0(
      "./",
      diretorio,
      "/proj_2050_1_CA/proj_2050_1_CA_Occurrence.",
      especie,
      ".grd"
    )
  )
projections.2050_1_CA <-
  subset(projections.2050_1_CA, c(names(projections.1)))
names(projections.2050_1_CA)

projections.2050_2_CA <-
  stack(
    paste0(
      "./",
      diretorio,
      "/proj_2050_2_CA/proj_2050_2_CA_Occurrence.",
      especie,
      ".grd"
    )
  )
projections.2050_2_CA <-
  subset(projections.2050_2_CA, c(names(projections.2)))
names(projections.2050_2_CA)

projections.2050_1_CS <-
  stack(
    paste0(
      "./",
      diretorio,
      "/proj_2050_1_CS/proj_2050_1_CS_Occurrence.",
      especie,
      ".grd"
    )
  )
projections.2050_1_CS <-
  subset(projections.2050_1_CS, c(names(projections.1)))
names(projections.2050_1_CS)
projections.2050_2_CS <-
  stack(
    paste0(
      "./",
      diretorio,
      "/proj_2050_2_CS/proj_2050_2_CS_Occurrence.",
      especie,
      ".grd"
    )
  )
projections.2050_2_CS <-
  subset(projections.2050_2_CS, c(names(projections.2)))
names(projections.2050_2_CS)

projections.2050_1_IP <-
  stack(
    paste0(
      "./",
      diretorio,
      "/proj_2050_1_IP/proj_2050_1_IP_Occurrence.",
      especie,
      ".grd"
    )
  )
projections.2050_1_IP <-
  subset(projections.2050_1_IP, c(names(projections.1)))
names(projections.2050_1_IP)
projections.2050_2_IP <-
  stack(
    paste0(
      "./",
      diretorio,
      "/proj_2050_2_IP/proj_2050_2_IP_Occurrence.",
      especie,
      ".grd"
    )
  )
projections.2050_2_IP <-
  subset(projections.2050_2_IP, c(names(projections.2)))
names(projections.2050_2_IP)


projections.all.ft1 <- stack(
  projections.2050_1_CA,
  projections.2050_1_CS,
  projections.2050_1_IP
)


projections.all.ft2 <- stack(
  projections.2050_2_CA,
  projections.2050_2_CS,
  projections.2050_2_IP
)

projections.all.ft<-stack(projections.all.ft1,projections.all.ft2)

#-------------------------------------------------------#          
# Consensus between Continuous Projections (future) ####           
#-----------------------------------------------------#


# Regression

fam.reg.f<-stack()
for (l in 1:length(RG)) {
  fam.reg.f<- stack(fam.reg.f, subset(projections.all.ft, grep(RG[l], names(projections.all.ft))))
}
fam.reg.f
fam.reg.m.f<-mean(fam.reg.f)
writeRaster(
  fam.reg.m.f,
  filename = paste0("./outputs/", especie, "_", "Regression - Future Climate.tif"),
  format = "GTiff",
  overwrite = TRUE
)


# Machine Learning

fam.mac.f<-stack()
for (l in 1:length(MC)) {
  fam.mac.f<- stack(fam.mac.f, subset(projections.all.ft, grep(MC[l], names(projections.all.ft))))
}
fam.mac.f
fam.mac.m.f<-mean(fam.mac.f)
writeRaster(
  fam.mac.m.f,
  filename = paste0("./outputs/", especie, "_", "Machine - Future Climate.tif"),
  format = "GTiff",
  overwrite = TRUE
)

ensemble_future_mean<-mean(fam.reg.m.f,fam.mac.m.f)/1000



writeRaster(
  ensemble_future_mean,
  filename = paste0("./outputs/", especie, "_", "Ensemble Future.tif"),
  formato = "GTiff",
  overwrite = TRUE
)

#--------------------------------------#
# Consensus maps: final binary ####
#------------------------------------#   


Convert.f<-th(ensemble_future_mean,s.m)
ensemble_future_bin <- Convert.f$pa.raster
writeRaster(
  ensemble_future_bin,
  filename = paste0(
    './outputs/',
    especie,"_",
    "Ensemble - Future Climate Binary.tif"
  ),
  format = "GTiff", overwrite = TRUE
)
#--------------------#          
# Move the files ####
#------------------#          

results<-list.files(
  "./outputs/",paste0(especie, "_"),
  full.names = TRUE
)

file.move((list.files(
  "./outputs/",paste0(especie, "_"),
  full.names = TRUE
)), (paste0("./outputs/", especie)), overwrite = TRUE)
#--------------------#          
# Time Compting ####
#------------------#    
sink("./outputs/tempo.txt", append = T)
print(especie)
print(Sys.time() - ini1)
sink()

}

#END
