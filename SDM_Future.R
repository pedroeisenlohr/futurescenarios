#####################################################
#### NICHE MODELS FORECASTED FOR FUTURE SCENARIOS ###
#####################################################

# Contact: Pedro V. Eisenlohr (pedro.eisenlohr@unemat.br)

###################### Acknowledgments ##########################
### Dr. Guarino Colli's team of Universidade de Brasília.
### Dr. Diogo Souza Bezerra Rocha (Instituto de Pesquisas Jardim Botânico/RJ).
### My students of Ecology Lab, mainly João Carlos Pires de Oliveira.



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

setwd(choose.dir()) #Define the folder "Habitat Suitability Models"
getwd()
dir() #Among folders, there must be "Environmental Layers" e "Shapefiles"




#----------------------------------------#
# INSTALL AND LIBRARY THE  PACKAGES ####
#--------------------------------------#

# install.packages("biomod2", dep=T)
# install.packages("car", dep=T)
# install.packages("colorRamps", dep=T)
# install.packages("doParallel", dep=T)
# install.packages("dplyr", dep=T)
# install.packages("maps", dep=T)
# install.packages("plotKML", dep=T)
# install.packages("rgdal", dep=T)
# install.packages("sdm", dep=T)
# install.packages("sqldf", dep=T)
# install.packages("testthat", dep=T)
# install.packages("usdm", dep=T)
# install.packages("FactoMineR", dep=T)
# install.packages("foreach", dep=T)
# install.packages("sdmvspecies",dep=T)
# install.packages("filesstrings",dep=T)


library(biomod2)
library(foreach)
library(raster)
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
if (dir.exists("outputs") == F) {
  dir.create("outputs")
}

### #Parallel processing
detectCores()
getDoParWorkers()
cl <- parallel::makeCluster(1, outfile =paste0("./outputs/", "log_models.log"))
#cl <- parallel::makeCluster(10, type = "MPI", outfile = "./outputs/joao.log")
registerDoParallel(cl)
getDoParWorkers()

# Memory
memory.limit(100000000)


#Download the Maxent
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

##MUDAR DIRETÓRIO DE ARQUIVOS TEMPORÁRIOS DO raster
## define the name of a temp directory where raster tmp files will be stored
# raster_tmp_dir <- "raster_tmp"
# if (dir.exists("raster_tmp") == F) {
#   ## create the directory (somente ao iniciar a modelagem):
#   dir.create(raster_tmp_dir,
#              showWarnings = F,
#              recursive = T)
# }
# 
# ## set raster options
# rasterOptions(tmpdir = raster_tmp_dir)
#Caso o PC desligue sem que você tenha concluído a modelagem,
#dê o comando acima antes de retomar a rotina.


#----------------------------------#
# IMPORTING AND CHECKING DATA ##
#--------------------------------#

# Importing climate data from present
# bio.crop <-
#   list.files(
#     "C:/Users/JBRJ/Documents/MEGA/Alunos/joao/env/Chelsa",
#     full.names = TRUE,
#     pattern = ".grd"
#   )

bio.crop <-
  list.files(
    "./Environmental layers/CHELSA/PCA",
    full.names = TRUE,
    pattern = ".tif$"
  )
#bio.crop <- list.files("./Environmental layers/CHELSA", full.names=TRUE, pattern=".grd")
bio.crop
bio.crop <- stack(bio.crop)
bio.crop<-rescale(bio.crop)
names(bio.crop)
res(bio.crop)

# Shapefiles
# neotrop <- readOGR("./Shapefiles/ShapeNeo/neotropic.shp")
# domains <-
#   readOGR("./Shapefiles/Shape_Dominios/Dominio_AbSaber.shp")


# Importing biotic data

spp<-read.table(file.choose(), header=T, sep=",")
# spp <- read.table("./testes.csv", header = T, sep = ";")
dim(spp)
head(spp, 10)

table(spp$sp)

especies <- unique(spp$sp)
especies

#------------------------#
#beginning of the modeling ####
#----------------------#

# Creating objects for models calibration

n.runs = 10 # number of RUNs
n.algo1 = 3 # number of algorithms
n.algo2 = 6 #numero de algorithms
n.conj.pa2 = 10 # set of pseudo-absences
env.selected = bio.crop
especie = especies[1] # To model without a loop, remove the '#' of this line and add it to the 'for', 'foreach' and '.packages'

#-------------------------#
#beginning of the loop####
#-----------------------#

# for(especie in especies[1]) { # For sequential loop (One species)
# foreach(especie = especies, # For parallel looping (Multiple Species)
            # .packages = c("raster", "biomod2", 'sp', "sdmvspecies", "filesstrings")) %dopar% {
# ini1 = Sys.time()
        occs <- spp[spp$sp == especie, c("lon", "lat")]
    
    # nome = strsplit(as.vector(especie), " ")
    # especie = paste(nome[[1]][1], nome[[1]][2], sep = ".")
    
    # Selecting spatially unique records #
    mask <- env.selected[[1]]
    cell <-
      cellFromXY(mask, occs[, 1:2]) # get the cell number for each point
    dup <- duplicated(cbind(occs[, 1:2], cell))
    occs <-
      occs[!dup, ]# select the records that are not duplicated
    dim(occs)
    
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
      paste0(
        './outputs/',
        especie,"_", ".csv"),
      row.names = F,
      sep = ";"
    )
})
    dim(occs)
    PA.number <- length(occs[, 1])
    PA.number #número de pontos de ocorrência espacialmente únicos
    
    diretorio = paste0("Occurrence.", especie)
    
    
    # Preparing for CTA, GBM e RF:
    sppBiomodData.PA.equal <- BIOMOD_FormatingData(
      resp.var = occurrence.resp,
      expl.var = env.selected,
      resp.xy = myRespXY,
      resp.name = diretorio,
      PA.nb.rep = n.conj.pa2,
      #número de datasets de pseudoausências
      PA.nb.absences = PA.number,
      #= número de pseudoausências = número de pontos espacialmente únicos
      PA.strategy = "disk",
      PA.dist.min = dist.min * 1000,
      PA.dist.max = dist.mean * 1000,
      na.rm = TRUE
      )
    sppBiomodData.PA.equal
    
    
    #Preparing for the other algorithms:
    sppBiomodData.PA.10000 <- BIOMOD_FormatingData(
      resp.var = occurrence.resp,
      expl.var = env.selected,
      resp.xy = myRespXY,
      resp.name = diretorio,
      PA.nb.rep = n.conj.pa2,
      PA.nb.absences = 1000,
      PA.strategy = "disk",
      PA.dist.min = dist.min * 1000,
      PA.dist.max = dist.mean * 1000,
      na.rm = TRUE
      )
    sppBiomodData.PA.10000

    myBiomodOption <-
      BIOMOD_ModelingOptions(MAXENT.Phillips = list(path_to_maxent.jar = jar))

    
    #---------------#
    # Modeling ####
    #-------------#
    
    # Partitioning test x train data:
    sppModelOut.PA.equal <- BIOMOD_Modeling(
      sppBiomodData.PA.equal,
      models = c("GBM", "CTA", "RF"),
      models.options = NULL,
      NbRunEval = n.runs,
      #número de repetições para cada algoritmo
      DataSplit = 70,
      #percentagem de pts para treino.
      Prevalence = 0.5,
      VarImport = 0,
      #caso queira avaliar a importância das variáveis, mudar para 10 ou 100 permutações
      models.eval.meth = c("TSS", "ROC"),
      SaveObj = TRUE,
      rescal.all.models = TRUE,
      do.full.models = FALSE,
      modeling.id = "spp_presente"
    )
    sppModelOut.PA.equal

    sppModelOut.PA.10000 <- BIOMOD_Modeling(
      sppBiomodData.PA.10000,
      models = c("MAXENT.Phillips", "GLM", "GAM", "ANN", "FDA", "MARS"),
      models.options = myBiomodOption,
      NbRunEval = n.runs,
      #número de repetições para cada algoritmo
      DataSplit = 70,
      #percentagem de pts para treino.
      Prevalence = 0.5,
      VarImport = 0,
      #caso queira avaliar a importância das variáveis, mudar para 10 ou 100 permutações
      models.eval.meth = c("TSS", "ROC"),
      SaveObj = TRUE,
      rescal.all.models = TRUE,
      do.full.models = FALSE,
      modeling.id = "spp_presente"
    )
    sppModelOut.PA.10000
    
    #---------------------------------#
    # EVALUATE MODELS USING BIOMOD2 ##
    #-------------------------------#
    
    # Please see:
    # http://www.cawcr.gov.au/projects/verification/#Methods_for_dichotomous_forecasts
    
    
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
    sdm.models1 <- c("GBM", "CTA", "RF") #3 models
    sdm.models1
    eval.methods1 <- c("TSS", "ROC") #2 evaluation methods
    eval.methods1
    
    
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
    write.table(summary.eval.equal.mean,"Models1_Evaluation_Mean.csv")
    
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
    
    
    sdm.models2 <-
      c("MAXENT.Phillips", "GLM", "GAM", "ANN", "FDA", "MARS") #7 models
    sdm.models2
    eval.methods2 <- c("TSS", "ROC") #2 evaluation methods
    eval.methods2
    
    
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
    write.table(summary.eval.10000.mean,"Models2_Evaluation_Mean.csv")
    
    
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
    # x1<-length(names(projections_1))
    x1<-length(na.omit(summary.eval.equal$Model))
    summary.eval.equal <-na.omit(summary.eval.equal)
    summary.eval.equal = summary.eval.equal[order(summary.eval.equal$PA),]
    summary.eval.equal = summary.eval.equal[order(summary.eval.equal$Run),]
    summary.eval.equal$ID = 1:x1


    summary.eval.equal$ID = 1:x1
    
    sel = summary.eval.equal[summary.eval.equal[, "TSS"] > 0.400,]
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
    # x2<-length(names(projections_2))
    x2<-length(na.omit(summary.eval.10000$Model))
    summary.eval.10000 <-na.omit(summary.eval.10000)
    summary.eval.10000 = summary.eval.10000[order(summary.eval.10000$PA),]
    summary.eval.10000 = summary.eval.10000[order(summary.eval.10000$Run),]
    summary.eval.10000$ID = 1:x2
    
    sel2 = summary.eval.10000[summary.eval.10000[, "TSS"] > 0.400,]
    sel2 <- na.omit(sel2)
    
    projections.2 = (subset(projections_2, sel2[, "ID"]))
    proj.select2 <- names(projections.2)

    #-----------------------------------------------#
    # Mean of the models by algorithm (Present) ####
    #---------------------------------------------#
try({
    projections.all1 <- stack(projections.1)#, projections_2)
    
    projections.RF.all <-
      subset(projections.all1, grep("RF", names(projections.all1)))
    projections.RF.mean <-
      mean(projections.RF.all) / 10
    writeRaster(
      projections.RF.mean,
      filename = paste0("./outputs/", especie, "_", "Current Climate_RF.tif"),
      formato = "GTiff",
      overwrite = TRUE
    )
})
    try({
    projections.GBM.all <-
      subset(projections.all1, grep("GBM", names(projections.all1)))
    projections.GBM.mean <-
      mean(projections.GBM.all) / 10
    writeRaster(
      projections.GBM.mean,
      filename = paste0("./outputs/", especie, "_", "Current Climate_GBM.tif"),
      formato = "GTiff",
      overwrite = TRUE
    )
})
    try({
    projections.CTA.all <-
      subset(projections.all1, grep("CTA", names(projections.all1)))
    projections.CTA.mean <-
      mean(projections.CTA.all) / 10
    writeRaster(
      projections.CTA.mean,
      filename = paste0("./outputs/", especie, "_", "Current Climate_CTA.tif"),
      formato = "GTiff",
      overwrite = TRUE
    )
})

  try({    
  projections.all2 <- stack(projections.2)
    projections.GLM.all <-
      subset(projections.all2, grep("GLM", names(projections.all2)))
    projections.GLM.mean <-
      mean(projections.GLM.all) / 10
    writeRaster(
      projections.GLM.mean,
      filename = paste0("./outputs/", especie, "_", "Current Climate_GLM.tif"),
      formato = "GTiff",
      overwrite = TRUE
    )
})
    try({ 
    projections.GAM.all <-
      subset(projections.all2, grep("GAM", names(projections.all2)))
    projections.GAM.mean <-
      mean(projections.GAM.all) / 10
    writeRaster(
      projections.GAM.mean,
      filename = paste0("./outputs/", especie, "_", "Current Climate_GAM.tif"),
      formato = "GTiff",
      overwrite = TRUE
    )
    })
    try({ 
    projections.ANN.all <-
      subset(projections.all2, grep("ANN", names(projections.all2)))
    projections.ANN.mean <-
      mean(projections.ANN.all) / 10
    writeRaster(
      projections.ANN.mean,
      filename = paste0("./outputs/", especie, "_", "Current Climate_ANN.tif"),
      formato = "GTiff",
      overwrite = TRUE
    )
    })
    # try({ 
    # projections.SRE.all <-
    #   subset(projections.all2, grep("SRE", names(projections.all2)))
    # projections.SRE.mean <-
    #   mean(projections.SRE.all) / 10
    # writeRaster(
    #   projections.SRE.mean,
    #   filename = paste0("./outputs/", especie, "_", "Current Climate_SRE.tif"),
    #   formato = "GTiff",
    #   overwrite = TRUE
    # )
    # })
    try({ 
    projections.MARS.all <-
      subset(projections.all2, grep("MARS", names(projections.all2)))
    projections.MARS.mean <-
      mean(projections.MARS.all) / 10
    writeRaster(
      projections.MARS.mean,
      filename = paste0("./outputs/", especie, "_", "Current Climate_MARS.tif"),
      formato = "GTiff",
      overwrite = TRUE
    )
    })
    try({ 
    projections.FDA.all <-
      subset(projections.all2, grep("FDA", names(projections.all2)))
    projections.FDA.mean <-
      mean(projections.FDA.all) / 10
    writeRaster(
      projections.FDA.mean,
      filename = paste0("./outputs/", especie, "_", "Current Climate_FDA.tif"),
      formato = "GTiff",
      overwrite = TRUE
    )
    })
    try({ 
    projections.MAXENT.Phillips.all <-
      subset(projections.all2, grep("MAXENT.Phillips", names(projections.all2)))
    projections.MAXENT.Phillips.mean <-
      mean(projections.MAXENT.Phillips.all) / 10
    writeRaster(
      projections.MAXENT.Phillips.mean,
      filename = paste0(
        "./outputs/",
        especie,
        "_",
        "Current Climate_MAXENT.Phillips.tif"
      ),
      formato = "GTiff",
      overwrite = TRUE
    )
    })

    #--------------------------------#
    # Ensemble - Current Climate ####
    #------------------------------#
# try({
    projections.all.mean <-
      mean(projections.1, projections.2) / 10# / dim(stack(projections_1, projections_2))[3]
    
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
    
    (scores_equal <- get_evaluations(sppModelOut.PA.equal))
          scores_ROC_equal <-
            as.numeric(scores_equal["ROC", "Cutoff", , , ])
          scores_ROC_equal[scores_ROC_equal=='-Inf']<-NA
          write.table(scores_ROC_equal, paste0("./outputs/",especie, "_", "scores_equal_.csv"))
          
          try({ 
            ##Scores GBM
            scores_ROC_GBM <-
              as.numeric(scores_equal["ROC", "Cutoff", "GBM", , ])
            scores_ROC_GBM[scores_ROC_GBM=="-Inf"]<-NA
            scores_ROC_GBM
            scores_ROC_GBM <- as.numeric(na.exclude(scores_ROC_GBM))
            th_GBM <- mean(scores_ROC_GBM) / 10
            th_GBM
            write.table(th_GBM, paste0("./outputs/",especie, "_", "scores_ROC_GBM_.csv"))
          })
          try({ 
            ##Scores CTA
            scores_ROC_CTA <-
              as.numeric(scores_equal["ROC", "Cutoff", "CTA", , ])
            scores_ROC_CTA[scores_ROC_CTA=="-Inf"]<-NA
            scores_ROC_CTA <- na.exclude(scores_ROC_CTA)
            th_CTA <- mean(scores_ROC_CTA) / 10
            th_CTA
            write.table(th_CTA, paste0("./outputs/",especie, "_", "scores_ROC_CTA_.csv"))
          })
          try({ 
            ##Scores RF
            scores_ROC_RF <-
              as.numeric(scores_equal["ROC", "Cutoff", "RF", , ])
            scores_ROC_RF <- na.exclude(scores_ROC_RF)
            scores_ROC_RF[scores_ROC_RF=="-Inf"]<-NA
            th_RF <- mean(scores_ROC_RF) / 10
            th_RF
            write.table(th_RF, paste0("./outputs/",especie, "_", "scores_ROC_RF_.csv"))
          })
          
          ## Evaluation Scores of the  Projections with PA.10000
          (scores_10000 <- get_evaluations(sppModelOut.PA.10000))
          scores_ROC_10000 <-
            as.numeric(scores_10000["ROC", "Cutoff", , , ])
          scores_ROC_10000[scores_ROC_10000=='-Inf']<-NA
          write.table(scores_ROC_10000, paste0("./outputs/",especie, "_", "scores_10000_.csv"))

          
          try({ 
            ##Scores GLM
            scores_ROC_GLM <-
              as.numeric(scores_10000["ROC", "Cutoff", "GLM", , ])
            scores_ROC_GLM[scores_ROC_GLM=="-Inf"]<-NA
            scores_ROC_GLM <- na.exclude(scores_ROC_GLM)
            th_GLM <- mean(scores_ROC_GLM) / 10
            th_GLM
            write.table(th_GLM, paste0("./outputs/",especie, "_", "scores_ROC_GLM_.csv"))
          })
          try({ 
            ##Scores GAM
            scores_ROC_GAM <-
              as.numeric(scores_10000["ROC", "Cutoff", "GAM", , ])
            scores_ROC_GAM[scores_ROC_GAM=="-Inf"]<-NA
            scores_ROC_GAM <- na.exclude(scores_ROC_GAM)
            th_GAM <- mean(scores_ROC_GAM) / 10
            th_GAM
            write.table(th_GAM, paste0("./outputs/",especie, "_", "scores_ROC_GAM_.csv"))
          })
          try({ 
            ##Scores ANN
            scores_ROC_ANN <-
              as.numeric(scores_10000["ROC", "Cutoff", "ANN", , ])
            scores_ROC_ANN[scores_ROC_ANN=="-Inf"]<-NA
            scores_ROC_ANN <- na.exclude(scores_ROC_ANN)
            th_ANN <- mean(scores_ROC_ANN) / 10
            th_ANN
            write.table(th_ANN, paste0("./outputs/",especie, "_", "scores_ROC_ANN_.csv"))
          })
          # ##Scores SRE
          # scores_ROC_SRE <-
          #   as.numeric(scores_10000["ROC", "Cutoff", "SRE", , ])
          # scores_ROC_SRE[scores_ROC_SRE=="-Inf"]<-NA
          # scores_ROC_SRE <- na.exclude(scores_ROC_SRE)
          # th_SRE <- mean(scores_ROC_SRE) / 10
          # th_SRE
          # write.table(th_SRE, paste0("./outputs/",especie, "_", "scores_ROC_SRE_.csv"))
          try({ 
            ##Scores FDA
            scores_ROC_FDA <-
              as.numeric(scores_10000["ROC", "Cutoff", "FDA", , ])
            scores_ROC_FDA[scores_ROC_FDA=="-Inf"]<-NA
            scores_ROC_FDA <- na.exclude(scores_ROC_FDA)
            th_FDA <- mean(scores_ROC_FDA) / 10
            th_FDA
            write.table(th_FDA, paste0("./outputs/",especie, "_", "scores_ROC_FDA_.csv"))
          })
          try({ 
            ##Scores MARS
            scores_ROC_MARS <-
              as.numeric(scores_10000["ROC", "Cutoff", "MARS", , ])
            scores_ROC_MARS[scores_ROC_MARS=="-Inf"]<-NA
            scores_ROC_MARS <- na.exclude(scores_ROC_MARS)
            th_MARS <- mean(scores_ROC_MARS) / 10
            th_MARS
            write.table(th_MARS, paste0("./outputs/",especie, "_", "scores_ROC_MARS_.csv"))
          })
          try({ 
            ##Scores MAXENT.Phillips
            scores_ROC_MAXENT.Phillips <-
              as.numeric(scores_10000["ROC", "Cutoff", "MAXENT.Phillips", , ])
            scores_ROC_MAXENT.Phillips[scores_ROC_MAXENT.Phillips=="-Inf"]<-NA
            scores_ROC_MAXENT.Phillips <-
              na.exclude(scores_ROC_MAXENT.Phillips)
            th_MAXENT.Phillips <-
              mean(scores_ROC_MAXENT.Phillips) / 10
            th_MAXENT.Phillips
            write.table(th_MAXENT.Phillips, paste0("./outputs/",especie, "_", "scores_ROC_MAXENT.Phillips_.csv"))
          })
          #Scores mean
          (th_mean<-mean(c(scores_ROC_10000,scores_ROC_equal), na.rm=T)/10)
          write.table(th_mean, paste0("./outputs/",especie, "_", "scores_ROC_mean.csv"))
          
#-------------------------------------------------------#
# Binary models by each algorithm (Current Climate) ####
#-----------------------------------------------------#
          try({ 
          projections.binary.RF <- BinaryTransformation(projections.RF.mean, th_RF) #Calcular th
          class(projections.binary.RF)
          summary(values(projections.binary.RF))
          writeRaster(
            projections.binary.RF,
            filename = paste0("./outputs/", especie, "_", "Current Climate_RF Binary.tif"),
            formato = "GTiff",
            overwrite = TRUE
          )
          })
          try({ 
          projections.binary.GBM <- BinaryTransformation(projections.GBM.mean, th_GBM) #Calcular th
          class(projections.binary.GBM)
          summary(values(projections.binary.GBM))
          writeRaster(
            projections.binary.GBM,
            filename = paste0("./outputs/", especie, "_", "Current Climate_GBM Binary.tif"),
            formato = "GTiff",
            overwrite = TRUE
          )
          })
          try({ 
          projections.binary.CTA <- BinaryTransformation(projections.CTA.mean, th_CTA) #Calcular th
          class(projections.binary.CTA)
          summary(values(projections.binary.CTA))
          writeRaster(
            projections.binary.CTA,
            filename = paste0("./outputs/", especie, "_", "Current Climate_CTA Binary.tif"),
            formato = "GTiff",
            overwrite = TRUE
          )
          })
          try({ 
          projections.binary.GLM <- BinaryTransformation(projections.GLM.mean, th_GLM) #Calcular th
          class(projections.binary.GLM)
          summary(values(projections.binary.GLM))
          writeRaster(
            projections.binary.GLM,
            filename = paste0("./outputs/", especie, "_", "Current Climate_GLM Binary.tif"),
            formato = "GTiff",
            overwrite = TRUE
          )
          })
          try({ 
          projections.binary.GAM <- BinaryTransformation(projections.GAM.mean, th_GAM) #Calcular th
          class(projections.binary.GAM)
          summary(values(projections.binary.GAM))
          writeRaster(
            projections.binary.GAM,
            filename = paste0("./outputs/", especie, "_", "Current Climate_GAM Binary.tif"),
            formato = "GTiff",
            overwrite = TRUE
          )
          })
          try({ 
          projections.binary.ANN <- BinaryTransformation(projections.ANN.mean, th_ANN) #Calcular th
          class(projections.binary.ANN)
          summary(values(projections.binary.ANN))
          writeRaster(
            projections.binary.ANN,
            filename = paste0("./outputs/", especie, "_", "Current Climate_ANN Binary.tif"),
            formato = "GTiff",
            overwrite = TRUE
          )
          })
          # try({ 
          # projections.binary.SRE <- BinaryTransformation(projections.SRE.mean, th_SRE) #Calcular th
          # class(projections.binary.SRE)
          # summary(values(projections.binary.SRE))
          # writeRaster(
          #   projections.binary.SRE,
          #   filename = paste0("./outputs/", especie, "_", "Current Climate_SRE Binary.tif"),
          #   formato = "GTiff",
          #   overwrite = TRUE
          # )
          # })
          try({ 
          projections.binary.MARS <- BinaryTransformation(projections.MARS.mean, th_MARS) #Calcular th
          class(projections.binary.MARS)
          summary(values(projections.binary.MARS))
          writeRaster(
            projections.binary.MARS,
            filename = paste0("./outputs/", especie, "_", "Current Climate_MARS Binary.tif"),
            formato = "GTiff",
            overwrite = TRUE
          )
          })
          try({ 
          projections.binary.FDA <- BinaryTransformation(projections.FDA.mean, th_FDA) #Calcular th
          class(projections.binary.FDA)
          summary(values(projections.binary.FDA))
          writeRaster(
            projections.binary.FDA,
            filename = paste0("./outputs/", especie, "_", "Current Climate_FDA Binary.tif"),
            formato = "GTiff",
            overwrite = TRUE
          )
          })
          try({ 
          projections.binary.MAXENT.Phillips <- BinaryTransformation(projections.MAXENT.Phillips.mean, th_MAXENT.Phillips) #Calcular th
          class(projections.binary.MAXENT.Phillips)
          summary(values(projections.binary.MAXENT.Phillips))
          writeRaster(
            projections.binary.MAXENT.Phillips,
            filename = paste0("./outputs/", especie, "_", "Current Climate_MAXENT.Phillips Binary.tif"),
            formato = "GTiff",
            overwrite = TRUE
          )
          })

          
          
#---------------------#          
# Ensenble Binary ####
#-------------------#          
          
          projections.binary.all <- BinaryTransformation(projections.all.mean, th_mean)
          writeRaster(
            projections.binary.all,
            filename = paste0("./outputs/", especie, "_","Ensemble Binary - Current Climate.tif"),
            format = "GTiff",
            overwrite = TRUE
          )
          
          # #--------------------#          
          # # Move the files ####
          # #------------------#          
          # 
          # results<-list.files(
          #   "./outputs/",paste0(especie, "_"),
          #   full.names = TRUE
          # )
          # 
          # file.move((list.files(
          #   "./outputs/",paste0(especie, "_"),
          #   full.names = TRUE
          # )), (paste0("./outputs/", especie)), overwrite = TRUE)
          # 
          # #--------------------#          
          # # Time Compting ####
          # #------------------#    
          # sink("./outputs/tempo.txt", append = T)
          # print(especie)
          # print(Sys.time() - ini)
          # sink()
          # 
          #   }
          #---------------------------------#
          ########## FUTURE #################
          #---------------------------------#
#--------------------------------------------------#          
# Loading raster layer with future projections ####
#------------------------------------------------#          
          ###GCM 1: CCSM4
          
          bio50.85_CC <-
            list.files(
              "./Environmental layers/CHELSA_Future/CCSM4",
              pattern = ".tif$",
              full.names = TRUE
            )
          bio50.85_CC
          bio50.85_CC <- stack(bio50.85_CC)
          # bio50.85_CC<-rescale(bio50.85_CC)
          environment50.85_CC <- rescale(bio50.85_CC)
          names(environment50.85_CC) #Atenção a esta sequência!
  
          ###GCM 2: CMCC_CM

          bio50.85_CM <-
            list.files(
              "./Environmental layers/CHELSA_Future/CMCC",
              pattern = ".tif$",
              full.names = TRUE
            )
          bio50.85_CM
          bio50.85_CM <- stack(bio50.85_CM)
          # bio50.85_CM<-rescale(bio50.85_CM)
           environment50.85_CM <- rescale(bio50.85_CM)
          names(environment50.85_CM) #Atenção a esta sequência!


          ###GCM 3: CSIRO_Mk3

          bio50.85_CS <-
            list.files(
              "./Environmental layers/CHELSA_Future/CSIRO",
              pattern = ".tif$",
              full.names = TRUE
            )
          bio50.85_CS
          bio50.85_CS <- stack(bio50.85_CS)
          environment50.85_CS <-rescale(bio50.85_CS)
          names(environment50.85_CS) #Atenção a esta sequência!


          ###GCM 4: GFDL_CM3

          bio50.85_GF <-
            list.files("./Environmental layers/CHELSA_Future/GFDL",
                       pattern = ".tif$",
                       full.names = TRUE)
          bio50.85_GF
          bio50.85_GF <- stack(bio50.85_GF)
          environment50.85_GF <- rescale(bio50.85_GF)
          names(environment50.85_GF) #Atenção a esta sequência!


          #GCM 5: HadGEM2

          bio50.85_HG <-
            list.files("./Environmental layers/CHELSA_Future/HadGEM2",
                       pattern = ".tif$",
                       full.names = TRUE)
          bio50.85_HG
          bio50.85_HG <- stack(bio50.85_HG)
          environment50.85_HG <- rescale(bio50.85_HG)
          names(environment50.85_HG) #Atenção a esta sequência!


          #GCM 6: MIROC-ESM

          bio50.85_MR <-
            list.files("./Environmental layers/CHELSA_Future/MIROC_ESM",
                       pattern = ".tif$",
                       full.names = TRUE)
          bio50.85_MR
          bio50.85_MR <- stack(bio50.85_MR)
          environment50.85_MR <- rescale(bio50.85_MR)
          names(environment50.85_MR) #Atenção a esta sequência!

          
          
#--------------------------------------------#
# Model projection for the future - 2050 ####
#------------------------------------------#
          
#-------------------------------------------------------#
# Listing all objects produced for future scenarios ####
#-----------------------------------------------------#
          scenario.list <- list(
            environment50.85_CC,
            environment50.85_CM,
            environment50.85_CS,
            environment50.85_GF,
            environment50.85_HG,
            environment50.85_MR
          )
          
          names(scenario.list) <- c(
            "env50.85_CC",
            "env50.85_CM",
            "env50.85_CS",
            "env50.85_GF",
            "env50.85_HG",
            "env50.85_MR"
          )
          
          # save.image()
          spp.projections.2050.rcp85_1_CC <- BIOMOD_Projection(
            modeling.output = sppModelOut.PA.equal,
            new.env = environment50.85_CC,
            proj.name = "2050.rcp85_1_CC",
            selected.models = proj.select1,
            output.format = ".grd"
          )
          spp.projections.2050.rcp85_2_CC <- BIOMOD_Projection(
            modeling.output = sppModelOut.PA.10000,
            new.env = environment50.85_CC,
            proj.name = "2050.rcp85_2_CC",
            selected.models = proj.select2,
            output.format = ".grd"
          )

          spp.projections.2050.rcp85_1_CM <- BIOMOD_Projection(
            modeling.output = sppModelOut.PA.equal,
            new.env = environment50.85_CM,
            proj.name = "2050.rcp85_1_CM",
            selected.models = proj.select1,
            output.format = ".grd"
          )
          spp.projections.2050.rcp85_2_CM <- BIOMOD_Projection(
            modeling.output = sppModelOut.PA.10000,
            new.env = environment50.85_CM,
            proj.name = "2050.rcp85_2_CM",
            selected.models = proj.select2,
            output.format = ".grd"
          )

          spp.projections.2050.rcp85_1_CS <- BIOMOD_Projection(
            modeling.output = sppModelOut.PA.equal,
            new.env = environment50.85_CS,
            proj.name = "2050.rcp85_1_CS",
            selected.models = proj.select1,
            output.format = ".grd"
          )
          spp.projections.2050.rcp85_2_CS <- BIOMOD_Projection(
            modeling.output = sppModelOut.PA.10000,
            new.env = environment50.85_CS,
            proj.name = "2050.rcp85_2_CS",
            selected.models = proj.select2,
            output.format = ".grd"
          )

          spp.projections.2050.rcp85_1_GF <- BIOMOD_Projection(
            modeling.output = sppModelOut.PA.equal,
            new.env = environment50.85_GF,
            proj.name = "2050.rcp85_1_GF",
            selected.models = proj.select1,
            output.format = ".grd"
          )
          spp.projections.2050.rcp85_2_GF <- BIOMOD_Projection(
            modeling.output = sppModelOut.PA.10000,
            new.env = environment50.85_GF,
            proj.name = "2050.rcp85_2_GF",
            selected.models = proj.select2,
            output.format = ".grd"
          )

          spp.projections.2050.rcp85_1_HG <- BIOMOD_Projection(
            modeling.output = sppModelOut.PA.equal,
            new.env = environment50.85_HG,
            proj.name = "2050.rcp85_1_HG",
            selected.models = proj.select1,
            output.format = ".grd"
          )
          spp.projections.2050.rcp85_2_HG <- BIOMOD_Projection(
            modeling.output = sppModelOut.PA.10000,
            new.env = environment50.85_HG,
            proj.name = "2050.rcp85_2_HG",
            selected.models = proj.select2,
            output.format = ".grd"
          )

          spp.projections.2050.rcp85_1_MR <- BIOMOD_Projection(
            modeling.output = sppModelOut.PA.equal,
            new.env = environment50.85_MR,
            proj.name = "2050.rcp85_1_MR",
            selected.models = proj.select1,
            output.format = ".grd"
          )
          spp.projections.2050.rcp85_2_MR <- BIOMOD_Projection(
            modeling.output = sppModelOut.PA.10000,
            new.env = environment50.85_MR,
            proj.name = "2050.rcp85_2_MR",
            selected.models = proj.select2,
            output.format = ".grd"
          )

        
#-----------------------#          
# Stack projections ####
#---------------------#          
          projections.2050.rcp85_1_CC <-
            stack(
              paste0(
                "./",
                diretorio,
                "/proj_2050.rcp85_1_CC/proj_2050.rcp85_1_CC_Occurrence.",
                especie,
                ".grd"
              )
            )
          projections.2050.rcp85_1_CC <-
            subset(projections.2050.rcp85_1_CC, c(names(projections.1)))
          names(projections.2050.rcp85_1_CC)
          
          projections.2050.rcp85_2_CC <-
            stack(
              paste0(
                "./",
                diretorio,
                "/proj_2050.rcp85_2_CC/proj_2050.rcp85_2_CC_Occurrence.",
                especie,
                ".grd"
              )
            )
          projections.2050.rcp85_2_CC <-
            subset(projections.2050.rcp85_2_CC, c(names(projections.2)))
          names(projections.2050.rcp85_2_CC)
          
          projections.2050.rcp85_1_CM <-
            stack(
              paste0(
                "./",
                diretorio,
                "/proj_2050.rcp85_1_CM/proj_2050.rcp85_1_CM_Occurrence.",
                especie,
                ".grd"
              )
            )
          projections.2050.rcp85_1_CM <-
            subset(projections.2050.rcp85_1_CM, c(names(projections.1)))
          names(projections.2050.rcp85_1_CM)
          projections.2050.rcp85_2_CM <-
            stack(
              paste0(
                "./",
                diretorio,
                "/proj_2050.rcp85_2_CM/proj_2050.rcp85_2_CM_Occurrence.",
                especie,
                ".grd"
              )
            )
          projections.2050.rcp85_2_CM <-
            subset(projections.2050.rcp85_2_CM, c(names(projections.2)))
          names(projections.2050.rcp85_2_CM)
          
          projections.2050.rcp85_1_CS <-
            stack(
              paste0(
                "./",
                diretorio,
                "/proj_2050.rcp85_1_CS/proj_2050.rcp85_1_CS_Occurrence.",
                especie,
                ".grd"
              )
            )
          projections.2050.rcp85_1_CS <-
            subset(projections.2050.rcp85_1_CS, c(names(projections.1)))
          names(projections.2050.rcp85_1_CS)
          projections.2050.rcp85_2_CS <-
            stack(
              paste0(
                "./",
                diretorio,
                "/proj_2050.rcp85_2_CS/proj_2050.rcp85_2_CS_Occurrence.",
                especie,
                ".grd"
              )
            )
          projections.2050.rcp85_2_CS <-
            subset(projections.2050.rcp85_2_CS, c(names(projections.2)))
          names(projections.2050.rcp85_2_CS)
          
          
          projections.2050.rcp85_1_GF <-
            stack(
              paste0(
                "./",
                diretorio,
                "/proj_2050.rcp85_1_GF/proj_2050.rcp85_1_GF_Occurrence.",
                especie,
                ".grd"
              )
            )
          projections.2050.rcp85_1_GF <-
            subset(projections.2050.rcp85_1_GF, c(names(projections.1)))
          names(projections.2050.rcp85_1_GF)
          projections.2050.rcp85_2_GF <-
            stack(
              paste0(
                "./",
                diretorio,
                "/proj_2050.rcp85_2_GF/proj_2050.rcp85_2_GF_Occurrence.",
                especie,
                ".grd"
              )
            )
          projections.2050.rcp85_2_GF <-
            subset(projections.2050.rcp85_2_GF, c(names(projections.2)))
          names(projections.2050.rcp85_2_GF)
          
          projections.2050.rcp85_1_HG <-
            stack(
              paste0(
                "./",
                diretorio,
                "/proj_2050.rcp85_1_HG/proj_2050.rcp85_1_HG_Occurrence.",
                especie,
                ".grd"
              )
            )
          projections.2050.rcp85_1_HG <-
            subset(projections.2050.rcp85_1_HG, c(names(projections.1)))
          names(projections.2050.rcp85_1_HG)
          projections.2050.rcp85_2_HG <-
            stack(
              paste0(
                "./",
                diretorio,
                "/proj_2050.rcp85_2_HG/proj_2050.rcp85_2_HG_Occurrence.",
                especie,
                ".grd"
              )
            )
          projections.2050.rcp85_2_HG <-
            subset(projections.2050.rcp85_2_HG, c(names(projections.2)))
          names(projections.2050.rcp85_2_HG)
          
          
          projections.2050.rcp85_1_MR <-
            stack(
              paste0(
                "./",
                diretorio,
                "/proj_2050.rcp85_1_MR/proj_2050.rcp85_1_MR_Occurrence.",
                especie,
                ".grd"
              )
            )
          projections.2050.rcp85_1_MR <-
            subset(projections.2050.rcp85_1_MR, c(names(projections.1)))
          names(projections.2050.rcp85_1_MR)
          projections.2050.rcp85_2_MR <-
            stack(
              paste0(
                "./",
                diretorio,
                "/proj_2050.rcp85_2_MR/proj_2050.rcp85_2_MR_Occurrence.",
                especie,
                ".grd"
              )
            )
          projections.2050.rcp85_2_MR <-
            subset(projections.2050.rcp85_2_MR, c(names(projections.2)))
          names(projections.2050.rcp85_2_MR)
          
          
          projections.all.ft1 <- stack(
            projections.2050.rcp85_1_CC,
            projections.2050.rcp85_1_CM,
            projections.2050.rcp85_1_CS,
            projections.2050.rcp85_1_GF,
            projections.2050.rcp85_1_HG,
            projections.2050.rcp85_1_MR
          )
          
          
          projections.all.ft2 <- stack(
            projections.2050.rcp85_2_CC,
            projections.2050.rcp85_2_CM,
            projections.2050.rcp85_2_CS,
            projections.2050.rcp85_2_GF,
            projections.2050.rcp85_2_HG,
            projections.2050.rcp85_2_MR
          )
          
          projections.all.ft<-stack(projections.all.ft1,projections.all.ft2)
#------------------------------------------------#          
### Mean of the models by algorithm (Future) ####
#----------------------------------------------#          
          try({ 
          projections.RF.all.2050.rcp85 <-
            subset(projections.all.ft1, grep("RF", names(projections.all.ft1)))
          projections.RF.mean.2050.rcp85 <-
            mean(projections.RF.all.2050.rcp85) / 10
          writeRaster(
            projections.RF.mean.2050.rcp85,
            filename = paste0("./outputs/", especie, "_", "Future Climate_RF.tif"),
            formato = "GTiff",
            overwrite = TRUE
          )
          })
          try({ 
          projections.GBM.all.2050.rcp85 <-
            subset(projections.all.ft1, grep("GBM", names(projections.all.ft1)))
          projections.GBM.mean.2050.rcp85 <-
            mean(projections.GBM.all.2050.rcp85) / 10
          writeRaster(
            projections.GBM.mean.2050.rcp85,
            filename = paste0("./outputs/", especie, "_", "Future Climate_GBM.tif"),
            formato = "GTiff",
            overwrite = TRUE
          )
          })
          
          try({           
          projections.CTA.all.2050.rcp85 <-
            subset(projections.all.ft1, grep("CTA", names(projections.all.ft1)))
          projections.CTA.mean.2050.rcp85 <-
            mean(projections.CTA.all.2050.rcp85) / 10
          writeRaster(
            projections.CTA.mean.2050.rcp85,
            filename = paste0("./outputs/", especie, "_", "Future Climate_CTA.tif"),
            formato = "GTiff",
            overwrite = TRUE
          )
          })
          
          try({ 
          projections.GLM.all.2050.rcp85 <-
            subset(projections.all.ft2, grep("GLM", names(projections.all.ft2)))
          projections.GLM.mean.2050.rcp85 <-
            mean(projections.GLM.all.2050.rcp85) / 10
          writeRaster(
            projections.GLM.mean.2050.rcp85,
            filename = paste0("./outputs/", especie, "_", "Future Climate_GLM.tif"),
            formato = "GTiff",
            overwrite = TRUE
          )
          })
          
          try({ 
          projections.GAM.all.2050.rcp85 <-
            subset(projections.all.ft2, grep("GAM", names(projections.all.ft2)))
          projections.GAM.mean.2050.rcp85 <-
            mean(projections.GAM.all.2050.rcp85) / 10
          writeRaster(
            projections.GAM.mean.2050.rcp85,
            filename = paste0("./outputs/", especie, "_", "Future Climate_GAM.tif"),
            formato = "GTiff",
            overwrite = TRUE
          )
          })
          
          try({ 
          projections.ANN.all.2050.rcp85 <-
            subset(projections.all.ft2, grep("ANN", names(projections.all.ft2)))
          projections.ANN.mean.2050.rcp85 <-
            mean(projections.ANN.all.2050.rcp85) / 10
          writeRaster(
            projections.ANN.mean.2050.rcp85,
            filename = paste0("./outputs/", especie, "_", "Future Climate_ANN.tif"),
            formato = "GTiff",
            overwrite = TRUE
          )
          })
          
          # try({
          # projections.SRE.all.2050.rcp85 <-
          #   subset(projections.all.ft2, grep("SRE", names(projections.all.ft2)))
          # projections.SRE.mean.2050.rcp85 <-
          #   mean(projections.SRE.all.2050.rcp85) / 10
          # writeRaster(
          #   projections.SRE.mean.2050.rcp85,
          #   filename = paste0("./outputs/", especie, "_", "Future Climate_SRE.tif"),
          #   formato = "GTiff",
          #   overwrite = TRUE
          # )
          # })
          
          try({ 
          projections.MARS.all.2050.rcp85 <-
            subset(projections.all.ft2, grep("MARS", names(projections.all.ft2)))
          projections.MARS.mean.2050.rcp85 <-
            mean(projections.MARS.all.2050.rcp85) / 10
          writeRaster(
            projections.MARS.mean.2050.rcp85,
            filename = paste0("./outputs/", especie, "_", "Future Climate_MARS.tif"),
            formato = "GTiff",
            overwrite = TRUE
          )
          })
          
          try({ 
          projections.FDA.all.2050.rcp85 <-
            subset(projections.all.ft2, grep("FDA", names(projections.all.ft2)))
          projections.FDA.mean.2050.rcp85 <-
            mean(projections.FDA.all.2050.rcp85) / 10
          writeRaster(
            projections.FDA.mean.2050.rcp85,
            filename = paste0("./outputs/", especie, "_", "Future Climate_FDA.tif"),
            formato = "GTiff",
            overwrite = TRUE
          )
          })
          
          try({ 
          projections.MAXENT.Phillips.all.2050.rcp85 <-
            subset(projections.all.ft2, grep("MAXENT.Phillips", names(projections.all.ft2)))
          projections.MAXENT.Phillips.mean.2050.rcp85 <-
            mean(projections.MAXENT.Phillips.all.2050.rcp85) / 10
          writeRaster(
            projections.MAXENT.Phillips.mean.2050.rcp85,
            filename = paste0("./outputs/", especie, "_", "Future Climate_MAXENT.Phillips.tif"),
            formato = "GTiff",
            overwrite = TRUE
          )
          })

#----------------------------------------------#
# Binary models by each algorithm (Future) ####
#--------------------------------------------#
          try({
            projections.binary.future.RF <- BinaryTransformation(projections.RF.mean.2050.rcp85, th_RF) #Calcular th
            class(projections.binary.future.RF)
            summary(values(projections.binary.future.RF))
            writeRaster(
              projections.binary.RF,
              filename = paste0("./outputs/", especie, "_", "Future Climate_RF Binary.tif"),
              formato = "GTiff",
              overwrite = TRUE
            )
          })
          
          try({
            projections.binary.future.GBM <- BinaryTransformation(projections.GBM.mean.2050.rcp85, th_GBM) #Calcular th
            class(projections.binary.future.GBM)
            summary(values(projections.binary.future.GBM))
            writeRaster(
              projections.binary.GBM,
              filename = paste0("./outputs/", especie, "_", "Future Climate_GBM Binary.tif"),
              formato = "GTiff",
              overwrite = TRUE
            )
          })
          
          try({
            projections.binary.future.CTA <- BinaryTransformation(projections.CTA.mean.2050.rcp85, th_CTA) #Calcular th
            class(projections.binary.future.CTA)
            summary(values(projections.binary.future.CTA))
            writeRaster(
              projections.binary.CTA,
              filename = paste0("./outputs/", especie, "_", "Future Climate_CTA Binary.tif"),
              formato = "GTiff",
              overwrite = TRUE
            )
          })
          
          try({
            projections.binary.future.GLM <- BinaryTransformation(projections.GLM.mean.2050.rcp85, th_GLM) #Calcular th
            class(projections.binary.future.GLM)
            summary(values(projections.binary.future.GLM))
            writeRaster(
              projections.binary.GLM,
              filename = paste0("./outputs/", especie, "_", "Future Climate_GLM Binary.tif"),
              formato = "GTiff",
              overwrite = TRUE
            )
          })
          
          try({
            projections.binary.future.GAM <- BinaryTransformation(projections.GAM.mean.2050.rcp85, th_GAM) #Calcular th
            class(projections.binary.future.GAM)
            summary(values(projections.binary.future.GAM))
            writeRaster(
              projections.binary.GAM,
              filename = paste0("./outputs/", especie, "_", "Future Climate_GAM Binary.tif"),
              formato = "GTiff",
              overwrite = TRUE
            )
          })
          
          try({
            projections.binary.future.ANN <- BinaryTransformation(projections.ANN.mean.2050.rcp85, th_ANN) #Calcular th
            class(projections.binary.future.ANN)
            summary(values(projections.binary.future.ANN))
            writeRaster(
              projections.binary.ANN,
              filename = paste0("./outputs/", especie, "_", "Future Climate_ANN Binary.tif"),
              formato = "GTiff",
              overwrite = TRUE
            )
          })
          
          # try({
          #   projections.binary.future.SRE <- BinaryTransformation(projections.SRE.mean.2050.rcp85, th_SRE) #Calcular th
          #   class(projections.binary.future.SRE)
          #   summary(values(projections.binary.future.SRE))
          #   writeRaster(
          #     projections.binary.SRE,
          #     filename = paste0("./outputs/", especie, "_", "Future Climate_SRE Binary.tif"),
          #     formato = "GTiff",
          #     overwrite = TRUE
          #   )
          # })
          
          try({
            projections.binary.future.MARS <- BinaryTransformation(projections.MARS.mean.2050.rcp85, th_MARS) #Calcular th
            class(projections.binary.future.MARS)
            summary(values(projections.binary.future.MARS))
            writeRaster(
              projections.binary.MARS,
              filename = paste0("./outputs/", especie, "_", "Future Climate_MARS Binary.tif"),
              formato = "GTiff",
              overwrite = TRUE
            )
          })
          
          try({
            projections.binary.future.FDA <- BinaryTransformation(projections.FDA.mean.2050.rcp85, th_FDA) #Calcular th
            class(projections.binary.future.FDA)
            summary(values(projections.binary.future.FDA))
            writeRaster(
              projections.binary.FDA,
              filename = paste0("./outputs/", especie, "_", "Future Climate_FDA Binary.tif"),
              formato = "GTiff",
              overwrite = TRUE
            )
          })
          
          try({
            projections.binary.future.MAXENT.Phillips <- BinaryTransformation(projections.MAXENT.Phillips.mean.2050.rcp85, th_MAXENT.Phillips) #Calcular th
            class(projections.binary.future.MAXENT.Phillips)
            summary(values(projections.binary.future.MAXENT.Phillips))
            writeRaster(
              projections.binary.MAXENT.Phillips,
              filename = paste0("./outputs/", especie, "_", "Future Climate_MAXENT.Phillips Binary.tif"),
              formato = "GTiff",
              overwrite = TRUE
            )
          })
          
          
          
#-------------------------------------------------------#          
# Consensus between Continuous Projections (future) ####           
#-----------------------------------------------------#
          # ensemble2050rcp8.5<-stack(projections.RF.all.2050.rcp85, projections.GBM.all.2050.rcp85, projections.CTA.all.2050.rcp85,
          #                      projections.GLM.all.2050.rcp85, projections.GAM.all.2050.rcp85, projections.ANN.all.2050.rcp85,
          #                      projections.SRE.all.2050.rcp85, projections.MARS.all.2050.rcp85, projections.FDA.all.2050.rcp85,
          #                      projections.MAXENT.Phillips.all.2050.rcp85)
        
          ensemble_2050rcp8.5_mean<-mean(projections.all.ft)/10
          
          
          
            writeRaster(
              ensemble_2050rcp8.5_mean,
            filename = paste0("./outputs/", especie, "_", "Ensemble Future 2050.tif"),
            formato = "GTiff",
            overwrite = TRUE
          )
        
        
#--------------------------------------#
# Consensus maps: final binary ####
#------------------------------------#          
          ensemble_2050rcp8.5_bin <-
            BinaryTransformation(ensemble_2050rcp8.5_mean, th_mean)
          writeRaster(
            ensemble_2050rcp8.5_bin,
            filename = paste0(
              './outputs/',
              especie,"_",
              "Ensemble - Final Future Climate Binary_2050_rcp8.5.tif"
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
          # Time Computing ####
          #------------------#    
          sink("./outputs/tempo.txt", append = T)
          print(especie)
          print(Sys.time() - ini1)
          sink()
          
        }

#END
