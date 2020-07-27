#install.packages(c("ape","phytools","geiger","fields","animation","phangorn","rgdal","sp"))


#install ImageMagick (http://www.imagemagick.org)
#for macs 10.9 and up, make sure to have X11 installed. you can get it here http://xquartz.macosforge.org/landing/
#setwd("C:/Users/jriver58/Dropbox/JR Postdoc/SceloporusBiogeo/SceloporusBiogeo/scelop.project")
setwd("C:/Users/jriver58/Dropbox/Postdoc/Share Emilia/JR scel biogeography")
setwd("C:/Users/julio/Dropbox/Postdoc/Share Emilia/JR scel biogeography")

library(ape)
library(phytools)
library(geiger)
source("source_functions.R")
library(fields)
library(animation)
library(phangorn)
library(rgdal)
library(sp)



#par(mfrow=c(2,2))
#getmap.jr(ma=c(2, 5, 13, 20), model= "PALEOMAP", xlim=c(-180,-45), ylim=c(20, 70))

#Set bounds for analysis
bounds <- list(sigsq = c(min = 0, max = 1000000), SE = c(0, 0.1), alpha = c(min = 0, max = 150), a = c(min = -1, max = 1), slope = c(min = -100, max = 100), lambda = c(min = 0, max = 1), kappa = c(min = 0, max = 1), delta = c(min = 0, max = 10), drift = c(min = -100, max = 100))

get(load("occurrences.RData"))

#ToKeep <- c("jarrovi", "clarkii", "virgatus", "magister", "graciosus", "scalaris")
#occu <- subset(occurrences, Species==ToKeep)

#newtree <- read.nexus("scelop_timecalib_2016.nex")

 
#To.Drop <- c("Gambelia_wislizenii",  "Liolaemus_darwinii",  "Uma_notata",  "Callisaurusum_draconoides_MVZunknown",      
#"Callisaurus_draconoides_MVZ265543", "Holbrookia_maculata",    "Cophosaurus_texanus" ,                      
#"Phrynosoma_asio",  "Phrynosoma_cornutum",     "Phrynosoma_modestum",   "Phrynosoma_goodei",   "Phrynosoma_platyrhinos",                    "Phrynosoma_ditmarsi",    "Phrynosoma_douglasi",  "Phrynosoma_hernandesi",  "Phrynosoma_obiculare",  "Phrynosoma_braconnieri",   
#"Phrynosoma_taurus",  "Phrynosoma_sherbrookei", "Phrynosoma_solare", "Phrynosoma_cerroense", "Phrynosoma_blainvillii", "Phrynosoma_coronatum",                      
#"Phrynosoma_mcallii", "Uta_stansburiana", "Urosaurus_ornatus" , "Petrosaurus_thalassinus","Sceloporus_teapensis", "Sceloporus_smithi",                         
#"Sceloporus_variabilis",  "Sceloporus_chrysostictus", "Sceloporus_parvus" , "Sceloporus_couchii", "Sceloporus_grandaevus",                     
#"Sceloporus_angustus", "Sceloporus_squamosus" , "Sceloporus_carinatus", "Sceloporus_siniferus_mvz236299" , "Sceloporus_siniferus_uwbm6653" ,            
#"Sceloporus_gadoviae" ,  "Sceloporus_maculosus" , "Sceloporus_pyrocephalus_utar53473", "Sceloporus_pyrocephalus_unknown","Sceloporus_nelsoni" ,                       
#"Sceloporus_jalapae","Sceloporus_ochoterenae",  "Sceloporus_subniger_rwb0686" , "Sceloporus_aeneus"  , "Sceloporus_subniger_rwb0645",               
#"Sceloporus_bicanthalis", "Sceloporus_aurantius" ,"Sceloporus_unicanthalis", "Sceloporus_scalaris_rwb06247" , "Sceloporus_scalaris_uwbm6589",                                  
#"Sceloporus_chaneyi", "Sceloporus_samcolemani_rwb06263" , "Sceloporus_samcolemani_jjw698", "Sceloporus_brownorum", "Sceloporus_arenicolus",                     
#"Sceloporus_vandenburgianus","Sceloporus_zosteromus_ADG74" ,  "Sceloporus_zosteromus_ADG49", "Sceloporus_lineatulus", "Sceloporus_magister_uniformis_mvz162077",   
#"Sceloporus_magister_cephaloflavus_uwbm7395", "Sceloporus_magister_uniformis_dgm474", "Sceloporus_magister_bimaculosus_dgm924",    
#"Sceloporus_licki", "Sceloporus_orcutti_uwbm7654", "Sceloporus_orcutti_rwm798", "Sceloporus_hunsakeri" ,  "Sceloporus_asper" ,                         
#"Sceloporus_macdougalli" ,  "Sceloporus_aureolus_rvt54" ,  "Sceloporus_aureolus_jac22409", "Sceloporus_mucronatus" , "Sceloporus_dugesii",                        
#"Sceloporus_minor", "Sceloporus_serrifer_utar39870", "Sceloporus_ornatus", "Sceloporus_cyanostictus" , "Sceloporus_poinsettii",                     
#"Sceloporus_cyanogenys" , "Sceloporus_lineolateralis", "Sceloporus_insignis_anmo1130", "Sceloporus_insignis", "Sceloporus_torquatus_uwbm6600",             
#"Sceloporus_torquatus_uogv2526",  "Sceloporus_bulleri", "Sceloporus_megalepidurus", "Sceloporus_megalepidurus_pictus", "Sceloporus_anahuacus",                      
#"Sceloporus_palaciosi", "Sceloporus_heterolepis",  "Sceloporus_shannonorum" ,  "Sceloporus_grammicus_microlepidotus",       
#"Sceloporus_grammicus", "Sceloporus_melanorhinus", "Sceloporus_horridus" , "Sceloporus_spinosus" , "Sceloporus_formosus_scitulus",              
#"Sceloporus_adleri", "Sceloporus_druckercolini", "Sceloporus_stejnegeri", "Sceloporus_subpictus",  "Sceloporus_cryptus",                        
#"Sceloporus_formosus_anmo1248", "Sceloporus_formosus_rvt76",  "Sceloporus_smaragdinus", "Sceloporus_malachiticus", "Sceloporus_taeniocnemis_mvz4213",           
#"Sceloporus_acanthinus_anmo1932" , "Sceloporus_internasalis",  "Sceloporus_cautus",   "Sceloporus_exsul" ,  "Sceloporus_olivaceus" ,                     
#"Sceloporus_woodi" , "Sceloporus_tristicus" , "Sceloporus_cowelsi" , "Sceloporus_consobrinus", "Sceloporus_virgatus" , "Sceloporus_occidentalis_uwbm6281" ,         
#"Sceloporus_occidentalis_mvz245697", "Sceloporus_edwardtaylori"  ,   "Sceloporus_merriami")


#newtree2 <- drop.tip(newtree, tip=To.Drop)


get(load("beastLeache.RData"))

#ToDrop <- c(
#"arenicolus",      "bicanthalis",     "cautus",          "couchii",         "cryptus",        
#"dugesii",         "edwardtaylori90", "formosus",        "gadovi",          "grammicus",      
#"grandaevus",      "heterolepis",     "horridus",        "hunsakeri",       "jalapae",               
#"licki",           "lineolateralis",  "macdougalli",     "maculosus",       "malachiticus",   
#"megalepidurus",   "melanorhinus",    "merriami",        "mucronatus",      "occidentalis",    "ochoterrane",    
#"olivaceus",       "orcutti",         "ornatus",         "palaciosi",       "parvus",          "pictus",         
#"poinsettii",      "pyrocephalus",    "siniferus",       "smithii",         "spinosus",       
#"stejnegeri",      "subpictus",       "taeniocnemis",    "torquatus",       "utiformis",      
#"vandenburgianus", "variabilis",      "woodi",           "zosteromus", "virgatus")

tree <- read.nexus("consensus_tree.nxs")
tree <- c(tree, tree, tree, tree, tree, tree, tree, tree, tree, tree,
tree, tree, tree, tree, tree, tree, tree, tree, tree, tree,
tree, tree, tree, tree, tree, tree, tree, tree, tree, tree,
tree, tree, tree, tree, tree, tree, tree, tree, tree, tree,
tree, tree, tree, tree, tree, tree, tree, tree, tree, tree,
tree, tree, tree, tree, tree, tree, tree, tree, tree, tree,
tree, tree, tree, tree, tree, tree, tree, tree, tree, tree,
tree, tree, tree, tree, tree, tree, tree, tree, tree, tree,
tree, tree, tree, tree, tree, tree, tree, tree, tree, tree,
tree, tree, tree, tree, tree, tree, tree, tree, tree, tree)
which_run <- sample(1:length(tree), 100)
ex_mytree <- tree[which_run]



#This run is with NO FOSSILS
#Runs may take a while, so do a few trial runs before you commit to many permutations.
trialest <- ppgm(occurrences = occurrences, fossils = F, trees = ex_mytree, fossils.edges = FALSE, model = "estimate", permut = 1, which.biovars = c(1, 4, 15), path = "scratch/p_", plot.TraitGram = T, plot.AnimatedMaps = F, plot.GeoRates = F)
trialest1 <- ppgm(occurrences = occurrences, fossils = F, trees = ex_mytree, fossils.edges = FALSE, model = "estimate", permut = 1, which.biovars = c(6), path = "scratch/p_", plot.TraitGram = T, plot.AnimatedMaps = F, plot.GeoRates = F)
#for one trait and one model, use code below to get the aicc mean for all 100 trees
mean(unlist(sapply(1:length(trialest1$model_min), function(trees) trialest1$model_min[[trees]][[1]][[1]]$fitted['aicc'])))

save(trialest, file="noFossilB1B4B15JR.Rdata")
save(trialest1, file="noFossilB6JR.Rdata")

ex_mytree <- lapply(1:100, function(x) trialest$model_min[[x]][[1]][[1]]$phy)

#Make table for results
models <- c("BM", "OU", "EB")
for(trees in 1:length(trialest$model_min)){
    for(traits in 1:length(trialest$model_min[[1]][[1]])){
      temp_min <- cbind(unlist(sapply(models, function(z) trialest$model_min[[trees]][[1]][[traits]]$fitted[[z]]['aicc'])))
      temp_max <- cbind(unlist(sapply(models, function(z) trialest$model_max[[trees]][[1]][[traits]]$fitted[[z]]['aicc'])))
      colnames(temp_min) <- colnames(temp_max) <- paste("aicc", trees, 1, traits, sep = "")
      if(trees == 1){
        print_table_min <- temp_min
        print_table_max <- temp_max
      }
      else {
        print_table_min <- cbind(print_table_min, temp_min)
        print_table_max <- cbind(print_table_max, temp_max)
      }
      rownames(print_table_min) <- rownames(print_table_max) <- models
    }
}
print_table_min <- print_table_min[-4,]
print_table_max <- print_table_max[-4,]

for(trees in 1:length(trialest1$model_min)){
      temp_min <- cbind(unlist(sapply(models, function(z) trialest1$model_min[[trees]][[1]][[1]]$fitted[[z]]['aicc'])))
      temp_max <- cbind(unlist(sapply(models, function(z) trialest1$model_max[[trees]][[1]][[1]]$fitted[[z]]['aicc'])))
      colnames(temp_min) <- colnames(temp_max) <- paste("aicc", trees, 1, 1, sep = "")
      if(trees == 1){
        print_table1_min <- temp_min
        print_table1_max <- temp_max
      }
      else {
        print_table1_min <- cbind(print_table1_min, temp_min)
        print_table1_max <- cbind(print_table1_max, temp_max)
      }
      rownames(print_table_min) <- rownames(print_table_max) <- models[-4]
}
#combine all phylo and permuts to get average aicc by trait
print(cbind(MAT = apply(print_table_min[, seq(1, dim(print_table_min)[2], 3)], 1, min), TS = apply(print_table_min[, seq(2, dim(print_table_min)[2],3)], 1, min), MTCM = apply(print_table1_min, 1, min), PS = apply(print_table_min[, seq(3, dim(print_table_min)[2], 3)], 1, min)), digits = 5)
print(cbind(MAT = apply(print_table_max[, seq(1, dim(print_table_max)[2], 3)], 1, min), TS = apply(print_table_max[, seq(2, dim(print_table_max)[2],3)], 1, min), MTCM = apply(print_table1_max, 1, min), PS = apply(print_table_max[, seq(3, dim(print_table_max)[2], 3)], 1, min)), digits = 5)





#NOW INCLUDING FOSSILS
#permuts only matter when including fossils!!!
permuts <- 1
load("fossils.RData")
load("paleoclimate.RData")
fossils$early_age <- ceiling(fossils$early_age)
fossils$late_age <- ceiling(fossils$late_age)
fossilRanges <- (apply(fossils[, 3:4], 1, function(x) {x[1] - x[2]}) + 1)

#get fossil table in proper format
manipulatedFosssils <- array(NA, dim = c(sum(fossilRanges), 4))
z <- 1
for(i in 1:length(fossils[, 1])){
  for(j in 1:fossilRanges[i]){
    manipulatedFosssils[z, ] <- unlist(c(fossils$late_age[i] + j - 1, fossils[i, c(5, 6, 2)]))
    z <- z + 1
  }
}
#Note: if supplying biovarFossils, supply all 19 or will throw an error
biovarFossils <- getBioclimVars(manipulatedFosssils[, 1:3], which.biovars = 1:19)

#find all the generic scelop, remove specific scelop
fossilsedges <- gsub(1, "", manipulatedFosssils[, 4])
fossilsedges[which(fossilsedges == "")] <- NA
fossilsedges[which(fossilsedges != T)] <- "undulatus"
#cut down fossil sample with subsample
ex_fossils <- biovarFossils[is.na(fossilsedges), ]
which.biovars = c(2, 3, 4, 6, 7, 8, 9, 11, 14, 15, 19)
ex_fossils <- as.data.frame(ex_fossils)

trialestWF <- ppgm(occurrences = occurrences, fossils = ex_fossils, trees = ex_mytree, fossils.edges = F, model = "estimate", permut = permuts, which.biovars = c(1, 4, 15), path = "scratch/q_", plot.TraitGram = TRUE, plot.AnimatedMaps = FALSE, plot.GeoRates = FALSE, bounds = list(alpha = c(0, 1)), control = list(niter = 20), only.biovars = TRUE, plot.BumpChart = FALSE)
trialestWF2 <- ppgm(occurrences = occurrences, fossils = ex_fossils, trees = ex_mytree, fossils.edges = F, model = "estimate", permut = permuts, which.biovars = c(6), path = "scratch/q_", plot.TraitGram = TRUE, plot.AnimatedMaps = FALSE, plot.GeoRates = FALSE, bounds = list(alpha = c(0, 1)), control = list(niter = 20), only.biovars = TRUE, plot.BumpChart = FALSE)
trialestWF3 <- ppgm(occurrences = occurrences, fossils = ex_fossils, trees = ex_mytree, fossils.edges = F, model = "estimate", permut = permuts, which.biovars = c(15), path = "scratch/q_", plot.TraitGram = TRUE, plot.AnimatedMaps = FALSE, plot.GeoRates = FALSE, bounds = list(alpha = c(0, 1)), control = list(niter = 20), only.biovars = TRUE, plot.BumpChart = FALSE)

save(trialestWF, file="FossilB1B4B15JR.Rdata")
save(trialestWF2, file="FossilB6JR.Rdata")

##########################################
##########################################
#### START HERE
##########################################
##########################################


#skip run, save time, and load data this is michelle work 
#trialestWF <- load("Fossils100treesB1B4B15.Rdata")
#trialestWF2 <- load("Fossils100treesB6.Rdata")
#trialest <- load("noFossilsResults100treesB6.Rdata")
#trialest1 <- load("noFoss100treesB1B4B15.Rdata")

#skip run, save time, and load data this is julios NEW work 
trialestWF <- get(load("fossils10treeBall.Rdata"))
trialestWF2 <- get(load("fossils10treeB6.Rdata"))
trialest1 <- get(load("nofossils10treeB6.Rdata"))
trialest <- get(load("nofossils10treeBall.Rdata"))


#make table for results with fossils
models <- c("BM", "OU", "EB")
for(trees in 1:length(trialestWF$model_min)){
  for(permut in 1:length(trialestWF$model_min[[trees]])){
    for(traits in 1:length(trialestWF$model_min[[trees]][[permut]])){
      temp_min <- cbind(unlist(sapply(models,function(z) trialestWF$model_min[[trees]][[permut]][[traits]]$fitted[[z]]['aicc'])))
      temp_max <- cbind(unlist(sapply(models,function(z) trialestWF$model_max[[trees]][[permut]][[traits]]$fitted[[z]]['aicc'])))
      colnames(temp_min) <- colnames(temp_max) <- paste("aicc", trees, permut, traits, sep = "")
      if(trees == 1){
        print_table_min <- temp_min
        print_table_max <- temp_max
      }
      else{
        print_table_min <- cbind(print_table_min, temp_min)
        print_table_max <- cbind(print_table_max, temp_max)
      }
      rownames(print_table_min) <- rownames(print_table_max) <- models
    }
  }
}

for(trees in 1:length(trialestWF2$model_min)){
  for(permut in 1:length(trialestWF2$model_min[[trees]])){
  temp_min <- cbind(unlist(sapply(models, function(z) trialestWF2$model_min[[trees]][[permut]][[1]]$fitted[[z]]['aicc'])))
  temp_max <- cbind(unlist(sapply(models, function(z) trialestWF2$model_max[[trees]][[permut]][[1]]$fitted[[z]]['aicc'])))
  colnames(temp_min) <- paste("aicc", trees, permut, 1, sep = "")
  colnames(temp_max) <- paste("aicc", trees, permut, 1, sep = "")
  if(trees == 1){
    print_table1_min <- temp_min
    print_table1_max <- temp_max
  } else {
    print_table1_min <- cbind(print_table1_min, temp_min)
    print_table1_max <- cbind(print_table1_max, temp_max)
  }
  rownames(print_table1_min) <- rownames(print_table1_max) <- models
  }
}


#combine all phylo and permuts to get average aicc by trait
print(cbind(MAT = apply(print_table_min[, seq(1, dim(print_table_min)[2], 3)], 1, min), TS = apply(print_table_min[, seq(2, dim(print_table_min)[2], 3)], 1, min), MTCM = apply(print_table1_min, 1, min), PS = apply(print_table_min[, seq(3, dim(print_table_min)[2], 3)], 1, min)), digits = 5)
print(cbind(MAT = apply(print_table_max[, seq(1, dim(print_table_max)[2], 3)], 1, min), TS = apply(print_table_max[, seq(2, dim(print_table_max)[2], 3)], 1, min), MTCM = apply(print_table1_max, 1, min), PS = apply(print_table_max[, seq(3, dim(print_table_max)[2], 3)], 1, min)), digits = 5)






#below are the average parameter estimates for the "best model" for the runs
parameters <- array(NA, dim = c(4, 4))
colnames(parameters) <- c("B1", "B4", "B6", "B15")
rownames(parameters) <- c("min", "max", "minWF", "maxWF")
parameters[1, c(1, 2, 4)] <- rowMeans(sapply(1:length(ex_mytree), function(trees) sapply(1:3, function(traits) trialest$model_min[[trees]][[1]][[traits]]$fitted$OU$alpha)), na.rm = T)
parameters[1, 3] <- mean(sapply(1:length(ex_mytree), function(trees) trialest1$model_min[[trees]][[1]][[1]]$fitted$OU$alpha), na.rm = T)
parameters[2, c(1, 2, 4)] <- rowMeans(sapply(1:length(ex_mytree), function(trees) sapply(1:3, function(traits) trialest$model_max[[trees]][[1]][[traits]]$fitted$OU$alpha)), na.rm = T)
parameters[2, 3] <- mean(sapply(1:length(ex_mytree), function(trees) trialest1$model_max[[trees]][[1]][[1]]$fitted$OU$alpha), na.rm = T)
parameters[3, c(1, 2, 4)] <- rowMeans(sapply(1:length(ex_mytree), function(trees) sapply(1:3, function(traits) trialestWF$model_min[[trees]][[1]][[traits]]$fitted$OU$alpha)), na.rm = T)
parameters[3, 3] <- mean(sapply(1:length(ex_mytree), function(trees) trialestWF2$model_min[[trees]][[1]][[1]]$fitted$OU$alpha), na.rm = T)
parameters[4, c(1, 2, 4)] <- rowMeans(sapply(1:length(ex_mytree), function(trees) sapply(1:3, function(traits) trialest$model_max[[trees]][[1]][[traits]]$fitted$lambda$lambda)))
parameters[4, 3] <- mean(sapply(1:length(ex_mytree), function(trees) trialestWF2$model_max[[trees]][[1]][[1]]$fitted$delta$delta), na.rm = T)
print(parameters)

#plot ppgmMESS
par(mar = c(2, 1, 1, 1))

#get enveoples 1,4,15 then 6 ... put appropriate order
cem_min <- cbind(trialestWF$cem[, 1], trialestWF$cem[, 2], trialestWF$cem[, 3], trialestWF$cem[, 4], trialestWF$cem[, 5], trialestWF2$cem[, 1], trialestWF$cem[, 6], trialestWF$cem[, 7], trialestWF$cem[, 8], trialestWF$cem[, 9], trialestWF$cem[, 10], trialestWF$cem[, 11], trialestWF$cem[, 12], trialestWF$cem[, 13], trialestWF$cem[, 14], trialestWF$cem[, 15], trialestWF$cem[, 16], trialestWF$cem[, 17], trialestWF$cem[, 18] )
cem_max <- cbind(trialestWF$cem[, 37], trialestWF$cem[, 38], trialestWF$cem[, 39], trialestWF$cem[, 40], trialestWF$cem[, 41], trialestWF2$cem[, 3], trialestWF$cem[, 42], trialestWF$cem[, 43], trialestWF$cem[, 44], trialestWF$cem[, 45], trialestWF$cem[, 46], trialestWF$cem[, 47], trialestWF$cem[, 48], trialestWF$cem[, 49], trialestWF$cem[, 50], trialestWF$cem[, 51], trialestWF$cem[, 52], trialestWF$cem[, 53], trialestWF$cem[, 54])
rownames(cem_min) <- rownames(cem_max) <- rownames(trialestWF$cem)



#extract one variable from the results of multivar run
relist1 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 3, 100))[,,1,x]), list)
relist2 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 3, 100))[,,2,x]), list)
relist3 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF2$node_est), dim = c(2, 52, 1, 100))[,,1,x]), list)
relist4 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 3, 100))[,,3,x]), list)

############use internal nodes to subset
#use "node" to indicate the node that you want
node=53
relist1 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 18, 100))[,,1,node]), list)
relist2 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 18, 100))[,,2,node]), list)
relist3 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 18, 100))[,,3,node]), list)
relist4 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 18, 100))[,,4,node]), list)
relist5 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 18, 100))[,,5,node]), list)
relist6 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF2$node_est), dim = c(2, 52, 1, 100))[,,1,node]), list)
relist7 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 18, 100))[,,6,node]), list)
relist8 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 18, 100))[,,7,node]), list)
relist9 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 18, 100))[,,8,node]), list)
relist10 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 18, 100))[,,9,node]), list)
relist11 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 18, 100))[,,10,node]), list)
relist12 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 18, 100))[,,11,node]), list)
relist13 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 18, 100))[,,12,node]), list)
relist14 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 18, 100))[,,13,node]), list)
relist15 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 18, 100))[,,14,node]), list)
relist16 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 18, 100))[,,15,node]), list)
relist17 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 18, 100))[,,16,node]), list)
relist18 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 18, 100))[,,17,node]), list)
relist19 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 18, 100))[,,18,node]), list)


## change = [,,1,22]) 22 specifically for different tips
#for single tree, jarrovii=48, clarkii=6, graciosus=18, undulatus=13, slivi=25, magister=2, virgatus=22
## 33=torquatus, 20=grammicus, 45=spinosus, 25=scalaris

#MESS, will produce plots in the working directory
mess <- ppgmMESS(cem_min, cem_max , est = list( relist2, relist6, relist7, relist8, relist9, relist11, relist14, relist19), tree = ex_mytree, fossils = ex_fossils, timeslice = c(1,2,3,4,5,6,7,8,9,10), which.biovars = c(2,6,7,8,9,11,14,19), which.plot = "all")

#this will give you a dataframe with suitable habitat values for 1 species.
#remember to change to object name for sp.mess so that it is the species you are looking at
#or the object will be overwritten.
tor.mess <- ppgmMESS.scores(cem_min, cem_max , est = list( relist2, relist6, relist7, relist8, relist9, relist11, relist14, relist19), tree = ex_mytree, fossils = ex_fossils, timeslice = c(1,2,3,4,5,6,7,8,9,10), which.biovars = c(2,6,7,8,9,11,14,19), which.plot = "all")

#remember to give the object sp.mess a new name by replacing the sp with the name of the species.
#next we will compare two different species in suitable habitat


compare.ppgmMESS(timeslice = c(1,2,3,4,5,6,7,8,9,10), which.biovars = c(2,6,7,8,9,11,14,19), which.plot = "all", sp.mess1=tor.mess, sp.mess2=sca.mess)





################################################################
################################################################
#comparingspecies values to fossil values for BIOCLIM variables

dat <- trialestWF$cem[,19:36]
dat6 <- trialestWF2$cem[,2]


################################bioclim 2
t2x <- max(dat[,2])
t2n <- min(dat[,2])

f2x <- max(ex_fossils[,5])
f2n <- min(ex_fossils[,5])

max2 <- rbind(t2x, f2x)
min2 <- rbind(t2n, f2n)

bio2 <- rbind(min2, max2)
bio2 <- as.data.frame(bio2)

t2 <- dat[,2]
f2 <- ex_fossils[,5]
tfdat2 <- cbind(t2, f2)
tfdat2 <- as.data.frame(tfdat2)
names(tfdat2) <- c("Taxa", "Fossils")
boxplot(tfdat2[,1:2], main="BIOCLIM 2", ylab="BIOCLIM values")

#################################bioclim 3
t3x <- max(dat[,3])
t3n <- min(dat[,3])

f3x <- max(ex_fossils[,6])
f3n <- min(ex_fossils[,6])

max3 <- rbind(t3x, f3x)
min3 <- rbind(t3n, f3n)

bio3 <- rbind(min3, max3)
bio3 <- as.data.frame(bio3)

t3 <- dat[,3]
f3 <- ex_fossils[,6]
tfdat3 <- cbind(t3, f3)
tfdat3 <- as.data.frame(tfdat3)
names(tfdat3) <- c("Taxa", "Fossils")
boxplot(tfdat3[,1:2], main="BIOCLIM 3", ylab="BIOCLIM values")


###############################bioclim 4
t4x <- max(dat[,4])
t4n <- min(dat[,4])

f4x <- max(ex_fossils[,7])
f4n <- min(ex_fossils[,7])

max4 <- rbind(t4x, f4x)
min4 <- rbind(t4n, f4n)

bio4 <- rbind(min4, max4)
bio4 <- as.data.frame(bio4)

t4 <- dat[,4]
f4 <- ex_fossils[,7]
tfdat4 <- cbind(t4, f4)
tfdat4 <- as.data.frame(tfdat4)
names(tfdat4) <- c("Taxa", "Fossils")
boxplot(tfdat4[,1:2], main="BIOCLIM 4", ylab="BIOCLIM values")


################################bioclim 6
t6x <- max(dat6)
t6n <- min(dat6)

f6x <- max(ex_fossils[,8])
f6n <- min(ex_fossils[,8])

max6 <- rbind(t6x, f6x)
min6 <- rbind(t6n, f6n)

bio6 <- rbind(min6, max6)
bio6 <- as.data.frame(bio6)


t6 <- dat6
f6 <- ex_fossils[,8]
tfdat6 <- cbind(t6, f6)
tfdat6 <- as.data.frame(tfdat6)
names(tfdat6) <- c("Taxa", "Fossils")
boxplot(tfdat6[,1:2], main="BIOCLIM 6", ylab="BIOCLIM values")


###############################bioclim 7
t7x <- max(dat[,6])
t7n <- min(dat[,6])

f7x <- max(ex_fossils[,10])
f7n <- min(ex_fossils[,10])

max7 <- rbind(t7x, f7x)
min7 <- rbind(t7n, f7n)

bio7 <- rbind(min7, max7)
bio7 <- as.data.frame(bio7)



t7 <- dat[,6]
f7 <- ex_fossils[,10]
tfdat7 <- cbind(t7, f7)
tfdat7 <- as.data.frame(tfdat7)
names(tfdat7) <- c("Taxa", "Fossils")
boxplot(tfdat7[,1:2], main="BIOCLIM 7", ylab="BIOCLIM values")


###############################bioclim 8
t8x <- max(dat[,7])
t8n <- min(dat[,7])

f8x <- max(ex_fossils[,11])
f8n <- min(ex_fossils[,11])

max8 <- rbind(t8x, f8x)
min8 <- rbind(t8n, f8n)

bio8 <- rbind(min8, max8)
bio8 <- as.data.frame(bio8)


t8 <- dat[,7]
f8 <- ex_fossils[,11]
tfdat8 <- cbind(t8, f8)
tfdat8 <- as.data.frame(tfdat8)
names(tfdat8) <- c("Taxa", "Fossils")
boxplot(tfdat8[,1:2], main="BIOCLIM 8", ylab="BIOCLIM values")


##############################bioclim 9
t9x <- max(dat[,8])
t9n <- min(dat[,8])

f9x <- max(ex_fossils[,12])
f9n <- min(ex_fossils[,12])

max9 <- rbind(t9x, f9x)
min9 <- rbind(t9n, f9n)

bio9 <- rbind(min9, max9)
bio9 <- as.data.frame(bio9)

t9 <- dat[,8]
f9 <- ex_fossils[,12]
tfdat9 <- cbind(t9, f9)
tfdat9 <- as.data.frame(tfdat9)
names(tfdat9) <- c("Taxa", "Fossils")
boxplot(tfdat9[,1:2], main="BIOCLIM 9", ylab="BIOCLIM values")


#################################bioclim 11
t11x <- max(dat[,10])
t11n <- min(dat[,10])

f11x <- max(ex_fossils[,14])
f11n <- min(ex_fossils[,14])

max11 <- rbind(t11x, f11x)
min11 <- rbind(t11n, f11n)

bio11 <- rbind(min11, max11)
bio11 <- as.data.frame(bio11)

t11 <- dat[,10]
f11 <- ex_fossils[,14]
tfdat11 <- cbind(t11, f11)
tfdat11 <- as.data.frame(tfdat11)
names(tfdat11) <- c("Taxa", "Fossils")
boxplot(tfdat11[,1:2], main="BIOCLIM 11", ylab="BIOCLIM values")

####################################bioclim 14
t14x <- max(dat[,13])
t14n <- min(dat[,13])

f14x <- max(ex_fossils[,17])
f14n <- min(ex_fossils[,17])

max14 <- rbind(t14x, f14x)
min14 <- rbind(t14n, f14n)

bio14 <- rbind(min14, max14)
bio14 <- as.data.frame(bio14)

t14 <- dat[,13]
f14 <- ex_fossils[,17]
tfdat14 <- cbind(t14, f14)
tfdat14 <- as.data.frame(tfdat14)
names(tfdat14) <- c("Taxa", "Fossils")
boxplot(tfdat14[,1:2], main="BIOCLIM 14", ylab="BIOCLIM values")


#####################################bioclim 15
t15x <- max(dat[,14])
t15n <- min(dat[,14])

f15x <- max(ex_fossils[,18])
f15n <- min(ex_fossils[,18])

max15 <- rbind(t15x, f15x)
min15 <- rbind(t15n, f15n)

bio15 <- rbind(min15, max15)
bio15 <- as.data.frame(bio15)

t15 <- dat[,14]
f15 <- ex_fossils[,18]
tfdat15 <- cbind(t15, f15)
tfdat15 <- as.data.frame(tfdat15)
names(tfdat15) <- c("Taxa", "Fossils")
boxplot(tfdat15[,1:2], main="BIOCLIM 15", ylab="BIOCLIM values")

####################################bioclim 19
t19x <- max(dat[,18])
t19n <- min(dat[,18])

f19x <- max(ex_fossils[,22])
f19n <- min(ex_fossils[,22])

max19 <- rbind(t19x, f19x)
min19 <- rbind(t19n, f19n)

bio19 <- rbind(min19, max19)
bio19 <- as.data.frame(bio19)

t19 <- dat[,18]
f19 <- ex_fossils[,22]
tfdat19 <- cbind(t19, f19)
tfdat19 <- as.data.frame(tfdat19)
names(tfdat19) <- c("Taxa", "Fossils")
boxplot(tfdat19[,1:2], main="BIOCLIM 19", ylab="BIOCLIM values")

############################################# putting stuff together
bio.all <- rbind(bio2, bio3)
bio.all <- rbind(bio.all, bio4)
bio.all <- rbind(bio.all, bio6)
bio.all <- rbind(bio.all, bio7)
bio.all <- rbind(bio.all, bio8)
bio.all <- rbind(bio.all, bio9)
bio.all <- rbind(bio.all, bio11)
bio.all <- rbind(bio.all, bio14)
bio.all <- rbind(bio.all, bio15)
bio.all <- rbind(bio.all, bio19)


par(mfrow=c(4,3))
barplot(height=bio.all$V1[1:4], names.arg=c("Species min", "Fossil min", "Species max", "Fossil max" ), las=1, ylim=c(0, 200), ylab="BIOCLIM value", main="BIOCLIM 2", space=c(0.1, 0.1, 0.5, 0.1))
barplot(height=bio.all$V1[5:8], names.arg=c("Species min", "Fossil min", "Species max", "Fossil max"), las=1, ylim=c(0, 100), ylab="BIOCLIM value", main="BIOCLIM 3", space=c(0.1, 0.1, 0.5, 0.1))
barplot(height=bio.all$V1[9:12], names.arg=c("Species min", "Fossil min", "Species max", "Fossil max" ), las=1, ylim=c(0, 11000), ylab="BIOCLIM value", main="BIOCLIM 4", space=c(0.1, 0.1, 0.5, 0.1))
barplot(height=bio.all$V1[13:16], names.arg=c("Species min", "Fossil min", "Species max", "Fossil max" ), las=1, ylim=c(-100, 500), ylab="BIOCLIM value", main="BIOCLIM 6", space=c(0.1, 0.1, 0.5, 0.1))
barplot(height=bio.all$V1[17:20], names.arg=c("Species min", "Fossil min", "Species max", "Fossil max" ), las=1, ylim=c(0, 500), ylab="BIOCLIM value", main="BIOCLIM 7", space=c(0.1, 0.1, 0.5, 0.1))
barplot(height=bio.all$V1[21:24], names.arg=c("Species min", "Fossil min", "Species max", "Fossil max" ), las=1, ylim=c(0, 400), ylab="BIOCLIM value", main="BIOCLIM 8", space=c(0.1, 0.1, 0.5, 0.1))
barplot(height=bio.all$V1[25:28], names.arg=c("Species min", "Fossil min", "Species max", "Fossil max" ), las=1, ylim=c(0, 400), ylab="BIOCLIM value", main="BIOCLIM 9", space=c(0.1, 0.1, 0.5, 0.1))
barplot(height=bio.all$V1[29:32], names.arg=c("Species min", "Fossil min", "Species max", "Fossil max" ), las=1, ylim=c(-100, 300), ylab="BIOCLIM value", main="BIOCLIM 11", space=c(0.1, 0.1, 0.5, 0.1))
barplot(height=bio.all$V1[33:36], names.arg=c("Species min", "Fossil min", "Species max", "Fossil max" ), las=1, ylim=c(0, 60), ylab="BIOCLIM value", main="BIOCLIM 14", space=c(0.1, 0.1, 0.5, 0.1))
barplot(height=bio.all$V1[37:40], names.arg=c("Species min", "Fossil min", "Species max", "Fossil max" ), las=1, ylim=c(0, 150), ylab="BIOCLIM value", main="BIOCLIM 15", space=c(0.1, 0.1, 0.5, 0.1))
barplot(height=bio.all$V1[41:44], names.arg=c("Species min", "Fossil min", "Species max", "Fossil max" ), las=1, ylim=c(0, 1250), ylab="BIOCLIM value", main="BIOCLIM 19", space=c(0.1, 0.1, 0.5, 0.1))

####################################################################
####################################################################









#PHYSIOLOGICAL model
load("data/paleoclimate.Rdata")
timeslice <- c(2, 5, 13, 20)
Tb <- c(28, 32, 35, 38)
par(mar = c(2, 1, 1, 1))
colorscheme <- c("lightgray", "red", "#990000", "#330000")

for(tb in 1:length(Tb)){
  for(p in 1:length(timeslice)){
    hr <- 6.12 + (0.74 * (paleoclimate[[(timeslice[p])]][, 11] / 10 - Tb[tb]))
    hr[hr < 4] <- 1
    hr[hr >= 4 & hr < 7] <- 2
    hr[hr >= 7 & hr < 10] <- 3
    hr[hr >= 10] <- 4
  
    spdata <- SpatialPoints(paleoclimate[[(timeslice[p] + 1)]][, 2:3])
    proj4string(spdata)  <- CRS("+init=epsg:4326")
    spdata <- spTransform(spdata, CRS("+init=epsg:26978"))
    if(sum(ex_fossils[, 1] == (timeslice[p] + 1)) != 0){
      spfossils <- SpatialPoints(ex_fossils[, 2:3])
      proj4string(spfossils)  <- CRS("+init=epsg:4326")
      spfossils <- spTransform(spfossils, CRS("+init=epsg:26978"))
      spfossils <- spfossils[ex_fossils[, 1] == (timeslice[p] + 1), ]
    }
    #pdf(paste("Physio", Tb[tb], "Time", timeslice[p], ".pdf", sep = ""), width = 80, height = 80, pointsize = 100, useDingbats = F)
    plot(spdata, cex = 1, xlab = "", ylab = "", axes = FALSE, pch = 16, col = colorscheme[hr])
    if(sum(ex_fossils[,1] == (timeslice[p] + 1)) != 0){
      points(spfossils, cex = 2, pch = 16, col = "black") #change cex to 4 to get large points if saving to pdf
    }
    #dev.off()
  }
}




###################################################
#################################
# start of animation

############use internal nodes to subset
relist1 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 3, 100))[,,1,41]), list)
relist2 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 3, 100))[,,2,41]), list)
relist3 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF2$node_est), dim = c(2, 52, 1, 100))[,,1,41]), list)
relist4 <- lapply(lapply(1:length(ex_mytree), function(x) array(unlist(trialestWF$node_est), dim = c(2, 52, 3, 100))[,,3,41]), list)



#MESS, will produce plots in the working directory

png(file="sampgraciosus%02d.png")

	for (i in ppgmMESS(cem_min, cem_max, est = list(relist1, relist2, relist3, relist4), tree = ex_mytree, fossils = ex_fossils, timeslice = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, 17,18,19, 20), which.biovars = c(1, 4, 6, 15), which.plot = "all")){
	plot.new()
	}
dev.off()


#####################################################
######################################################
# Ancestral reconstruction nad phylognetic signal for Rivera 2019 paper
#using Bioclim 4 only to show ancestral reconstruction for most imporant variable

library(ape)
library(phytools)
library(geiger)
require(plyr)
require(data.table)


#get the tree
tree2 <- beastLeache

#for physignal

# get bioclim data
dat <- trialestWF[1] 
dat <- as.data.frame(dat)

dat2 <- trialestWF2[1]
dat2 <- as.data.frame(dat2)

alldat <- cbind(dat, dat2)



			

df <- list()			
for (i in 1:1000) {

xoo <- alldat[,2]     #this is what you have to change
names(xoo) <- rownames(alldat)
df[[i]] <- phylosig(tree2[[i]], x=xoo, method="lambda")
write.csv(df, file="physig.csv")

physig.1 <- read.csv("physig.csv")
physig.1 <- t(physig.1)			
			
df.new <- physig.1[seq(2, nrow(physig.1), 2), ]

std <- function(x) sd(x)/sqrt(length(x))
se <- std(df.new)

print(se)
print(mean(df.new))  }
			
			
	

#old loop	
#for (i in 1:57) {
#
#print(
#dat <- phylosig(tree2[[1]], alldat[,i], nsim=1000, method="lambda")
#
#     )
#	             }

				 
				 

#strong phylogenetic signal in bioclim 19 max only
			
### Ancestral reconstruction
#fig for ms	
b4 <- alldat[,4]
names(b4)=rownames(alldat)

b22 <- alldat[,22]
names(b22)=rownames(alldat)

b40 <- alldat[,40]
names(b40)=rownames(alldat)

par(mfrow=c(1,3))
contMap(tree=tree2, b4, sig=2, fsize=1)	
contMap(tree=tree2, b22, sig=2, fsize=1)
contMap(tree=tree2, b40, sig=2, fsize=1)


## fig for supp

#going to make trees using 2 Bioclim variables at a time
## for 2 and 3
b2 <- alldat$cem.2Min
names(b2)=rownames(alldat)

b3 <- alldat$cem.2Mean
names(b3)=rownames(alldat)

b4 <- alldat$cem.2Max
names(b4)=rownames(alldat)

b5 <- alldat$cem.3Min
names(b5)=rownames(alldat)

b6 <- alldat$cem.3Mean
names(b6)=rownames(alldat)

b7 <- alldat$cem.3Max
names(b7)=rownames(alldat)


par(mfrow=c(2,3))
contMap(tree=tree2, b2, sig=2, fsize=.75)	
contMap(tree=tree2, b3, sig=2, fsize=.75)
contMap(tree=tree2, b4, sig=2, fsize=.75)
contMap(tree=tree2, b5, sig=2, fsize=.75)	
contMap(tree=tree2, b6, sig=2, fsize=.75)
contMap(tree=tree2, b7, sig=2, fsize=.75)

###### for 6 and 7
b2 <- alldat$cem.6Min
names(b2)=rownames(alldat)

b3 <- alldat$cem.6Mean
names(b3)=rownames(alldat)

b4 <- alldat$cem.6Max
names(b4)=rownames(alldat)

b5 <- alldat$cem.7Min
names(b5)=rownames(alldat)

b6 <- alldat$cem.7Mean
names(b6)=rownames(alldat)

b7 <- alldat$cem.7Max
names(b7)=rownames(alldat)


par(mfrow=c(2,3))
contMap(tree=tree2, b2, sig=2, fsize=.75)	
contMap(tree=tree2, b3, sig=2, fsize=.75)
contMap(tree=tree2, b4, sig=2, fsize=.75)
contMap(tree=tree2, b5, sig=2, fsize=.75)	
contMap(tree=tree2, b6, sig=2, fsize=.75)
contMap(tree=tree2, b7, sig=2, fsize=.75)

###### for 8 and 9
b2 <- alldat$cem.8Min
names(b2)=rownames(alldat)

b3 <- alldat$cem.8Mean
names(b3)=rownames(alldat)

b4 <- alldat$cem.8Max
names(b4)=rownames(alldat)

b5 <- alldat$cem.9Min
names(b5)=rownames(alldat)

b6 <- alldat$cem.9Mean
names(b6)=rownames(alldat)

b7 <- alldat$cem.9Max
names(b7)=rownames(alldat)


par(mfrow=c(2,3))
contMap(tree=tree2, b2, sig=2, fsize=.75)	
contMap(tree=tree2, b3, sig=2, fsize=.75)
contMap(tree=tree2, b4, sig=2, fsize=.75)
contMap(tree=tree2, b5, sig=2, fsize=.75)	
contMap(tree=tree2, b6, sig=2, fsize=.75)
contMap(tree=tree2, b7, sig=2, fsize=.75)

###### for 11 and 14
b2 <- alldat$cem.11Min
names(b2)=rownames(alldat)

b3 <- alldat$cem.11Mean
names(b3)=rownames(alldat)

b4 <- alldat$cem.11Max
names(b4)=rownames(alldat)

b5 <- alldat$cem.14Min
names(b5)=rownames(alldat)

b6 <- alldat$cem.14Mean
names(b6)=rownames(alldat)

b7 <- alldat$cem.14Max
names(b7)=rownames(alldat)


par(mfrow=c(2,3))
contMap(tree=tree2, b2, sig=2, fsize=.75)	
contMap(tree=tree2, b3, sig=2, fsize=.75)
contMap(tree=tree2, b4, sig=2, fsize=.75)
contMap(tree=tree2, b5, sig=2, fsize=.75)	
contMap(tree=tree2, b6, sig=2, fsize=.75)
contMap(tree=tree2, b7, sig=2, fsize=.75)



###### for 15 and 19
b2 <- alldat$cem.15Min
names(b2)=rownames(alldat)

b3 <- alldat$cem.15Mean
names(b3)=rownames(alldat)

b4 <- alldat$cem.15Max
names(b4)=rownames(alldat)

b5 <- alldat$cem.19Min
names(b5)=rownames(alldat)

b6 <- alldat$cem.19Mean
names(b6)=rownames(alldat)

b7 <- alldat$cem.19Max
names(b7)=rownames(alldat)


par(mfrow=c(2,3))
contMap(tree=tree2, b2, sig=2, fsize=.75)	
contMap(tree=tree2, b3, sig=2, fsize=.75)
contMap(tree=tree2, b4, sig=2, fsize=.75)
contMap(tree=tree2, b5, sig=2, fsize=.75)	
contMap(tree=tree2, b6, sig=2, fsize=.75)
contMap(tree=tree2, b7, sig=2, fsize=.75)