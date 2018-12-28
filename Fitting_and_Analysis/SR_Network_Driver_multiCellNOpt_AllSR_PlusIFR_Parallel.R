# Driver script for multiCellNopt package - SR dataset

rm(list=ls()) # clear all variables
cat("\014") # clear screen
if (length(dev.list()>0)) {dev.off()} # clear figures

# ===== Automated optimisation pipeline ===== #

NrRoundsNeed  <- 5
NrInRoundRep  <- 3
TimeRoundRep  <- 3600 # in seconds
MultiModPlot  <- TRUE

# IndivCpd <- 4 # BFA (Brefidin A)
# IndivCpd <- 8 # CDDP (Cisplatin)
# IndivCpd <- 10 # CYA (Cyclosporin A)
# IndivCpd <- 20 # MYT (Mitomycin C)
# IndivCpd <- 29 # THAP (Thapsigargin)
# IndivCpd <- 31 # TUN (Tunicamycin)
IndivCpd <- NULL # or NULL

RunIdx <- 7 # Index of optimisation run (to be include in file names)

# =========================================== #

# install.packages('devtools')
library(devtools)
# install("multiCellNOpt-master/")
library(multiCellNOpt)
# source("https://bioconductor.org/biocLite.R")
# biocLite("MEIGOR")
library(MEIGOR)

library(doParallel)
argsJob=commandArgs(trailingOnly = TRUE)
repIndex1 <- as.numeric(argsJob[1])

# install_github("saezlab/CellNOptR/packages/CellNOptR")

All_SR_Models_Multi  <- list()
All_SR_Results_Multi <- list()

# Compound_List = c("CDDP","ETO","MYT","PCM") 
Compound_List = t(read.table(file = "Compound_List.csv",header = F,sep = "\n",stringsAsFactors = F))
if (!is.null(IndivCpd)) {
  Compound_List <- Compound_List[IndivCpd]
}

for (Round_ID in 1:NrRoundsNeed) {
    
  All_SR_Results <- list()
  All_SR_Models <- list()

  for (counter in 1:length(Compound_List)) {
  # counter <- repIndex1
  # counter <- 1
  
    print(paste0("Running optimisation for compound: ",counter,"/",length(Compound_List), " - ",Compound_List[counter]))
    
    # ModelFile = "SR_LogicODE_NoCrosstalk.sif"
    # ModelFile = "SR_LogicODE_NoCrosstalk_minimal.sif"
    # ModelFile = "SR_LogicODE_NoCrosstalk_minimal_withAbsorbNode.sif"
    # ModelFile = "SR_LogicODE_NoCrosstalk_minimal_withCT.sif"
    # ModelFile = "SR_LogicODE_NoCrosstalk_minimal_withCT_withApoptosisAND.sif"
    # ModelFile = "SR_LogicODE_NoCrosstalk_minimal_withCT_withApoptosisOR.sif"
    
    # ModelFile = "SR_LogicODE_withCrosstalk_minimal_ApoptosisAND.sif"
    # ModelFile = "SR_LogicODE_withCrosstalk_minimal_ApoptosisOR.sif"
    ModelFile = "SR_LogicODE_withCrosstalk_minimal_ApoptosisANDOR.sif"
    
    # ModelFile = "SR_LogicODE_NoCrosstalk_minimal_withAbsorbNode_withCT.sif"
    # ModelFile = "SR_LogicODE_NoCrosstalk_LukasCommentNr1.sif"
    # ModelFile = "SR_LogicODE_NoCrosstalk_withAbsorbNode.sif"
    # ModelFile = "SR_LogicODE_NoCrosstalk_withAbsorbNode_LukasCommentNr1.sif"
    # DataFile = paste0("20180712_leiden_AllSR/MIDAS_output/SR_MIDAS_",Compound_List[counter],"_AllSR_Corrected_NormPerRepl_NoNormCellCt_v4.csv",sep="")
    # DataFile = paste0("20180712_leiden_AllSR/MIDAS_output/SR_MIDAS_",Compound_List[counter],"_AllSR_Corrected_NormPerRepl_NoNormCellCt_v4_NormLog10Doses.csv",sep="")
    DataFile = paste0("20180712_leiden_AllSR/MIDAS_output/SR_MIDAS_",Compound_List[counter],"_AllSR_Corrected_NormPerRepl_NoNormCellCt_v5_NormLog10Doses_Apoptosis.csv",sep="")

    All_SR_Models[[counter]] = logicODEModel$new(SIFfile = ModelFile, exps = DataFile,initCond = 0)
    All_SR_Models[[counter]]$addControlExperiment(basalCue = 0,basalSignal = 0)
    # All_SR_Models[[counter]]$addControlExperiment(basalCue = 0,basalSignal = 0,forNonMeasuredNodes = TRUE)
    # All_SR_Models[[counter]]$preprocessing(inhibANDExpansion = TRUE)
    
    # pdf(paste0("PlotData_",ModelFile,"_",Compound_List[counter],"_InRoundRep",NrInRoundRep,"_Time",TimeRoundRep,".pdf"))
    # pdf(paste0("PlotData_",DataFile,"_",Compound_List[counter],"_InRoundRep",NrInRoundRep,"_Time",TimeRoundRep,".pdf"))
    # pdf(paste0("PlotData_",Compound_List[counter],"_PlusApoptosis.pdf"))
    # All_SR_Models[[counter]]$plotData()
    # dev.off()
    # All_SR_Models[[counter]]$plotData()

    pdf(paste0("PlotModel_",ModelFile,".pdf"))
    All_SR_Models[[counter]]$plotModel()
    dev.off()
    # All_SR_Models[[counter]]$plotModel()
    
    All_SR_Models[[counter]]$initODE_parameters(opt_n = FALSE)
    # All_SR_Models[[counter]]$initODE_parameters(opt_n = TRUE)
    # All_SR_Models[[counter]]$ode_parameters
    
    # Assign parameter bounds
    # LB_n <- 1
    LB_n <- 3
    LB_k <- 0
    LB_tau <- 0
    # UB_n <- 10
    UB_n <- 3
    UB_k <- 10
    # UB_k <- 1
    UB_tau <- 10
    # UB_tau <- 1

    All_SR_Models[[counter]]$ode_parameters$LB[All_SR_Models[[counter]]$ode_parameters$index_n]   <- LB_n
    All_SR_Models[[counter]]$ode_parameters$LB[All_SR_Models[[counter]]$ode_parameters$index_k]   <- LB_k
    All_SR_Models[[counter]]$ode_parameters$LB[All_SR_Models[[counter]]$ode_parameters$index_tau] <- LB_tau
    
    All_SR_Models[[counter]]$ode_parameters$UB[All_SR_Models[[counter]]$ode_parameters$index_n]   <- UB_n
    All_SR_Models[[counter]]$ode_parameters$UB[All_SR_Models[[counter]]$ode_parameters$index_k]   <- UB_k
    All_SR_Models[[counter]]$ode_parameters$UB[All_SR_Models[[counter]]$ode_parameters$index_tau] <- UB_tau
    
    All_SR_Models[[counter]]$transfer_function <- 4
    
    All_SR_Models[[counter]]$objectiveFunction <- All_SR_Models[[counter]]$getDefaultLS(SSpenalty_fac = 10,SScontrolPenalty_fac = 10)
    # All_SR_Models[[counter]]$objectiveFunction <- M$getDefaultLS(SSpenalty_fac = 0,SScontrolPenalty_fac = 0)
    All_SR_Models[[counter]]$objectiveFunction <- All_SR_Models[[counter]]$getDefaultLS(lambda_tau = 0.001, lambda_k = 0.001)
    
    for (RoundRep in 1:NrInRoundRep) {
      print("======================================")
      print("--------------------------------------")
      print(paste("Optimising Exp:",toString(counter),"- Round:",  toString(Round_ID),"- Rep:", toString(RoundRep)))
      print("--------------------------------------")
      print("======================================")
      All_SR_Models[[counter]]$fit(maxTime = TimeRoundRep,nRun = 1, nCores = 1)
    }
    
    if(length(dev.list())>0){dev.off()}
    
    pdf(paste0("PlotFit_",ModelFile,"_",Compound_List[counter],"_InRoundRep",NrInRoundRep,"_Time",TimeRoundRep,"_OnlyMeas","Run",RunIdx,".pdf"))
    All_SR_Models[[counter]]$plotFit(measuredNodesOnly = TRUE)
    dev.off()
    pdf(paste0("PlotFit_",ModelFile,"_",Compound_List[counter],"_InRoundRep",NrInRoundRep,"_Time",TimeRoundRep,"_AllNodes","Run",RunIdx,".pdf"))
    All_SR_Models[[counter]]$plotFit(measuredNodesOnly = FALSE)
    dev.off()
    All_SR_Results[[counter]] <- All_SR_Models[[counter]]$ode_parameters
  
  }
  
  All_SR_Models_Multi[[Round_ID]]  <- All_SR_Models
  All_SR_Results_Multi[[Round_ID]] <- All_SR_Results
  
  save(All_SR_Models_Multi,file=paste0("Results_SingleModels_",ModelFile,"_InRoundRep",NrInRoundRep,"_Time",TimeRoundRep,"_Run",RunIdx,"_",Compound_List[counter],".RData"))
  save(All_SR_Results_Multi,file=paste0("Results_SingleResults_",ModelFile,"_InRoundRep",NrInRoundRep,"_Time",TimeRoundRep,"_Run",RunIdx,"_",Compound_List[counter],".RData"))

}       

# for (counter in 1:length(All_SR_Models)) {
#   All_SR_Models[[counter]]$ode_parameters <- All_SR_Results[[counter]]
# }

# For all models from all compounds

# if (!is.null(IndivCpd)) {
#   save(All_SR_Models_Multi,file=paste0("Results_AllModelsMulti_",ModelFile,"_InRoundRep",NrInRoundRep,"_Time",TimeRoundRep,"_Run",RunIdx,"_",Compound_List,".RData"))
#   save(All_SR_Results_Multi,file=paste0("Results_AllResultsMulti_",ModelFile,"_InRoundRep",NrInRoundRep,"_Time",TimeRoundRep,"_Run",RunIdx,"_",Compound_List,".RData"))
# } else {
#   save(All_SR_Models_Multi,file=paste0("Results_AllModelsMulti_",ModelFile,"_InRoundRep",NrInRoundRep,"_Time",TimeRoundRep,"_Run",RunIdx,".RData"))
#   save(All_SR_Results_Multi,file=paste0("Results_AllResultsMulti_",ModelFile,"_InRoundRep",NrInRoundRep,"_Time",TimeRoundRep,"_Run",RunIdx,".RData"))
# }



# # === Multi experiment plotting === #
# 
# if (MultiModPlot) {
# 
#   lm = multiLogicODEModel$new(SIFfile="SR_Network_HandCurated_CellNOpt_v4.txt",
#                               exps=c("MIDAS_Datasets/SR_MIDAS_CDDP_TimeCourse.csv",
#                                      "MIDAS_Datasets/SR_MIDAS_ETO_TimeCourse.csv",
#                                      "MIDAS_Datasets/SR_MIDAS_MYT_TimeCourse.csv",
#                                      "MIDAS_Datasets/SR_MIDAS_PCM_TimeCourse.csv"))
#   
#   Compound_List = t(read.table(file = "Compound_List.csv",header = F,sep = "\n",stringsAsFactors = F))
#   
#   MIDAS_filenames <- NULL
#   for (counter in 1:length(Compound_List)) {
#     # MIDAS_filenames <- c(MIDAS_filenames,paste0("MIDAS_Datasets/SR_MIDAS_",Compound_List[counter],"_TimeCourse.csv"))
#     MIDAS_filenames <- c(MIDAS_filenames,paste0("MIDAS_Datasets/SR_MIDAS_",Compound_List[counter],"_TimeCourse_withApoptosis.csv"))
#   }
#   
#   # lm = multiLogicODEModel$new(SIFfile="SR_Network_HandCurated_CellNOpt_v4.txt",
#   lm = multiLogicODEModel$new(SIFfile="SR_Network_HandCurated_WithApoptosis.txt",
#                                                           exps=MIDAS_filenames)
# 
#   lm$addControlExperiment(basalCue = 0,basalSignal = 0)
#   
#   lm$ode_parameters$set_odeParameters(All_SR_Results_Multi[[1]])  
# 
#   lm$plotFit()
#   
#   # pdf("PlotFit_SR_MultiCellNOpt.pdf")
#   pdf("PlotFit_SR_MultiCellNOpt_withApoptosis.pdf")
#   # lm$plotFit(observablesOnly = FALSE)
#   lm$plotFit(observablesOnly = TRUE)
#   dev.off()
# 
# }

# # Plotting Data
# Compound_List = t(read.table(file = "Compound_List.csv",header = F,sep = "\n",stringsAsFactors = F))
# 
# for (counter in 1:length(Compound_List)) {
#   
#   print(paste0("Running optimisation for compound: ",counter,"/",length(Compound_List), " - ",Compound_List[counter]))
#   
#   ModelFile = "SR_LogicODE_NoCrosstalk.sif"
#   DataFile = paste0("20180712_leiden_AllSR/MIDAS_output/SR_MIDAS_",Compound_List[counter],"_AllSR_Corrected_NormPerRepl_NoNormCellCt_v4.csv",sep="")
#   
#   All_SR_Models[[counter]] = logicODEModel$new(SIFfile = ModelFile, exps = DataFile,initCond = 0)
#   All_SR_Models[[counter]]$addControlExperiment(basalCue = 0,basalSignal = 0)
#   # All_SR_Models[[counter]]$addControlExperiment(basalCue = 0,basalSignal = 0,forNonMeasuredNodes = TRUE)
#   # All_SR_Models[[counter]]$preprocessing(inhibANDExpansion = TRUE)
#   
#   pdf(paste0("PlotData_",Compound_List[counter],"_PreCheck.pdf"))
#   # pdf(paste0("PlotData_",DataFile,"_",Compound_List[counter],"_InRoundRep",NrInRoundRep,"_Time",TimeRoundRep,".pdf"))
#   All_SR_Models[[counter]]$plotData()
#   dev.off()
#   # All_SR_Models[[counter]]$plotData()
# }

# --- End of the script --- #
