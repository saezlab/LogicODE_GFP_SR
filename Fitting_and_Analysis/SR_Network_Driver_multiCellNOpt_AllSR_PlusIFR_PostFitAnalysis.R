# Driver script for multiCellNopt package - SR dataset

# ====================================== #
# == CLEAN WORKSPACE & LOAD LIBRARIES == #
# ====================================== #

rm(list=ls()) # clear all variables
cat("\014") # clear screen
if (length(dev.list()>0)) {dev.off()} # clear figures

library(ggplot2)
library(reshape2)

# ====================================== #
# ========= PRE-LOAD DATASETS ========== #
# ====================================== #


# Load dosing information (AllDoses)
load("20180712_leiden_AllSR/StressResponse_Variables_AllSR_Corrected_withSD_CellCt_v4.RData")
rm(list=setdiff(ls(),"AllDoses"))

# DMSO and HEP have only 2 doses (to be investigated)
AllDoses[[14]] <- AllDoses[[14]][2]
AllDoses[[17]] <- AllDoses[[17]][2]

# Fix Swapped DMSO/DMEM doses (Idx 13 and 14)
TempDMSOdoses <- AllDoses[[14]]; TempDMEMDoses <- AllDoses[[13]]
AllDoses[[13]] <- TempDMSOdoses; AllDoses[[14]] <- TempDMEMDoses
# Replace empty dose with very low numbers
# AllDoses[[13]][1] <- 1e-10;
# AllDoses[[17]][1] <- 1e-10;

# Load the list of compounds
Compound_List = t(read.table(file = "Compound_List.csv",header = F,sep = "\n",stringsAsFactors = F))

# Load information on BestRun
BestFitTable <- as.data.frame(read.table("AllFitCost_Table_MCP_101218.tsv",sep = "\t",header=T,stringsAsFactors = F))
BestRun <- as.numeric(substr(BestFitTable$BestRun,0,1))

# ====================================== #
# ======== LONG TERM SIMULATION ======== #
# ====================================== #

# Long term simulation (72hr [3 days] -> e.g. 240hr [10 days])
SimulateTime <- 240

for (CounterCpd in 1:length(Compound_List)) {
  
  print(paste0("Long-term simulation for ",CounterCpd,"/",length(Compound_List)," : ",Compound_List[CounterCpd]))
  
  if (CounterCpd==13) {CounterCpd=14} else if (CounterCpd==14) {CounterCpd=13} # Fix wrong indices of DMSO and DMEM in result filess
  load(paste0("Results_GFP_101218/Results_SingleModels_SR_LogicODE_withCrosstalk_minimal_ApoptosisANDOR.sif_InRoundRep3_Time3600_Run7_",Compound_List[CounterCpd],".RData"))
  if (CounterCpd==13) {CounterCpd=14} else if (CounterCpd==14) {CounterCpd=13} # Return the indices
  BestModelIdx <- BestRun[CounterCpd]
  ModelCpd <- All_SR_Models_Multi[[BestModelIdx]][[CounterCpd]]
  SimCpd <- ModelCpd$simulate(timeSignals = 1:SimulateTime)
  SimCpdArray <- array(data = NA,dim = c(nrow(SimCpd$signals[[1]]),length(SimCpd$signals),ncol(SimCpd$signals[[1]])),dimnames = list(c('d0',paste0('d',AllDoses[[CounterCpd]])),SimCpd$timepoints,colnames(SimCpd$signals[[1]])))
  for (counter in 1:length(SimCpd$signals)) {
    for (counter2 in 1:ncol(SimCpd$signals[[1]])) {
      SimCpdArray[,counter,counter2] <- SimCpd$signals[[counter]][,counter2]
    }
  }
  
  pdf(paste0("LongSim_",Compound_List[CounterCpd],"_",SimulateTime,"hr.pdf"))
  
  for (counter2 in 1:dim(SimCpdArray)[3]) {
    SimCpdGFP <- SimCpdArray[,,counter2]
    SimCpdMeltedGFP <- melt(SimCpdGFP); colnames(SimCpdMeltedGFP) <- c('Doses','Time','Value')
    
    print(ggplot() +
      geom_line(data=SimCpdMeltedGFP,aes(x=Time,y=Value ,color=Doses),linetype="solid",size=1,alpha=0.5)+
      ggtitle(paste0("Long-term simulation for ",colnames(SimCpd$signals[[1]])[counter2],": ",Compound_List[CounterCpd]))+
      scale_y_continuous(limits=c(0,1))+
      geom_vline(xintercept=72,linetype="dotted",color="gray50") +
      theme_light()) +
      xlab('Time (hr)')
  }
  dev.off()
}

# ====================================== #
# === PARAMETER SENSITIVITY ANALYSIS === #
# ====================================== #

# Initial a variable to extract the Top3 parameters
MaxChange <- as.data.frame(matrix(NA,6,length(Compound_List)))
colnames(MaxChange) <- Compound_List
rownames(MaxChange) <- c("Param1st","Change1st","Param2nd","Change2nd","Param3rd","Change3rd")

for (CounterCpd in 1:length(Compound_List)) {
  # for (CounterCpd in 13:14) {
  
  print(paste0("Mapping SA results ",CounterCpd,"/",length(Compound_List)," : ",Compound_List[CounterCpd]))
  # CounterCpd <- 7
  
  # Load fitted model  + Correct mismatched of compound mapping for DMSO and DMEM
  # if (CounterCpd==13) {
  #   load(paste0("Results_GFP_101218/Results_SingleModels_SR_LogicODE_withCrosstalk_minimal_ApoptosisANDOR.sif_InRoundRep3_Time3600_Run7_",Compound_List[14],".RData"))
  #   BestModelIdx <- BestRun[14]
  #   ModelOrig <- All_SR_Models_Multi[[BestModelIdx]][[14]]
  #   CounterCpd_Dose <- 14
  # } else if (CounterCpd==14) {
  #   load(paste0("Results_GFP_101218/Results_SingleModels_SR_LogicODE_withCrosstalk_minimal_ApoptosisANDOR.sif_InRoundRep3_Time3600_Run7_",Compound_List[13],".RData"))
  #   BestModelIdx <- BestRun[13]
  #   ModelOrig <- All_SR_Models_Multi[[BestModelIdx]][[13]]
  #   CounterCpd_Dose <- 13
  # } else {
  
  if (CounterCpd==13) {CounterCpd=14} else if (CounterCpd==14) {CounterCpd=13} # Fix wrong indices of DMSO and DMEM in result filess
  load(paste0("Results_GFP_101218/Results_SingleModels_SR_LogicODE_withCrosstalk_minimal_ApoptosisANDOR.sif_InRoundRep3_Time3600_Run7_",Compound_List[CounterCpd],".RData"))
  if (CounterCpd==13) {CounterCpd=14} else if (CounterCpd==14) {CounterCpd=13} # Fix wrong indices of DMSO and DMEM in result filess
  BestModelIdx <- BestRun[CounterCpd]
  ModelOrig <- All_SR_Models_Multi[[BestModelIdx]][[CounterCpd]]
  CounterCpd_Dose <- CounterCpd
  # }
  ParamOrig <- ModelOrig$ode_parameters$parValues
  
  # Simulate original value
  SimOrig <- ModelOrig$simulate() 
  # SimOrigArray <- array(data = NA,dim = c(nrow(SimOrig$signals[[1]]),length(SimOrig$signals),ncol(SimOrig$signals[[1]])),dimnames = list(paste0('d',1:nrow(SimOrig$signals[[1]])),SimOrig$timepoints,colnames(SimOrig$signals[[1]])))
  SimOrigArray <- array(data = NA,dim = c(nrow(SimOrig$signals[[1]]),length(SimOrig$signals),ncol(SimOrig$signals[[1]])),dimnames = list(c('d0',paste0('d',AllDoses[[CounterCpd_Dose]])),SimOrig$timepoints,colnames(SimOrig$signals[[1]])))
  for (counter in 1:length(SimOrig$signals)) {
    for (counter2 in 1:ncol(SimOrig$signals[[1]])) {
      SimOrigArray[,counter,counter2] <- SimOrig$signals[[counter]][,counter2]
    }
  }
  SimOrigApop <- SimOrigArray[,,which(colnames(SimOrig$signals[[1]])=="Apoptosis")]
  
  # Local sensitivity analysis at 1% or 10%
  # PertSensiv <- 0.01 # 1%
  PertSensiv <- 0.1 # 10%
  
  # Extract k and tau parameters
  Idx_k_tau <- c(ModelOrig$ode_parameters$index_k,ModelOrig$ode_parameters$index_tau)
  AvgAbsChange <- rep(NA,length(Idx_k_tau))
  
  pdf(paste0("SA_TimeCouseApoptosis_",Compound_List[CounterCpd],"_",PertSensiv*100,"percent.pdf"))
  
  for (counter_PertParam in 1:length(Idx_k_tau)) {
    
    SelectedParamIdx <- Idx_k_tau[counter_PertParam]
    # Idx_k[37] == "DDIT3_k_Apoptosis"
    
    # Increase parameter value
    ModelPertUp <- ModelOrig
    ModelPertUp$ode_parameters$parValues <- ParamOrig
    
    ModelPertUp$ode_parameters$parValues[SelectedParamIdx] <- ModelPertUp$ode_parameters$parValues[SelectedParamIdx] + 
      PertSensiv*ModelPertUp$ode_parameters$parValues[SelectedParamIdx]
    
    # print(ModelPertUp$ode_parameters$parValues[SelectedParamIdx])
    
    SimPertUp <- ModelPertUp$simulate()
    # SimPertUpArray <- array(data = NA,dim = c(nrow(SimPertUp$signals[[1]]),length(SimPertUp$signals),ncol(SimPertUp$signals[[1]])),dimnames = list(paste0('d',1:nrow(SimPertUp$signals[[1]])),SimPertUp$timepoints,colnames(SimPertUp$signals[[1]])))
    SimPertUpArray <- array(data = NA,dim = c(nrow(SimPertUp$signals[[1]]),length(SimPertUp$signals),ncol(SimPertUp$signals[[1]])),dimnames = list(c('d0',paste0('d',AllDoses[[CounterCpd_Dose]])),SimPertUp$timepoints,colnames(SimPertUp$signals[[1]])))
    for (counter in 1:length(SimPertUp$signals)) {
      for (counter2 in 1:ncol(SimPertUp$signals[[1]])) {
        SimPertUpArray[,counter,counter2] <- SimPertUp$signals[[counter]][,counter2]
      }
    }
    SimPertUpApop <- SimPertUpArray[,,which(colnames(SimPertUp$signals[[1]])=="Apoptosis")]
    
    # Decrease parameter value
    ModelPertDn <- ModelOrig
    
    ModelPertDn$ode_parameters$parValues <- ParamOrig
    
    ModelPertDn$ode_parameters$parValues[SelectedParamIdx] <- ModelPertDn$ode_parameters$parValues[SelectedParamIdx] - 
      PertSensiv*ModelPertDn$ode_parameters$parValues[SelectedParamIdx]
    
    # print(ModelPertDn$ode_parameters$parValues[SelectedParamIdx])
    
    SimPertDn <- ModelPertDn$simulate()
    # SimPertDnArray <- array(data = NA,dim = c(nrow(SimPertDn$signals[[1]]),length(SimPertDn$signals),ncol(SimPertDn$signals[[1]])),dimnames = list(paste0('d',1:nrow(SimPertDn$signals[[1]])),SimPertDn$timepoints,colnames(SimPertDn$signals[[1]])))
    SimPertDnArray <- array(data = NA,dim = c(nrow(SimPertDn$signals[[1]]),length(SimPertDn$signals),ncol(SimPertDn$signals[[1]])),dimnames = list(c('d0',paste0('d',AllDoses[[CounterCpd_Dose]])),SimPertDn$timepoints,colnames(SimPertDn$signals[[1]])))
    for (counter in 1:length(SimPertDn$signals)) {
      for (counter2 in 1:ncol(SimPertDn$signals[[1]])) {
        SimPertDnArray[,counter,counter2] <- SimPertDn$signals[[counter]][,counter2]
      }
    }
    SimPertDnApop <- SimPertDnArray[,,which(colnames(SimPertDn$signals[[1]])=="Apoptosis")]
    
    # Down-sampling of time?
    # SimApopTimeDS <- SimApop[,seq(1,ncol(SimApop),length=50)]
    
    # Calculate avarage absolute change per dose
    AvgAbsChange[counter_PertParam] <- mean(rowSums(abs(SimOrigApop-SimPertUpApop) + abs(SimOrigApop-SimPertDnApop)))
    
    # Start plotting
    library(reshape2)
    # SimPlot <- array(data = NA,dim = c(nrow(SimOrigApop),ncol(SimOrigApop),3),dimnames = list(rownames(SimOrigApop),colnames(SimOrigApop),c("Orig","Up","Dn")))
    # SimPlot[,,1] <- SimOrigApop; SimPlot[,,2] <- SimPertUpApop; SimPlot[,,3] <- SimPertDnApop
    # SimPlotMelted <- melt(SimPlot); colnames(SimPlotMelted) <- c('Doses','Time','Type','Value')
    
    SimOrigMelted <- melt(SimOrigApop); colnames(SimOrigMelted) <- c('Doses','Time','Value')
    SimPertUpMelted <- melt(SimPertUpApop); colnames(SimPertUpMelted) <- c('Doses','Time','Value')
    SimPertDnMelted <- melt(SimPertDnApop); colnames(SimPertDnMelted) <- c('Doses','Time','Value')
    
    library(ggplot2)
    print(ggplot() +
            # geom_point(data=SimOrigMelted,aes(x=Time,y=Value ,color=Doses),shape=1,size=0.1,alpha=0.5)+
            geom_line(data=SimOrigMelted,aes(x=Time,y=Value ,color=Doses),linetype="solid",size=1,alpha=0.5)+
            geom_point(data=SimPertUpMelted,aes(x=Time,y=Value ,color=Doses),shape=2,size=0.5,alpha=0.5)+
            geom_point(data=SimPertDnMelted,aes(x=Time,y=Value ,color=Doses),shape=6,size=0.5,alpha=0.5)+
            # geom_line(data=SimPertUpMelted,aes(x=Time,y=Value ,color=Doses),linetype="solid",size=1,alpha=0.5)+
            # geom_line(data=SimPertDnMelted,aes(x=Time,y=Value ,color=Doses),linetype="solid",size=1,alpha=0.5)+
            
            # geom_ribbon(aes(ymax = SimPertUpApop,
            #                  ymin = SimPertDnApop,
            #                  fill = Doses),
            #              alpha = 0.1) +
            # geom_ribbon(data = aes(ymax = meanInt + sdInt,
            #                        ymin = meanInt - sdInt,
            #                        fill = dose_uM),
            #             alpha = 0.1) +
          ggtitle(paste0("SimApoptosis: ",Compound_List[CounterCpd]," ; ",ModelOrig$ode_parameters$parNames[SelectedParamIdx], "=",round(ParamOrig[SelectedParamIdx],digit=6)," +/- ",PertSensiv*100,"%"))+
            scale_y_continuous(limits=c(0,1))+
            theme_light()
          
          # theme_classic()
    )
    rm(ModelPertDn);rm(ModelPertUp);rm(SimPertUp);rm(SimPertDn)
  }
  # print(g)
  dev.off()
  
  names(AvgAbsChange) <- names(ParamOrig)[Idx_k_tau]
  AvgAbsChange <- as.data.frame(sort(AvgAbsChange,decreasing = T))
  AvgAbsChange <- cbind(AvgAbsChange,rownames(AvgAbsChange))
  colnames(AvgAbsChange) <- c("Change","Param")
  AvgAbsChange <- AvgAbsChange[order(AvgAbsChange[,1],decreasing = T),]
  AvgAbsChange$Param <- factor(AvgAbsChange$Param,levels=AvgAbsChange$Param)
  
  MaxChange[1,CounterCpd] <- as.character(AvgAbsChange[1,2])
  MaxChange[2,CounterCpd] <- AvgAbsChange[1,1]
  MaxChange[3,CounterCpd] <- as.character(AvgAbsChange[2,2])
  MaxChange[4,CounterCpd] <- AvgAbsChange[2,1]
  MaxChange[5,CounterCpd] <- as.character(AvgAbsChange[3,2])
  MaxChange[6,CounterCpd] <- AvgAbsChange[3,1]
  
  pdf(paste0("SA_AvgAccumApoptosis_",Compound_List[CounterCpd],"_",PertSensiv*100,"percent.pdf"))
  
  print(ggplot(data=AvgAbsChange,aes(x=Param,y=Change)) + 
          geom_bar(stat="identity",aes(fill=Change)) +
          # scale_fill_brewer(palette="Spectral") + 
          # theme(axis.text.x=element_text(angle=90,hjust=1)) +
          labs(title = paste0("Accumulated %changes of Apoptosis per dose: ",Compound_List[CounterCpd]), subtitle = paste0("Perturbed parameters ",PertSensiv*100,"%")) +
          # xlab("time (days)") +
          # ylab("ng/ml") +
          theme_light()+
          theme(plot.title = element_text(size=15)) +
          theme(plot.subtitle = element_text(size=12)) +
          theme(axis.text.x = element_text(colour="grey20",size=7.5,angle=90,hjust=1,vjust=0,face="plain"),
                axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
                axis.title.x = element_blank(),
                axis.title.y = element_blank()))
  dev.off()
  
}

save(MaxChange,file = "MaxChangeApoptosis.RData")

# Add compound names at the end of the top parameter names
for (counter in 1:ncol(MaxChange)) {
  MaxChange[c(1,3,5),counter] <- paste0(MaxChange[c(1,3,5),counter],"_",colnames(MaxChange)[counter])
}

MaxChangePlot <- as.data.frame(t(cbind(MaxChange[c(1,2),],MaxChange[c(3,4),],MaxChange[c(5,6),])),stringsAsFactors = F)
MaxChangePlot[,2] <- as.numeric(MaxChangePlot[,2])
MaxChangePlot <- MaxChangePlot[order(MaxChangePlot[,2],decreasing = T),]
MaxChangePlot <- MaxChangePlot[-which(MaxChangePlot[,2]==0),] # Remove all entries with zero value
MaxChangePlot$Param1st <- factor(MaxChangePlot$Param1st,levels=MaxChangePlot$Param1st)
colnames(MaxChangePlot) <- c("Param","Change")

pdf(paste0("SA_MaxChangeApop_Top3_",PertSensiv*100,"percent.pdf"))

print(ggplot(data=MaxChangePlot,aes(x=Param,y=Change)) + 
        geom_bar(stat="identity",aes(fill=Change)) +
        # scale_fill_brewer(palette="Spectral") + 
        # theme(axis.text.x=element_text(angle=90,hjust=1)) +
        labs(title = paste0("Top 3 param for accumulated %changes of Apoptosis per dose"), subtitle = paste0("Perturbed parameters ",PertSensiv*100,"%")) +
        # xlab("time (days)") +
        # ylab("ng/ml") +
        theme_light()+
        theme(plot.title = element_text(size=15)) +
        theme(plot.subtitle = element_text(size=12)) +
        theme(axis.text.x = element_text(colour="grey20",size=6,angle=90,hjust=1,vjust=0,face="plain"),
              axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
              axis.title.x = element_blank(),
              axis.title.y = element_blank()))
dev.off()


# ====================================== #
# ========= KNOCK-OUT ANALYSIS ========= #
# ====================================== #

# Focus on two most important reactions: DDIT3_k_Apoptosis, TP53_k_Apoptosis & nodes tau_DDIT3, tau_TP53
# KO_IntAct <- c('DDIT3_k_Apoptosis','TP53_k_Apoptosis')
# KO_Node <- c('tau_DDIT3','tau_TP53')
# KO_All <- c(KO_IntAct,KO_Node)

# Initial a variable to extract the Top3 parameters
MaxKOChange <- as.data.frame(matrix(NA,6,length(Compound_List)))
colnames(MaxKOChange) <- Compound_List
rownames(MaxKOChange) <- c("Param1st","Change1st","Param2nd","Change2nd","Param3rd","Change3rd")

for (CounterCpd in 1:length(Compound_List)) {
# for (CounterCpd in 13:length(Compound_List)) {
    
  print(paste0("Mapping KO results ",CounterCpd,"/",length(Compound_List)," : ",Compound_List[CounterCpd]))

  if (CounterCpd==13) {CounterCpd=14} else if (CounterCpd==14) {CounterCpd=13} # Fix wrong indices of DMSO and DMEM in result filess
  load(paste0("Results_GFP_101218/Results_SingleModels_SR_LogicODE_withCrosstalk_minimal_ApoptosisANDOR.sif_InRoundRep3_Time3600_Run7_",Compound_List[CounterCpd],".RData"))
  if (CounterCpd==13) {CounterCpd=14} else if (CounterCpd==14) {CounterCpd=13} # Fix wrong indices of DMSO and DMEM in result filess
  BestModelIdx <- BestRun[CounterCpd]
  ModelOrig <- All_SR_Models_Multi[[BestModelIdx]][[CounterCpd]]
  CounterCpd_Dose <- CounterCpd
  ParamOrig <- ModelOrig$ode_parameters$parValues
  
  # Simulate original value
  SimOrig <- ModelOrig$simulate() 
  # SimOrigArray <- array(data = NA,dim = c(nrow(SimOrig$signals[[1]]),length(SimOrig$signals),ncol(SimOrig$signals[[1]])),dimnames = list(paste0('d',1:nrow(SimOrig$signals[[1]])),SimOrig$timepoints,colnames(SimOrig$signals[[1]])))
  SimOrigArray <- array(data = NA,dim = c(nrow(SimOrig$signals[[1]]),length(SimOrig$signals),ncol(SimOrig$signals[[1]])),dimnames = list(c('d0',paste0('d',AllDoses[[CounterCpd_Dose]])),SimOrig$timepoints,colnames(SimOrig$signals[[1]])))
  for (counter in 1:length(SimOrig$signals)) {
    for (counter2 in 1:ncol(SimOrig$signals[[1]])) {
      SimOrigArray[,counter,counter2] <- SimOrig$signals[[counter]][,counter2]
    }
  }
  SimOrigApop <- SimOrigArray[,,which(colnames(SimOrig$signals[[1]])=="Apoptosis")]
  
  # In silico knock-out experiment
  
  IdxKO <- c(ModelOrig$ode_parameters$index_k,ModelOrig$ode_parameters$index_tau)
  
  # IdxKO <- NULL
  # for (counter in 1:length(KO_All)) {
  #   IdxKO <- c(IdxKO,which(KO_All[counter]==ModelOrig$ode_parameters$parNames))
  # }
  
  
  AvgAbsChange <- rep(NA,length(IdxKO))
  
  pdf(paste0("KO_TimeCouseApoptosis_",Compound_List[CounterCpd],".pdf"))
  
  for (counter_KO in 1:length(IdxKO)) {
    
    SelectedParamIdx <- IdxKO[counter_KO]
    # Idx_k[37] == "DDIT3_k_Apoptosis"
    
    # Increase parameter value
    ModelKO <- ModelOrig
    ModelKO$ode_parameters$parValues <- ParamOrig
    
    ModelKO$ode_parameters$parValues[SelectedParamIdx] <- 1e-10
    # ModelKO$ode_parameters$parValues[SelectedParamIdx] <- 0
    
    # print(ModelKO$ode_parameters$parValues[SelectedParamIdx])
    
    SimKO <- ModelKO$simulate()
    # SimKOArray <- array(data = NA,dim = c(nrow(SimKO$signals[[1]]),length(SimKO$signals),ncol(SimKO$signals[[1]])),dimnames = list(paste0('d',1:nrow(SimKO$signals[[1]])),SimKO$timepoints,colnames(SimKO$signals[[1]])))
    SimKOArray <- array(data = NA,dim = c(nrow(SimKO$signals[[1]]),length(SimKO$signals),ncol(SimKO$signals[[1]])),dimnames = list(c('d0',paste0('d',AllDoses[[CounterCpd_Dose]])),SimKO$timepoints,colnames(SimKO$signals[[1]])))
    for (counter in 1:length(SimKO$signals)) {
      for (counter2 in 1:ncol(SimKO$signals[[1]])) {
        SimKOArray[,counter,counter2] <- SimKO$signals[[counter]][,counter2]
      }
    }
    SimKOApop <- SimKOArray[,,which(colnames(SimKO$signals[[1]])=="Apoptosis")]
    
    # Calculate avarage absolute change per dose
    AvgAbsChange[counter_KO] <- mean(rowSums(abs(SimOrigApop-SimKOApop)))
    
    # Start plotting
    library(reshape2)
    # SimPlot <- array(data = NA,dim = c(nrow(SimOrigApop),ncol(SimOrigApop),3),dimnames = list(rownames(SimOrigApop),colnames(SimOrigApop),c("Orig","Up","Dn")))
    # SimPlot[,,1] <- SimOrigApop; SimPlot[,,2] <- SimPertUpApop; SimPlot[,,3] <- SimPertDnApop
    # SimPlotMelted <- melt(SimPlot); colnames(SimPlotMelted) <- c('Doses','Time','Type','Value')
    
    SimOrigMelted <- melt(SimOrigApop); colnames(SimOrigMelted) <- c('Doses','Time','Value')
    SimKOMelted <- melt(SimKOApop); colnames(SimKOMelted) <- c('Doses','Time','Value')

    library(ggplot2)
    print(ggplot() +
            # geom_point(data=SimOrigMelted,aes(x=Time,y=Value ,color=Doses),shape=1,size=0.1,alpha=0.5)+
            geom_line(data=SimOrigMelted,aes(x=Time,y=Value ,color=Doses),linetype="solid",size=1,alpha=0.5)+
            geom_point(data=SimKOMelted,aes(x=Time,y=Value ,color=Doses),shape=6,size=0.5,alpha=0.5)+
            # geom_point(data=SimPertDnMelted,aes(x=Time,y=Value ,color=Doses),shape=6,size=0.5,alpha=0.5)+
            # geom_line(data=SimPertUpMelted,aes(x=Time,y=Value ,color=Doses),linetype="solid",size=1,alpha=0.5)+
            # geom_line(data=SimPertDnMelted,aes(x=Time,y=Value ,color=Doses),linetype="solid",size=1,alpha=0.5)+
            
            # geom_ribbon(aes(ymax = SimPertUpApop,
            #                  ymin = SimPertDnApop,
            #                  fill = Doses),
            #              alpha = 0.1) +
            # geom_ribbon(data = aes(ymax = meanInt + sdInt,
            #                        ymin = meanInt - sdInt,
            #                        fill = dose_uM),
            #             alpha = 0.1) +
          ggtitle(paste0("SimKO: ",Compound_List[CounterCpd]," ; ",ModelOrig$ode_parameters$parNames[SelectedParamIdx]))+
            scale_y_continuous(limits=c(0,1))+
            theme_light()
          
          # theme_classic()
    )
    rm(ModelKO);rm(SimKO)
  }
  # print(g)
  dev.off()
  
  names(AvgAbsChange) <- names(ParamOrig)[IdxKO]
  AvgAbsChange <- as.data.frame(sort(AvgAbsChange,decreasing = T))
  AvgAbsChange <- cbind(AvgAbsChange,rownames(AvgAbsChange))
  colnames(AvgAbsChange) <- c("Change","Param")
  AvgAbsChange <- AvgAbsChange[order(AvgAbsChange[,1],decreasing = T),]
  AvgAbsChange$Param <- factor(AvgAbsChange$Param,levels=AvgAbsChange$Param)

  MaxKOChange[1,CounterCpd] <- as.character(AvgAbsChange[1,2])
  MaxKOChange[2,CounterCpd] <- AvgAbsChange[1,1]
  MaxKOChange[3,CounterCpd] <- as.character(AvgAbsChange[2,2])
  MaxKOChange[4,CounterCpd] <- AvgAbsChange[2,1]
  MaxKOChange[5,CounterCpd] <- as.character(AvgAbsChange[3,2])
  MaxKOChange[6,CounterCpd] <- AvgAbsChange[3,1]
  
  pdf(paste0("KO_AvgAccumApoptosis_",Compound_List[CounterCpd],".pdf"))
  
  print(ggplot(data=AvgAbsChange,aes(x=Param,y=Change)) + 
          geom_bar(stat="identity",aes(fill=Change)) +
          # scale_fill_brewer(palette="Spectral") + 
          # theme(axis.text.x=element_text(angle=90,hjust=1)) +
          labs(title = paste0("Accumulated %changes of Apoptosis per dose: ",Compound_List[CounterCpd]), subtitle = paste0("KO parameters")) +
          # xlab("time (days)") +
          # ylab("ng/ml") +
          theme_light()+
          theme(plot.title = element_text(size=15)) +
          theme(plot.subtitle = element_text(size=12)) +
          theme(axis.text.x = element_text(colour="grey20",size=7.5,angle=90,hjust=1,vjust=0,face="plain"),
                axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
                axis.title.x = element_blank(),
                axis.title.y = element_blank()))
  dev.off()
  
}

save(MaxKOChange,file = "MaxKOChangeApoptosis.RData")

# Add compound names at the end of the top parameter names
for (counter in 1:ncol(MaxKOChange)) {
  MaxKOChange[c(1,3,5),counter] <- paste0(MaxKOChange[c(1,3,5),counter],"_",colnames(MaxKOChange)[counter])
}

MaxKOChangePlot <- as.data.frame(t(cbind(MaxKOChange[c(1,2),],MaxKOChange[c(3,4),],MaxKOChange[c(5,6),])),stringsAsFactors = F)
MaxKOChangePlot[,2] <- as.numeric(MaxKOChangePlot[,2])
MaxKOChangePlot <- MaxKOChangePlot[order(MaxKOChangePlot[,2],decreasing = T),]
MaxKOChangePlot <- MaxKOChangePlot[-which(MaxKOChangePlot[,2]==0),] # Remove all entries with zero value
MaxKOChangePlot$Param1st <- factor(MaxKOChangePlot$Param1st,levels=MaxKOChangePlot$Param1st)
colnames(MaxKOChangePlot) <- c("Param","Change")

pdf(paste0("KO_MaxChangeApop_Top3.pdf"))

print(ggplot(data=MaxKOChangePlot,aes(x=Param,y=Change)) +
        geom_bar(stat="identity",aes(fill=Change)) +
        # scale_fill_brewer(palette="Spectral") +
        # theme(axis.text.x=element_text(angle=90,hjust=1)) +
        labs(title = paste0("Top 3 param for accumulated %changes of Apoptosis per dose"), subtitle = paste0("KO parameters")) +
        # xlab("time (days)") +
        # ylab("ng/ml") +
        theme_light()+
        theme(plot.title = element_text(size=15)) +
        theme(plot.subtitle = element_text(size=12)) +
        theme(axis.text.x = element_text(colour="grey20",size=6,angle=90,hjust=1,vjust=0,face="plain"),
              axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"),
              axis.title.x = element_blank(),
              axis.title.y = element_blank()))
dev.off()


# =========================================== #

# --- End of the script --- #
