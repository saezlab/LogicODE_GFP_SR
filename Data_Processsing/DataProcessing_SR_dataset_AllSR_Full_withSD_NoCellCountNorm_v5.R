# Data Processing pipeline for (new) Lukas' stress response dataset
# Full DDR pathway 

rm(list=ls());cat("\014");if(length(dev.list())>0){dev.off()}

# install.packages("dplyr")
library(dplyr)

# List all compounds
AllRpts <- c("BTG2","MDM2","P21","TP53","HMOX1","NRF2","SRXN1","ATF4","BIP","CHOP","XBP1","A20","ICAM1")
AllSignalVariablesGFP <- c("Cytoplasm_Intensity_IntegratedIntensity_image_GFP", # BTG2
                           "Nuclei_Intensity_MeanIntensity_image_GFP", # MDM2
                           "Nuclei_Intensity_MeanIntensity_image_GFP", # P21
                           "Nuclei_Intensity_MeanIntensity_image_GFP", # TP53
                           "Cytoplasm_Intensity_IntegratedIntensity_image_GFP", # HMOX1
                           "Nuclei_Intensity_IntegratedIntensity_image_GFP", # NRF2
                           "Cytoplasm_Intensity_IntegratedIntensity_image_GFP", # SRXN1
                           "Nuclei_Intensity_MeanIntensity_image_GFP", # ATF4
                           "Cytoplasm_Intensity_IntegratedIntensity_image_GFP", # BIP
                           "Nuclei_Intensity_MeanIntensity_image_GFP", # CHOP
                           "Nuclei_Intensity_MeanIntensity_image_GFP", # XBP1
                           "Cytoplasm_Intensity_IntegratedIntensity_image_GFP", # A20
                           "Cytoplasm_Intensity_IntegratedIntensity_image_GFP") # ICAM1

AllExtractedData <- list()
# AllCellCountData <- list()

# for (counter_rpt in 1:length(AllRpts)) {
for (counter_rpt in 1:11) {
  
  AllFilesGFP <- list.files(pattern=paste0(AllRpts[counter_rpt],"_GFP"))
  AllFilesAnVPI <- list.files(pattern=paste0(AllRpts[counter_rpt],"_AnVPI"))
  
  RawGFP <- NULL
  RawAnVPI <- NULL
  for (counter_file in 1:length(AllFilesGFP)) {
    print(paste0("Reading File: ", AllFilesGFP[counter_file]))
    RawGFP <- rbind(RawGFP,read.table(file = AllFilesGFP[counter_file],header = T,row.names = 1,sep = "\t",stringsAsFactors = F))
    print(paste0("Reading File: ", AllFilesAnVPI[counter_file]))
    RawAnVPI <- rbind(RawAnVPI,read.table(file = AllFilesAnVPI[counter_file],header = T,row.names = 1,sep = "\t",stringsAsFactors = F))
  }
  
  # Filter only relevant columns and rows
  # RawGFP <- RawGFP %>% select(one_of(c("treatment","dose_uM","timeID","replID","timeAfterExposure","variable","value")))
  # RawGFP <- RawGFP %>% filter(variable=="imageCountParentObj" | variable==AllSignalVariablesGFP[counter_rpt])
  # RawAnVPI <- RawAnVPI %>% select(one_of(c("treatment","dose_uM","timeID","replID","timeAfterExposure","variable","value")))
  # RawAnVPI <- RawAnVPI %>% filter(variable=="imageCountParentObj" | variable=="count_AnV_masked_primaryID_AreaShape_Area.DIV.Nuclei_AreaShape_Area_larger_0.1_" | variable=="count_PI_masked_primaryID_AreaShape_Area.DIV.Nuclei_AreaShape_Area_larger_0.1_")
  
  # Non dplyr version
  IdxColGFP <- which(colnames(RawGFP)=="treatment" | colnames(RawGFP)=="dose_uM" | colnames(RawGFP)=="timeID" | colnames(RawGFP)=="replID" | colnames(RawGFP)=="timeAfterExposure" | colnames(RawGFP)=="variable" | colnames(RawGFP)=="value")
  IdxRowGFP <- which(RawGFP$variable=="imageCountParentObj" | RawGFP$variable==AllSignalVariablesGFP[counter_rpt])
  RawGFP <- RawGFP[IdxRowGFP,IdxColGFP]
  
  IdxColAnVPI <- which(colnames(RawAnVPI)=="treatment" | colnames(RawAnVPI)=="dose_uM" | colnames(RawAnVPI)=="timeID" | colnames(RawAnVPI)=="replID" | colnames(RawAnVPI)=="timeAfterExposure" | colnames(RawAnVPI)=="variable" | colnames(RawAnVPI)=="value")
  IdxRowAnVPI <- which(RawAnVPI$variable=="imageCountParentObj" | RawAnVPI$variable=="count_AnV_masked_primaryID_AreaShape_Area.DIV.Nuclei_AreaShape_Area_larger_0.1_" | RawAnVPI$variable=="count_PI_masked_primaryID_AreaShape_Area.DIV.Nuclei_AreaShape_Area_larger_0.1_")
  RawAnVPI <- RawAnVPI[IdxRowAnVPI,IdxColAnVPI] 
  
  # ===== MANUAL CORRECTION OF INCONSISTENT FILE FORMATS ===== #
  
  # == Correct files that has the extra column "threshholdValues" == #
  # 1. "20180712_LeidenUniv_MDM2_GFP_R3_Clean.txt"
  # 2. "20180712_LeidenUniv_P21_GFP_R3_Clean.txt"
  # 3. "20180712_LeidenUniv_P21_GFP_R4_Clean.txt"
  # 4. "20180725_LeidenUniv_SRXN1_GFP_R2_Clean.txt"
  # 5. "20180725_LeidenUniv_SRXN1_GFP_R3_Clean.txt"
  # temp <- read.table(file = AllFilesGFP[counter_file],header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
  # temp <- temp[-which(colnames(temp)=="threshholdValues")]
  # write.table(x = temp,file = AllFilesGFP[counter_file],quote = F,sep = "\t")
  
  # == Correct files that has missing column "cell_line" == #
  # 1. "20180725_LeidenUniv_NRF2_GFP_R2_Clean_fitted.txt"
  # temp <- read.table(file = "20180725_LeidenUniv_NRF2_GFP_R2_Clean_fitted.txt",header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
  # temp2 <- read.table(file = "20180725_LeidenUniv_NRF2_GFP_R1_Clean.txt",header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
  # temp <- cbind(temp,rep(NA,nrow(temp)))
  # colnames(temp)[ncol(temp)] <- "cell_line"
  # temp <- temp[-which(colnames(temp)=="temp")]
  # temp <- temp[,c(1,2,3,4,10,5,6,8,7,9)]
  # write.table(x = temp,file = "20180725_LeidenUniv_NRF2_GFP_R2_Clean_fitted.txt",quote = F,sep = "\t")
  
  # == Correct files that do not have row number == #
  # 1. "20180728_LeidenUniv_BIP_GFP_R1_Clean.txt"
  # 2. "20180914_LeidenUniv_BIP_AnVPI_R1_Clean.txt"
  # RawGFP <- rbind(RawGFP,read.table(file = AllFilesGFP[counter_file],header = T,sep = "\t",stringsAsFactors = F))
  # write.table(x = RawGFP,file = AllFilesGFP[counter_file],quote = F,sep = "\t")
  # RawAnVPI <- rbind(RawAnVPI,read.table(file = AllFilesAnVPI[counter_file],header = T,sep = "\t",stringsAsFactors = F))
  # write.table(x = RawAnVPI,file = AllFilesAnVPI[counter_file],quote = F,sep = "\t")
  
  # == Correct files that have extra NA rows that couldn't be loaded == #
  # 1. "20180914_LeidenUniv_BIP_AnVPI_R4_Clean.txt"
  # Use Excel to read and remove NA columns at the bottom
  
  # ============================================== #
  
  AllCompounds <- sort(unique(RawGFP$treatment))
  # AllTreated <- AllCompounds[-c(which(AllCompounds=="DMSO "),which(AllCompounds=="HEP "))]
  # AllTreated <- AllCompounds[-c(which(AllCompounds=="DMSO "),which(AllCompounds=="HEP "),which(AllCompounds=="DMEM"))]
  AllTreated <- AllCompounds # Now collect everything
  AllTimeID <- sort(unique(RawGFP$timeID))
  AllreplID <- sort(unique(RawGFP$replID))
  
  # Processing/Labelling DoseID
  # RawGFP_CellCt <- RawGFP %>% filter(variable=="imageCountParentObj") %>% mutate(doseID=NA)
  
  RawGFP_CellCt <- RawGFP[which(RawGFP$variable=="imageCountParentObj"),]
  RawGFP_CellCt <- cbind(RawGFP_CellCt,rep(NA,nrow(RawGFP_CellCt)))
  colnames(RawGFP_CellCt)[ncol(RawGFP_CellCt)] <- "doseID"
  
  # for (counter in 1:length(unique(RawGFP_CellCt$treatment))) {
  #   print(paste0("Mapping GFP treatment: ",counter,"/",length(unique(RawGFP_CellCt$treatment))))
  #   RawGFP_Tx <- RawGFP_CellCt %>% filter(treatment==unique(RawGFP_CellCt$treatment)[counter])
  #   for (counter2 in 1:length(sort(unique(RawGFP_Tx$dose_uM)))) {
  #     RawGFP_CellCt$doseID[RawGFP_CellCt$treatment==unique(RawGFP_CellCt$treatment)[counter] &
  #                            RawGFP_CellCt$dose_uM==sort(unique(RawGFP_Tx$dose_uM))[counter2]] <- counter2
  #   }
  # }
  # 
  # # Real doses
  # AvgCellCtPerDose_GFP_RealDose <- RawGFP %>% filter(variable=="imageCountParentObj") %>% group_by(treatment,dose_uM) %>% summarise(meanVal=mean(value),stdVal=sd(value))
  # # library(ggplot2)
  # # pdf(paste0(AllRpts[counter_rpt],"_CellCount_GFP_withSD_RealDose.pdf"))
  # # print(ggplot(data = AvgCellCtPerDose_GFP_RealDose, aes(x=dose_uM,y=meanVal))
  # #       # + geom_pointrange(aes(ymin=meanVal-stdVal, ymax=meanVal+stdVal),size=0.2)
  # #       # + geom_smooth(method="loess", se=FALSE)
  # #       + geom_ribbon(aes(ymin=meanVal-stdVal, ymax=meanVal+stdVal),size=0.2,alpha=0.2)
  # #       + geom_point(aes(col=dose_uM),size=1)
  # #       + ylim(c(min(AvgCellCtPerDose_GFP_RealDose$meanVal),max(AvgCellCtPerDose_GFP_RealDose$meanVal)))
  # #       + facet_wrap( ~ treatment, scales="free")
  # #       + theme(axis.text.x=element_blank())
  # #       + ggtitle(paste0("Number of cell counts ",AllRpts[counter_rpt]," GFP"))
  # #       + theme(legend.position = "right",legend.title=element_text(size=10),legend.text=element_text(size=8)))
  # # dev.off()
  # 
  # # DoseID
  # # GFP
  # AvgCellCtPerDose_GFP <- RawGFP_CellCt %>% group_by(treatment,doseID) %>% summarise(meanVal=mean(value),stdVal=sd(value))
  # # pdf(paste0(AllRpts[counter_rpt],"_CellCount_GFP_withSD.pdf"))
  # # print(ggplot(data = AvgCellCtPerDose_GFP, aes(x=doseID,y=meanVal))
  # #       # + geom_pointrange(aes(ymin=meanVal-stdVal, ymax=meanVal+stdVal),size=0.2)
  # #       # + geom_smooth(method="loess", se=FALSE)
  # #       + geom_ribbon(aes(ymin=meanVal-stdVal, ymax=meanVal+stdVal),size=0.2,alpha=0.2)
  # #       + geom_point(aes(col=doseID),size=1)
  # #       + ylim(c(min(AvgCellCtPerDose_GFP$meanVal),max(AvgCellCtPerDose_GFP$meanVal)))
  # #       + facet_wrap( ~ treatment, scales="free")
  # #       + theme(axis.text.x=element_blank())
  # #       + ggtitle(paste0("Number of cell counts ",AllRpts[counter_rpt]," GFP"))
  # #       + theme(legend.position = "right",legend.title=element_text(size=10),legend.text=element_text(size=8)))
  # # dev.off()
  
  # # AnVPI
  # RawAnVPI_CellCt <- RawAnVPI %>% filter(variable=="imageCountParentObj") %>% mutate(doseID=NA)
  
  RawAnVPI_CellCt <- RawAnVPI[which(RawAnVPI$variable=="imageCountParentObj"),]
  RawAnVPI_CellCt <- cbind(RawAnVPI_CellCt,rep(NA,nrow(RawAnVPI_CellCt)))
  colnames(RawAnVPI_CellCt)[ncol(RawAnVPI_CellCt)] <- "doseID"
  
  # for (counter in 1:length(unique(RawAnVPI_CellCt$treatment))) {
  #   print(paste0("Mapping AnV Treatment: ",counter,"/",length(unique(RawAnVPI_CellCt$treatment))))
  #   RawAnVPI_Tx <- RawAnVPI_CellCt %>% filter(treatment==unique(RawAnVPI_CellCt$treatment)[counter])
  #   for (counter2 in 1:length(sort(unique(RawAnVPI_Tx$dose_uM)))) {
  #     RawAnVPI_CellCt$doseID[RawAnVPI_CellCt$treatment==unique(RawAnVPI_CellCt$treatment)[counter] &
  #                            RawAnVPI_CellCt$dose_uM==sort(unique(RawAnVPI_Tx$dose_uM))[counter2]] <- counter2
  #   }
  # }
  # AvgCellCtPerDose_AnVPI <- RawAnVPI_CellCt %>% group_by(treatment,doseID) %>% summarise(meanVal=mean(value),stdVal=sd(value))
  # pdf(paste0(AllRpts[counter_rpt],"_CellCount_AnVPI_withSD.pdf"))
  # print(ggplot(data = AvgCellCtPerDose_AnVPI, aes(x=doseID,y=meanVal))
  #       # + geom_pointrange(aes(ymin=meanVal-stdVal, ymax=meanVal+stdVal),size=0.2)
  #       # + geom_smooth(method="loess", se=FALSE)
  #       + geom_ribbon(aes(ymin=meanVal-stdVal, ymax=meanVal+stdVal),size=0.2,alpha=0.2)
  #       + geom_point(aes(col=doseID),size=1)
  #       + ylim(c(min(AvgCellCtPerDose_AnVPI$meanVal),max(AvgCellCtPerDose_AnVPI$meanVal)))
  #       + facet_wrap( ~ treatment, scales="free")
  #       + theme(axis.text.x=element_blank())
  #       + ggtitle(paste0("Number of cell counts ",AllRpts[counter_rpt]," AnVPI"))
  #       + theme(legend.position = "right",legend.title=element_text(size=10),legend.text=element_text(size=8)))
  # dev.off()
  
  
  # AllCellCountData[[counter_rpt]] <-
  #   list(Reporter = AllRpts[counter_rpt],
  #        AvgCellCtPerDose_GFP_RealDose = AvgCellCtPerDose_GFP_RealDose,
  #        AvgCellCtPerDose_GFP = AvgCellCtPerDose_GFP,
  #        AvgCellCtPerDose_AnVPI = AvgCellCtPerDose_AnVPI)
  
  # }
  
  # Prepare a list to store data from DMSO/DMEM (solvent) and from different treatments
  DMSOdoses <- sort(unique(RawGFP$dose_uM[which(RawGFP$treatment=="DMSO ")]))
  # AvgDMSO <- list(signal = matrix(NA,length(DMSOdoses),length(AllTimeID)),
  #                 AnV = matrix(NA,length(DMSOdoses),length(AllTimeID)),
  #                 PI = matrix(NA,length(DMSOdoses),length(AllTimeID)))
  AvgDMSO <- list(signal = array(NA,c(length(DMSOdoses),length(AllTimeID),length(AllreplID))),
                  AnV = array(NA,c(length(DMSOdoses),length(AllTimeID),length(AllreplID))),
                  PI = array(NA,c(length(DMSOdoses),length(AllTimeID),length(AllreplID))))
  for (counter_DMSO in 1:length(AvgDMSO)) {
    rownames(AvgDMSO[[counter_DMSO]]) <- paste0('d',DMSOdoses)
    colnames(AvgDMSO[[counter_DMSO]]) <- paste0('t',AllTimeID)
  }
  StdDMSO <- AvgDMSO
  
  DMEMdoses <- sort(unique(RawGFP$dose_uM[which(RawGFP$treatment=="DMEM")]))
  # AvgDMEM <- list(signal = matrix(NA,length(DMEMdoses),length(AllTimeID)),
  #                 AnV = matrix(NA,length(DMEMdoses),length(AllTimeID)),
  #                 PI = matrix(NA,length(DMEMdoses),length(AllTimeID)))
  AvgDMEM <- list(signal = array(NA,c(length(DMEMdoses),length(AllTimeID),length(AllreplID))),
                  AnV = array(NA,c(length(DMEMdoses),length(AllTimeID),length(AllreplID))),
                  PI = array(NA,c(length(DMEMdoses),length(AllTimeID),length(AllreplID))))
  for (counter_DMEM in 1:length(AvgDMEM)) {
    rownames(AvgDMEM[[counter_DMEM]]) <- paste0('d',DMEMdoses)
    colnames(AvgDMEM[[counter_DMEM]]) <- paste0('t',AllTimeID)
  }
  StdDMEM <- AvgDMEM
  
  AvgTreatment <- list()
  AllDoses <- list()
  for (counter_Tx in 1:length(AllTreated)) {
    TxDoses <- sort(unique(RawGFP$dose_uM[which(RawGFP$treatment==AllTreated[counter_Tx])]))
    AvgTreatment[[counter_Tx]] <- list(name = if (substr(AllTreated[counter_Tx],nchar(AllTreated[counter_Tx]),nchar(AllTreated[counter_Tx]))==" ") {
      substr(AllTreated[counter_Tx],1,nchar(AllTreated[counter_Tx])-1)
    } else { AllTreated[counter_Tx] },
    # signal = matrix(NA,length(TxDoses),length(AllTimeID)),
    # AnV = matrix(NA,length(TxDoses),length(AllTimeID)),
    # PI = matrix(NA,length(TxDoses),length(AllTimeID)))
    signal = array(NA,c(length(TxDoses),length(AllTimeID),length(AllreplID))),
    AnV = array(NA,c(length(TxDoses),length(AllTimeID),length(AllreplID))),
    PI = array(NA,c(length(TxDoses),length(AllTimeID),length(AllreplID))),
    CellCtGFP = array(NA,c(length(TxDoses),length(AllTimeID),length(AllreplID))),
    CellCtAnVPI = array(NA,c(length(TxDoses),length(AllTimeID),length(AllreplID))))
    for (counter_Dx in 2:length(AvgTreatment[[counter_Tx]])) {
      rownames(AvgTreatment[[counter_Tx]][[counter_Dx]]) <- paste0('d',TxDoses)
      colnames(AvgTreatment[[counter_Tx]][[counter_Dx]]) <- paste0('t',AllTimeID)
    }
    AllDoses[[counter_Tx]] <- TxDoses
  }
  StdTreatment <- AvgTreatment
  
  # Calculate average time from replicates and extract mean signals of DMSO (solvent control) and all compounds
  AllAvgTime <- rep(NA,length(AllTimeID))
  
  for (counter in 1:length(AllTimeID)) {
    print(" ")
    print(paste0("Mapping ",AllRpts[counter_rpt]," Data from timeID: ",counter,"/",length(AllTimeID)))
    print(" ")
    
    CurrentReplID <- unique(RawGFP$replID[which(RawGFP$timeID==AllTimeID[counter])])
    TimeTick <- rep(NA,length(CurrentReplID))
    for (counter2 in 1:length(TimeTick)) {
      TimeTick[counter2] <- unique(RawGFP$timeAfterExposure[which(RawGFP$timeID==AllTimeID[counter] & RawGFP$replID==CurrentReplID[counter2])])
    }
    AllAvgTime[counter] <- mean(TimeTick)
    
    for (counter_DMSO in 1:length(DMSOdoses)) {
      for (counter2_DMSO in 1:length(AllreplID)) {
        IdxSignal <- which(RawGFP$variable==AllSignalVariablesGFP[counter_rpt] & RawGFP$timeID==AllTimeID[counter] & RawGFP$dose_uM==DMSOdoses[counter_DMSO] & RawGFP$treatment=="DMSO " & RawGFP$replID==AllreplID[counter2_DMSO])
        IdxAnv <- which(RawAnVPI$variable=="count_AnV_masked_primaryID_AreaShape_Area.DIV.Nuclei_AreaShape_Area_larger_0.1_" & RawAnVPI$timeID==AllTimeID[counter] & RawAnVPI$dose_uM==DMSOdoses[counter_DMSO] & RawAnVPI$treatment=="DMSO " & RawAnVPI$replID==AllreplID[counter2_DMSO])
        IdxPI <- which(RawAnVPI$variable=="count_PI_masked_primaryID_AreaShape_Area.DIV.Nuclei_AreaShape_Area_larger_0.1_" & RawAnVPI$timeID==AllTimeID[counter] & RawAnVPI$dose_uM==DMSOdoses[counter_DMSO] & RawAnVPI$treatment=="DMSO " & RawAnVPI$replID==AllreplID[counter2_DMSO])
        AvgDMSO$signal[counter_DMSO,counter,counter2_DMSO] <- mean(RawGFP$value[IdxSignal])
        AvgDMSO$AnV[counter_DMSO,counter,counter2_DMSO] <- mean(RawAnVPI$value[IdxAnv])
        AvgDMSO$PI[counter_DMSO,counter,counter2_DMSO] <- mean(RawAnVPI$value[IdxPI])
        StdDMSO$signal[counter_DMSO,counter,counter2_DMSO] <- sd(RawGFP$value[IdxSignal])
        StdDMSO$AnV[counter_DMSO,counter,counter2_DMSO] <- sd(RawAnVPI$value[IdxAnv])
        StdDMSO$PI[counter_DMSO,counter,counter2_DMSO] <- sd(RawAnVPI$value[IdxPI])
      }
    }
    
    
    for (counter_DMEM in 1:length(DMEMdoses)) {
      for (counter2_DMEM in 1:length(AllreplID)) {
        IdxSignal <- which(RawGFP$variable==AllSignalVariablesGFP[counter_rpt] & RawGFP$timeID==AllTimeID[counter] & RawGFP$dose_uM==DMEMdoses[counter_DMEM] & RawGFP$treatment=="DMEM" & RawGFP$replID==AllreplID[counter2_DMEM])
        IdxAnv <- which(RawAnVPI$variable=="count_AnV_masked_primaryID_AreaShape_Area.DIV.Nuclei_AreaShape_Area_larger_0.1_" & RawAnVPI$timeID==AllTimeID[counter] & RawAnVPI$dose_uM==DMEMdoses[counter_DMEM] & RawAnVPI$treatment=="DMEM" & RawAnVPI$treatment=="DMEM" & RawAnVPI$replID==AllreplID[counter2_DMEM])
        IdxPI <- which(RawAnVPI$variable=="count_PI_masked_primaryID_AreaShape_Area.DIV.Nuclei_AreaShape_Area_larger_0.1_" & RawAnVPI$timeID==AllTimeID[counter] & RawAnVPI$dose_uM==DMEMdoses[counter_DMEM] & RawAnVPI$treatment=="DMEM" & RawAnVPI$treatment=="DMEM" & RawAnVPI$replID==AllreplID[counter2_DMEM])
        AvgDMEM$signal[counter_DMEM,counter,counter2_DMEM] <- mean(RawGFP$value[IdxSignal])
        AvgDMEM$AnV[counter_DMEM,counter,counter2_DMEM] <- mean(RawAnVPI$value[IdxAnv])
        AvgDMEM$PI[counter_DMEM,counter,counter2_DMEM] <- mean(RawAnVPI$value[IdxPI])
        StdDMEM$signal[counter_DMEM,counter,counter2_DMEM] <- sd(RawGFP$value[IdxSignal])
        StdDMEM$AnV[counter_DMEM,counter,counter2_DMEM] <- sd(RawAnVPI$value[IdxAnv])
        StdDMEM$PI[counter_DMEM,counter,counter2_DMEM] <- sd(RawAnVPI$value[IdxPI])
      }
    }
    
    for (counter4 in 1:length(AvgTreatment)) {
      print(paste0("Mapping ",AllRpts[counter_rpt]," Data from compounds: ",counter4,"/",length(AvgTreatment)))
      
      for (counter5 in 1:nrow(AvgTreatment[[counter4]]$signal)) {
        
        for (counter6 in 1:length(AllreplID)) { # Also separate replicates
          
          # IdxSignal <- which(RawGFP$variable==AllSignalVariablesGFP[counter_rpt] & RawGFP$timeID==AllTimeID[counter] & RawGFP$dose_uM==substr(rownames(AvgTreatment[[counter4]]$signal)[counter5],2,nchar(rownames(AvgTreatment[[counter4]]$signal)[counter5])) & RawGFP$treatment==AllTreated[counter4])
          # IdxAnv <- which(RawAnVPI$variable=="count_AnV_masked_primaryID_AreaShape_Area.DIV.Nuclei_AreaShape_Area_larger_0.1_" & RawAnVPI$timeID==AllTimeID[counter] & RawAnVPI$dose_uM==substr(rownames(AvgTreatment[[counter4]]$signal)[counter5],2,nchar(rownames(AvgTreatment[[counter4]]$signal)[counter5])) & RawAnVPI$treatment==AllTreated[counter4])
          # IdxPI <- which(RawAnVPI$variable=="count_PI_masked_primaryID_AreaShape_Area.DIV.Nuclei_AreaShape_Area_larger_0.1_" & RawAnVPI$timeID==AllTimeID[counter] & RawAnVPI$dose_uM==substr(rownames(AvgTreatment[[counter4]]$signal)[counter5],2,nchar(rownames(AvgTreatment[[counter4]]$signal)[counter5])) & RawAnVPI$treatment==AllTreated[counter4])
          IdxSignal <- which(RawGFP$variable==AllSignalVariablesGFP[counter_rpt] & RawGFP$timeID==AllTimeID[counter] & RawGFP$dose_uM==substr(rownames(AvgTreatment[[counter4]]$signal)[counter5],2,nchar(rownames(AvgTreatment[[counter4]]$signal)[counter5])) & RawGFP$treatment==AllTreated[counter4] & RawGFP$replID==AllreplID[counter6])
          IdxAnv <- which(RawAnVPI$variable=="count_AnV_masked_primaryID_AreaShape_Area.DIV.Nuclei_AreaShape_Area_larger_0.1_" & RawAnVPI$timeID==AllTimeID[counter] & RawAnVPI$dose_uM==substr(rownames(AvgTreatment[[counter4]]$signal)[counter5],2,nchar(rownames(AvgTreatment[[counter4]]$signal)[counter5])) & RawAnVPI$treatment==AllTreated[counter4] & RawAnVPI$replID==AllreplID[counter6])
          IdxPI <- which(RawAnVPI$variable=="count_PI_masked_primaryID_AreaShape_Area.DIV.Nuclei_AreaShape_Area_larger_0.1_" & RawAnVPI$timeID==AllTimeID[counter] & RawAnVPI$dose_uM==substr(rownames(AvgTreatment[[counter4]]$signal)[counter5],2,nchar(rownames(AvgTreatment[[counter4]]$signal)[counter5])) & RawAnVPI$treatment==AllTreated[counter4] & RawAnVPI$replID==AllreplID[counter6])
          IdxCellCtGFP <- which(RawGFP_CellCt$timeID==AllTimeID[counter] & RawGFP_CellCt$dose_uM==substr(rownames(AvgTreatment[[counter4]]$signal)[counter5],2,nchar(rownames(AvgTreatment[[counter4]]$signal)[counter5])) & RawGFP_CellCt$treatment==AllTreated[counter4] & RawGFP_CellCt$replID==AllreplID[counter6])
          IdxCellCtAnVPI <- which(RawAnVPI_CellCt$timeID==AllTimeID[counter] & RawAnVPI_CellCt$dose_uM==substr(rownames(AvgTreatment[[counter4]]$signal)[counter5],2,nchar(rownames(AvgTreatment[[counter4]]$signal)[counter5])) & RawAnVPI_CellCt$treatment==AllTreated[counter4] & RawAnVPI_CellCt$replID==AllreplID[counter6])
          # AvgTreatment[[counter4]]$signal[counter5,counter] <- mean(RawGFP$value[IdxSignal])
          # AvgTreatment[[counter4]]$AnV[counter5,counter] <- mean(RawAnVPI$value[IdxAnv])
          # AvgTreatment[[counter4]]$PI[counter5,counter] <- mean(RawAnVPI$value[IdxPI])
          # StdTreatment[[counter4]]$signal[counter5,counter] <- sd(RawGFP$value[IdxSignal])
          # StdTreatment[[counter4]]$AnV[counter5,counter] <- sd(RawAnVPI$value[IdxAnv])
          # StdTreatment[[counter4]]$PI[counter5,counter] <- sd(RawAnVPI$value[IdxPI])
          AvgTreatment[[counter4]]$signal[counter5,counter,counter6] <- mean(RawGFP$value[IdxSignal])
          AvgTreatment[[counter4]]$AnV[counter5,counter,counter6] <- mean(RawAnVPI$value[IdxAnv])
          AvgTreatment[[counter4]]$PI[counter5,counter,counter6] <- mean(RawAnVPI$value[IdxPI])
          AvgTreatment[[counter4]]$CellCtGFP[counter5,counter,counter6] <- mean(RawGFP_CellCt$value[IdxCellCtGFP])
          AvgTreatment[[counter4]]$CellCtAnVPI[counter5,counter,counter6] <- mean(RawAnVPI_CellCt$value[IdxCellCtAnVPI])
          
          StdTreatment[[counter4]]$signal[counter5,counter,counter6] <- sd(RawGFP$value[IdxSignal])
          StdTreatment[[counter4]]$AnV[counter5,counter,counter6] <- sd(RawAnVPI$value[IdxAnv])
          StdTreatment[[counter4]]$PI[counter5,counter,counter6] <- sd(RawAnVPI$value[IdxPI])
          StdTreatment[[counter4]]$CellCtGFP[counter5,counter,counter6] <- sd(RawGFP_CellCt$value[IdxCellCtGFP])
          StdTreatment[[counter4]]$CellCtAnVPI[counter5,counter,counter6] <- sd(RawAnVPI_CellCt$value[IdxCellCtAnVPI])
          
        }
      }
    }
    
  }
  
  # =======================================
  # Background subtraction with DMSO signal
  # =======================================
  
  # Calculate average background signal from multiple doses of DMSO/DMEM for each timeID
  
  # AvgDMSO_signal <- rep(NA,ncol(AvgDMSO$signal))
  # AvgDMSO_AnV <- rep(NA,ncol(AvgDMSO$AnV))
  # AvgDMSO_PI <- rep(NA,ncol(AvgDMSO$PI))
  AvgDMSO_signal <- matrix(NA,length(AllreplID),ncol(AvgDMSO$signal))
  AvgDMSO_AnV <- matrix(NA,length(AllreplID),ncol(AvgDMSO$AnV))
  AvgDMSO_PI <- matrix(NA,length(AllreplID),ncol(AvgDMSO$PI))
  for (counter_AvgDMSO in 1:length(AllreplID)) {
    for (counter2_AvgDMSO in 1:ncol(AvgDMSO_signal)) {
      # AvgDMSO_signal[counter_AvgDMSO,counter2_AvgDMSO] <- mean(AvgDMSO$signal[,counter_AvgDMSO])
      # AvgDMSO_AnV[counter_AvgDMSO,counter2_AvgDMSO] <- mean(AvgDMSO$AnV[,counter_AvgDMSO])
      # AvgDMSO_PI[counter_AvgDMSO,counter2_AvgDMSO] <- mean(AvgDMSO$PI[,counter_AvgDMSO])
      AvgDMSO_signal[counter_AvgDMSO,counter2_AvgDMSO] <- mean(AvgDMSO$signal[,counter2_AvgDMSO,counter_AvgDMSO],na.rm=T)
      AvgDMSO_AnV[counter_AvgDMSO,counter2_AvgDMSO] <- mean(AvgDMSO$AnV[,counter2_AvgDMSO,counter_AvgDMSO],na.rm=T)
      AvgDMSO_PI[counter_AvgDMSO,counter2_AvgDMSO] <- mean(AvgDMSO$PI[,counter2_AvgDMSO,counter_AvgDMSO],na.rm=T)
    }
  }
  
  # AvgDMEM_signal <- rep(NA,ncol(AvgDMEM$signal))
  # AvgDMEM_AnV <- rep(NA,ncol(AvgDMEM$AnV))
  # AvgDMEM_PI <- rep(NA,ncol(AvgDMEM$PI))
  AvgDMEM_signal <- matrix(NA,length(AllreplID),ncol(AvgDMEM$signal))
  AvgDMEM_AnV <- matrix(NA,length(AllreplID),ncol(AvgDMEM$AnV))
  AvgDMEM_PI <- matrix(NA,length(AllreplID),ncol(AvgDMEM$PI))
  for (counter_AvgDMEM in 1:length(AllreplID)) {
    for (counter2_AvgDMEM in 1:ncol(AvgDMEM_signal)) {
      # AvgDMEM_signal[counter_AvgDMEM,counter2_AvgDMEM] <- mean(AvgDMEM$signal[,counter_AvgDMEM])
      # AvgDMEM_AnV[counter_AvgDMEM,counter2_AvgDMEM] <- mean(AvgDMEM$AnV[,counter_AvgDMEM])
      # AvgDMEM_PI[counter_AvgDMEM,counter2_AvgDMEM] <- mean(AvgDMEM$PI[,counter_AvgDMEM])
      AvgDMEM_signal[counter_AvgDMEM,counter2_AvgDMEM] <- mean(AvgDMEM$signal[,counter2_AvgDMEM,counter_AvgDMEM],na.rm=T)
      AvgDMEM_AnV[counter_AvgDMEM,counter2_AvgDMEM] <- mean(AvgDMEM$AnV[,counter2_AvgDMEM,counter_AvgDMEM],na.rm=T)
      AvgDMEM_PI[counter_AvgDMEM,counter2_AvgDMEM] <- mean(AvgDMEM$PI[,counter2_AvgDMEM,counter_AvgDMEM],na.rm=T)
    }
  }
  
  
  # Subtract DMSO background by timeID
  
  AvgTreatment_BgSub <- AvgTreatment
  
  # for (counter_BgSub in 1:length(AllTimeID)) {
  #   for (counter2_BgSub in 1:length(AvgTreatment_BgSub)) {
  #     if(AvgTreatment_BgSub[[counter2_BgSub]]$name!="CDDP ") { # For any other molecule apart from Cisplatin (CDDP) -> DMSO subtraction
  #       AvgTreatment_BgSub[[counter2_BgSub]]$signal[,counter_BgSub] <- AvgTreatment[[counter2_BgSub]]$signal[,counter_BgSub] - AvgDMSO_signal[counter_BgSub]
  #       AvgTreatment_BgSub[[counter2_BgSub]]$AnV[,counter_BgSub] <- AvgTreatment[[counter2_BgSub]]$AnV[,counter_BgSub] - AvgDMSO_AnV[counter_BgSub]
  #       AvgTreatment_BgSub[[counter2_BgSub]]$PI[,counter_BgSub] <- AvgTreatment[[counter2_BgSub]]$PI[,counter_BgSub] - AvgDMSO_PI[counter_BgSub]
  #     } else { # For cisplatin -> DMEM subtraction
  #       AvgTreatment_BgSub[[counter2_BgSub]]$signal[,counter_BgSub] <- AvgTreatment[[counter2_BgSub]]$signal[,counter_BgSub] - AvgDMEM_signal[counter_BgSub]
  #       AvgTreatment_BgSub[[counter2_BgSub]]$AnV[,counter_BgSub] <- AvgTreatment[[counter2_BgSub]]$AnV[,counter_BgSub] - AvgDMEM_AnV[counter_BgSub]
  #       AvgTreatment_BgSub[[counter2_BgSub]]$PI[,counter_BgSub] <- AvgTreatment[[counter2_BgSub]]$PI[,counter_BgSub] - AvgDMEM_PI[counter_BgSub]
  #     }
  #   }
  # }
  
  for (counter_BgSub in 1:length(AllTimeID)) {
    for (counter2_BgSub in 1:length(AvgTreatment_BgSub)) {
      for(counter3_BgSub in 1:length(AllreplID)) {
        if (AvgTreatment_BgSub[[counter2_BgSub]]$name!="CDDP" & AvgTreatment_BgSub[[counter2_BgSub]]$name!="DMEM") { # For any other molecule apart from Cisplatin (CDDP) -> DMSO subtraction
          AvgTreatment_BgSub[[counter2_BgSub]]$signal[,counter_BgSub,counter3_BgSub] <- AvgTreatment[[counter2_BgSub]]$signal[,counter_BgSub,counter3_BgSub] - AvgDMSO_signal[counter3_BgSub,counter_BgSub]
          AvgTreatment_BgSub[[counter2_BgSub]]$AnV[,counter_BgSub,counter3_BgSub] <- AvgTreatment[[counter2_BgSub]]$AnV[,counter_BgSub,counter3_BgSub] - AvgDMSO_AnV[counter3_BgSub,counter_BgSub]
          AvgTreatment_BgSub[[counter2_BgSub]]$PI[,counter_BgSub,counter3_BgSub] <- AvgTreatment[[counter2_BgSub]]$PI[,counter_BgSub,counter3_BgSub] - AvgDMSO_PI[counter3_BgSub,counter_BgSub]
        } else { # For cisplatin or DMEM itself -> DMEM subtraction
          AvgTreatment_BgSub[[counter2_BgSub]]$signal[,counter_BgSub,counter3_BgSub] <- AvgTreatment[[counter2_BgSub]]$signal[,counter_BgSub,counter3_BgSub] - AvgDMEM_signal[counter3_BgSub,counter_BgSub]
          AvgTreatment_BgSub[[counter2_BgSub]]$AnV[,counter_BgSub,counter3_BgSub] <- AvgTreatment[[counter2_BgSub]]$AnV[,counter_BgSub,counter3_BgSub] - AvgDMEM_AnV[counter3_BgSub,counter_BgSub]
          AvgTreatment_BgSub[[counter2_BgSub]]$PI[,counter_BgSub,counter3_BgSub] <- AvgTreatment[[counter2_BgSub]]$PI[,counter_BgSub,counter3_BgSub] - AvgDMEM_PI[counter3_BgSub,counter_BgSub]
        }
      }
    }
  }
  
  
  # =======================================
  # (!!NEED A DISCUSSION!!) Set all negative values to be zero as signals are below DMSO (!!)
  # =======================================
  
  AvgTreatment_BgSub_NonZero <- AvgTreatment_BgSub
  
  # for (counter_NonZero in 1:length(AvgTreatment_BgSub_NonZero)) {
  #   for (counter2_NonZero in 2:length(AvgTreatment_BgSub_NonZero[[counter_NonZero]])) {
  #     AvgTreatment_BgSub_NonZero[[counter_NonZero]][[counter2_NonZero]][AvgTreatment_BgSub_NonZero[[counter_NonZero]][[counter2_NonZero]]<=0] <- 0
  #   }
  # }
  
  for (counter_NonZero in 1:length(AvgTreatment_BgSub_NonZero)) {
    for (counter2_NonZero in 2:length(AvgTreatment_BgSub_NonZero[[counter_NonZero]])) {
      # for (counter3_NonZero in 1:length(AllreplID)) {
      AvgTreatment_BgSub_NonZero[[counter_NonZero]][[counter2_NonZero]][AvgTreatment_BgSub_NonZero[[counter_NonZero]][[counter2_NonZero]]<=0] <- 0
      # }
    }
  }
  
  
  # rm(list = ls(pattern = "counter"))
  
  # Push extracted results from comprehensive list
  
  AllExtractedData[[counter_rpt]] <-
    list(Reporter = AllRpts[counter_rpt],
         AllCompounds = AllCompounds,
         AllTreated = AllTreated,
         AllTimeID = AllTimeID,
         AllreplID = AllreplID,
         AllDoses = AllDoses,
         DMSOdoses = DMSOdoses,
         AvgDMSO = AvgDMSO,
         StdDMSO = StdDMSO,
         DMEMdoses = DMEMdoses,
         AvgDMEM = AvgDMEM,
         StdDMEM = StdDMEM,
         AvgTreatment = AvgTreatment,
         StdTreatment = StdTreatment,
         AllAvgTime = AllAvgTime,
         AvgDMSO_signal = AvgDMSO_signal,
         AvgDMSO_AnV = AvgDMSO_AnV,
         AvgDMSO_PI = AvgDMSO_PI,
         AvgDMEM_signal = AvgDMEM_signal,
         AvgDMEM_AnV = AvgDMEM_AnV,
         AvgDMEM_PI = AvgDMEM_PI,
         AvgTreatment_BgSub = AvgTreatment_BgSub,
         AvgTreatment_BgSub_NonZero = AvgTreatment_BgSub_NonZero)
}

# save.image(file = "StressResponse_Variables_DDR_Corrected_withSD.RData")
# save.image(file = "StressResponse_Variables_AllSR_Corrected_withSD.RData")
# save.image(file = "StressResponse_Variables_AllSR_Corrected_withSD_CellCt_v2.RData")
save.image(file = "StressResponse_Variables_AllSR_Corrected_withSD_CellCt_v4.RData")

# ===========================================================
# ===========================================================
# ===========================================================

# Plot raw data

# rm(list=ls());cat("\014");if(length(dev.list())>0){dev.off()}

load("StressResponse_Variables_AllSR_Corrected_withSD_CellCt_v4.RData")

AllRpts <- c("BTG2","MDM2","P21","TP53","HMOX1","NRF2","SRXN1","ATF4","BIP","CHOP","XBP1","A20","ICAM1")

for (counter_rpt in 1:length(AllRpts)) {
  
  # Plot raw intensity
  library(reshape2)
  library(ggplot2)
  library(RColorBrewer)
  
  MeltedData_Signal <- NULL; MeltedSD_Signal <- NULL
  MeltedData_AnV <- NULL; MeltedSD_AnV <- NULL
  MeltedData_PI <- NULL; MeltedSD_PI <- NULL
  MeltedData_CellCtGFP <- NULL; MeltedSD_CellCtGFP <- NULL
  MeltedData_CellCtAnVPI <- NULL; MeltedSD_CellCtAnVPI <- NULL
  
  for (counter in 1:length(AllExtractedData[[counter_rpt]]$AllTreated)) {
    
    print(paste0("Processing ",AllRpts[counter_rpt]," compound: ",counter,"/",length(AllExtractedData[[counter_rpt]]$AllTreated)))
    
    # === New === # 
    
    TempExtractedData <- AllExtractedData
    
    # Normalise dataset -- REMOVE
    for (counter_DatSig in 1:length(AllExtractedData[[counter_rpt]]$AllreplID)) {
      # TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,,counter_DatSig] <- AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,,counter_DatSig]
      # TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,,counter_DatSig] <- AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,,counter_DatSig]
      # TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,,counter_DatSig] <- AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,,counter_DatSig]
      # TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$CellCtGFP[,,counter_DatSig] <- AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$CellCtGFP[,,counter_DatSig]
      # TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$CellCtAnVPI[,,counter_DatSig] <- AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$CellCtAnVPI[,,counter_DatSig]
      TempExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$signal[,,counter_DatSig] <- AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$signal[,,counter_DatSig]
      TempExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$AnV[,,counter_DatSig] <- AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$AnV[,,counter_DatSig]
      TempExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$PI[,,counter_DatSig] <- AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$PI[,,counter_DatSig]
      TempExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$CellCtGFP[,,counter_DatSig] <- AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$CellCtGFP[,,counter_DatSig]
      TempExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$CellCtAnVPI[,,counter_DatSig] <- AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$CellCtAnVPI[,,counter_DatSig]
    }
    # MeanGFPSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal,c(1,2),mean,na.rm=T)
    # StdGFPSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal,c(1,2),sd,na.rm=T)
    # MeanAnVSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV,c(1,2),mean,na.rm=T)
    # StdAnVSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV,c(1,2),sd,na.rm=T)
    # MeanPISig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI,c(1,2),mean,na.rm=T)
    # StdPISig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI,c(1,2),sd,na.rm=T)
    # MeanCellCtGFPSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$CellCtGFP,c(1,2),mean,na.rm=T)
    # StdCellCtGFPSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$CellCtGFP,c(1,2),sd,na.rm=T)
    # MeanCellCtAnVPISig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$CellCtAnVPI,c(1,2),mean,na.rm=T)
    # StdCellCtAnVPISig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$CellCtAnVPI,c(1,2),sd,na.rm=T)
    MeanGFPSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$signal,c(1,2),mean,na.rm=T)
    StdGFPSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$signal,c(1,2),sd,na.rm=T)
    MeanAnVSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$AnV,c(1,2),mean,na.rm=T)
    StdAnVSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$AnV,c(1,2),sd,na.rm=T)
    MeanPISig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$PI,c(1,2),mean,na.rm=T)
    StdPISig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$PI,c(1,2),sd,na.rm=T)
    MeanCellCtGFPSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$CellCtGFP,c(1,2),mean,na.rm=T)
    StdCellCtGFPSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$CellCtGFP,c(1,2),sd,na.rm=T)
    MeanCellCtAnVPISig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$CellCtAnVPI,c(1,2),mean,na.rm=T)
    StdCellCtAnVPISig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$CellCtAnVPI,c(1,2),sd,na.rm=T)
    
    # =========== #
    
    for (counter2 in 1:length(AllExtractedData[[counter_rpt]]$AllAvgTime)) {
      # print(paste0("Processing TimePoint: ",counter2,"/",length(AllAvgTime)))
      
      # === Data === #
      # DataCurrentTime_Signal <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$name,ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,counter2,]))),
      #                                 rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,counter2,]))),
      #                                 melt(MeanGFPSig[,counter2]),
      #                                 melt(StdGFPSig[,counter2]),
      #                                 as.character(1:ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,counter2,]))))
      DataCurrentTime_Signal <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$name,ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$signal[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$signal[,counter2,]))),
                                      rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$signal[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$signal[,counter2,]))),
                                      melt(MeanGFPSig[,counter2]),
                                      melt(StdGFPSig[,counter2]),
                                      as.character(1:ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$signal[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$signal[,counter2,]))))
      
      # DataCurrentTime_AnV <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$name,ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,counter2,]))),
      #                              rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,counter2,]))),
      #                              melt(MeanAnVSig[,counter2]),
      #                              melt(StdAnVSig[,counter2]),
      #                              as.character(1:ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,counter2,]))))
      DataCurrentTime_AnV <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$name,ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$AnV[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$AnV[,counter2,]))),
                                   rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$AnV[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$AnV[,counter2,]))),
                                   melt(MeanAnVSig[,counter2]),
                                   melt(StdAnVSig[,counter2]),
                                   as.character(1:ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$AnV[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$AnV[,counter2,]))))
      
      # DataCurrentTime_PI <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$name,ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,counter2,]))),
      #                             rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,counter2,]))),
      #                             melt(MeanPISig[,counter2]),
      #                             melt(StdPISig[,counter2]),
      #                             as.character(1:ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,counter2,]))))
      DataCurrentTime_PI <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$name,ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$PI[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$PI[,counter2,]))),
                                  rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$PI[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$PI[,counter2,]))),
                                  melt(MeanPISig[,counter2]),
                                  melt(StdPISig[,counter2]),
                                  as.character(1:ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$PI[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$PI[,counter2,]))))
      
      DataCurrentTime_CellCtGFP <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$name,ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$CellCtGFP[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$CellCtGFP[,counter2,]))),
                                         rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$CellCtGFP[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$CellCtGFP[,counter2,]))),
                                         melt(MeanCellCtGFPSig[,counter2]),
                                         melt(StdCellCtGFPSig[,counter2]),
                                         as.character(1:ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$CellCtGFP[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$CellCtGFP[,counter2,]))))
      
      DataCurrentTime_CellCtAnVPI <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$name,ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$CellCtAnVPI[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$CellCtAnVPI[,counter2,]))),
                                           rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$CellCtAnVPI[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$CellCtAnVPI[,counter2,]))),
                                           melt(MeanCellCtAnVPISig[,counter2]),
                                           melt(StdCellCtAnVPISig[,counter2]),
                                           as.character(1:ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$CellCtAnVPI[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment[[counter]]$CellCtAnVPI[,counter2,]))))
      
      
      MeltedData_Signal <- rbind(MeltedData_Signal,DataCurrentTime_Signal)
      MeltedData_AnV <- rbind(MeltedData_AnV,DataCurrentTime_AnV)
      MeltedData_PI <- rbind(MeltedData_PI,DataCurrentTime_PI)
      MeltedData_CellCtGFP <- rbind(MeltedData_CellCtGFP,DataCurrentTime_CellCtGFP)
      MeltedData_CellCtAnVPI <- rbind(MeltedData_CellCtAnVPI,DataCurrentTime_CellCtAnVPI)
      
    }
  }
  
  colnames(MeltedData_Signal) <- c("treatment","time","signal","std","doseNr")
  colnames(MeltedData_AnV) <- c("treatment","time","signal","std","doseNr")
  colnames(MeltedData_PI) <- c("treatment","time","signal","std","doseNr")
  colnames(MeltedData_CellCtGFP) <- c("treatment","time","signal","std","doseNr")
  colnames(MeltedData_CellCtAnVPI) <- c("treatment","time","signal","std","doseNr")
  
  # 25.09.18 --- Introduce removing zero entries here
  
  # MeltedData_Signal$signal[which(MeltedData_Signal$signal<0)] <- 0
  # MeltedData_Signal$std[which(MeltedData_Signal$signal<0)] <- 0
  # MeltedData_AnV$signal[which(MeltedData_AnV$signal<0)] <- 0
  # MeltedData_AnV$std[which(MeltedData_AnV$signal<0)] <- 0
  # MeltedData_PI$signal[which(MeltedData_PI$signal<0)] <- 0
  # MeltedData_PI$std[which(MeltedData_PI$signal<0)] <- 0
  
  
  # Remove normalisation step
  NormMeltedData_Signal <- MeltedData_Signal; NormMeltedData_Signal$signal <- NormMeltedData_Signal$signal
  NormMeltedData_AnV <- MeltedData_AnV; NormMeltedData_AnV$signal <- NormMeltedData_AnV$signal
  NormMeltedData_PI <- MeltedData_PI; NormMeltedData_PI$signal <- NormMeltedData_PI$signal
  NormMeltedData_CellCtGFP <- MeltedData_CellCtGFP; NormMeltedData_CellCtGFP$signal <- NormMeltedData_CellCtGFP$signal
  NormMeltedData_CellCtAnVPI <- MeltedData_CellCtAnVPI; NormMeltedData_CellCtAnVPI$signal <- NormMeltedData_CellCtAnVPI$signal
  
  # Plot profiles of signal, AnV and PI
  pdf(paste0(AllRpts[counter_rpt],"_Signal_Raw_withSD_PerRepl.pdf"))
  print(ggplot(data = NormMeltedData_Signal, aes(x=time,y=signal))
        + geom_pointrange(aes(ymin=signal-std, ymax=signal+std),size=0.1,alpha=0.2)
        # + geom_smooth(method="loess", se=FALSE)
        # + geom_ribbon(aes(ymin=signal-std, ymax=signal+std),size=0.2,alpha=0.2)
        + geom_point(aes(col=doseNr),size=1)
        # + ylim(c(-1e-10,1)) 
        + facet_wrap( ~ treatment, scales="free")
        + theme(axis.text.x=element_blank())
        + ggtitle(paste0("Raw Data ",AllRpts[counter_rpt]," GFP signal"))
        + theme(legend.position = "right",legend.title=element_text(size=10),legend.text=element_text(size=8)))
  dev.off()
  
  pdf(paste0(AllRpts[counter_rpt],"_AnV_Raw_withSD_PerRepl.pdf"))
  print(ggplot(data = NormMeltedData_AnV, aes(x=time,y=signal))
        + geom_pointrange(aes(ymin=signal-std, ymax=signal+std),size=0.1,alpha=0.2)
        # + geom_smooth(method="loess", se=FALSE)
        # + geom_ribbon(aes(ymin=signal-std, ymax=signal+std),size=0.2,alpha=0.2)
        + geom_point(aes(col=doseNr),size=1)
        + ylim(c(-0.25,1.1))
        + facet_wrap( ~ treatment, scales="free")
        + theme(axis.text.x=element_blank())
        + ggtitle(paste0("Raw Data ",AllRpts[counter_rpt]," AnV"))
        + theme(legend.position = "right",legend.title=element_text(size=10),legend.text=element_text(size=8)))
  dev.off()
  
  pdf(paste0(AllRpts[counter_rpt],"_PI_Raw_withSD_PerRepl.pdf"))
  print(ggplot(data = NormMeltedData_PI, aes(x=time,y=signal))
        + geom_pointrange(aes(ymin=signal-std, ymax=signal+std),size=0.1,alpha=0.2)
        # + geom_smooth(method="loess", se=FALSE)
        # + geom_ribbon(aes(ymin=signal-std, ymax=signal+std),size=0.2,alpha=0.2)
        + geom_point(aes(col=doseNr),size=1)
        + ylim(c(-0.25,1))
        + facet_wrap( ~ treatment, scales="free")
        + theme(axis.text.x=element_blank())
        + ggtitle(paste0("Raw Data ",AllRpts[counter_rpt]," PI"))
        + theme(legend.position = "right",legend.title=element_text(size=10),legend.text=element_text(size=8)))
  dev.off()
  
  # Plot cell counts GFP
  pdf(paste0(AllRpts[counter_rpt],"_CellCtGFP_Raw_withSD_PerRepl.pdf"))
  print(ggplot(data = NormMeltedData_CellCtGFP, aes(x=time,y=signal))
        + geom_pointrange(aes(ymin=signal-std, ymax=signal+std),size=0.1,alpha=0.2)
        # + geom_smooth(method="loess", se=FALSE)
        # + geom_ribbon(aes(ymin=signal-std, ymax=signal+std),size=0.2,alpha=0.2)
        + geom_point(aes(col=doseNr),size=1)
        # + ylim(c(-0.25,1))
        + facet_wrap( ~ treatment, scales="free")
        + theme(axis.text.x=element_blank())
        + ggtitle(paste0("Raw Data ",AllRpts[counter_rpt]," CellCtGFP"))
        + theme(legend.position = "right",legend.title=element_text(size=10),legend.text=element_text(size=8)))
  dev.off()
  
  # Plot cell counts AnVPI
  pdf(paste0(AllRpts[counter_rpt],"_CellCtAnVPI_Raw_withSD_PerRepl.pdf"))
  print(ggplot(data = NormMeltedData_CellCtAnVPI, aes(x=time,y=signal))
        + geom_pointrange(aes(ymin=signal-std, ymax=signal+std),size=0.1,alpha=0.2)
        # + geom_smooth(method="loess", se=FALSE)
        # + geom_ribbon(aes(ymin=signal-std, ymax=signal+std),size=0.2,alpha=0.2)
        + geom_point(aes(col=doseNr),size=1)
        # + ylim(c(-0.25,1))
        + facet_wrap( ~ treatment, scales="free")
        + theme(axis.text.x=element_blank())
        + ggtitle(paste0("Raw Data ",AllRpts[counter_rpt]," CellCtAnVPI"))
        + theme(legend.position = "right",legend.title=element_text(size=10),legend.text=element_text(size=8)))
  dev.off()
  
  
  
}


# ===============================================
# ===============================================
# ===============================================

# Calculate and plot average AnV, PI, CellCt across reporters for each compound

# !!! Issue with non-overlapping time and timeID from each reporter

# rm(list=ls());cat("\014");if(length(dev.list())>0){dev.off()}
# 
# load("StressResponse_Variables_AllSR_Corrected_withSD_CellCt_v4.RData")
# 
# # AllCompounds # Loaded from RData
# AllRpts <- c("BTG2","MDM2","P21","TP53","HMOX1","NRF2","SRXN1","ATF4","BIP","CHOP","XBP1","A20","ICAM1")
# # Measures <- c("AnV","PI","CtGFP","CtAnVPI")
# # 
# # ExtractedPhenotype <- array(data = NA,dim = c(length(AllCompounds),length(AllRpts),length(Measures)),dimnames = list(AllCompounds,AllRpts,Measures))
# # 
# # for (counter in 1:length(AllCompounds)) {
# #   for (counter2 in 1:length(AllRpts)) {
# #     ExtractedPhenotype[counter,counter2,1] <- AllExtractedData[[counter2]]$AvgTreatment_BgSub[[counter]]$AnV
# #     ExtractedPhenotype[counter,counter2,2] <- AllExtractedData[[counter2]]$AvgTreatment_BgSub[[counter]]$PI
# #     ExtractedPhenotype[counter,counter2,3] <- AllExtractedData[[counter2]]$AvgTreatment_BgSub[[counter]]$CellCtGFP
# #     ExtractedPhenotype[counter,counter2,4] <- AllExtractedData[[counter2]]$AvgTreatment_BgSub[[counter]]$CellCtAnVPI
# #   }
# # }
# 
# MeltedData_AnV <- NULL; MeltedSD_AnV <- NULL
# MeltedData_PI <- NULL; MeltedSD_PI <- NULL
# 
# ExtractedPhenotype <- list()
# 
# for (counter_rpt in 1:length(AllRpts)) {
#   
#   MaxValueByCompound <- list("signal"=matrix(NA,length(AllExtractedData[[counter_rpt]]$AllTreated),length(AllExtractedData[[counter_rpt]]$AllreplID)),
#                              "AnV"=matrix(NA,length(AllExtractedData[[counter_rpt]]$AllTreated),length(AllExtractedData[[counter_rpt]]$AllreplID)),
#                              "PI"=matrix(NA,length(AllExtractedData[[counter_rpt]]$AllTreated),length(AllExtractedData[[counter_rpt]]$AllreplID)))
#   # for (counter in 1:length(AllExtractedData[[counter_rpt]]$AllTreated)) {
#   #     MaxValueByCompound$signal[counter] <- max(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal,na.rm = T)
#   #     MaxValueByCompound$AnV[counter] <- max(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV,na.rm = T)
#   #     MaxValueByCompound$PI[counter] <- max(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI,na.rm = T)
#   # }
#   for (counter in 1:length(AllExtractedData[[counter_rpt]]$AllTreated)) {
#     for (counter2 in 1:length(AllExtractedData[[counter_rpt]]$AllreplID)) {
#       MaxValueByCompound$signal[counter,counter2] <- max(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,,counter2],na.rm = T)
#       MaxValueByCompound$AnV[counter,counter2] <- max(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,,counter2],na.rm = T)
#       MaxValueByCompound$PI[counter,counter2] <- max(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,,counter2],na.rm = T)
#     }
#   }
#   
#   # Missing data returned -Inf -> turn to NA
#   MaxValueByCompound$signal[which(MaxValueByCompound$signal==-Inf)] <- NA
#   MaxValueByCompound$AnV[which(MaxValueByCompound$AnV==-Inf)] <- NA
#   MaxValueByCompound$PI[which(MaxValueByCompound$PI==-Inf)] <- NA
#   
#   for (counter in 1:length(AllExtractedData[[counter_rpt]]$AllTreated)) {
#     # for (counter in 1:8) {
#     print(paste0("Processing ",AllRpts[counter_rpt]," compound: ",counter,"/",length(AllExtractedData[[counter_rpt]]$AllTreated)))
#     
#     # === New === # 
#     
#     TempExtractedData <- AllExtractedData
#     
#     # Normalise dataset
#     for (counter_DatSig in 1:length(AllExtractedData[[counter_rpt]]$AllreplID)) {
#       # TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,,counter_DatSig] <- AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,,counter_DatSig]/max(MaxValueByCompound$signal[,counter_DatSig],na.rm = T)
#       # TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,,counter_DatSig] <- AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,,counter_DatSig]/max(MaxValueByCompound$signal[,counter_DatSig],na.rm = T)
#       # TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,,counter_DatSig] <- AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,,counter_DatSig]/max(MaxValueByCompound$signal[,counter_DatSig],na.rm = T)
#       # TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,,counter_DatSig] <- AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,,counter_DatSig]/max(MaxValueByCompound$signal[,counter_DatSig],na.rm = T)
#       TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,,counter_DatSig] <- AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,,counter_DatSig]/max(MaxValueByCompound$signal[,counter_DatSig],na.rm = T)
#       TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,,counter_DatSig] <- AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,,counter_DatSig]/max(MaxValueByCompound$signal[,counter_DatSig],na.rm = T)
#     }
#     # MeanGFPSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal,c(1,2),mean,na.rm=T)
#     # StdGFPSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal,c(1,2),sd,na.rm=T)
#     # MeanAnVSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV,c(1,2),mean,na.rm=T)
#     # StdAnVSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV,c(1,2),sd,na.rm=T)
#     # MeanPISig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI,c(1,2),mean,na.rm=T)
#     # StdPISig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI,c(1,2),sd,na.rm=T)
#     # MeanGFPSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal,c(1,2),mean,na.rm=T)
#     # StdGFPSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal,c(1,2),sd,na.rm=T)
#     MeanAnVSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV,c(1,2),mean,na.rm=T)
#     StdAnVSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV,c(1,2),sd,na.rm=T)
#     MeanPISig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI,c(1,2),mean,na.rm=T)
#     StdPISig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI,c(1,2),sd,na.rm=T)
#     
#     # NormMeanGFPSig <- MeanGFPSig/max(MeanGFPSig,na.rm = T)
#     # NormStdGFPSig <- StdGFPSig/max(MeanGFPSig,na.rm = T)
#     # NormMeanAnVSig <- MeanAnVSig/max(MeanAnVSig,na.rm = T)
#     # NormStdAnVSig <- StdAnVSig/max(MeanAnVSig,na.rm = T)
#     # NormMeanPISig <- MeanPISig/max(MeanPISig,na.rm = T)
#     # NormStdPISig <- StdPISig/max(MeanPISig,na.rm = T)
#     
#     # # Missing data returned -Inf -> turn to NA
#     # NormMeanGFPSig[which(NormMeanGFPSig==-Inf)] <- NA
#     # NormStdGFPSig[which(NormStdGFPSig==-Inf)] <- NA
#     # NormMeanAnVSig[which(NormMeanAnVSig==-Inf)] <- NA
#     # NormStdAnVSig[which(NormStdAnVSig==-Inf)] <- NA
#     # NormMeanPISig[which(NormMeanPISig==-Inf)] <- NA
#     # NormStdPISig[which(NormStdPISig==-Inf)] <- NA
# 
#   }
#   ExtractedPhenotype[[counter_rpt]] <- list(
#     Reporter = AllRpts[counter_rpt],
#     MeanAnVSig = MeanAnVSig,
#     MeanPISig = MeanPISig)
# }
# 
# AvgExtractedPhenotype <- ExtractedPhenotype[[1]]
# AvgExtractedPhenotype$Reporter <- "AllRpts"
# AvgExtractedPhenotype$MeanAnVSig <- array(0,dim(AvgExtractedPhenotype$MeanAnVSig)) # Reset values to zero
# AvgExtractedPhenotype$MeanPISig <-  array(0,dim(AvgExtractedPhenotype$MeanPISig)) # Reset values to zero
# 
# for (counter in 1:length(ExtractedPhenotype)) {
#   AvgExtractedPhenotype$MeanAnVSig <- AvgExtractedPhenotype$MeanAnVSig + ExtractedPhenotype[[counter]]$MeanAnVSig
#   AvgExtractedPhenotype$MeanPISig <- AvgExtractedPhenotype$MeanPISig + ExtractedPhenotype[[counter]]$MeanPISig
# }
# AvgExtractedPhenotype$MeanAnVSig <- AvgExtractedPhenotype$MeanAnVSig/length(ExtractedPhenotype)
# AvgExtractedPhenotype$MeanPISig <- AvgExtractedPhenotype$MeanPISig/length(ExtractedPhenotype)
# 
# 
#     
#     # =========== #
#     
#     for (counter2 in 1:length(AllExtractedData[[counter_rpt]]$AllAvgTime)) {
#       # print(paste0("Processing TimePoint: ",counter2,"/",length(AllAvgTime)))
#       
#       # === Data === #
#       
#       # DataCurrentTime_Signal <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$name,length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,counter2])),
#       #                                 rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,counter2])),
#       #                                 melt(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,counter2]/max(MaxValueByCompound$signal)),
#       #                                 as.character(1:length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,counter2])))
#       
#       # DataCurrentTime_Signal <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$name,length(rowMeans(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,counter2,]))),
#       #                                rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],length(rowMeans(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,counter2,]))),
#       # DataCurrentTime_Signal <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$name,ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,counter2,]))),
#       #                                 rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,counter2,]))),
#       #                                 melt(MeanGFPSig[,counter2]),
#       #                                 melt(StdGFPSig[,counter2]),
#       #                                 #                                as.character(1:length(rowMeans(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,counter2,]))))
#       #                                 as.character(1:ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,counter2,]))))
#       
#       # DataCurrentTime_AnV <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$name,length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,counter2])),
#       #                              rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,counter2])),
#       #                              melt(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,counter2]/max(MaxValueByCompound$AnV)),
#       #                              as.character(1:length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,counter2])))
#       
#       # DataCurrentTime_AnV <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$name,length(rowMeans(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,counter2,]))),
#       #                                 rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],length(rowMeans(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,counter2,]))),
#       DataCurrentTime_AnV <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$name,ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,counter2,]))),
#                                    rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,counter2,]))),
#                                    melt(MeanAnVSig[,counter2]),
#                                    melt(StdAnVSig[,counter2]),
#                                    #                              as.character(1:length(rowMeans(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,counter2,]))))
#                                    as.character(1:ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,counter2,]))))
#       
#       # DataCurrentTime_PI <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$name,length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,counter2])),
#       #                             rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,counter2])),
#       #                             melt(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,counter2]/max(MaxValueByCompound$PI)),
#       #                             as.character(1:length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,counter2])))
#       
#       # DataCurrentTime_PI <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$name,length(rowMeans(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,counter2,]))),
#       #                              rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],length(rowMeans(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,counter2,]))),
#       DataCurrentTime_PI <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$name,ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,counter2,]))),
#                                   rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,counter2,]))),
#                                   melt(MeanPISig[,counter2]),
#                                   melt(StdPISig[,counter2]),
#                                   #                             as.character(1:length(rowMeans(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,counter2,]))))
#                                   as.character(1:ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,counter2,]))))
#       
#       # MeltedData_Signal <- rbind(MeltedData_Signal,DataCurrentTime_Signal)
#       MeltedData_AnV <- rbind(MeltedData_AnV,DataCurrentTime_AnV)
#       MeltedData_PI <- rbind(MeltedData_PI,DataCurrentTime_PI)
#       
#       # === SD === #
#       
#       # SDCurrentTime_Signal <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$name,length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,counter2])),
#       #                                 rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,counter2])),
#       #                                 melt(AllExtractedData[[counter_rpt]]$StdTreatment[[counter]]$signal[,counter2]/max(MaxValueByCompound$signal)),
#       #                                 as.character(1:length(AllExtractedData[[counter_rpt]]$StdTreatment[[counter]]$signal[,counter2])))
#       # 
#       # SDCurrentTime_AnV <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$name,length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,counter2])),
#       #                              rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,counter2])),
#       #                              melt(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,counter2]/max(MaxValueByCompound$AnV)),
#       #                              as.character(1:length(AllExtractedData[[counter_rpt]]$StdTreatment[[counter]]$AnV[,counter2])))
#       # 
#       # SDCurrentTime_PI <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$name,length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,counter2])),
#       #                             rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,counter2])),
#       #                             melt(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,counter2]/max(MaxValueByCompound$PI)),
#       #                             as.character(1:length(AllExtractedData[[counter_rpt]]$StdTreatment[[counter]]$PI[,counter2])))
#       # 
#       # MeltedSD_Signal <- rbind(MeltedSD_Signal,SDCurrentTime_Signal)
#       # MeltedSD_AnV <- rbind(MeltedSD_AnV,SDCurrentTime_AnV)
#       # MeltedSD_PI <- rbind(MeltedSD_PI,SDCurrentTime_PI)
#       
#     }
#   # }
#   
#   # colnames(MeltedData_Signal) <- c("treatment","time","signal","doseNr"); colnames(MeltedSD_Signal) <- colnames(MeltedData_Signal)
#   # colnames(MeltedData_AnV) <- c("treatment","time","signal","doseNr"); colnames(MeltedSD_AnV) <- colnames(MeltedData_AnV)
#   # colnames(MeltedData_PI) <- c("treatment","time","signal","doseNr"); colnames(MeltedSD_PI) <- colnames(MeltedData_PI)
#   
#   # colnames(MeltedData_Signal) <- c("treatment","time","signal","std","doseNr")
#   colnames(MeltedData_AnV) <- c("treatment","time","signal","std","doseNr")
#   colnames(MeltedData_PI) <- c("treatment","time","signal","std","doseNr")
#   
#   # 25.09.18 --- Introduce removing zero entries here
#   
#   # MeltedData_Signal$signal[which(MeltedData_Signal$signal<0)] <- 0
#   # MeltedData_Signal$std[which(MeltedData_Signal$signal<0)] <- 0
#   MeltedData_AnV$signal[which(MeltedData_AnV$signal<0)] <- 0
#   MeltedData_AnV$std[which(MeltedData_AnV$signal<0)] <- 0
#   MeltedData_PI$signal[which(MeltedData_PI$signal<0)] <- 0
#   MeltedData_PI$std[which(MeltedData_PI$signal<0)] <- 0
#   
#   # NormMeltedData_Signal <- MeltedData_Signal; NormMeltedData_Signal$signal <- NormMeltedData_Signal$signal/max(NormMeltedData_Signal$signal,na.rm = T)
#   NormMeltedData_AnV <- MeltedData_AnV; NormMeltedData_AnV$signal <- NormMeltedData_AnV$signal/max(NormMeltedData_AnV$signal,na.rm = T)
#   NormMeltedData_PI <- MeltedData_PI; NormMeltedData_PI$signal <- NormMeltedData_PI$signal/max(NormMeltedData_PI$signal,na.rm = T)
#   
# # }

# ===========================================================
# ===========================================================
# ===========================================================

# # [Debugging] Retrieve back mapped information
# rm(list=ls());cat("\014");if(length(dev.list())>0){dev.off()}
# # load("StressResponse_Variables_DDR_Corrected_withSD.RData")
# # load("StressResponse_Variables_AllSR_Corrected_withSD.RData")
# load("StressResponse_Variables_AllSR_Corrected_withSD_CellCt.RData")
# 
# # counter_rpt=1
# # 
# # Reporter = AllExtractedData[[counter_rpt]]$Reporter
# # AllCompounds = AllExtractedData[[counter_rpt]]$AllCompounds
# # AllTreated = AllExtractedData[[counter_rpt]]$AllTreated
# # AllTimeID = AllExtractedData[[counter_rpt]]$AllTimeID
# # AllreplID = AllExtractedData[[counter_rpt]]$AllreplID
# # AllDoses = AllExtractedData[[counter_rpt]]$AllDoses
# # DMSOdoses = AllExtractedData[[counter_rpt]]$DMSOdoses
# # AvgDMSO = AllExtractedData[[counter_rpt]]$AvgDMSO
# # StdDMSO = AllExtractedData[[counter_rpt]]$StdDMSO
# # DMEMdoses = AllExtractedData[[counter_rpt]]$DMEMdoses
# # AvgDMEM = AllExtractedData[[counter_rpt]]$AvgDMEM
# # StdDMEM = AllExtractedData[[counter_rpt]]$StdDMEM
# # AvgTreatment = AllExtractedData[[counter_rpt]]$AvgTreatment
# # StdTreatment = AllExtractedData[[counter_rpt]]$StdTreatment
# # AllAvgTime = AllExtractedData[[counter_rpt]]$AllAvgTime
# # AvgDMSO_signal = AllExtractedData[[counter_rpt]]$AvgDMSO_signal
# # AvgDMSO_AnV = AllExtractedData[[counter_rpt]]$AvgDMSO_AnV
# # AvgDMSO_PI = AllExtractedData[[counter_rpt]]$AvgDMSO_PI
# # AvgDMEM_signal = AllExtractedData[[counter_rpt]]$AvgDMEM_signal
# # AvgDMEM_AnV = AllExtractedData[[counter_rpt]]$AvgDMEM_AnV
# # AvgDMEM_PI = AllExtractedData[[counter_rpt]]$AvgDMEM_PI
# # AvgTreatment_BgSub = AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub
# # AvgTreatment_BgSub_NonZero = AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero

# =========================================================== 
# ===========================================================
# ===========================================================
# 
# # Load results from BQ cluster
# AllRpts <- c("BTG2","MDM2","P21","TP53","HMOX1","NRF2","SRXN1","ATF4","BIP","CHOP","XBP1")
# 
# AllExtractedData_Parallel <- list()
# 
# for (counter_rpt in 1:11) {
# #for (counter_rpt in 4:11) {
#     load(paste0("SR_Mapped_Full/SR_Mapped_Parallel_",counter_rpt,"_",AllRpts[counter_rpt],".RData"))
#   AllExtractedData_Parallel[[counter_rpt]] <- AllExtractedData[[counter_rpt]]
# }
# 
# AllExtractedData <- AllExtractedData_Parallel
# 
# rm(AllDoses,AllExtractedData_Parallel,AvgDMEM,AvgDMEM_AnV,AvgDMEM_PI,AvgDMEM_signal,AvgDMSO,AvgDMSO_AnV,AvgDMSO_PI,AvgDMSO_signal,
#    AvgTreatment,AvgTreatment_BgSub,AvgTreatment_BgSub_NonZero,RawAnVPI,RawGFP,StdDMEM,StdDMSO,StdTreatment,
#    AllAvgTime,AllCompounds,AllRpts,AllFilesAnVPI,AllFilesGFP,AllreplID,AllSignalVariablesGFP,AllTimeID,AllTreated,argsJob,
#    CurrentReplID,DMEMdoses,DMSOdoses,IdxAnv,IdxColAnVPI,IdxColGFP,IdxPI,IdxRowAnVPI,IdxRowGFP,IdxSignal,repIndex,TimeTick,TxDoses)
# rm(list = ls(pattern = "counter")) 
# 
# save.image(file = "StressResponse_Variables_AllSR_Corrected_withSD_CellCt_BQcluster.RData")


# ================================== #
# === Normalisation and plotting === #
# ================================== #

rm(list=ls()); cat("\014")
# load("StressResponse_Variables_DDR_Corrected_withSD.RData")
# load("StressResponse_Variables_AllSR_Corrected_withSD.RData")
# load("StressResponse_Variables_AllSR_Corrected_withSD_CellCt.RData")
# load("StressResponse_Variables_AllSR_Corrected_withSD_CellCt_PIfilled_v2.RData")

load("StressResponse_Variables_AllSR_Corrected_withSD_CellCt_v4.RData")
# load("StressResponse_Variables_AllSR_Corrected_withSD_CellCt_v4_BkUp_DataMapped.RData")

# for (counter_rm in 12:13) { # Remove "DMEM no-TNF" experiment of A20 and ICAM1 from the list of MIDAS mapping
#   AllExtractedData[[counter_rm]]$AllCompounds <- AllExtractedData[[counter_rm]]$AllCompounds[-14]
#   AllExtractedData[[counter_rm]]$AllTreated <- AllExtractedData[[counter_rm]]$AllTreated[-14]
#   AllExtractedData[[counter_rm]]$AllDoses[[14]] <- NULL 
#   AllExtractedData[[counter_rm]]$AvgTreatment[[14]] <- NULL
#   AllExtractedData[[counter_rm]]$StdTreatment[[14]] <- NULL
#   AllExtractedData[[counter_rm]]$AvgTreatment_BgSub[[14]] <- NULL
#   AllExtractedData[[counter_rm]]$AvgTreatment_BgSub_NonZero[[14]] <- NULL
# }

# rm(list = ls(pattern = "counter")) 

# AllRpts <- c("BTG2","MDM2","P21","TP53")
AllRpts <- c("BTG2","MDM2","P21","TP53","HMOX1","NRF2","SRXN1","ATF4","BIP","CHOP","XBP1","A20","ICAM1")

AllExtractedDataMaxNorm <- list()
# AllExtractedSDMaxNorm <- list()

# 25.09.18 --- Extract data from AllExtractedData$AvgTreatment_BgSub instead of AllExtractedData$AvgTreatment_BgSub_NonZero
# To remove zero after normalisation and calculating to the mean of all replicates

# Normalise by cell count

# for (counter_rpt in 1:length(AllRpts)) {
#   
#   for (counter_tx in 1:length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub)) {
# 
#     print(paste0("Mapping ",counter_rpt,"/",length(AllRpts)," - Treatment - ",counter_tx,"/",length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub)))
#        
#     for (counter_repl in 1:(dim(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter_tx]]$signal)[3])) {
# 
#       AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter_tx]]$signal[,,counter_repl] <-
#         AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter_tx]]$signal[,,counter_repl] / AllExtractedData[[counter_rpt]]$AvgTreatment[[counter_tx]]$CellCtGFP[,,counter_repl]
# 
#       AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter_tx]]$AnV[,,counter_repl] <-
#         AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter_tx]]$AnV[,,counter_repl] / AllExtractedData[[counter_rpt]]$AvgTreatment[[counter_tx]]$CellCtAnVPI[,,counter_repl]
#       
#       AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter_tx]]$PI[,,counter_repl] <-
#         AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter_tx]]$PI[,,counter_repl] / AllExtractedData[[counter_rpt]]$AvgTreatment[[counter_tx]]$CellCtAnVPI[,,counter_repl]
#       
#     }
#      
#   }
#   
# }



for (counter_rpt in 1:length(AllRpts)) {
  # for (counter_rpt in 4:length(AllRpts)) {
  
  # Identify the strongest signal across all conditions & Normalise
  # MaxValueByCompound <- list("signal"=rep(NA,length(AllExtractedData[[counter_rpt]]$AllTreated)),
  #                            "AnV"=rep(NA,length(AllExtractedData[[counter_rpt]]$AllTreated)),
  #                            "PI"=rep(NA,length(AllExtractedData[[counter_rpt]]$AllTreated)))
  MaxValueByCompound <- list("signal"=matrix(NA,length(AllExtractedData[[counter_rpt]]$AllTreated),length(AllExtractedData[[counter_rpt]]$AllreplID)),
                             "AnV"=matrix(NA,length(AllExtractedData[[counter_rpt]]$AllTreated),length(AllExtractedData[[counter_rpt]]$AllreplID)),
                             "PI"=matrix(NA,length(AllExtractedData[[counter_rpt]]$AllTreated),length(AllExtractedData[[counter_rpt]]$AllreplID)))
  # for (counter in 1:length(AllExtractedData[[counter_rpt]]$AllTreated)) {
  #     MaxValueByCompound$signal[counter] <- max(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal,na.rm = T)
  #     MaxValueByCompound$AnV[counter] <- max(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV,na.rm = T)
  #     MaxValueByCompound$PI[counter] <- max(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI,na.rm = T)
  # }
  for (counter in 1:length(AllExtractedData[[counter_rpt]]$AllTreated)) {
    for (counter2 in 1:length(AllExtractedData[[counter_rpt]]$AllreplID)) {
      MaxValueByCompound$signal[counter,counter2] <- max(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,,counter2],na.rm = T)
      MaxValueByCompound$AnV[counter,counter2] <- max(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,,counter2],na.rm = T)
      MaxValueByCompound$PI[counter,counter2] <- max(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,,counter2],na.rm = T)
    }
  }
  
  # Missing data returned -Inf -> turn to NA
  MaxValueByCompound$signal[which(MaxValueByCompound$signal==-Inf)] <- NA
  MaxValueByCompound$AnV[which(MaxValueByCompound$AnV==-Inf)] <- NA
  MaxValueByCompound$PI[which(MaxValueByCompound$PI==-Inf)] <- NA
  
  # Plot raw intensity
  library(reshape2)
  library(ggplot2)
  library(RColorBrewer)
  
  MeltedData_Signal <- NULL; MeltedSD_Signal <- NULL
  MeltedData_AnV <- NULL; MeltedSD_AnV <- NULL
  MeltedData_PI <- NULL; MeltedSD_PI <- NULL
  
  for (counter in 1:length(AllExtractedData[[counter_rpt]]$AllTreated)) {
    # for (counter in 1:8) {
    print(paste0("Processing ",AllRpts[counter_rpt]," compound: ",counter,"/",length(AllExtractedData[[counter_rpt]]$AllTreated)))
    
    # === New === # 
    
    TempExtractedData <- AllExtractedData
    
    # Normalise dataset
    for (counter_DatSig in 1:length(AllExtractedData[[counter_rpt]]$AllreplID)) {
      # TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,,counter_DatSig] <- AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,,counter_DatSig]/max(MaxValueByCompound$signal[,counter_DatSig],na.rm = T)
      # TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,,counter_DatSig] <- AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,,counter_DatSig]/max(MaxValueByCompound$signal[,counter_DatSig],na.rm = T)
      # TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,,counter_DatSig] <- AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,,counter_DatSig]/max(MaxValueByCompound$signal[,counter_DatSig],na.rm = T)
      TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,,counter_DatSig] <- AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,,counter_DatSig]/max(MaxValueByCompound$signal[,counter_DatSig],na.rm = T)
      TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,,counter_DatSig] <- AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,,counter_DatSig]/max(MaxValueByCompound$AnV[,counter_DatSig],na.rm = T)
      TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,,counter_DatSig] <- AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,,counter_DatSig]/max(MaxValueByCompound$PI[,counter_DatSig],na.rm = T)
    }
    # MeanGFPSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal,c(1,2),mean,na.rm=T)
    # StdGFPSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal,c(1,2),sd,na.rm=T)
    # MeanAnVSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV,c(1,2),mean,na.rm=T)
    # StdAnVSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV,c(1,2),sd,na.rm=T)
    # MeanPISig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI,c(1,2),mean,na.rm=T)
    # StdPISig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI,c(1,2),sd,na.rm=T)
    MeanGFPSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal,c(1,2),mean,na.rm=T)
    StdGFPSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal,c(1,2),sd,na.rm=T)
    MeanAnVSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV,c(1,2),mean,na.rm=T)
    StdAnVSig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV,c(1,2),sd,na.rm=T)
    MeanPISig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI,c(1,2),mean,na.rm=T)
    StdPISig <- apply(TempExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI,c(1,2),sd,na.rm=T)
    
    # NormMeanGFPSig <- MeanGFPSig/max(MeanGFPSig,na.rm = T)
    # NormStdGFPSig <- StdGFPSig/max(MeanGFPSig,na.rm = T)
    # NormMeanAnVSig <- MeanAnVSig/max(MeanAnVSig,na.rm = T)
    # NormStdAnVSig <- StdAnVSig/max(MeanAnVSig,na.rm = T)
    # NormMeanPISig <- MeanPISig/max(MeanPISig,na.rm = T)
    # NormStdPISig <- StdPISig/max(MeanPISig,na.rm = T)
    
    # # Missing data returned -Inf -> turn to NA
    # NormMeanGFPSig[which(NormMeanGFPSig==-Inf)] <- NA
    # NormStdGFPSig[which(NormStdGFPSig==-Inf)] <- NA
    # NormMeanAnVSig[which(NormMeanAnVSig==-Inf)] <- NA
    # NormStdAnVSig[which(NormStdAnVSig==-Inf)] <- NA
    # NormMeanPISig[which(NormMeanPISig==-Inf)] <- NA
    # NormStdPISig[which(NormStdPISig==-Inf)] <- NA
    
    # =========== #
    
    for (counter2 in 1:length(AllExtractedData[[counter_rpt]]$AllAvgTime)) {
      # print(paste0("Processing TimePoint: ",counter2,"/",length(AllAvgTime)))
      
      # === Data === #
      
      # DataCurrentTime_Signal <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$name,length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,counter2])),
      #                                 rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,counter2])),
      #                                 melt(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,counter2]/max(MaxValueByCompound$signal)),
      #                                 as.character(1:length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,counter2])))
      
      # DataCurrentTime_Signal <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$name,length(rowMeans(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,counter2,]))),
      #                                rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],length(rowMeans(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,counter2,]))),
      DataCurrentTime_Signal <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$name,ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,counter2,]))),
                                      rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,counter2,]))),
                                      melt(MeanGFPSig[,counter2]),
                                      melt(StdGFPSig[,counter2]),
                                      #                                as.character(1:length(rowMeans(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,counter2,]))))
                                      as.character(1:ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$signal[,counter2,]))))
      
      # DataCurrentTime_AnV <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$name,length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,counter2])),
      #                              rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,counter2])),
      #                              melt(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,counter2]/max(MaxValueByCompound$AnV)),
      #                              as.character(1:length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,counter2])))
      
      # DataCurrentTime_AnV <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$name,length(rowMeans(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,counter2,]))),
      #                                 rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],length(rowMeans(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,counter2,]))),
      DataCurrentTime_AnV <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$name,ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,counter2,]))),
                                   rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,counter2,]))),
                                   melt(MeanAnVSig[,counter2]),
                                   melt(StdAnVSig[,counter2]),
                                   #                              as.character(1:length(rowMeans(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,counter2,]))))
                                   as.character(1:ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$AnV[,counter2,]))))
      
      # DataCurrentTime_PI <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$name,length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,counter2])),
      #                             rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,counter2])),
      #                             melt(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,counter2]/max(MaxValueByCompound$PI)),
      #                             as.character(1:length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,counter2])))
      
      # DataCurrentTime_PI <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$name,length(rowMeans(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,counter2,]))),
      #                              rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],length(rowMeans(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,counter2,]))),
      DataCurrentTime_PI <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$name,ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,counter2,]))),
                                  rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,counter2,]))),
                                  melt(MeanPISig[,counter2]),
                                  melt(StdPISig[,counter2]),
                                  #                             as.character(1:length(rowMeans(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,counter2,]))))
                                  as.character(1:ifelse(test = is.null(nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,counter2,])),yes = 1,no = nrow(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub[[counter]]$PI[,counter2,]))))
      
      MeltedData_Signal <- rbind(MeltedData_Signal,DataCurrentTime_Signal)
      MeltedData_AnV <- rbind(MeltedData_AnV,DataCurrentTime_AnV)
      MeltedData_PI <- rbind(MeltedData_PI,DataCurrentTime_PI)
      
      # === SD === #
      
      # SDCurrentTime_Signal <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$name,length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,counter2])),
      #                                 rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$signal[,counter2])),
      #                                 melt(AllExtractedData[[counter_rpt]]$StdTreatment[[counter]]$signal[,counter2]/max(MaxValueByCompound$signal)),
      #                                 as.character(1:length(AllExtractedData[[counter_rpt]]$StdTreatment[[counter]]$signal[,counter2])))
      # 
      # SDCurrentTime_AnV <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$name,length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,counter2])),
      #                              rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,counter2])),
      #                              melt(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$AnV[,counter2]/max(MaxValueByCompound$AnV)),
      #                              as.character(1:length(AllExtractedData[[counter_rpt]]$StdTreatment[[counter]]$AnV[,counter2])))
      # 
      # SDCurrentTime_PI <- cbind(rep(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$name,length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,counter2])),
      #                             rep(AllExtractedData[[counter_rpt]]$AllAvgTime[counter2],length(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,counter2])),
      #                             melt(AllExtractedData[[counter_rpt]]$AvgTreatment_BgSub_NonZero[[counter]]$PI[,counter2]/max(MaxValueByCompound$PI)),
      #                             as.character(1:length(AllExtractedData[[counter_rpt]]$StdTreatment[[counter]]$PI[,counter2])))
      # 
      # MeltedSD_Signal <- rbind(MeltedSD_Signal,SDCurrentTime_Signal)
      # MeltedSD_AnV <- rbind(MeltedSD_AnV,SDCurrentTime_AnV)
      # MeltedSD_PI <- rbind(MeltedSD_PI,SDCurrentTime_PI)
      
    }
  }
  
  # colnames(MeltedData_Signal) <- c("treatment","time","signal","doseNr"); colnames(MeltedSD_Signal) <- colnames(MeltedData_Signal)
  # colnames(MeltedData_AnV) <- c("treatment","time","signal","doseNr"); colnames(MeltedSD_AnV) <- colnames(MeltedData_AnV)
  # colnames(MeltedData_PI) <- c("treatment","time","signal","doseNr"); colnames(MeltedSD_PI) <- colnames(MeltedData_PI)
  
  colnames(MeltedData_Signal) <- c("treatment","time","signal","std","doseNr")
  colnames(MeltedData_AnV) <- c("treatment","time","signal","std","doseNr")
  colnames(MeltedData_PI) <- c("treatment","time","signal","std","doseNr")
  
  # 25.09.18 --- Introduce removing zero entries here
  
  MeltedData_Signal$signal[which(MeltedData_Signal$signal<0)] <- 0
  MeltedData_Signal$std[which(MeltedData_Signal$signal<0)] <- 0
  MeltedData_AnV$signal[which(MeltedData_AnV$signal<0)] <- 0
  MeltedData_AnV$std[which(MeltedData_AnV$signal<0)] <- 0
  MeltedData_PI$signal[which(MeltedData_PI$signal<0)] <- 0
  MeltedData_PI$std[which(MeltedData_PI$signal<0)] <- 0
  
  NormMeltedData_Signal <- MeltedData_Signal; NormMeltedData_Signal$signal <- NormMeltedData_Signal$signal/max(NormMeltedData_Signal$signal,na.rm = T)
  NormMeltedData_AnV <- MeltedData_AnV; NormMeltedData_AnV$signal <- NormMeltedData_AnV$signal/max(NormMeltedData_AnV$signal,na.rm = T)
  NormMeltedData_PI <- MeltedData_PI; NormMeltedData_PI$signal <- NormMeltedData_PI$signal/max(NormMeltedData_PI$signal,na.rm = T)
  
  # Plot profiles of signal, AnV and PI
  pdf(paste0(AllRpts[counter_rpt],"_Signal_Norm_BgSub_withSD_PerRepl_v4_NoNormCellCt_Rep2.pdf"))
  print(ggplot(data = NormMeltedData_Signal, aes(x=time,y=signal))
        + geom_pointrange(aes(ymin=signal-std, ymax=signal+std),size=0.1,alpha=0.2)
        # + geom_smooth(method="loess", se=FALSE)
        # + geom_ribbon(aes(ymin=signal-std, ymax=signal+std),size=0.2,alpha=0.2)
        + geom_point(aes(col=doseNr),size=1)
        + ylim(c(-1e-10,1)) 
        + facet_wrap( ~ treatment, scales="free")
        + theme(axis.text.x=element_blank())
        + ggtitle(paste0("Normalised Bg-subtracted ",AllRpts[counter_rpt]," signal"))
        + theme(legend.position = "right",legend.title=element_text(size=10),legend.text=element_text(size=8)))
  dev.off()
  
  pdf(paste0(AllRpts[counter_rpt],"_AnV_Norm_BgSub_withSD_PerRepl_v4_NoNormCellCt_Rep2.pdf"))
  print(ggplot(data = NormMeltedData_AnV, aes(x=time,y=signal))
        + geom_pointrange(aes(ymin=signal-std, ymax=signal+std),size=0.1,alpha=0.2)
        # + geom_smooth(method="loess", se=FALSE)
        # + geom_ribbon(aes(ymin=signal-std, ymax=signal+std),size=0.2,alpha=0.2)
        + geom_point(aes(col=doseNr),size=1)
        + ylim(c(-1e-10,1))
        + facet_wrap( ~ treatment, scales="free")
        + theme(axis.text.x=element_blank())
        + ggtitle(paste0("Normalised Bg-subtracted ",AllRpts[counter_rpt]," AnV"))
        + theme(legend.position = "right",legend.title=element_text(size=10),legend.text=element_text(size=8)))
  dev.off()
  
  pdf(paste0(AllRpts[counter_rpt],"_PI_Norm_BgSub_withSD_PerRepl_v4_NoNormCellCt_Rep2.pdf"))
  print(ggplot(data = NormMeltedData_PI, aes(x=time,y=signal))
        + geom_pointrange(aes(ymin=signal-std, ymax=signal+std),size=0.1,alpha=0.2)
        # + geom_smooth(method="loess", se=FALSE)
        # + geom_ribbon(aes(ymin=signal-std, ymax=signal+std),size=0.2,alpha=0.2)
        + geom_point(aes(col=doseNr),size=1)
        + ylim(c(-1e-10,1))
        + facet_wrap( ~ treatment, scales="free")
        + theme(axis.text.x=element_blank())
        + ggtitle(paste0("Normalised Bg-subtracted ",AllRpts[counter_rpt]," PI"))
        + theme(legend.position = "right",legend.title=element_text(size=10),legend.text=element_text(size=8)))
  dev.off()
  
  AllExtractedDataMaxNorm[[counter_rpt]] <- 
    list(Reporter = AllRpts[counter_rpt],
         MaxValueByCompound = MaxValueByCompound,
         MeanGFPSig = MeanGFPSig,
         StdGFPSig = StdGFPSig,
         MeanAnVSig = MeanAnVSig,
         StdAnVSig = StdAnVSig,
         MeanPISig = MeanPISig,
         StdPISig = StdPISig,
         MeltedData_Signal = MeltedData_Signal,
         MeltedData_AnV = MeltedData_AnV,
         MeltedData_PI = MeltedData_PI,
         NormMeltedData_Signal = NormMeltedData_Signal,
         NormMeltedData_AnV = NormMeltedData_AnV,
         NormMeltedData_PI = NormMeltedData_PI)
  
  # AllExtractedSDMaxNorm[[counter_rpt]] <- 
  #   list(Reporter = AllRpts[counter_rpt],
  #        MeltedSD_Signal = MeltedSD_Signal,
  #        MeltedSD_AnV = MeltedSD_AnV,
  #        MeltedSD_PI = MeltedSD_PI)
}

# save.image(file = "StressResponse_Variables_DDR_Corrected_withSD.RData")
# save.image(file = "StressResponse_Variables_AllSR_Corrected_withSD.RData")
# save.image(file = "StressResponse_Variables_AllSR_Corrected_withSD_CellCt.RData")
# save.image(file = "StressResponse_Variables_AllSR_Corrected_withSD_CellCt_PIfilled_v2.RData")
# save.image(file = "StressResponse_Variables_AllSR_Corrected_withSD_CellCt_BQcluster_4to11.RData")

save.image(file = "StressResponse_Variables_AllSR_Corrected_withSD_CellCt_v4.RData")
# save.image(file = "StressResponse_Variables_AllSR_Corrected_withSD_CellCt_v4_Rep2.RData")

# ==================================== #
# Preparation of MIDAS file 
# ==================================== #

# Case study 2: Take all time course data from all reps from the dose with highest response [NEED DISCUSSION!!]

# Compound_List <- AllTreated
# 
# # Remove space at the end of compound name (if any)
# for (counter in 1:length(Compound_List)) {
#   Compound_List[counter] <- if (substr(Compound_List[counter],nchar(Compound_List[counter]),nchar(Compound_List[counter]))==" ") {
#     substr(Compound_List[counter],1,nchar(Compound_List[counter])-1)
#   } else { Compound_List[counter] }
#   if (grepl(pattern = " ",x = Compound_List[counter],fixed = T)) {
#     Compound_List[counter] <- gsub(x = Compound_List[counter],pattern = " ",replacement = "_",fixed = T)
#   }
#   if (grepl(pattern = ".",x = Compound_List[counter],fixed = T)) {
#     Compound_List[counter] <- gsub(x = Compound_List[counter],pattern = ".",replacement = "",fixed = T)
#   }
# }

# write.table(x = Compound_List,file = "Compound_List.csv",sep="\n",quote=FALSE,row.names = FALSE,col.names = F)

# CellLine_List <- c("BTG2","MDM2","P21","TP53")

# Compound_TimeCourse_Data <- list()
# Compound_TimeCourse_Apoptosis <- list()
# Selected_Dose_MaxSignal <- matrix(NA,2,length(Compound_List))
# 
# for (counter_CL in 1:length(CellLine_List)) {
#   
#   for (counter in 1:length(Compound_List)) {
#     
#     SelectedDoseLine <- which(max(AllExtractedData[[counter_CL]]$AvgTreatment_BgSub_NonZero[[counter]]$signal)==AllExtractedData[[counter_CL]]$AvgTreatment_BgSub_NonZero[[counter]]$signal,arr.ind = T)[1]
#     Selected_Dose_MaxSignal[1,counter] <- rownames(AvgTreatment_BgSub_NonZero[[counter]]$signal)[SelectedDoseLine]
#     Selected_Dose_MaxSignal[2,counter] <- SelectedDoseLine
#     
#     SR_Compound_List_AllRep <- AvgTreatment_BgSub_NonZero[[counter]]$signal[SelectedDoseLine,]/max(MaxValueByCompound$signal)
#     Unique_TimeCourse <- AllAvgTime
#     SR_Prep_MIDAS_TimeCourse <- as.data.frame(matrix(NA,length(CellLine_List),length(Unique_TimeCourse)))
#     colnames(SR_Prep_MIDAS_TimeCourse) <- Unique_TimeCourse
#     rownames(SR_Prep_MIDAS_TimeCourse) <- CellLine_List
#     
#     # for (counter2 in 1:length(CellLine_List)) {
#     #   
#     #   Current_CellLine_Entries <- SR_Compound_List_AllRep
#     #   
#     #   for (counter3 in 1:nrow(Current_CellLine_Entries)) {
#     #     Col_TimePoint_Index <- which(Current_CellLine_Entries$timeAfterExposure[counter3]==Unique_TimeCourse)
#     #     SR_Prep_MIDAS_TimeCourse[counter2,Col_TimePoint_Index] <- Current_CellLine_Entries$value[counter3]  
#     #   }
#     # }
#     # SR_Prep_MIDAS_TimeCourse <- SR_Compound_List_AllRep
#     
#     
#     Compound_TimeCourse_Data[[counter]] <- SR_Prep_MIDAS_TimeCourse
#     # Compound_TimeCourse_Data[[counter]] <- SR_Prep_MIDAS_TimeCourse
#     
#     
#     # --- Isolated script for Apoptosis alone --- #
#     
#     SR_Compound_List_AllRep_Apoptosis <- AvgTreatment_BgSub_NonZero[[counter]]$AnV[SelectedDoseLine,]/max(MaxValueByCompound$AnV)
#     Unique_TimeCourse_Apoptosis <- AllAvgTime
#     SR_Prep_MIDAS_TimeCourse_Apoptosis <- as.data.frame(matrix(NA,length(CellLine_List),length(Unique_TimeCourse_Apoptosis)))
#     colnames(SR_Prep_MIDAS_TimeCourse_Apoptosis) <- Unique_TimeCourse_Apoptosis
#     rownames(SR_Prep_MIDAS_TimeCourse_Apoptosis) <- "Apoptosis"
#     
#     SR_Prep_MIDAS_TimeCourse_Apoptosis <- SR_Compound_List_AllRep_Apoptosis
#     Compound_TimeCourse_Apoptosis[[counter]] <- SR_Prep_MIDAS_TimeCourse_Apoptosis
#     
#     # --- Isolated script for Apoptosis alone --- #
#     
#   }
#   
#   
# }
# print(Selected_Dose_MaxSignal)

# ===== Setting up MIDAS file ===== #

rm(list=ls()); cat("\014")
# load("StressResponse_Variables_DDR_Corrected_withSD.RData")
# load("StressResponse_Variables_AllSR_Corrected_withSD.RData")
# load("StressResponse_Variables_AllSR_Corrected_withSD_CellCt.RData")
# load("StressResponse_Variables_AllSR_Corrected_withSD_CellCt_PIfilled_v2.RData")
load("StressResponse_Variables_AllSR_Corrected_withSD_CellCt_v4.RData")

# for (counter_rm in 12:13) { # Remove "DMEM no-TNF" experiment of A20 and ICAM1 from the list of MIDAS mapping
#   AllExtractedData[[counter_rm]]$AllCompounds <- AllExtractedData[[counter_rm]]$AllCompounds[-14]
#   AllExtractedData[[counter_rm]]$AllTreated <- AllExtractedData[[counter_rm]]$AllTreated[-14]
#   AllExtractedData[[counter_rm]]$AllDoses[[14]] <- NULL 
#   AllExtractedData[[counter_rm]]$AvgTreatment[[14]] <- NULL
#   AllExtractedData[[counter_rm]]$StdTreatment[[14]] <- NULL
#   AllExtractedData[[counter_rm]]$AvgTreatment_BgSub[[14]] <- NULL
#   AllExtractedData[[counter_rm]]$AvgTreatment_BgSub_NonZero[[14]] <- NULL
# }

rm(list = ls(pattern = "counter")) 

CellLine_List <- c("BTG2","MDM2","P21","TP53","HMOX1","NRF2","SRXN1","ATF4","BIP","CHOP","XBP1","A20","ICAM1")

# Loop set-up

dir.create("MIDAS_output")

SR_MIDAS_Comprehensive <- list()

for (counter_CL in 1:length(CellLine_List)) {
  # for (counter_CL in 4:length(CellLine_List)) {
  
  # !! TEMP extraction!! -> better to extract doses since the beginning
  # AllDoses <- list()
  # for (counter_rpt in 1:length(AllExtractedData[[counter_CL]]$AllTreated)){
  #   AllDoses[[counter_rpt]] <- as.numeric(substring(rownames(AllExtractedData[[counter_CL]]$AvgTreatment[[counter_rpt]]$signal),2))
  # }  
  
  AllDoses <- AllExtractedData[[counter_CL]]$AllDoses
  
  # Normalise all Doses for logic modelling
  AllNormDoses <- list()
  for (counter_doses in 1:length(AllDoses)) {
    
    # Take Log10 and do linear normalisation of doses
    DosesCurrent <- AllDoses[[counter_doses]]
    DosesCurrent <- c(DosesCurrent[1]*0.1,DosesCurrent) # Temporary add 10% of lowest dose for scaling
    Log10Doses <- log10(DosesCurrent)
    NormLog10Doses <- (Log10Doses - min(Log10Doses))/(max(Log10Doses)-min(Log10Doses))
    NormLog10Doses <- NormLog10Doses[-1] # Remove the added temporary dose for scaling
    AllNormDoses[[counter_doses]] <- NormLog10Doses
    
    # AllNormDoses[[counter_doses]] <- AllDoses[[counter_doses]]/max(AllDoses[[counter_doses]]) # Norm to max
  }
  
  # Extract compound list
  Compound_List <- AllExtractedData[[counter_CL]]$AllTreated
  
  # Remove space at the end of compound name (if any)
  for (counter in 1:length(Compound_List)) {
    Compound_List[counter] <- if (substr(Compound_List[counter],nchar(Compound_List[counter]),nchar(Compound_List[counter]))==" ") {
      substr(Compound_List[counter],1,nchar(Compound_List[counter])-1)
    } else { Compound_List[counter] }
    if (grepl(pattern = " ",x = Compound_List[counter],fixed = T)) {
      Compound_List[counter] <- gsub(x = Compound_List[counter],pattern = " ",replacement = "_",fixed = T)
    }
    if (grepl(pattern = ".",x = Compound_List[counter],fixed = T)) {
      Compound_List[counter] <- gsub(x = Compound_List[counter],pattern = ".",replacement = "",fixed = T)
    }
  }
  
  
  # Convert common gene names to (official) gene symbols
  CellLine_List <- toupper(CellLine_List)
  for (counter in 1:length(CellLine_List)) {
    if (CellLine_List[counter]=="P21"){
      CellLine_List[counter] <- "CDKN1A"}
    else if (CellLine_List[counter]=="NRF2"){
      CellLine_List[counter] <- "NFE2L2"}
    else if (CellLine_List[counter]=="BIP"){
      CellLine_List[counter] <- "HSPA5"}
    else if (CellLine_List[counter]=="CHOP"){
      CellLine_List[counter] <- "DDIT3"}
    else if (CellLine_List[counter]=="A20"){
      CellLine_List[counter] <- "TNFAIP3"}
    else if (CellLine_List[counter]=="IKBA"){
      CellLine_List[counter] <- "NFKBIA"}
  }
  
  
  #   } else if (CellLine_List_GeneSymbol[counter]=="CHOP"){
  #     CellLine_List_GeneSymbol[counter] <- "DDIT3"
  # }
  
  # CellLine_List_GeneSymbol <- "TP53"
  # CellLine_List_GeneSymbol <- c("TP53","Apoptosis")
  # CellLine_List_GeneSymbol <- CellLine_List
  
  SR_MIDAS_Compound <- list()
  
  for (counter in 1:length(Compound_List)) {
    
    print(paste0("Mapping MIDAS for ",CellLine_List[counter_CL]," - Compound: ",counter,"/",length(Compound_List)))
    
    NrCond      <- length(AllNormDoses[[counter]]) # Number of experimental condition (perturbation/doses)
    NrTimePoint <- length(AllExtractedData[[counter_CL]]$AllAvgTime) # Number of control condition (all time course)
    NrExpSet    <- 1 # Number of experimental set (future compatible with multiCellNOpt)
    NrInput     <- 1 # Number of input (now put one at-a-time for drug vs no-drug)
    NrReadOut   <- length(CellLine_List)+1 # length(CellLine_List)
    # NrReadOut   <- 2 # length(CellLine_List) + Apoptosis
    
    FinalTimeCourse <- AllExtractedData[[counter_CL]]$AllAvgTime
    
    
    SR_MIDAS_TimeCourse <- as.data.frame(matrix(NA,NrCond*NrTimePoint,NrExpSet+NrInput+(NrReadOut*2)))
    ColNames_MIDAS_TimeCourse <- rep(NA,dim(SR_MIDAS_TimeCourse)[2])
    ColNames_MIDAS_TimeCourse[1] <- "TR:mock:CellLine"
    ColNames_MIDAS_TimeCourse[2] <- "TR:Drug"
    
    CellLine_List_Apop <- c(CellLine_List,"Apoptosis")
    
    for (counter2 in 1:NrReadOut) {
      ColNames_MIDAS_TimeCourse[2+counter2] <- paste("DA:",CellLine_List_Apop[counter2],sep="")
      ColNames_MIDAS_TimeCourse[2+NrReadOut+counter2] <- paste("DV:",CellLine_List_Apop[counter2],sep="")
    }  
    
    
    
    colnames(SR_MIDAS_TimeCourse) <- ColNames_MIDAS_TimeCourse
    
    # Current_SR_TimeCourse <- Compound_TimeCourse_Data[[counter]]
    # Current_SR_TimeCourse <- rbind(Compound_TimeCourse_Data[[counter]],Compound_TimeCourse_Apoptosis[[counter]]) # with Apoptosis
    
    SR_MIDAS_TimeCourse[,1] <- 1 # All for each cell line
    # SR_MIDAS_TimeCourse[,2] <- 1 # All condition with drugs
    
    # for (counter3 in 1:ncol(Current_SR_TimeCourse)) {
    # for (counter3 in 1:length(Current_SR_TimeCourse)) {
    
    for (counter_do in 1:NrCond) {
      for (counter3 in 1:length(AllExtractedData[[counter_CL]]$AllAvgTime)) {
        SR_MIDAS_TimeCourse[((counter_do-1)*length(AllExtractedData[[counter_CL]]$AllAvgTime))+counter3,2] <- AllNormDoses[[counter]][counter_do]
        SR_MIDAS_TimeCourse[((counter_do-1)*length(AllExtractedData[[counter_CL]]$AllAvgTime))+counter3,3:(2+NrReadOut)] <- rep(FinalTimeCourse[counter3],NrReadOut)
        
        # if (length(SR_MIDAS_TimeCourse[counter3,(2+NrReadOut+1):dim(SR_MIDAS_TimeCourse)[2]])==1) {
        #   SR_MIDAS_TimeCourse[counter3,(2+NrReadOut+1):dim(SR_MIDAS_TimeCourse)[2]] <- Current_SR_TimeCourse[counter3]
        # } else {
        #   SR_MIDAS_TimeCourse[counter3,(2+NrReadOut+1):dim(SR_MIDAS_TimeCourse)[2]] <- Current_SR_TimeCourse[,counter3]
        # }
        
        if (length(AllExtractedDataMaxNorm[[counter_CL]]$MeltedData_Signal$signal[which(
          AllExtractedDataMaxNorm[[counter_CL]]$MeltedData_Signal$treatment==Compound_List[counter] &
          AllExtractedDataMaxNorm[[counter_CL]]$MeltedData_Signal$doseNr==counter_do &
          AllExtractedDataMaxNorm[[counter_CL]]$MeltedData_Signal$time==FinalTimeCourse[counter3])])>0) {
          
          SR_MIDAS_TimeCourse[((counter_do-1)*length(AllExtractedData[[counter_CL]]$AllAvgTime))+counter3,
                              (2+NrReadOut)+counter_CL] <-
            AllExtractedDataMaxNorm[[counter_CL]]$MeltedData_Signal$signal[which(
              AllExtractedDataMaxNorm[[counter_CL]]$MeltedData_Signal$treatment==Compound_List[counter] &
                AllExtractedDataMaxNorm[[counter_CL]]$MeltedData_Signal$doseNr==counter_do &
                AllExtractedDataMaxNorm[[counter_CL]]$MeltedData_Signal$time==FinalTimeCourse[counter3])]
          
          # Map Apoptosis

          SR_MIDAS_TimeCourse[((counter_do-1)*length(AllExtractedData[[counter_CL]]$AllAvgTime))+counter3,
                              (2+NrReadOut)+length(CellLine_List)+1] <-
            
            # Check first if the value is NA?
            ifelse(test = is.na(SR_MIDAS_TimeCourse[((counter_do-1)*length(AllExtractedData[[counter_CL]]$AllAvgTime))+counter3,
                                          (2+NrReadOut)+length(CellLine_List)+1]),
                   yes = AllExtractedDataMaxNorm[[counter_CL]]$MeltedData_AnV$signal[which(
                     AllExtractedDataMaxNorm[[counter_CL]]$MeltedData_AnV$treatment==Compound_List[counter] &
                       AllExtractedDataMaxNorm[[counter_CL]]$MeltedData_AnV$doseNr==counter_do &
                       AllExtractedDataMaxNorm[[counter_CL]]$MeltedData_AnV$time==FinalTimeCourse[counter3])],
                   no = (SR_MIDAS_TimeCourse[((counter_do-1)*length(AllExtractedData[[counter_CL]]$AllAvgTime))+counter3,
                                            (2+NrReadOut)+length(CellLine_List)+1] +
                     AllExtractedDataMaxNorm[[counter_CL]]$MeltedData_AnV$signal[which(
                       AllExtractedDataMaxNorm[[counter_CL]]$MeltedData_AnV$treatment==Compound_List[counter] &
                         AllExtractedDataMaxNorm[[counter_CL]]$MeltedData_AnV$doseNr==counter_do &
                         AllExtractedDataMaxNorm[[counter_CL]]$MeltedData_AnV$time==FinalTimeCourse[counter3])])
                  )
            
            
          
          
        }
      }
    }
    
    SR_MIDAS_Compound[[counter]] <- SR_MIDAS_TimeCourse
    
  }
  
  SR_MIDAS_Comprehensive[[counter_CL]] <- SR_MIDAS_Compound
  
  # Control_Condition <- matrix(0,1,ncol(SR_MIDAS_TimeCourse))
  # Control_Condition[,1] <- 1 # Only one cell line
  # Control_Condition[2,2] <- 1 # Basal condition at zero
  
  
  
  # if (counter==1) {
  #   SR_multiCellNOpt <- SR_MIDAS_TimeCourse_Write2File
  # } else {
  #   SR_MIDAS_TimeCourse_Write2File[,1] <- counter
  #   SR_multiCellNOpt <- rbind(SR_multiCellNOpt,SR_MIDAS_TimeCourse_Write2File)
  # }
  
}

# Combind all the readouts from all reporters
for (counter in 1:length(Compound_List)) {
  # for (counter in 4:length(Compound_List)) {
  SR_MIDAS_ToWrite <- NULL
  for (counter2 in 1:length(SR_MIDAS_Comprehensive)) {
    SR_MIDAS_ToWrite <- rbind(SR_MIDAS_ToWrite,SR_MIDAS_Comprehensive[[counter2]][[counter]])
  }
  
  # Add control condition
  Control_Condition <- rep(0,ncol(SR_MIDAS_ToWrite))
  Control_Condition[1] <- 1 # Only one cell line
  SR_MIDAS_Write2File <- rbind(Control_Condition,as.matrix(SR_MIDAS_ToWrite))
  
  # write.csv(x = SR_MIDAS_Write2File,file = paste0("MIDAS_output/SR_MIDAS_", Compound_List[counter] ,".csv",sep = ""),quote=FALSE,row.names = FALSE)
  # write.csv(x = SR_MIDAS_Write2File,file = paste0("MIDAS_output/SR_MIDAS_", Compound_List[counter] ,"_DDR.csv",sep = ""),quote=FALSE,row.names = FALSE)
  # write.csv(x = SR_MIDAS_Write2File,file = paste0("MIDAS_output/SR_MIDAS_", Compound_List[counter] ,"_AllSR.csv",sep = ""),quote=FALSE,row.names = FALSE)
  # write.csv(x = SR_MIDAS_Write2File,file = paste0("MIDAS_output/SR_MIDAS_", Compound_List[counter] ,"_AllSR_Corrected.csv",sep = ""),quote=FALSE,row.names = FALSE)
  # write.csv(x = SR_MIDAS_Write2File,file = paste0("MIDAS_output/SR_MIDAS_", Compound_List[counter] ,"_AllSR_Corrected_NormPerRepl_NoNormCellCt.csv",sep = ""),quote=FALSE,row.names = FALSE)
  # write.csv(x = SR_MIDAS_Write2File,file = paste0("MIDAS_output/SR_MIDAS_", Compound_List[counter] ,"_AllSR_Corrected_NormPerRepl_NoNormCellCt_v4.csv",sep = ""),quote=FALSE,row.names = FALSE)
  # write.csv(x = SR_MIDAS_Write2File,file = paste0("MIDAS_output/SR_MIDAS_", Compound_List[counter] ,"_AllSR_Corrected_NormPerRepl_NoNormCellCt_v4_NormLog10Doses.csv",sep = ""),quote=FALSE,row.names = FALSE)
  write.csv(x = SR_MIDAS_Write2File,file = paste0("MIDAS_output/SR_MIDAS_", Compound_List[counter] ,"_AllSR_Corrected_NormPerRepl_NoNormCellCt_v5_NormLog10Doses_Apoptosis.csv",sep = ""),quote=FALSE,row.names = FALSE)
  # write.csv(x = SR_MIDAS_Write2File,file = paste0("MIDAS_output/SR_MIDAS_", Compound_List[counter] ,"_withApoptosis.csv",sep = ""),quote=FALSE,row.names = FALSE)
  
}

# --- End of the script --- #
