# Fitting and parameter analyses for LogicODE models of GFP dataset

# clear workspace
rm(list=ls());cat("\014");if(length(dev.list())>0){dev.off()}

# Load compound names
Cpd <- t(read.table(file = "Compound_List.csv",header = F,sep = "\n",stringsAsFactors = F))

# FileNameTag - M,A,C
# DateTag <- "291118"
# DateTag <- c("031218","031218")
# DateTag <- c("291118","291118","031218","031218")
# DateTag <- c("291118","291118","031218","051218")
# DateTag <- c("291118","031218")
# DateTag <- "071218"
DateTag <- "101218"
# DateTag <- "030119"
# ModelTag <- c('M','MA')
# ModelTag <- c('MC','MAC')
# ModelTag <- c('M','MA','MC','MAC')
# ModelTag <- c('M','MC')
ModelTag <- 'MCP'
# ModelTag <- 'MCP'
# ModelSIF <- c("Results_SingleModels_SR_LogicODE_NoCrosstalk_minimal.sif",
#               "Results_SingleModels_SR_LogicODE_NoCrosstalk_minimal_withAbsorbNode.sif")
# ModelSIF <- c("Results_SingleModels_SR_LogicODE_withCrosstalk_minimal.sif",
#               "Results_SingleModels_SR_LogicODE_withCrosstalk_minimal_withAbsorbNode.sif")
# ModelSIF <- c("Results_SingleModels_SR_LogicODE_NoCrosstalk_minimal.sif",
#               "Results_SingleModels_SR_LogicODE_NoCrosstalk_minimal_withAbsorbNode.sif",
#               "Results_SingleModels_SR_LogicODE_withCrosstalk_minimal.sif",
#               "Results_SingleModels_SR_LogicODE_withCrosstalk_minimal_withAbsorbNode.sif")
# ModelSIF <- c("Results_SingleModels_SR_LogicODE_NoCrosstalk_minimal.sif",
#               "Results_SingleModels_SR_LogicODE_withCrosstalk_minimal.sif")
# ModelSIF <- "Results_SingleModels_SR_LogicODE_withCrosstalk_minimal_ApoptosisAND.sif"
ModelSIF <- "Results_SingleModels_SR_LogicODE_withCrosstalk_minimal_ApoptosisANDOR.sif"
# RunTag <- c("3","3")
# RunTag <- c("4","4")
# RunTag <- c("3","3","4","4")
# RunTag <- c("3","3","4","5")
# RunTag <- c("3","4")
# RunTag <- 6
RunTag <- 7
# RunTag <- 8
# TimeTag <- c("1800","1800","1800","3600")
# TimeTag <- c("1800","1800")
# TimeTag <- 3600
TimeTag <- 3600
InRoundRep <- 3
# InRoundRep <- 5

# Extract models from individual compounds 
# AllModels_M <- list()
# AllModels_MA <- list()
AllModels <- list()
CombinedFitCost <- list()
NoResultCpd <- NULL

for (counter_mod in 1:length(ModelTag)) {
  
  for (counter in 1:length(Cpd)) {
    print(paste0("Mapping Model ",ModelTag[counter_mod]," Compound: ",counter,"/",length(Cpd)))
    if (file.exists(paste0("Results_GFP_",DateTag[counter_mod],"/",ModelSIF[counter_mod],"_InRoundRep",InRoundRep[counter_mod],"_Time",TimeTag[counter_mod],"_Run",RunTag[counter_mod],"_",Cpd[counter],".RData"))) {
      if (counter==13) {counter=14} else if (counter==14) {counter=13} # Fix wrong indices of DMSO and DMEM in result files
      load(paste0("Results_GFP_",DateTag[counter_mod],"/",ModelSIF[counter_mod],"_InRoundRep",InRoundRep[counter_mod],"_Time",TimeTag[counter_mod],"_Run",RunTag[counter_mod],"_",Cpd[counter],".RData"))
      if (counter==13) {counter=14} else if (counter==14) {counter=13} # Return the indices
      AllModels[[counter]] <- All_SR_Models_Multi  
      rm(All_SR_Models_Multi)
    } else {
      print(paste0("The compound ",Cpd[counter]," has no modelling result!"))
      AllModels[[counter]] <- NULL
      NoResultCpd <- c(NoResultCpd,counter)
    }
  }
  
  # if (!is.null(NoResultCpd)) {
  #   Cpd <- Cpd[-NoResultCpd]
  # }
  
  print("==============================")
  
  # ============================= # 
  # Extract fitting cost and plot
  # ============================= # 
  
  AllFitCost <- matrix(NA,length(Cpd),4)
  
  AllIdxNonOutlier <- list()
  
  library(outliers)
  
  for (counter in 1:length(Cpd)) {
    
    if (!is.null(AllModels[[counter]])) {
      
      # Identify best cost
      TempCost <- rep(NA,length(AllModels[[counter]]))
      for (counter2 in 1:length(AllModels[[counter]])) {
        TempCost[counter2] <- AllModels[[counter]][[counter2]][[counter]]$fitResults[[1]]$fbest
      }
      IdxBestCost <- which(TempCost==min(TempCost))
      
      # Get non-outlier indices
      if (length(which(!outlier(TempCost,logical = T)))>0) {
        IdxNonOutlier <- which(!outlier(TempCost,logical = T)) 
      } else {
        IdxNonOutlier <- 1:length(TempCost)
      }
      
      # If outlier is indeed the better one(s) comparing to the mean cost -> take only the best outlier one(s)
      if (length(IdxNonOutlier)<length(1:length(AllModels[[counter]]))) {
        if (mean(TempCost[IdxNonOutlier]) > mean(TempCost[setdiff(1:length(AllModels[[counter]]),IdxNonOutlier)])) {
          IdxNonOutlier <- setdiff(1:length(AllModels[[counter]]),IdxNonOutlier)
        }
      }
      AllIdxNonOutlier[[counter]] <- IdxNonOutlier
      
      # remove outlier (or not?)
      if (length(rm.outlier(TempCost))>0) {
        TempCost <- TempCost[IdxNonOutlier]
      }
      
      AllFitCost[counter,1] <- mean(TempCost,na.rm = T)
      AllFitCost[counter,2] <- sd(TempCost,na.rm = T)
      AllFitCost[counter,3] <- paste(IdxBestCost,collapse = "_")
      AllFitCost[counter,4] <- paste(IdxNonOutlier,collapse = "_")
      
    }
  }
  
  AllFitCost_Plot <- data.frame(matrix(data = NA,nrow = length(Cpd),ncol = 5))
  AllFitCost_Plot[,1] <- t(Cpd)
  AllFitCost_Plot[,2] <- as.numeric(AllFitCost[,1])
  AllFitCost_Plot[,3] <- as.numeric(AllFitCost[,2])
  AllFitCost_Plot[,4] <- AllFitCost[,3]
  AllFitCost_Plot[,5] <- AllFitCost[,4]
  colnames(AllFitCost_Plot) <- c("Compound","MeanCost","SDCost","BestRun","GoodRun")
  
  # Convert NA from SD value of single ran to be zero
  AllFitCost_Plot$SDCost[which(is.na(AllFitCost_Plot$SDCost))] <- 0
  
  # Write all Fitcosts per model
  write.table(x = AllFitCost_Plot,file = paste0("AllFitCost_Table_",ModelTag[counter_mod],"_",DateTag[counter_mod],".tsv"),quote = F,sep = "\t")
  
  CombinedFitCost[[counter_mod]] <- AllFitCost_Plot
  
  library(reshape2)
  Melted_AllFitCost_Mean <- melt(AllFitCost_Plot[,c(1,2)],id.var="Compound")
  colnames(Melted_AllFitCost_Mean) <- c("Compound","Model","FitCost")
  Melted_AllFitCost_SD <- melt(AllFitCost_Plot[,c(1,3)],id.var="Compound")
  colnames(Melted_AllFitCost_SD) <- c("Compound","SD","FitCost")
  
  library(ggplot2)
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  col_gg <- gg_color_hue(1)
  
  # pdf(paste0("PlotFitCost_All_withOutlier_",ModelTag[counter_mod],"_",DateTag[counter_mod],".pdf"))
  pdf(paste0("PlotFitCost_All_withoutOutlier_",ModelTag[counter_mod],"_",DateTag[counter_mod],".pdf"))
  print(ggplot() +
          # geom_point(data = Melted_AllFitCost_Mean,aes(x=Compound,y=FitCost,col=Model),size=3,alpha=0.7) +
          geom_point(data = Melted_AllFitCost_Mean,aes(x=Compound,y=FitCost,col=Model),size=3,alpha=0.7) +
          geom_errorbar(data = AllFitCost_Plot,aes(x=Compound,ymin=MeanCost-SDCost,ymax=MeanCost+SDCost),width=0.2,alpha=1,col=col_gg[1]) +
          # geom_errorbar(data = AllFitCost,aes(x=Compound,ymin=MinAbs-SD_MinAbs,ymax=MinAbs+SD_MinAbs),width=0.2,alpha=1,col=col_gg[2]) +
          # labs(title = paste0("Fitting Cost - ",ModelTag[counter_mod]), subtitle ="All data (with outliers)") +
          labs(title = paste0("Fitting Cost - ",ModelTag[counter_mod]), subtitle ="Best data (without outliers)") +
          xlab("Compounds") +
          ylab("FitCost") +
          theme_light()+
          
          theme(text = element_text(size=15),
                axis.text.x = element_text(angle=90, vjust=1))
  )
  dev.off()
  
  
  # ================================= # 
  # Extract parameter values and plot # 
  # ================================= # 
  
  ParamName <- AllModels[[1]][[1]][[1]]$ode_parameters$parNames
  ParamValue <- AllModels[[1]][[1]][[1]]$ode_parameters$parValues
  NrParam <- length(ParamName)
  Idx_k <- which(grepl(x = ParamName,pattern = "_k_",fixed = T))
  Idx_tau <- which(grepl(x = ParamName,pattern = "tau_",fixed = T))
  
  # ===========================================
  # k parameter
  k_all <- array(NA,dim = c(length(AllModels[[1]]),length(Idx_k),length(Cpd)))
  
  for (counter in 1:length(AllModels[[1]])) {
    for (counter2 in 1:length(Cpd)) {
      if (counter %in% AllIdxNonOutlier[[counter2]]) { # get only good ones
        k_all[counter, ,counter2] <- AllModels[[counter2]][[counter]][[counter2]]$ode_parameters$parValues[Idx_k]
      }
    }
  }
  
  k_all[k_all<1e-6] <- 1e-6 # lower signal to 1e-6
  # k_all[k_all<1e-8] <- 1e-8 # lower signal to 1e-8
  
  k_all_log10 <- log10(k_all)
  # k_all_log10[k_all_log10==-Inf] <- 0
  
  k_mean <- colMeans(k_all_log10,na.rm = T)
  k_mean <- t(k_mean)
  k_sd <- NULL
  for (counter in 1:length(Cpd)) {
    k_sd_current <- apply(k_all[ , , counter],2,sd,na.rm=T)
    k_sd <- rbind(k_sd,k_sd_current)
  }
  rownames(k_sd) <- NULL
  k_sd <- round(k_sd,digits = 2)
  
  colnames(k_mean) <- ParamName[Idx_k]
  rownames(k_mean) <- Cpd
  colnames(k_sd) <- ParamName[Idx_k]
  rownames(k_sd) <- Cpd
  
  # Remove rows with no data
  IdxNaN <- which(is.nan(rowSums(k_mean)))
  if (length(IdxNaN)>0) {
    k_mean <- k_mean[-IdxNaN,]
    k_sd <- k_sd[-IdxNaN,]
  }
  
  # library(pheatmap)
  # pheatmap(k_mean,clustering_distance_cols = "correlation",clustering_distance_rows = "correlation")
  # pheatmap(k_mean,clustering_distance_cols = "correlation")
  # pheatmap(k_mean,clustering_distance_rows = "correlation")
  
  library(gplots)
  library(RColorBrewer)
  #MyColours <- brewer.pal(11,"RdYlGn")
  MyColours_k <- colorRampPalette(brewer.pal(9,"Blues"))(10)
  MyColours_k[1] <- "#808080"
  # MyColours_k[60:100] <- MyColours_k[100]
  par(oma = c(6, 0, 0, 0))
  pdf(paste0("SR_Allk_Log10_WithSD_",ModelTag[counter_mod],"_",DateTag[counter_mod],".pdf"),height = 7)
  heatmap.2(x = k_mean,cellnote = k_sd,notecex = 0.2, notecol="black",
            labRow = as.vector(rownames(k_mean)),labCol = colnames(k_mean),
            cexRow = 0.7,cexCol = 0.5,main = paste0("All 'k' parameters - ",ModelTag[counter_mod]),
            col=MyColours_k,tracecol=NA)
  dev.off()
  
  
  # ===========================================
  # tau parameter
  tau_all <- array(NA,dim = c(length(AllModels[[1]]),length(Idx_tau),length(Cpd)))
  
  for (counter in 1:length(AllModels[[1]])) {
    for (counter2 in 1:length(Cpd)) {
      if (counter %in% AllIdxNonOutlier[[counter2]]) { # get only good ones
        tau_all[counter, ,counter2] <- AllModels[[counter2]][[counter]][[counter2]]$ode_parameters$parValues[Idx_tau]
      }
    }
  }
  
  tau_all[tau_all<1e-6] <- 1e-6 # lower signal to 1e-8
  # tau_all[tau_all<1e-8] <- 1e-8 # lower signal to 1e-8
  
  tau_all_log10 <- log10(tau_all)
  # tau_all_log10[tau_all_log10==-Inf] <- 0
  # tau_all_log10 <- tau_all
  
  tau_mean <- colMeans(tau_all_log10,na.rm = T)
  tau_mean <- t(tau_mean)
  tau_sd <- NULL
  for (counter in 1:length(Cpd)) {
    tau_sd_current <- apply(tau_all[ , , counter],2,sd,na.rm=T)
    tau_sd <- rbind(tau_sd,tau_sd_current)
  }
  rownames(tau_sd) <- NULL
  tau_sd <- round(tau_sd,digits = 2)
  
  colnames(tau_mean) <- ParamName[Idx_tau]
  rownames(tau_mean) <- Cpd
  colnames(tau_sd) <- ParamName[Idx_tau]
  rownames(tau_sd) <- Cpd
  
  # Remove rows with no data
  IdxNaNtau <- which(is.nan(rowSums(tau_mean)))
  if (length(IdxNaNtau)>0) {
    tau_mean <- tau_mean[-IdxNaNtau,]
    tau_sd <- tau_sd[-IdxNaNtau,]
  }
  
  # library(pheatmap)
  # pheatmap(tau_mean,clustering_distance_cols = "correlation",clustering_distance_rows = "correlation")
  # pheatmap(tau_mean,clustering_distance_cols = "correlation")
  # pheatmap(tau_mean,clustering_distance_rows = "correlation")
  
  library(gplots)
  library(RColorBrewer)
  #MyColours <- brewer.pal(11,"RdYlGn")
  MyColours_tau <- colorRampPalette(brewer.pal(9,"Greens"))(10)
  MyColours_tau[1] <- "#808080"
  par(oma = c(6, 0, 0, 0))
  pdf(paste0("SR_Alltau_Log10_WithSD_",ModelTag[counter_mod],"_",DateTag[counter_mod],".pdf"),height = 7)
  heatmap.2(x = tau_mean,cellnote = tau_sd,notecex = 0.2, notecol="black",
            labRow = as.vector(rownames(tau_mean)),labCol = colnames(tau_mean),
            cexRow = 0.7,cexCol = 0.5,main = paste0("All 'tau' parameters - ",ModelTag[counter_mod]),
            col=MyColours_tau,tracecol=NA)
  dev.off()
}

# ===================================
# Extract k and tau parameter for machine learning tasks in Weka

# install.packages('foreign')
library(foreign)

k_tau_mean <- cbind(k_mean,tau_mean)
k_tau_mean_class <- cbind(k_tau_mean,
                          c("DILI","NegCt","DILI","UPR","NegCt","DIKI","OSR","DIKI",
                            "DILI","DIKI","OSR","DILI","Solvent","Solvent","DDR","NegCt",
                            "Solvent","DILI","NegCt","DIKI","NegCt","NegCt","DILI","DILI",
                            "DILI","DILI","DILI","DILI","UPR","DILI","UPR","NegCt","DILI","NegCt"))
colnames(k_tau_mean_class)[ncol(k_tau_mean_class)] <- "Class"
k_tau_mean_class <- as.data.frame(k_tau_mean_class)
write.arff(x = k_tau_mean_class,file = "GFP_to_ML.arff")

k_tau_mean_class_compact <- cbind(k_tau_mean,
                                  c("DILI","NegCt","DILI","PosCt","NegCt","DIKI","PosCt","DIKI",
                                    "DILI","DIKI","PosCt","DILI","NegCt","NegCt","PosCt","NegCt",
                                    "NegCt","DILI","NegCt","DIKI","NegCt","NegCt","DILI","DILI",
                                    "DILI","DILI","DILI","DILI","PosCt","DILI","PosCt","NegCt","DILI","NegCt"))
colnames(k_tau_mean_class_compact)[ncol(k_tau_mean_class_compact)] <- "Class"
k_tau_mean_class_compact <- as.data.frame(k_tau_mean_class_compact)
write.arff(x = k_tau_mean_class_compact,file = "GFP_to_ML_compact.arff")

# sub-select only DILI and DIKI
k_tau_mean_class_DILI_DIKI <- k_tau_mean_class[c(which(k_tau_mean_class$Class=="DILI"),which(k_tau_mean_class$Class=="DIKI")),]
k_tau_mean_class_DILI_DIKI <- as.data.frame(k_tau_mean_class_DILI_DIKI)
write.arff(x = k_tau_mean_class_DILI_DIKI,file = "GFP_to_ML_only_DILI_DIKI.arff")


# ===================================
# Plot comparison of all FitCosts from all (best) model variants

CombinedFitCost_Plot <- data.frame(matrix(data = NA,nrow = length(Cpd),ncol = 1+length(CombinedFitCost)*2))
CombinedFitCost_Plot[,1] <- t(Cpd)
colnames(CombinedFitCost_Plot)[1] <- "Compound"

for (counter_mod in 1:length(CombinedFitCost)) {
  CombinedFitCost_Plot[,1+(counter_mod*2-1)] <- CombinedFitCost[[counter_mod]][,2]
  CombinedFitCost_Plot[,1+(counter_mod*2-1)+1] <- CombinedFitCost[[counter_mod]][,3]
  colnames(CombinedFitCost_Plot)[1+(counter_mod*2-1)] <- ModelTag[counter_mod]
  colnames(CombinedFitCost_Plot)[1+(counter_mod*2-1)+1] <- paste0(ModelTag[counter_mod],"_sd")
}

library(reshape2)
Melted_CombinedFitCost_Mean <- melt(CombinedFitCost_Plot[,c(1,2*(1:length(CombinedFitCost)))],id.var="Compound")
colnames(Melted_CombinedFitCost_Mean) <- c("Compound","Model","FitCost")
Melted_CombinedFitCost_SD <- melt(CombinedFitCost_Plot[,c(1,2*(1:length(CombinedFitCost))+1)],id.var="Compound")
colnames(Melted_CombinedFitCost_SD) <- c("Compound","SD","FitCost")

library(ggplot2)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
col_gg <- gg_color_hue(length(CombinedFitCost))

# pdf(paste0("PlotFitCost_All_withOutlier_",ModelTag[counter_mod],"_",DateTag[counter_mod],".pdf"))
# pdf(paste0("CombinedFitCost_All_withoutOutlier_",paste(ModelTag, collapse="_"),".pdf"))
pdf(paste0("CombinedFitCost_All_withoutOutlier_",paste(ModelTag, collapse="_"),"_LongMAC.pdf"))
print(ggplot() +
        # geom_point(data = Melted_AllFitCost_Mean,aes(x=Compound,y=FitCost,col=Model),size=3,alpha=0.7) +
        geom_point(data = Melted_CombinedFitCost_Mean,aes(x=Compound,y=FitCost,col=Model),size=3,alpha=0.7) +
        geom_errorbar(data = CombinedFitCost_Plot,aes(x=Compound,ymin=M-M_sd,ymax=M+M_sd),width=0.2,alpha=1,col=col_gg[1]) +
        # geom_errorbar(data = CombinedFitCost_Plot,aes(x=Compound,ymin=MA-MA_sd,ymax=MA+MA_sd),width=0.2,alpha=1,col=col_gg[2]) +
        geom_errorbar(data = CombinedFitCost_Plot,aes(x=Compound,ymin=MC-MC_sd,ymax=MC+MC_sd),width=0.2,alpha=1,col=col_gg[3]) +
        # geom_errorbar(data = CombinedFitCost_Plot,aes(x=Compound,ymin=MAC-MAC_sd,ymax=MAC+MAC_sd),width=0.2,alpha=1,col=col_gg[4]) +
        # labs(title = paste0("Fitting Cost - ",ModelTag[counter_mod]), subtitle ="All data (with outliers)") +
        labs(title = paste0("Combined Fitting Cost - ",paste(ModelTag, collapse="_")), subtitle ="Best data (without outliers)") +
        xlab("Compounds") +
        ylab("FitCost") +
        theme_light()+
        
        theme(text = element_text(size=15),
              axis.text.x = element_text(angle=90, vjust=1))
)
dev.off()


# --- End of the script --- # 
