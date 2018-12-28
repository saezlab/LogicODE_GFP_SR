# Crosstalk Identification with OminpathR

# clear workspace
rm(list=ls());cat("\014");if(length(dev.list())>0){dev.off()}

options(stringsAsFactors = F)

# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("saezlab/omnipathR")  

library(OmnipathR)

Nodes_MinModel <- c("KEAP1","NFE2L2","HMOX1","SRXN1","TP53","MDM2","CDKN1A","BTG2","NFKB1","TNFAIP3","ICAM1","HSPA5","ATF4","DDIT3","XBP1","EIF2A")
Nodes_FullModel <- c("KEAP1","NFE2L2","HMOX1","SRXN1","TP53BP1","TP53","MDM2","CDKN1A","BTG2","NFKB1","NFKBIA","RELA","TNFAIP3","ICAM1","HSPA5","ATF4","DDIT3","XBP1","EIF2A","ATF6","ERN1","EIF2AK3","PPP1R15A")
# Note: "PPP1R15A" = GADD34, "EIF2AK3" = PERK, "ERN1" = IRE1, "DDIT3" = CHOP, "HSPA5" = BIP, "TNFAIP3" = A20, "NFKBIA" = IkBa, "NFE2L2" = Nrf2, "CDKN1A" = p21

# InvestigateModel <- Nodes_MinModel; FileTag <- "MinModel"
InvestigateModel <- Nodes_FullModel; FileTag <- "FullModel"

# interactions = import_Omnipath_Interactions(filter_databases=c("SignaLink3","PhosphoSite", "Signor")) # filter for 3 databases
interactions = import_Omnipath_Interactions() # no filter
OPI_g = interaction_graph(interactions = interactions)

AllNodes_OP <- sort(unique(c(interactions$source_genesymbol,interactions$target_genesymbol)))
Idx_InOP <- is.element(InvestigateModel,AllNodes_OP)
InvestigateModel <- InvestigateModel[Idx_InOP]


# MinMod_Comb <- t(combn(Nodes_MinModel,2))

# install.packages('gtools')
library('gtools')
MinMod_Comb <- permutations(n=length(InvestigateModel),r=2,v=InvestigateModel,repeats.allowed = F)

MinMod_CT <- list()
Path0Idx <- NULL; Path1Idx <- NULL; Path2Idx <- NULL; Path3Idx <- NULL

for (counter in 1:nrow(MinMod_Comb)) {
  print(paste0(counter,"/",nrow(MinMod_Comb)))
  CurrentCT <- printPath_es(shortest_paths(OPI_g,from = MinMod_Comb[counter,1],to = MinMod_Comb[counter,2], output = 'epath')$epath[[1]],OPI_g)
  MinMod_CT[[counter]] <- list("IntAct" = paste0(MinMod_Comb[counter,1],"_",MinMod_Comb[counter,2]),"Length_CT" = length(CurrentCT$source), "Shortest_Path" = CurrentCT)
  if (length(CurrentCT$source)==0) {Path0Idx <- c(Path0Idx,counter)}
  if (length(CurrentCT$source)==1) {Path1Idx <- c(Path1Idx,counter)}
  if (length(CurrentCT$source)==2) {Path2Idx <- c(Path2Idx,counter)}
  if (length(CurrentCT$source)==3) {Path3Idx <- c(Path3Idx,counter)}
}

sink()

sink(paste0("Crosstalk_Report_",FileTag,".txt"))
cat("=============================")
cat("\n")
cat(paste0("Direct interaction: n=",length(Path1Idx)))
cat("\n")
cat("=============================")
cat("\n")
cat("\n")
for (counter in 1:length(Path1Idx)) {
  cat(MinMod_CT[[Path1Idx[counter]]]$IntAct)
  cat("\n")
  for (counter2 in 1:nrow(MinMod_CT[[Path1Idx[counter]]]$Shortest_Path)) {
    cat(unlist(MinMod_CT[[Path1Idx[counter]]]$Shortest_Path[counter2,]))
    cat("\n")
  }
  cat("\n")
}
cat("=============================")
cat("\n")
cat(paste0("1 neighborhood interaction: n=",length(Path2Idx)))
cat("\n")
cat("=============================")
cat("\n")
cat("\n")
for (counter in 1:length(Path2Idx)) {
  cat(MinMod_CT[[Path2Idx[counter]]]$IntAct)
  cat("\n")
  for (counter2 in 1:nrow(MinMod_CT[[Path2Idx[counter]]]$Shortest_Path)) {
    cat(unlist(MinMod_CT[[Path2Idx[counter]]]$Shortest_Path[counter2,]))
    cat("\n")
  }
  cat("\n")
}
cat("=============================")
cat("\n")
cat(paste0("2 neighborhood interaction: n=",length(Path3Idx)))
cat("\n")
cat("=============================")
cat("\n")
cat("\n")
for (counter in 1:length(Path3Idx)) {
  cat(MinMod_CT[[Path3Idx[counter]]]$IntAct)
  cat("\n")
  for (counter2 in 1:nrow(MinMod_CT[[Path3Idx[counter]]]$Shortest_Path)) {
    cat(unlist(MinMod_CT[[Path3Idx[counter]]]$Shortest_Path[counter2,]))
    cat("\n")
  }
  cat("\n")
}

sink()


# --- End of the script --- # 
