library(dplyr)

merge_maxEvidenceLevel_toFisherDF <- function(Cooccur_FreqDF,drugDF){
  # drugDF <- read.csv(DrugDF_path,row.names = NULL,header = TRUE )
  drugvsLevel <- drugDF[c("Classified_Drug_Name","level")]
  drugvsLevel$level <- factor(drugDF$level,levels=c("B3","A3","B2","A2","B1","A1"),ordered=TRUE)
  drugvsLevel <- na.omit(drugvsLevel)
  drugvsLevel <- unique(drugvsLevel)
  #drugvsLevel <- aggregate(drugvsLevel$level, by=list(drugvsLevel$Classified_Drug_Name), FUN=max)
  drugvsLevel <- aggregate(level ~ Classified_Drug_Name, data = drugvsLevel, max)
  
  Cooccur_FreqDF <- merge(Cooccur_FreqDF, drugvsLevel, by.x = "Drug2" ,by.y = "Classified_Drug_Name")
  Cooccur_FreqDF <- Cooccur_FreqDF %>% relocate(level2 = level, .after = Drug2) 
  
  Cooccur_FreqDF <- merge(Cooccur_FreqDF, drugvsLevel, by.x = "Drug1" ,by.y = "Classified_Drug_Name")
  Cooccur_FreqDF <- Cooccur_FreqDF %>% relocate(level1 = level, .after = Drug1)
  return(Cooccur_FreqDF)
}
get_DrugComb_FreqDF <- function(drugDF){
  
  sample_size <- length(unique(drugDF$Sample_ID))
  
  drugDF_pre <- drugDF[c("Sample_ID","Classified_Drug_Name")] 
  drugDF_pre <- unique(drugDF_pre)
  drugDF_pre <- na.omit(drugDF_pre)
  #drugDF$Drug_Occur <- 1
  
  V <- crossprod(table(drugDF_pre[1:2]))
  diag(V) <- 0
  V[lower.tri(V)]<-0
  CooccurDF<- as.data.frame(V)
  CooccurDF$Drug <- rownames(CooccurDF)
  Cooccur_FreqDF <- melt(CooccurDF,id.vars = c("Drug"))
  #Cooccur_FreqDF<-Cooccur_FreqDF[Cooccur_FreqDF$value!=0,]
  
  names(Cooccur_FreqDF) <- c( "Drug1","Drug2","Matched_Count")
  Cooccur_FreqDF<-Cooccur_FreqDF[Cooccur_FreqDF$Drug1!=Cooccur_FreqDF$Drug2,]
  Cooccur_FreqDF<-Cooccur_FreqDF[order(-Cooccur_FreqDF$Matched_Count),]
  Cooccur_FreqDF$NonMatched_Count <- sample_size - Cooccur_FreqDF$Matched_Count
  
  Cooccur_FreqDF$Percentage <- round(Cooccur_FreqDF$Matched_Count /sample_size,digits = 4)*100
  Cooccur_FreqDF$Drug_comb <- paste(Cooccur_FreqDF$Drug1,"+",Cooccur_FreqDF$Drug2,sep=" ")
  Cooccur_FreqDF <- merge_maxEvidenceLevel_toFisherDF(Cooccur_FreqDF,drugDF)


  return(Cooccur_FreqDF)
}
# target_DrugFreqDF <- get_DrugComb_FreqDF(cHL_drugDF)
# control_DrugFreqDF <- get_DrugComb_FreqDF(PMBL_drugDF)
# Cellline_DrugFreqDF <- get_DrugComb_FreqDF(Cellline_drugDF)

FisherTestDF_generate <- function(target_DrugFreqDF,control_DrugFreqDF,Freq_cutoff,testType){
  if(missing(testType)) 
    testType <- "greater"
  control_Count_DrugDF <- control_DrugFreqDF %>%  dplyr::select(Drug_comb,Matched_Count,NonMatched_Count) %>% 
    dplyr::rename(Side_Matched_Count = Matched_Count,
           Side_NonMatched_Count = NonMatched_Count) %>% distinct()
  FisherTestDF <- merge(target_DrugFreqDF,control_Count_DrugDF,by = "Drug_comb",all.x = T) 
  ### set Side_Matched_Count Na as 0 and Side_NonMatched_Count as max
  FisherTestDF$Side_Matched_Count[is.na(FisherTestDF$Side_Matched_Count)] <- 0
  FisherTestDF$Side_NonMatched_Count[is.na(FisherTestDF$Side_NonMatched_Count)] <- max(FisherTestDF$Side_Matched_Count+FisherTestDF$Side_NonMatched_Count,na.rm=TRUE)
  FisherTestDF <- FisherTestDF %>% relocate(Percentage, .after = last_col())   ### reorder the columns put Percentage column to the last.
  FisherTestDF <- FisherTestDF[FisherTestDF$Percentage>= Freq_cutoff,]  ###  Freq cutoff > or >=
  
  fisherTestResult <- apply(FisherTestDF, 1,
                            function(x) {
                              tbl <- matrix(as.numeric(c(x["Matched_Count"], x["NonMatched_Count"], x["Side_Matched_Count"], x["Side_NonMatched_Count"])), ncol=2, byrow=F)
                              fisherTestResult <- fisher.test(tbl, alternative=testType) #### oneside test vs two side test, two.sided/greater/less
                              c(p.value=fisherTestResult$p.value , oddsRatio=as.numeric(fisherTestResult$estimate))
                            })
  FisherTestResultDF<-as.data.frame(t(fisherTestResult))
  FisherTestDF<-cbind(FisherTestDF,FisherTestResultDF)
  
  FisherTestDF$adjust_p.value <- p.adjust(FisherTestDF$p.value ,method="BH")
  
  FisherTestDF$p.value <- round(FisherTestDF$p.value,digits=16)
  FisherTestDF$adjust_p.value <- round(FisherTestDF$adjust_p.value,digits=16)
  
  # FisherTestDF <- merge_maxEvidenceLevel_toFisherDF(FisherTestDF,drugDF)
  FisherTestDF<-FisherTestDF[order(FisherTestDF$adjust_p.value),]
}
# FisherTestDF <- FisherTestDF_generate(target_DrugFreqDF,control_DrugFreqDF,20)


FisherTestDF_appendOriginalDrugs <- function(FisherTestDF,drugDF){
  # drug_DF <- read.csv(drugDF_Path,row.names = NULL,header = TRUE,check.names=F )
  
  Onkopus_DrugClassDF <-drugDF[c("Origin_Drug_Name","Classified_Drug_Name")]
  Onkopus_DrugClassDF <-na.omit(unique(Onkopus_DrugClassDF))
  
  #  DrugTraceTBAggre <-Onkopus_DrugClassDF %>% 
  #    group_by(Classified_Drug_Name) %>% 
  #    summarise(Origin_Drug_Name = paste(unique(Origin_Drug_Name), collapse = ', '))    ### group by column and concatenate unique strings.
  DrugTraceTBAggre <-aggregate(Origin_Drug_Name ~Classified_Drug_Name, Onkopus_DrugClassDF, toString)  ### group by column and concatenate unique strings.
  
  Primary_FCTB_WithSourceDrug <-merge(FisherTestDF,DrugTraceTBAggre,by.x = "Drug2",by.y = "Classified_Drug_Name" , all.x = T)   ## Drug2 vs Origin_Drug_Name.x 
  Primary_FCTB_WithSourceDrug<-merge(Primary_FCTB_WithSourceDrug,DrugTraceTBAggre,by.x = "Drug1",by.y = "Classified_Drug_Name" , all.x = T)   ## Drug1 vs Origin_Drug_Name.y
  
  Primary_FCTB_WithSourceDrug <- Primary_FCTB_WithSourceDrug %>% relocate(Drug_comb,Drug1,level1,Drug2,level2)
  Primary_FCTB_WithSourceDrug <- dplyr::rename(Primary_FCTB_WithSourceDrug, Origin_Drug1 = Origin_Drug_Name.y, Origin_Drug2 = Origin_Drug_Name.x )
  Primary_FCTB_WithSourceDrug <- Primary_FCTB_WithSourceDrug %>% relocate(Origin_Drug2, .after = Origin_Drug1)
  return(Primary_FCTB_WithSourceDrug)
}
# Primary_FisherComparativeTB <- FisherTestDF_appendOriginalDrugs(FisherTestDF,cHL_drugDF)


Significant_DurgCombs_Extract<- function(Primary_FisherComparativeTB,Cellline_drugDF){
  colnames<-colnames(Primary_FisherComparativeTB)
  TargetAdjustPvalue_Index<- grep("adjust_p.value",colnames)[1]
  TargetAdjustPvalue<-colnames[TargetAdjustPvalue_Index]
  
  Significant_DurgCombsTB <- Primary_FisherComparativeTB
  
  Cellline_DrugFreqDF <- get_DrugComb_FreqDF(Cellline_drugDF)
  CelllineConfirmed_DrugCombs<- intersect(Significant_DurgCombsTB$Drug_comb, Cellline_DrugFreqDF$Drug_comb)
  TargetCellline_comparativeDrugCombDF<- Cellline_DrugFreqDF[Cellline_DrugFreqDF$Drug_comb %in% CelllineConfirmed_DrugCombs,]
 
  TargetCellline_comparativeDrugCombDF<-TargetCellline_comparativeDrugCombDF[c("Drug_comb","Percentage","Matched_Count","NonMatched_Count")]
  colnames(TargetCellline_comparativeDrugCombDF) <- c("Drug_comb","CelllineConfirmed(%)","Cellline_MatchCount","Cellline_UnmatchCount")
  Significant_DurgCombsTB<-merge(Significant_DurgCombsTB,TargetCellline_comparativeDrugCombDF,,by.x ="Drug_comb", by.y = "Drug_comb" ,all.x = T)
  Significant_DurgCombsTB$Cellline_UnmatchCount[is.na(Significant_DurgCombsTB$Cellline_UnmatchCount)] <- max(Significant_DurgCombsTB$Cellline_MatchCount+Significant_DurgCombsTB$Cellline_UnmatchCount,na.rm=TRUE)
  Significant_DurgCombsTB[is.na(Significant_DurgCombsTB)] <-0
  Cellline_drugIDpre_DF <- Cellline_drugDF %>% dplyr::select(Sample_ID,Classified_Drug_Name) %>% distinct()
  Cellline_drugCombIDDF <- Cellline_drugIDpre_DF %>%
    inner_join(Cellline_drugIDpre_DF, by = "Sample_ID") %>%
    filter(Classified_Drug_Name.x != Classified_Drug_Name.y) %>%
    mutate(Drug_comb = paste(pmin(Classified_Drug_Name.x, Classified_Drug_Name.y), pmax(Classified_Drug_Name.x, Classified_Drug_Name.y), sep = " + ")) %>%
    distinct(Sample_ID, Drug_comb) %>%
    dplyr::select(Drug_comb,Sample_ID) %>%
    group_by(Drug_comb) %>%
    summarise(CelllineSample_IDs = paste(Sample_ID, collapse = ", "))
  Significant_CelllineID_DurgCombsTB <- Significant_DurgCombsTB %>% left_join(Cellline_drugCombIDDF,by = "Drug_comb")
  #write.csv(Significant_DurgCombsTB,Output_Path,row.names = FALSE,na=" ")
  return(Significant_CelllineID_DurgCombsTB)
}
# Significant_DurgCombsTB<-Significant_DurgCombs_Extract(Primary_FisherComparativeTB,Cellline_drugDF)


DrugComb_analysis <- function(target_DrugDF,control_DrugDF,Cellline_DrugDF,Freq_cutoff,testType){
  if(missing(testType)) 
    testType <- "greater"
  target_DrugFreqDF <- get_DrugComb_FreqDF(target_DrugDF)
  control_DrugFreqDF <- get_DrugComb_FreqDF(control_DrugDF)
  
  FisherTestDF <- FisherTestDF_generate(target_DrugFreqDF,control_DrugFreqDF,Freq_cutoff,testType)
  Primary_FisherComparativeTB <- FisherTestDF_appendOriginalDrugs(FisherTestDF,target_DrugDF)
  if(missing(Cellline_DrugDF)){
    return(Primary_FisherComparativeTB)
  } else {  Cellline_DrugFreqDF <- get_DrugComb_FreqDF(Cellline_DrugDF)
  Significant_DurgCombsTB<-Significant_DurgCombs_Extract(Primary_FisherComparativeTB,Cellline_DrugDF)
  return(Significant_DurgCombsTB)
  }
}


# test <- DrugComb_analysis(cHL_drugDF,PMBL_drugDF,Cellline_drugDF,Freq_cutoff = 20)
