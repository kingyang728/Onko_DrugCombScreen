suppressPackageStartupMessages(library(maftools, quietly=TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(VariantAnnotation, quietly=TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(GenomicFeatures, quietly=TRUE, warn.conflicts = FALSE))
# suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19, quietly=TRUE, warn.conflicts = FALSE))
# suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38, quietly=TRUE, warn.conflicts = FALSE))
# Combined_DrugDB_path <- "Data/Combined_DrugDB.csv"
Comb_DrugClassified_DB_path <- "Data/Combined_DrugClassifiedDB.csv"
test_primary_file <- "test_data/Primary_Her2_mutationSNVCNV_DF.csv"
test_comparison_file <- "test_data/Comparison_Normallike_BRCA_mutation_DF.csv"
test_cellline_file <- "test_data/Cellline_CCLEBRCA_Mutation_DF.csv"

# Combined_DrugDB <- read.csv(Combined_DrugDB_path,row.names = NULL,header = TRUE,check.names = FALSE)
DrugClassified_DB <- read.csv(Comb_DrugClassified_DB_path,row.names = NULL,header = TRUE,check.names = FALSE)
test_Target_data <- read.csv(test_primary_file, header = TRUE, sep = ",")
test_Comparison_data <- read.csv(test_comparison_file, header = TRUE, sep = ",")
test_Cellline_data <- read.csv(test_cellline_file, header = TRUE, sep = ",")

# DLBCLc1Significant_DurgCombID_Path <- "Data/DLBCLc1_SignificantDrugCombsIDDF.csv"
# DLBCLc2Significant_DurgCombID_Path <- "Data/DLBCLc2_SignificantDrugCombsIDDF.csv"
# DLBCLc3Significant_DurgCombID_Path <- "Data/DLBCLc3_SignificantDrugCombsIDDF.csv"
# DLBCLc4Significant_DurgCombID_Path <- "Data/DLBCLc4_SignificantDrugCombsIDDF.csv"
# DLBCLc5Significant_DurgCombID_Path <- "Data/DLBCLc5_SignificantDrugCombsIDDF.csv"
# 
# cHL_Significant_DurgCombID_Path<- "Data/cHL_SignificantDrugCombsIDDF.csv"
# PMBL_Significant_DurgCombID_Path<- "Data/PMBL_SignificantDrugCombsIDDF.csv"
# DLBCL_Significant_DurgCombID_Path<- "Data/DLBCL_SignificantDrugCombsIDDF.csv"
# 
# DLBCLc1Significant_DurgCombIDDF <- read.csv(DLBCLc1Significant_DurgCombID_Path,row.names = NULL,header = TRUE,check.names = FALSE)
# DLBCLc2Significant_DurgCombIDDF <- read.csv(DLBCLc2Significant_DurgCombID_Path,row.names = NULL,header = TRUE,check.names = FALSE)
# DLBCLc3Significant_DurgCombIDDF <- read.csv(DLBCLc3Significant_DurgCombID_Path,row.names = NULL,header = TRUE,check.names = FALSE)
# DLBCLc4Significant_DurgCombIDDF <- read.csv(DLBCLc4Significant_DurgCombID_Path,row.names = NULL,header = TRUE,check.names = FALSE)
# DLBCLc5Significant_DurgCombIDDF <- read.csv(DLBCLc5Significant_DurgCombID_Path,row.names = NULL,header = TRUE,check.names = FALSE)
# 
# cHL_Significant_DurgCombIDDF <- read.csv(cHL_Significant_DurgCombID_Path,row.names = NULL,header = TRUE,check.names = FALSE)
# PMBL_Significant_DurgCombIDDF <- read.csv(PMBL_Significant_DurgCombID_Path,row.names = NULL,header = TRUE,check.names = FALSE)
# DLBCL_Significant_DurgCombIDDF <- read.csv(DLBCL_Significant_DurgCombID_Path,row.names = NULL,header = TRUE,check.names = FALSE)


### TXDB generation
GFFTXDBDataload <- function(gff_File){
  txdbName <- tools::file_path_sans_ext(gff_File)
  txDBPath <-paste(txdbName,"txdb.sqlite",sep = "_")

  if(file.exists(txDBPath)  ){
    txdb <<- loadDb(txDBPath)
    seqlevelsStyle(txdb) <- "UCSC"
    txDBPath<<-txDBPath
  } else if(!file.exists(txDBPath)){
    txdb <<- makeTxDbFromGFF(gff_File,organism="Homo sapiens")
    seqlevelsStyle(txdb) <- "UCSC"
    txdbName <- tools::file_path_sans_ext(gff_File)
    txDBPath <<-paste(txdbName,"txdb.sqlite",sep = "_")
    saveDb(txdb, file=txDBPath)

  } else if(! exists("txDBPath") ){
    txdb <<- makeTxDbFromGFF(gff_File,organism="Homo sapiens")
    seqlevelsStyle(txdb) <- "UCSC"
    txdbName <- tools::file_path_sans_ext(gff_File)
    txDBPath <<-paste(txdbName,"txdb.sqlite",sep = "_")
    saveDb(txdb, file=txDBPath)
  }
  return(txDBPath)
}
# ###Saving and Loading a TxDb Object
# 
# 
# TxDB_generate_from_URI <- function(GFFURI,save_path){
#   gff_File<-download_GFF_fromURI(GFFURI,save_path)
#   txdb <<- makeTxDbFromGFF(gff_File,organism="Homo sapiens")
#   seqlevelsStyle(txdb) <- "UCSC"
#   txDBtempPath<-dirname(gff_File)
#   txdbName <- tools::file_path_sans_ext(gff_File)
#   txDBPath <<-paste(txdbName,"txdb.sqlite",sep = "_")
#   saveDb(txdb, file=txDBPath)
#   #  txdb <- loadDb(txDBPath)
#   return(txDBPath)
# }
# gff_File = "Data/Homo_sapiens.GRCh38.110.gff3.gz"
# txDBPath = "Data/Homo_sapiens.GRCh38.110.gff3_txdb.sqlite"
# txDBPath<-GFFTXDBDataload(gff_File)

#############
VCF_to_SNV <- function(Input_Vcf,txDBPath){
  #txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  #txDBPath ="gff_annotation.sqlite"
  txdb <- loadDb(txDBPath)
  vcf <- readVcf(Input_Vcf)
  print("t3")
  seqlevelsStyle(vcf) <- "UCSC"
  seqlevelsStyle(txdb) <- "UCSC"
  seqlevels(vcf,pruning.mode="coarse")<-intersect(seqlevels(vcf), seqlevels(txdb))
  
  coding <- predictCoding(vcf, txdb, seqSource=Hsapiens)
  colCoding <-mcols(coding)
  colCoding <- as.data.frame(colCoding)
  
  Protein_Change <- paste(colCoding$REFAA,colCoding$PROTEINLOC,colCoding$VARAA,sep="")
  #replace item last char as '=' which amino acid code did not changed.
  Protein_Change[colCoding$REFAA==colCoding$VARAA] <- paste(colCoding$REFAA,colCoding$PROTEINLOC,"=",sep="")[colCoding$REFAA==colCoding$VARAA]
  Hugo_Symbol <- colCoding$GENEID
  Variant_Classification<-colCoding$CONSEQUENCE
  
  #calculate insertion or deletion by length difference of varAllele and REF
  Indel_status<- lengths(colCoding$varAllele)-lengths(colCoding$REF)
  SNVTable<- data.frame(Hugo_Symbol,Variant_Classification,Protein_Change,Indel_status)
  
  ####convert variant_classification
  print("t4")
  ###remove silent mutation
  SNVTable<-SNVTable[!(SNVTable[,"Variant_Classification"]=="synonymous"),]
  
  #replace consequence name
  SNVTable$Variant_Classification<- as.character(SNVTable$Variant_Classification)
  SNVTable$Variant_Classification[SNVTable$Variant_Classification=="nonsense"] <-  "Nonsense_Mutation"
  SNVTable$Variant_Classification[SNVTable$Variant_Classification=="nonsynonymous"] <-  "Missense_Mutation"
  SNVTable$Variant_Classification[SNVTable$Indel_status > 0 & SNVTable$Variant_Classification=="frameshift"] <- "Frame_Shift_Ins"
  SNVTable$Variant_Classification[SNVTable$Indel_status < 0 & SNVTable$Variant_Classification=="frameshift"] <- "Frame_Shift_Del"
  SNVTable$Variant_Classification[SNVTable$Indel_status > 0 & SNVTable$Variant_Classification=="missense"] <- "In_Frame_Ins"
  SNVTable$Variant_Classification[SNVTable$Indel_status < 0 & SNVTable$Variant_Classification=="missense"] <- "In_Frame_Del"
  SNVTable$Indel_status <- NULL
  
  ##################################
  SNVTable = unique(SNVTable)
  rownames(SNVTable)<-c(1:nrow(SNVTable))
  return(SNVTable)
}
getSNVfromMaf <- function(mafInput) {
  maftable<-mafInput@data
  
  
  ######### unify synonyms############
  if("HGVSp_Short" %in% names(maftable)&"Protein_Change" %in% names(maftable)){
    print("Problem: duplicate protein change column ")
  } else if("HGVSp_Short" %in% names(maftable)&!("Protein_Change" %in% names(maftable))) {
    protein_Change_Index <- which(names(maftable)=="HGVSp_Short")
    names(maftable)[protein_Change_Index] <- "Protein_Change"
  }
  ##########################
  SNV<-maftable
  
  #  SNV<-maftable[maftable$Variant_Type!="SNP",]    #remove SNP
  SNV<-as.data.frame(SNV)                                     #convert data.table to data.frame
  SNV                = SNV[SNV$Variant_Classification!='Silent',]
  SNV$Protein_Change = gsub("p.","",SNV[,"Protein_Change"])
  SNV                = SNV[,c("Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification","Protein_Change")]
  SNV                = SNV[which(SNV[,"Protein_Change"] != ""),]
  colnames(SNV)[colnames(SNV) == "Tumor_Sample_Barcode"] <- "Sample_ID"
  SNV                = unique(SNV)
  return(SNV)
}

MAF_to_SNVTable<-function(MafInputFile,txDBPath){
  mafInput = read.maf(MafInputFile)
  mafFields<-getFields(mafInput)
  if("HGVSp_Short" %in% mafFields |"Protein_Change" %in% mafFields  ){
    SNVTable<-getSNVfromMaf(mafInput)
    SNVTable = unique(SNVTable)
    return(SNVTable)
  } else {
    maftable<-mafInput@data
    Chromosome <- sub("chr","", maftable$Chromosome )
    #maftable[,"Start_Position"]
    
    Coordinate_input_Table<-data.frame(Chromosome,maftable$Start_Position,maftable$End_Position,maftable$Reference_Allele,
                                       maftable$Tumor_Seq_Allele2,maftable$Strand)
    Coordinate_input_Table <-unique(Coordinate_input_Table)
    colnames(Coordinate_input_Table) <- c("Chromosome",	"start","end","refAllele","varAllele","strand")
    txdb <- loadDb(txDBPath)
    if("strand" %in%  colnames(Coordinate_input_Table)){
      Coordinate_Granges <- GRanges(seqnames = Rle(Coordinate_input_Table$Chromosome),
                                    ranges = IRanges(start = Coordinate_input_Table$start, end = Coordinate_input_Table$end),
                                    strand = Coordinate_input_Table$strand,
                                    REF = Coordinate_input_Table$refAllele)
    } else {
      Coordinate_Granges <- GRanges(seqnames = Rle(Coordinate_input_Table$Chromosome),
                                    ranges = IRanges(start = Coordinate_input_Table$start, end = Coordinate_input_Table$end),
                                    REF = Coordinate_input_Table$refAllele)
    }
    seqlevelsStyle(txdb) <- "UCSC"
    seqlevelsStyle(Coordinate_Granges) <- "UCSC"
    
    VarAllele<-DNAStringSet(Coordinate_input_Table$varAllele)
    
    coding <- predictCoding(Coordinate_Granges, txdb, Hsapiens,VarAllele)
    colCoding <-mcols(coding)
    colCoding <- as.data.frame(colCoding)
    Protein_Change <- paste(colCoding$REFAA,colCoding$PROTEINLOC,colCoding$VARAA,sep="")
    #replace item last char as '=' which amino acid code did not changed.
    Protein_Change[colCoding$REFAA==colCoding$VARAA] <- paste(colCoding$REFAA,colCoding$PROTEINLOC,"=",sep="")[colCoding$REFAA==colCoding$VARAA]
    
    SNVTable<- data.frame(maftable$Tumor_Sample_Barcode,maftable$Hugo_Symbol,maftable$Variant_Classification,Protein_Change)
    colnames(SNV)[colnames(SNV) == "Tumor_Sample_Barcode"] <- "Sample_ID"
    SNVTable = unique(SNVTable)
    return(SNVTable)
  }
}



#### Get Drug Table
Drug_predict<- function(Mutation_DF,Drug_DB,patient_cancerType){
  Drug_DF <- merge(Mutation_DF, Drug_DB, by.x="Hugo_Symbol",by.y = "Gene")
  #Tumor_Solid_Status <- T
  Drug_DF$level <- ""
  
  ### set drug level based on the patint cancertype. Cancer type like any patient cancer type with general cancer and Unspecified Cancer can be set as A
  if(grepl(patient_cancerType, "Lymphoma",ignore.case = T)){
    if(length(Drug_DF[grepl(patient_cancerType, Drug_DF$Cancer,ignore.case = T), ]$level)!=0) {
      Drug_DF[grepl(patient_cancerType, Drug_DF$Cancer,ignore.case = T), ]$level <- paste0("A",Drug_DF[grepl(patient_cancerType, Drug_DF$Cancer,ignore.case = T), ]$DBlevel)
    }
    if(length(Drug_DF[!(grepl(patient_cancerType, Drug_DF$Cancer,ignore.case = T) ) , ]$level)!=0) {
      Drug_DF[!(grepl(patient_cancerType, Drug_DF$Cancer,ignore.case = T) ) , ]$level <- paste0("B",Drug_DF[!(grepl(patient_cancerType, Drug_DF$Cancer,ignore.case = T) ), ]$DBlevel)
    }
  } else if(grepl(patient_cancerType, "Leukemias|Blood",ignore.case = T)){
    if(length(Drug_DF[grepl(patient_cancerType, Drug_DF$Cancer,ignore.case = T), ]$level)!=0) {
      Drug_DF[grepl(patient_cancerType, Drug_DF$Cancer,ignore.case = T), ]$level <- paste0("A",Drug_DF[grepl(patient_cancerType, Drug_DF$Cancer,ignore.case = T), ]$DBlevel)
    }
    if(length(Drug_DF[!(grepl(patient_cancerType, Drug_DF$Cancer,ignore.case = T)), ]$level)!=0){
      Drug_DF[!(grepl(patient_cancerType, Drug_DF$Cancer,ignore.case = T)), ]$level <- paste0("B",Drug_DF[!(grepl(patient_cancerType, Drug_DF$Cancer,ignore.case = T)), ]$DBlevel)
    }
    
  } else {
    if(length(Drug_DF[grepl(patient_cancerType, Drug_DF$Cancer,ignore.case = T)|grepl("Unspecified Cancer|General Cancer",Drug_DF$Cancer,ignore.case = T), ]$level)!=0){
      Drug_DF[grepl(patient_cancerType, Drug_DF$Cancer,ignore.case = T)|grepl("Unspecified Cancer|General Cancer",Drug_DF$Cancer,ignore.case = T), ]$level <- paste0("A",Drug_DF[grepl(patient_cancerType, Drug_DF$Cancer,ignore.case = T)|grepl("Unspecified Cancer|General Cancer",Drug_DF$Cancer,ignore.case = T), ]$DBlevel)
    }
    if(length(Drug_DF[!(grepl(patient_cancerType, Drug_DF$Cancer,ignore.case = T)|grepl("Unspecified Cancer|General Cancer",Drug_DF$Cancer,ignore.case = T)), ]$level)!=0){
      Drug_DF[!(grepl(patient_cancerType, Drug_DF$Cancer,ignore.case = T)|grepl("Unspecified Cancer|General Cancer",Drug_DF$Cancer,ignore.case = T)), ]$level <- paste0("B",Drug_DF[!(grepl(patient_cancerType, Drug_DF$Cancer,ignore.case = T)|grepl("Unspecified Cancer|General Cancer",Drug_DF$Cancer,ignore.case = T)), ]$DBlevel)
    }
    
  }
  Drug_DF$DBlevel <- NULL
  
  if(!("Protein_Change" %in% colnames(Drug_DF))) Drug_DF$Protein_Change <- ""
  
  Drug_DF[Drug_DF$Protein_Change=="",]$Protein_Change <- as.character(Drug_DF[Drug_DF$Protein_Change=="",]$Variant_Classification)
  Drug_DF <- Drug_DF %>% dplyr::rename("Pat Var" = Protein_Change, Gene = Hugo_Symbol)
  #Drug_DF <- Drug_DF %>% relocate(level, .before = Gene)
  if("Sample_ID" %in% colnames(Drug_DF)){
    Drug_DF <- Drug_DF %>% relocate(any_of(c("Sample_ID","level", "Gene","Pat Var")))
  } else {
    Drug_DF <- Drug_DF %>% relocate(any_of(c("level", "Gene","Pat Var")))
  }
  
  # Drug_DF$Variant_Classification <- NULL
  Drug_DF<- Drug_DF[!grepl("A4|B4", Drug_DF$level,ignore.case = T),]
  Drug_DF<- unique(Drug_DF)
  Drug_DF
  
}

FilterDrug_DF <- function(Drug_DF){
  
  if(nrow(Drug_DF[grepl("\\w-\\w",Drug_DF$Pat.Var),])==0){
    fusionDrug_DF <- Drug_DF[grepl("\\w-\\w",Drug_DF$Pat.Var),]
  } else {
    fusionDrug_DFtemp1 <- Drug_DF[grepl("rearrangement|fusion",Drug_DF$Known.Var,ignore.case = T) & grepl("\\w-\\w",Drug_DF$Pat.Var), ]    ### GDKD/OnkoKB fusions
    fusionDrug_DFtemp2 <- Drug_DF %>% 
      rowwise() %>% 
      mutate(Matches = grepl(paste(unlist(strsplit(strsplit(Known.Var, " ")[[1]][1], split = "-")),collapse = "|"),Pat.Var) & grepl("[A-Z][A-Z0-9]{2,}-[A-Z][A-Z0-9]{2,}", strsplit(Known.Var, " ")[[1]][1]))   ### deal with CIVIC DB
    fusionDrug_DFtemp2 <- fusionDrug_DFtemp2[fusionDrug_DFtemp2$Matches,]
    fusionDrug_DFtemp2$Matches <- NULL    ##### CIVIC fusions
    fusionDrug_DF <- rbind(fusionDrug_DFtemp1,fusionDrug_DFtemp2)
  }
  
  CNVDrug_DF <- Drug_DF[(grepl("amplification|expression",Drug_DF$Known.Var,ignore.case = T) & grepl("amplification",Drug_DF$Pat.Var,ignore.case = T)  & !grepl("underexpression",Drug_DF$Known.Var,ignore.case = T))|
                          grepl("deletion|del",Drug_DF$Known.Var,ignore.case = T) & grepl("deletion|lof",Drug_DF$Pat.Var,ignore.case = T), ]
  # SNVDrug_DF <- Drug_DF[(!grepl("amplification|deletion|del|rearrangement|fusion|expression",Drug_DF$Known.Var,ignore.case = T) | 
  #                          (grepl("splice",Drug_DF$Pat.Var,ignore.case = T) & grepl("splice|loss|lof",Drug_DF$Known.Var,ignore.case = T)) |
  #                          (!(grepl("splice",Drug_DF$Pat.Var,ignore.case = T) & grepl("splice|loss|lof",Drug_DF$Known.Var,ignore.case = T)) & grepl("any",Drug_DF$Known.Var,ignore.case = T))) & !grepl("amplification|deletion|\\w-\\w",Drug_DF$Pat.Var), ]
  
  tmpSNVDrug_DF <- Drug_DF[!grepl("amplification|rearrangement|fusion|expression",Drug_DF$Known.Var,ignore.case = T) & !grepl("amplification|deletion|\\w-\\w",Drug_DF$Pat.Var),]
  tmpSNVDrug_DF <- tmpSNVDrug_DF[!(grepl("deletion|del",tmpSNVDrug_DF$Known.Var,ignore.case = T) & !grepl("any",tmpSNVDrug_DF$Known.Var,ignore.case = T)) ,]
  tmpSNVDrug_DF <- tmpSNVDrug_DF[grepl("splice",tmpSNVDrug_DF$Pat.Var,ignore.case = T) & grepl("loss|lof|splice|Truncating Mutations|fs",tmpSNVDrug_DF$Known.Var,ignore.case = T) | !grepl("splice",tmpSNVDrug_DF$Pat.Var,ignore.case = T),]
  SNVDrug_DF <- tmpSNVDrug_DF[grepl("frame",tmpSNVDrug_DF$Variant_Classification,ignore.case = T) & grepl("loss|lof|Truncating Mutations|any|fs|indel",tmpSNVDrug_DF$Known.Var,ignore.case = T) |grepl("in_frame_ins",tmpSNVDrug_DF$Variant_Classification,ignore.case = T) & grepl("ins",tmpSNVDrug_DF$Known.Var,ignore.case = T) | 
                                grepl("in_frame_del",tmpSNVDrug_DF$Variant_Classification,ignore.case = T) & grepl("del",tmpSNVDrug_DF$Known.Var,ignore.case = T) | !grepl("frame",tmpSNVDrug_DF$Variant_Classification,ignore.case = T),]
  
  Filtered_DrugDF <- do.call("rbind", list(SNVDrug_DF, CNVDrug_DF, fusionDrug_DF))
  Filtered_DrugDF$Variant_Classification <- NULL
  Filtered_DrugDF
}


Get_MTBreporter_DF <- function(Drug_DB,patient_cancerType,Mutation_DF,Filter_status){
  # Mutation_DF <- Read_CSV_reform(inp_file,inp_type)
  Drug_DF <- Drug_predict(Mutation_DF,Drug_DB,patient_cancerType)
  Drug_DF <- Drug_DF %>% dplyr::rename(Pat.Var = "Pat Var", Known.Var = "Known Var")
  Drug_DF <- Drug_DF %>%
    # |         Pat.Var        |  PatVar.INTPN    |
    # |:----------------------:|:----------------:|
    # |      In_Frame_Ins      | ins(exact match) |
    # |      In_Frame_Del      | del(exact match) |
    # |     Frame_Shift_Ins    |       loss       |
    # |     Frame_Shift_Del    |       loss       |
    # |       Splice_site      |       loss       |
    # |      amplification     |       gain       |
    # |        deletion        |       loss       |
    # |    Missense_Mutation   |     mutation     |
    # |    Nonsense_Mutation   |       loss       |
    # |    Nonstop_Mutation    |    exact match   |
    # | Translation_Start_Site |    exact match   |
  dplyr::mutate(PatVar.INTPN = case_when(
      grepl("In_Frame_Ins|inframe_ins", Variant_Classification, ignore.case = TRUE) ~ "ins",
      grepl("In_Frame_Del|inframe_deletion", Variant_Classification, ignore.case = TRUE) ~ "del",
      grepl("frameshift|Frame_Shift_Ins|Frame_Shift_Del|Splice_site|nonsense|deletion|stop_gain|splice_donor_variant", Variant_Classification, ignore.case = TRUE) ~ "loss",
      grepl("amplification", Variant_Classification, ignore.case = TRUE) ~ "gain",
      grepl("Missense|coding_sequence_variant|protein_altering_variant|incomplete_terminal_codon_variant|Start_Codon_SNP", Variant_Classification, ignore.case = TRUE) ~ "mutation",
      grepl("Nonstop_Mutation|Translation_Start_Site", Variant_Classification, ignore.case = TRUE) ~ "exact match",
      TRUE ~ as.character("NA") # default case to NA
    )) %>% 
    # |         Known.Var         |   knownVar.INTPN       |
    # |:-------------------------:|:----------------------:|
    # |          “(GoF)”          |          gain          |
    # |          “(LoF)”          |          loss          |
    # |          “splice”         |          loss          |
    # |          “delins”         |          loss          |
    # |           “ins”           |           ins          |
    # |           “del”           |           del          |
    # |          “indel”          |          loss          |
    # |            “fs”           |          loss          |
    # |         “deletion”        |          loss          |
    # |      “amplification”      |          gain          |
    # |          “dup”            |          gain          |
    # |            mut            |        mutation        |
    # |            any            |        mutation        |
    # |   loss/loss-of-function   |          loss          |
    # |         “mutation”        |        mutation        |
    # |       “^expression”       |          gain          |
    # |      “Overexpression”     |          gain          |
    # |     “Underexpression”     |          loss          |
    # |    Truncating Mutations   |          loss          |
    # |    Oncogenic Mutations    |        mutation        |
    # |        “FRAMESHIFT”       |          loss          |
    # |       “FRAME SHIFT”       |          loss          |
    # |     Exon 17 mutations     | mutation (exact match) |
    # |      Exon 19 Deletion     |   loss  (exact match)  |
    # | EXON 14 SKIPPING MUTATION | mutation (exact match) |
    # |  ^[A-Z][0-9]+[A-Z]        | mutation (exact match) |
    # |         "wild"            |         wildtype       |
    # | "[A-Za-z0-9_]+-[A-Za-z0-9_]|fusion|rearrangement" | fusion |
  dplyr::mutate(knownVar.INTPN = case_when(
      grepl("\\(GoF\\)", Known.Var, ignore.case = TRUE) ~ "gain",
      grepl("\\(LoF\\)", Known.Var, ignore.case = TRUE) ~ "loss",
      grepl("splice|delins|indel|fs|deletion|Underexpression|FRAMESHIFT|FRAME SHIFT|Truncating Mutations|loss|loss-of-function|VIII|LOH|Methylation|Biallelic Inactivation", Known.Var, ignore.case = TRUE) ~ "loss",
      grepl("amplification|expression|Overexpression|dup|PHOSPHORYLATION", Known.Var, ignore.case = TRUE) ~ "gain",
      grepl("([a-zA-Z0-9_]+ins[a-zA-Z0-9]*)|insertion|ITD", Known.Var, perl = TRUE, ignore.case = TRUE)  ~ "ins",
      grepl("([a-zA-Z0-9_]+del[a-zA-Z0-9]*)", Known.Var, perl = TRUE, ignore.case = TRUE) ~ "del",
      Known.Var %in% c("Exon 17 mutations", "EXON 14 SKIPPING MUTATION") ~ "exact match",
      Known.Var == "Exon 19 Deletion" ~ "exact match",
      grepl("[A-Z][0-9]+\\*", Known.Var, ignore.case = TRUE) ~ "loss",  ## capture the stop codon
      grepl("^[A-Z][0-9]+[A-Z]*", Known.Var, ignore.case = TRUE) ~ "mutation",
      grepl("mut|any|mutation|Oncogenic Mutations|ALTERATION", Known.Var, ignore.case = TRUE) ~ "mutation",
      grepl("wild", Known.Var, ignore.case = TRUE) ~ "wildtype",
      grepl("[A-Za-z0-9_]+-[A-Za-z0-9_]|fusion|rearrangement", Known.Var, ignore.case = TRUE) ~ "fusion",
      TRUE ~ as.character("mutation") # default case to NA
    )) %>% 
    dplyr::filter(knownVar.INTPN != "fusion" & knownVar.INTPN != "wildtype") %>%  
    dplyr::filter(
      (PatVar.INTPN == "gain" & knownVar.INTPN == "gain") | 
        (PatVar.INTPN == "loss" & knownVar.INTPN == "loss") | 
        (PatVar.INTPN == "ins" & knownVar.INTPN == "ins") | 
        (PatVar.INTPN == "del" & knownVar.INTPN == "del") |
        (knownVar.INTPN %in% c("mutation", "exact match"))
    )
  
  Drug_DF$Match_Sign <-"NONE"
  if(dim(Drug_DF[mapply(grepl, gsub("\\*", "\\\\*", Drug_DF$Pat.Var), gsub("\\*", "\\\\*", Drug_DF$Known.Var),MoreArgs = list(ignore.case = TRUE)),])[1]!=0){
    Drug_DF[mapply(grepl, gsub("\\*", "\\\\*", Drug_DF$Pat.Var), gsub("\\*", "\\\\*", Drug_DF$Known.Var),MoreArgs = list(ignore.case = TRUE)),]$Match_Sign <- "Exact"}
  
  if(dim(Drug_DF[mapply(grepl, "any|mut",  Drug_DF$Known.Var,MoreArgs = list(ignore.case = TRUE)),])[1]!=0) {
    Drug_DF[mapply(grepl, "any|mut",  Drug_DF$Known.Var,MoreArgs = list(ignore.case = TRUE)),]$Match_Sign <- "AnyMut"}
  Drug_DF <- Drug_DF[order(Drug_DF$Sample_ID),]  
  Drug_DF <-Drug_DF[!(Drug_DF$Drugs == "NA"),]
  
  Drug_DF <- Drug_DF %>% dplyr:: filter(( (knownVar.INTPN == "exact match"|PatVar.INTPN =="exact match") & Match_Sign == "Exact")|
                                          (PatVar.INTPN == "ins" & knownVar.INTPN == "ins" & Match_Sign == "Exact")|
                                          (PatVar.INTPN == "del" & knownVar.INTPN == "del" & Match_Sign == "Exact") |
                                          (PatVar.INTPN == "gain" & knownVar.INTPN == "gain") | 
                                          (PatVar.INTPN == "loss" & knownVar.INTPN == "loss") | 
                                           knownVar.INTPN == "mutation"
  ) %>% dplyr::select((-c(Variant_Classification,PatVar.INTPN, knownVar.INTPN)))
  return(Drug_DF)
  # if(missing(Filter_status)){
  #   FilteredDrug_DF<-FilterDrug_DF(Drug_DF)
  #   return(FilteredDrug_DF)
  # } else if(Filter_status){
  #   FilteredDrug_DF<-FilterDrug_DF(Drug_DF)
  #   return(FilteredDrug_DF)
  # }else{
  #   return(Drug_DF)
  # }
}

Get_DrugPredictsType_DF <- function(Drugs_DF,PredictType){
  # PredictsTypeDrug_DF <- NULL
  if(grepl("Positive|Negative|All",PredictType,ignore.case = T)){
    if(grepl("Positive",PredictType,ignore.case = T)){
      PredictsTypeDrug_DF <- Drugs_DF[grepl("^sensitivity$|^response$|^sensitivity/response$",Drugs_DF$Predicts,ignore.case = T),]  ### be careful for the case
    } else if(grepl("Negative",PredictType,ignore.case = T)){
      PredictsTypeDrug_DF <- Drugs_DF[grepl("^Adverse Response$|^resistance$|^no response$",Drugs_DF$Predicts,ignore.case = T),]
    } else if(grepl("All",PredictType,ignore.case = T)){
      PredictsTypeDrug_DF <- Drugs_DF
      
    }
    
  } else{
    print("Wrong predict type input")
  }
  return(PredictsTypeDrug_DF)
  
}
####

PlotDF_Create <- function(Significant_DurgCombsTB){
  FisherTestPlotDF <- Significant_DurgCombsTB %>%  mutate(
    log2oddsRatio = log2(oddsRatio),
    P.value = Significant_DurgCombsTB$p.value,
    AdjP.value = Significant_DurgCombsTB$adjust_p.value) %>% 
    dplyr::select(Drug_comb, Drug1, level1, Drug2, level2, Percentage, `CelllineConfirmed(%)`,log2oddsRatio,P.value,AdjP.value) 
   
  return(FisherTestPlotDF)
}

#FisherTestDF<-PlotDF_Create(cHL_Significant_DurgCombIDDF)

#### CirclePlot for multiple -----------
# CirclePlot_Preprocess <- function(...){
#   Significant_DurgCombID_PathList <- c(as.list(environment()), list(...))
#   #Cancer_type_name <- strsplit(basename(cHL_drug_Path), "_|-")[[1]][1]
#   ListNames <- lapply(Significant_DurgCombID_PathList,function(x) strsplit(basename(x), "_|-")[[1]][1])
#   Titles <- paste(ListNames,"Drug-combs chordDiagram")
#   names(Significant_DurgCombID_PathList) <- ListNames
#   Significant_DurgCombID_DFlist <- lapply(Significant_DurgCombID_PathList,read.csv,row.names = NULL,header = TRUE,check.names = FALSE)
#   #maxCelllineSupportList<-lapply(Significant_DurgCombID_DFlist,function(x) max(x$`CelllineConfirmed(%)`))
#   #Significant_DurgCombID_DFlist <- lapply(Significant_DurgCombID_DFlist, function(x) x[x$`CelllineConfirmed(%)`>=min(85,max(x$`CelllineConfirmed(%)`)),])  ### Cellline support level Celllinecutoff means percentage
#   CirclePlot_DFlist <- lapply(Significant_DurgCombID_DFlist, function(x) dplyr::select(x,Drug1,Drug2,Percentage))
#   Drugname_list <- lapply(CirclePlot_DFlist, function(x) union(x$Drug1,x$Drug2))
#   Drugnames<- Reduce(union, Drugname_list) 
#   DrugnameDF<-data.frame(Drugname = Drugnames,Colorindex=rainbow(length(Drugnames)))
#   
#   #gridColorlist <- lapply(Drugname_list, function(x) DrugnameDF[DrugnameDF$Drugname %in% x,]$Colorindex)
#   DrugNameDFlist <- lapply(Drugname_list, function(x) DrugnameDF[DrugnameDF$Drugname %in% x,])
#   
#   CirclePlot_dataList<-list(CirclePlot_DFlist=CirclePlot_DFlist,DrugNameDFlist=DrugNameDFlist)
#   return(CirclePlot_dataList)
#   
# }
# ----------------

CirclePlot_DF_Create <- function(DrugComb_Analysis_Table){
  CirclePlot_DF <- DrugComb_Analysis_Table %>% dplyr::select(Drug1,Drug2,Percentage) 
  Drugnames <- union(CirclePlot_DF$Drug1,CirclePlot_DF$Drug2)
  DrugnameDF<-data.frame(Drugname = Drugnames,Colorindex=rainbow(length(Drugnames)))
  CirclePlot_DFList<-list(CirclePlot_DF=CirclePlot_DF,DrugnameDF=DrugnameDF)
  return(CirclePlot_DFList)
  
}
#Significant_DurgCombID_PathList <- list(cHL_Significant_DurgCombID_Path,PMBL_Significant_DurgCombID_Path,DLBCL_Significant_DurgCombID_Path,DLBCLc1Significant_DurgCombID_Path,DLBCLc2Significant_DurgCombID_Path,DLBCLc3Significant_DurgCombID_Path,DLBCLc4Significant_DurgCombID_Path,DLBCLc5Significant_DurgCombID_Path)
# CirclePlot_dataList<-CirclePlot_Preprocess(cHL_Significant_DurgCombID_Path,PMBL_Significant_DurgCombID_Path,DLBCL_Significant_DurgCombID_Path,DLBCLc1Significant_DurgCombID_Path,DLBCLc2Significant_DurgCombID_Path,DLBCLc3Significant_DurgCombID_Path,DLBCLc4Significant_DurgCombID_Path,DLBCLc5Significant_DurgCombID_Path)


#### onkopus annotated drug data

