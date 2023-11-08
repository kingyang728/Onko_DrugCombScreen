library(ggrepel)
library(dplyr)
library(reshape2)
library(ggalluvial)
library(circlize)
library(ComplexHeatmap)
##################################
# REACTIVE VOLCANO PLOT FUNCTION #
##################################

## Volcano plot ----------------
plotVolcano <- function(data, 
                        pvalue_col, 
                        pvalue_thresh, 
                        logfc_thresh,
                        de_vec,
                        show_logfc_thresh,
                        show_pvalue_thresh,
                        highlight_drugcombs = NULL,
                        x_label,
                        y_label,
                        xlim,
                        ylim,
                        specific_category) {
  # check that columns exist
  # FIXME don't include drug_comb as required
  if (!all(c("log2oddsRatio", pvalue_col) %in% colnames(data))) {
    stop("provided column names do not match dataset")
  }
  
  # convert pval to -log10(pval)
  data <- mutate(data,
                 log_pval = -log10(data[[pvalue_col]]))
  
  # Reorder de_vec factors so that TRUE is first
  de_vec <- factor(de_vec, levels=c("TRUE", "FALSE"))
  data$Significant <- de_vec

  # build base of plot
  # volcano <- ggplot(data, aes(x = .data[["log2oddsRatio"]], y = log_pval)) + geom_point(alpha = .6, aes(color = de_vec))
  
  data$delabel <- NA
  data$delabel[data$Significant != "FALSE"] <- data$Drug_comb[data$Significant != "FALSE"]
  data$`level1+level2` <-  paste(data$level1,data$level2,sep="+")
  data$`level1+level2` <- factor(data$`level1+level2`)
  
  volcano<-ggplot(data, aes(x=log2oddsRatio, y=log_pval, size = Percentage, shape = `level1+level2`,col=Significant, label=delabel)) +
    geom_point() + 
    theme_minimal() +
    geom_text_repel() +
    #scale_color_manual(values=colorVector) +
    # scale_color_manual(values=c("black", "red")) +
    ggtitle(paste(specific_category,"Volcano plot")) + 
    scale_shape_manual(values=1:nlevels(data$`level1+level2`)) +
    scale_size_continuous(range = c(2, 5))

  # if( all(unique(data$Significant) == c("NO","UP"))){
  #   testType = "greater"
  #   colorVector = c("black", "red")
  #   #colorVector = c("black", "red")
  # } else if(all(unique(data$Significant) == c("DOWN","NO"))){
  #   testType = "less"
  #   colorVector = c("blue","black")
  #   # colorVector = c("black", "red")
  # } else if(all(unique(data$Significant) == c("DOWN","NO","UP"))){
  #   testType = "two.sided"
  #   colorVector = c("blue","black", "red")
  #   # colorVector = c("blue","black", "red")
  # } else { testType = "greater"
  # colorVector = c("black", "red")}
  
  if (any(!is.null(c(xlim, ylim)))) {
    volcano <- volcano + coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE)
  }

  # if show_pvalue_thresh = true add vline layer
  if (show_pvalue_thresh) {
    volcano <- volcano + 
      geom_hline(yintercept = -log10(pvalue_thresh), linetype = "dashed", col = "grey", size = 1)
  }
  
  # if show_logfc_thresh = true add hline layer
  if (show_logfc_thresh) {
    volcano <- volcano + 
      geom_vline(xintercept = c(logfc_thresh, -logfc_thresh), linetype = "dashed", col = "grey", size = 1)
  }
  
  # if highlight_drugcombs isn't null, add to plot
  # add error handling for missing gene_Col
  if (!is.null(highlight_drugcombs)) {
    # create vector of gene ids for label
    labels <- data[data[["Drug_comb"]] %in% highlight_drugcombs,]
    volcano <- suppressWarnings(volcano +
        ggrepel::geom_label_repel(data = data[data[["Drug_comb"]] %in% highlight_drugcombs,], 
                                  aes(label = labels[["Drug_comb"]]))) 
  }
  
  volcanoPlot <- volcano +
    labs(x = x_label, y = y_label)+
    theme_classic(base_size = 12)
  # display plot
  volcanoPlot
}
## Heatmap ----------------
plotHeatmap <- function(data,
                        Cellline_cutoff,
                        col_high,
                        col_mid,
                        col_low
                        ){
  data[data$`CelllineConfirmed(%)`<=Cellline_cutoff,]$Percentage <- data[data$`CelllineConfirmed(%)`<=Cellline_cutoff,]$Percentage*(0)
  data <- data[data$Percentage!=0,]
  
  Significant_DurgCombMatDF <- dcast(data, Drug2~Drug1, value.var = 'Percentage', fill = 0,fun.aggregate = mean)
  rownames(Significant_DurgCombMatDF) <- Significant_DurgCombMatDF$Drug2
  Significant_DurgCombMatDF$Drug2 = NULL
  Significant_DurgCombMat<-as.matrix(Significant_DurgCombMatDF)
  
  col_fun = colorRamp2(c(0,20, 30, 40), c("white",col_low,col_mid,col_high))
  ht<-Heatmap(Significant_DurgCombMat,name = "Percent",col = col_fun, 
              column_title = "Drug 1", 
              row_title = "Drug 2",
              layer_fun = function(j, i, x, y, width, height, fill) {
                v = pindex(Significant_DurgCombMat, i, j)
                l = v != 0
                grid.text(sprintf("%.1f", v[l]), x[l], y[l], gp = gpar(fontsize = 12))
                #grid.text(sprintf("%.1f", pindex(Significant_DurgCombMat, i, j)), x, y, gp = gpar(fontsize = 12))
                
              })
  
  ht
  draw(ht)
}
## Circle Plot ----------------
plotCircle <- function(CirclePlot_DF,
                       DrugNameDF,
                       title_name
                        ){
  chordDiagram(CirclePlot_DF,annotationTrack = "grid",order=DrugNameDF$Drugname, grid.col = DrugNameDF$Colorindex,
               transparency = 0.5,
               big.gap = 10, small.gap = 1)
  circos.track(track.index = 1, panel.fun = function(CirclePlot_DF, DrugNameDF)  {
    xlim = get.cell.meta.data("xlim")
    xplot = get.cell.meta.data("xplot")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    
    if(abs(xplot[2] - xplot[1]) < 10) {
      circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
                  niceFacing = TRUE, adj = c(0, 1),cex = 0.7)
    } else {
      circos.text(mean(xlim), ylim[1], sector.name,  facing = "outside",niceFacing = TRUE,adj = c(0.5, 4),cex = 0.7)
    }
    # circos.text(mean(xlim), ylim[1], sector.name,  facing = "outside",niceFacing = TRUE,adj = c(0.5, 5),cex = 0.7)
  }, 
  bg.border = NA) 
  title(title_name, cex = 0.6)
}

## UpsetPlot ----------------
plotUpset <- function(drugDF,
                      Predict_response,
                      title_name
                      ){
  
  if(missing(Predict_response)){
    Predict_response <- "Positive"
    drugDF<- drugDF[grepl("^sensitivity$|^sensitivity/response$|^response$",drugDF$Predicts,ignore.case = T),]  # filter response drugs 
  } else if(Predict_response == "all"){
    drugDF <- drugDF
  } else if(Predict_response == "Positive"){
    drugDF<- drugDF[grepl("^sensitivity$|^sensitivity/response$|^response$",drugDF$Predicts,ignore.case = T),]  # filter response drugs 
  } else if(Predict_response == "Negative"){
    drugDF<- drugDF[grepl("^Adverse Response$|^resistance$|^no response$",drugDF$Predicts,ignore.case = T),]   # filter response drugs 
  }
  # drugDF<- drugDF[drugDF$Predicts %in% c("sensitivity","response","sensitivity/response"),]  # filter response drugs 
  drugDF <- drugDF[c("Sample_ID","Classified_Drug_Name")] 
  drugDF <- unique(drugDF)
  drugDF <- na.omit(drugDF)
  drugDF$Drug_Occur <- 1
  reshaped_drugDF<-dcast(drugDF, Sample_ID ~ Classified_Drug_Name, value.var = 'Drug_Occur', fill = 0,fun.aggregate = sum)
  
  
  frequency_drugDF<-unique(drugDF)
  frequency_drugDF <- as.data.frame(table(frequency_drugDF$Classified_Drug_Name))
  colnames(frequency_drugDF)<-c("drug","Freq")
  frequency_drugDF$Percentage <- round(frequency_drugDF$Freq /length(unique(drugDF$Sample_ID)),digits = 4)*100
  frequency_drugDF <- frequency_drugDF[order(-frequency_drugDF$Freq),]
  top_frequency_drugDF<-head(frequency_drugDF, n = 30)
  topFrequencyDrugs<-as.character(top_frequency_drugDF$drug)
  
  reshaped_drugDF<-reshaped_drugDF[,c("Sample_ID",topFrequencyDrugs)]
  m = make_comb_mat(reshaped_drugDF,top_n_sets = 15, mode = "intersect")
  m_combination <- m[comb_degree(m) == 2]
  
  p<-UpSet(m_combination,comb_order = order(-comb_size(m_combination)), column_title = paste0(title_name," ",Predict_response," UpsetPlot"))
  p
  
}
## Alluvial ----------------
plotAlluvial <- function(drugDF,
                         classified_drugName,
                         Alluvial_title){
  drugDF$level <- factor(drugDF$level,levels=c("B3","A3","B2","A2","B1","A1"),ordered=TRUE)
  specificDrugDF <- drugDF[drugDF$Classified_Drug_Name == classified_drugName, ]
  specificDrugDF <- na.omit(specificDrugDF)
  specificDrugDF <- unique(specificDrugDF[c("level","Classified_Drug_Name","Gene","Pat.Var","Sample_ID")]) 

  p <- ggplot(specificDrugDF, aes(
  axis1 = Classified_Drug_Name, axis2 = Gene,axis3= Pat.Var, axis4 = Sample_ID))+
  geom_alluvium(aes(fill = level),
                width = 0, knot.pos = 0, reverse = FALSE)+
  geom_stratum(width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),reverse = FALSE)+
  scale_x_discrete(limits = c("Drug", "Gene","Pat.Var","Sample_ID"), expand = c(.05, .05)) +
  ggtitle(paste(Alluvial_title,"alluvial plot")) 
  p
  #print(p)
}

## Barplot ----------------
plotfrequencyBarplot <- function(drugDF,
                                 Barplot_title
){
  
  drugDF <- drugDF[c("Sample_ID","Classified_Drug_Name")] 
  frequency_drugDF<-unique(drugDF)
  frequency_drugDF <- as.data.frame(table(frequency_drugDF$Classified_Drug_Name))
  colnames(frequency_drugDF)<-c("drug","Freq")
  frequency_drugDF$Percentage <- round(frequency_drugDF$Freq /length(unique(drugDF$Sample_ID)),digits = 4)*100
  frequency_drugDF <- frequency_drugDF[order(-frequency_drugDF$Freq),]
  frequency_drugDF$Drugs_frequency <- paste0(frequency_drugDF$drug,"(",frequency_drugDF$Percentage,")")
  #,fill=Drugs_frequency
  p<-ggplot(data=frequency_drugDF, aes(x=reorder(drug,Percentage), y=Percentage,order=Percentage)) +
    geom_text(aes(label=Percentage), hjust=-0.3, size=2.5) + 
    geom_bar(color="blue",stat="identity")
  p <- p + coord_flip() + ggtitle(paste(Barplot_title,"Single Drug Frequency")) + xlab("Classified_Drugs") + ylab("Percentage") 
  # +   guides(fill=guide_legend(nrow=32,byrow=TRUE))
  p
  
}