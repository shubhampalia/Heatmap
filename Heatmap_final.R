library(shiny)
library(pheatmap)

library(readxl)

#library(tidyverse)
load("data.rda")
load("expression_browser.dra")



# ee1=read_xlsx("count_matrix.xlsx", sheet=1)
# ee2=read_xlsx("count_matrix.xlsx", sheet=1)
# apc= lapply(1:6, function(x)
#     read_xlsx("APC.xlsx", sheet=x))
#tt = read_xls("Metadata.xls",sheet=1,col_names = F)
#save(tt, ee1, ee2, apc, file= "data.rda")
#master_data<- excel_sheets("APC.xlsx") %>% map_df(~read_xlsx("APC.xlsx", .x))

b<- excel_sheets("APC.xlsx") 
ch<- "All Comparisions"

ui<- navbarPage(
  title="DIFFERENTIALLY EXPRESSED GENES",
  
  
  
  
  tabPanel("Heatmap",
           sidebarLayout(
             sidebarPanel(
               radioButtons("pvalue", "Choose The Probability Value", choices = c(0.01, 0.05)),
               selectInput("sheet", 
                           "You are currently observing the comparision between :", 
                           choices = c(excel_sheets("APC.xlsx"), "All Comparisions"))    #excel sheets are selected and are given as choices to the users
             ),
             mainPanel(
               print("This heatmap shows the Differencially Expressed Genes with respect to Log2 Fold Change values on the x- axis
                           "),
               plotOutput("heatmap"),
               textOutput("sheetname")
             )
           )
  )
  
  
)






server<- function(input, output){
  
  #the following function just renders/ displays the line within the quoted text in print("") statement
  output$sheetname<- ({
    
    renderText({
      input$sheet
    })
  })
  
  output$heatmap<- renderPlot({
    
    #b will contain all the excel sheets according to their order, calling b[1] will print name of the first sheet of APC.xlsx
    c<- input$sheet #reactive functions like input$sheet can not be directly be compared to any other variable
    for (i in 1:6){
      if (b[i]==c)
      {
        #apc[[i]] will call the i-th sheet of acp file
        mm1  = match(apc[[i]]$Name[apc[[i]]$`FDR p-value`< input$pvalue], ee1$Name)     
        mm = match( paste0("gene:",apc[[i]]$Name[apc[[i]]$`FDR p-value`< input$pvalue]), ee1$Name)

        for (j in c(1:length(mm))){
          if (is.na(mm[j]==TRUE)){
            mm[j]=mm1[j]
          }
        }
        
        colids = grep("RPKM",colnames(ee1))
        newdf = as.matrix( ee1[mm,colids])
        coln = colnames(ee1)[colids]
        rownames(newdf) = apc[[i]]$Name[apc[[i]]$`FDR p-value`< input$pvalue]
        colnames(newdf) = sapply( coln, function(x) strsplit(x,"_sequence")[[1]][[1]])
        goodid = setdiff( rownames(newdf), rownames(newdf)[is.na(newdf[,1])])
        newdf = newdf[which(rownames(newdf) %in% goodid),]
        dict = unlist(tt[,2])
        names(dict) = sapply( unlist(tt[,1]), function(x) strsplit(x,"_sequence")[[1]][[1]] )
        colnames(newdf) = dict[colnames(newdf)]
        
        
        
        pheatmap::pheatmap(log10(newdf+0.1),fontsize_row = 2, cluster_cols = FALSE, cluster_rows = T,
                           annotation_names_col = TRUE, annotation_names_row = TRUE,
                           # ((min(k, na.rm = T): max(k, na.rm = T))),
                           labels_col = c("endosperm 1 ", "endosperm 2", "endosperm 3",
                                          "pericarp 1", "pericaprp 2", "pericarp 3",
                                          "embryo 1", "embryo 2", "embryo 3",
                                          "aleurone 1", "aleurone 2", "aleurone 3"), 
                           fontsize_col = 14, angle_col = 90, treeheight_row = 0, treeheight_col = 0,
                           main = "HEATMAP")
        
      }
    }
    
    
    if (c==ch){
      #master_data<- excel_sheets("APC.xlsx") %>% map_df(~read_xlsx("APC.xlsx", .x)) 
      #apc=master_data
      #newdf= masterdf
      
      mm12  = match(apc[[i]]$Name[apc[[i]]$`FDR p-value`< input$pvalue], ee1$Name)     
      mm2 = match( paste0("gene:",apc[[i]]$Name[apc[[i]]$`FDR p-value`< input$pvalue]), ee1$Name)
      
      for (j in c(1:length(mm2))){
        if (is.na(mm2[j]==TRUE)){
          mm2[j]=mm12[j]
        }
      }
      
      mm2 = match( paste0("gene:",master_data$Name[master_data$`FDR p-value`< input$pvalue]), ee1$Name)     
      colids2 = grep("RPKM",colnames(ee1))
      masterdf = as.matrix( ee1[mm2,colids2])
      coln2 = colnames(ee1)[colids2]
      rownames(masterdf) = master_data$Name[master_data$`FDR p-value`< input$pvalue]
      colnames(masterdf) = sapply( coln2, function(x) strsplit(x,"_sequence")[[1]][[1]])
      goodid2 = setdiff( rownames(masterdf), rownames(masterdf)[is.na(masterdf[,1])])
      masterdf = masterdf[which(rownames(masterdf) %in% goodid2),]
      dict2 = unlist(tt[,2])
      names(dict2) = sapply( unlist(tt[,1]), function(x) strsplit(x,"_sequence")[[1]][[1]] )
      colnames(masterdf) = dict2[colnames(masterdf)]
      pheatmap::pheatmap(log10(masterdf+0.1),fontsize_row = 2, cluster_cols = FALSE, cluster_rows = T,
                         annotation_names_col = TRUE, annotation_names_row = TRUE,
                         fontsize_col = 14, angle_col = 90, treeheight_row = 0, treeheight_col = 0,
                         main = "HEATMAP", labels_col = c("endosperm 1 ", "endosperm 2", "endosperm 3",
                                                          "pericarp 1", "pericaprp 2", "pericarp 3",
                                                          "embryo 1", "embryo 2", "embryo 3",
                                                          "aleurone 1", "aleurone 2", "aleurone 3"))
    }
  })
  
}

shinyApp(ui=ui, server=server)