# library(shiny)
# library(pheatmap)
# library(UpSetR)
# library(readxl)
load("data.rda")
load("expression_browser.dra")


#ee1=read_xlsx("count_matrix.xlsx", sheet=1)
#ee2=read_xlsx("count_matrix.xlsx", sheet=1)
# apc= lapply(1:6, function(x)
 #    read_xlsx("APC.xlsx", sheet=x))
# tt = read_xls("Metadata.xls",sheet=1,col_names = F)
#save(tt, ee1, ee2, apc, file= "data.rda")



#apc[[i]] will call the i-th sheet of acp file
#apc[[1]]$`FDR p-value`< 0.05------ it will result in TRUE if this condition satisfies
#apc[[1]]$Name[TRUE]-----it will print all values
#apc[[1]]$Name[apc[[1]]$`FDR p-value`< 0.05------- will only print values for which the condition is true
# as every name in ee excel file starts like "gene : <name>" apc sheets are also modified accordinly 
         #to perform the match function
#match function returns the position of the 2 vectors where their values/ char are matched.
mm = match( paste0("gene:",apc[[1]]$Name[apc[[1]]$`FDR p-value`< 0.05]), ee1$Name)     #input$pvalue=0.05


#colids will contain the column numbers containing the "RKPM" word in ee1 file
#its output------ [1] 22 26 30 34 38 42 46 50 54 58 62 66
colids = grep("RPKM",colnames(ee1))


#a new dataframe is defined below
#this dataframe contains 12 columnns(i.e lenght(colids)=12) and length(mm) number of rows
#as.matrix(ee1[mm, colids])------- dataframe having the values of ee1 with respective to columns and rows
newdf = as.matrix( ee1[mm,colids])

#coln contains the names of all columns of file count_matrix having RKPM
coln = colnames(ee1)[colids]

#the following lines allot the respective names to columns and rows
rownames(newdf) = apc[[1]]$Name[apc[[1]]$`FDR p-value`< 0.05]
colnames(newdf) = sapply( coln, function(x) strsplit(x,"_sequence")[[1]][[1]])

#the following code checks if the rows are unique
goodid = setdiff( rownames(newdf), rownames(newdf)[is.na(newdf[,1])])
newdf = newdf[which(rownames(newdf) %in% goodid),]

dict = unlist(tt[,2])
names(dict) = sapply( unlist(tt[,1]), function(x) strsplit(x,"_sequence")[[1]][[1]] )
colnames(newdf) = dict[colnames(newdf)]
k<-apc[[1]][5]



pheatmap::pheatmap(log10(newdf+0.1),fontsize_row = 2, cluster_cols = FALSE, cluster_rows = T,
                   annotation_names_col = TRUE, annotation_names_row = TRUE,
                   # ((min(k, na.rm = T): max(k, na.rm = T))),
                   labels_col = c("endosperm 1 ", "endosperm 2", "endosperm 3",
                   "pericarp 1", "pericaprp 2", "pericarp 3",
                   "embryo 1", "embryo 2", "embryo 3",
                   "aleurone 1", "aleurone 2", "aleurone 3"), 
                   fontsize_col = 14, angle_col = 90, treeheight_row = 0, treeheight_col = 0,
                   main = "HEATMAP")
      
    
