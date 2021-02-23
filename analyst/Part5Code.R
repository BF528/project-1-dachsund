#Imports
library(readr)
library(dplyr)

#Sets working directory for meta, just in case I do something stupid later
setwd("/project/bf528/project_1/doc")
meta <- read.csv("proj_metadata.csv") #Read the metadata

#Set Working directory for import/export datasets
setwd("/projectnb/bf528/users/dachshund/project_1/analyst") #gotta store things here

#Read files
step4.4 <- read.csv("step4_4.csv") #Data from step 4.4 (Standard data)    
step4.5 <- read.csv("step4_5.csv") #Data from step 4.5 (Because we have a biologist)
temp <- read.csv("project_1_step_5_6.csv")

#Change row names to first column value
rownames(step4.4) <- step4.4[,1]
step4.4 <- step4.4[,-1]
colnames(step4.4) <- stringr::str_extract(names(step4.4), "^[^_]+(?=_)")#trims names
#Same, but for biologist-exclusive data
rownames(step4.5) <- step4.5[,1]
step4.5 <- step4.5[,-1]
colnames(step4.5) <- stringr::str_extract(names(step4.5), "^[^_]+(?=_)")#trims names
#Change row names to first column value since meta just like that yo
rownames(meta) <- meta$geo_accession


#5.1
#Hierarchal clustering function
hc <- function(mtx){
  #Clean data and scale it
  mtx <- na.omit(mtx)
  mtx <- scale(mtx)
  mtx <- t(mtx)
  
  # Dissimilarity matrix
  d <- dist(mtx[,], method = "euclidean")
  # Hierarchical clustering
  hc1 <- hclust(d, method = "ward.D" )
}

#Run clustering and plot
cluster <- hc(step4.4)
plot(cluster , cex = 0.6, hang = -1)


#5.2
#Cut dendogram into 2 clusters 
list <- cutree(cluster, k=2)
print(as.data.frame(table(list))) #gets the number of values in each cut of the cluster
#1- 58
#2- 76

#5.3
#Create a heatmap
#clus determines if pt is subtype C3 from metadata or not
clus <- function(s){
  l2 <- rep("blue", length(s))
  for(i in 1:length(s)) {
    if(meta[colnames(s)[i], "SixSubtypesClassification"] == "C3"){l2[i] <- "red"}
  }
  l2
}
  
#Creates heatmap
heatmap <- heatmap(data.matrix(step4.4),  
                   #col = hsv(seq(0.7,1,length.out = 12),1, 1),
                   margins=c(5,10),
                   ColSideColors = clus(step4.4),
)

#5.4
#Ranking Welch T test
output_function <- function(dataframe){
  output <- data.frame( t=as.double(),
                        p_value=as.double(),
                        adj_p=as.double())
  ttest1 <- dataframe[which(meta$SixSubtypesClassification == "C3")]
  ttest2 <- dataframe[which(meta$SixSubtypesClassification != "C3")]
  
  for(name in row.names(dataframe)){
    temp <- t.test(ttest1[name,], ttest2[name,])
    output[name ,"t"] <- temp$statistic
    output[name ,"p_value"] <- temp$p.value
    output[name ,"adj_p"] <- 5.00
  }
  temp <- p.adjust(output$p_value,method="fdr")
  for(i in 1:nrow(output)){
    output[i, "adj_p"] <- temp[i]
  }
  return(output)
}

ttest <- output_function(step4.4)

#Separates by 5.2's determined cluster1/2 and then tells how many adj_p < 0.05
#ttest1 <- output_function(step4.4[which(meta$SixSubtypesClassification == "C3")])
sum(ttest$adj_p < 0.05) #Run to print out # of adj_p<0.05 to the console
#1297 outputted

#Writes to csv
write.csv(output, "project_1_step5_5_2.csv",
          quote = FALSE,
          row.names = TRUE)

#Step 5.6
#Run ranking t-test via( output_function()) on both clusters and combines outputs into 1 data frame
output5.6 <- output_function(step4.5)
sum(output5.6$adj_p < 0.05)
#22729 outputted

#Writes to csv
write.csv(output5.6, "project_1_step_5_6_2.csv",
          quote = FALSE,
          row.names = TRUE)
