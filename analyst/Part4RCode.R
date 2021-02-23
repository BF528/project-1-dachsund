library(readr)
setwd("/projectnb/bf528/users/dachshund/project_1")
data <- read_csv("programmer/ComBat_adjusted_data.csv")

setwd("/projectnb/bf528/users/dachshund/project_1/analyst") #gotta store things here

data <- as.data.frame(data)#makes data a data.frame if it isn't already
rownames(data) <- data[,1]
data <- data[,-1]
#Filter 1
#Apply filter 1, gene expression in >=20% of pts @ > log2(15)
apply1 <- function(f1){
  f1_threshold <- log2(15) #threshold value; constant
  filter1 <- f1[
    apply(
      f1,
      1, #Indicates by rows
      function(y) 
        sum( #tallies the number of cells that:
          y > f1_threshold) #have cell values over threshold value
      / length(y)) #divides by the number of cells in row to get the percentage of cells that exceed threshold
    >= 0.2, #Only keep is %exceeds is greater than 0.2
  ]
}

filter1 <- apply1(data)

#filter2
#determine which  set is going into this function

apply2 <- function(f2){
  #Finds 
  dof = length(f2)-1 #set degrees of freedom
  f2.rv <- apply(data[,-1], MARGIN=1, FUN=var, na.rm = TRUE) #find all row variance using data rather than filter1 findings
  f2.rv.median <- median(f2.rv) #gets the median row variance
  
  #Create function to test for variance
  test <- function(row){(var(row)-f2.rv.median)^2/f2.rv.median}
  filter2.1 <- f2[
    apply
    (
      f2,#skip first row since it is the column names
      1,#specify row
      FUN = test ) < qchisq(0.01/2, df = (dof), lower.tail = F) #checks if variance is on one end
    ,]
  
  filter2.2 <- f2[
    apply(
      f2,#skip first row since it is the column names
      1,#specify row
      FUN = test) > qchisq(1-0.01/2, df=dof, lower.tail=F) #checks for variance on the other end
    ,]
  
  filter2.1.1 <- filter2.1[!(rownames(filter2.1) %in% rownames(filter2.2)),] #get values unique to upper end
  filter2.2.1 <- filter2.2[!(rownames(filter2.2) %in% rownames(filter2.1)),] #get values unique to lower end
  #we want values that apply to one filter OR the other, not both
  filter2 <- rbind(filter2.1.1, filter2.2.1) #merge filtered2.2/1 values together 
}
filter2 <- apply2(filter1)
write.csv(filter2, "step4_5.csv", 
          quote = FALSE,
          row.names = TRUE)



#Filter 3
apply3 <- function(f3){
  cov.thresh <- 0.186 #Threshold
  #create coefficient of variation function
  cv <- function(row){sd(row)/mean(row)}
  f3[apply(
    f3[,-1], #skip first row since it is the column names
    MARGIN=1, #specify row
    FUN=cv)> cov.thresh #coefficient of variance is sd/mean. make sure it is above threshold
    ,]}
filter3 <- apply3(filter2) #run for filter 3
#Write out file for step 4.4
write.csv(filter3, "step4_4.csv",
          quote = FALSE,
          row.names = TRUE)
