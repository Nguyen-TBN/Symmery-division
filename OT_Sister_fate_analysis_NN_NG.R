library(ggplot2)
#choose the location
options(stringsAsFactors=FALSE)
setwd("D:/2 - Nguyen-folder/manuscript_CRC met paper/Figure 2-new/Organoid-Tracker analysis/raw csv files")
temp = list.files(pattern="*.csv")#this is to have them in the same file
myfiles = lapply(temp, read.csv, header = FALSE)
temp

###PLEASE check your data is similar to the other data, after loading they have a spare row 

#set the place for output data
setwd("D:/2 - Nguyen-folder/manuscript_CRC met paper/Figure 2-new/Organoid-Tracker analysis/Sister fate analysis")
source("Function_sister_fate_NG.R")
##Fill in the file you want to plot (file number)
##For comparing the diferences between two sister cells -> only check at the sisters not the mothers

for (n in 1:length(myfiles)){
  OR_name <- substr(temp[[n]],1,21)
  dt <- data.frame(myfiles[[n]])
  dt <- data.frame(dt[-1,]) #be careful this maybe not needed
  
  ##Get the cell identity out
  cell.name <- unlist(dt[1, substring(dt[1,],1,4) == "Cell"])
  
  ##get STAR red data
  chan561 <- "Channel 4"
  red.data <- dt[3:nrow(dt), dt[2,] == chan561]
  colnames(red.data) <- paste(substr(cell.name,6,nchar(cell.name)),"_C4",sep="")
  #convert data in each table into numeric (they were character)
  red.data[] <- lapply(red.data, function(x) {as.numeric(as.character(x))})
  
  ##------# MAKING LINEAGES TREE
  mom_list <- gsub("\\(daughter of ", "", dt[1,grep("daughter", dt[1,])])#get the mother
  mom_list <- gsub("\\)", "", mom_list)
  dau_list <- unlist(dt[1,grep("daughter", dt[1,])-1])#get the daughter
  dau_list <- substr(dau_list,6,nchar(dau_list))
  
  #Get mother-daughter relationships and put in list
  mds <- list() 
  for (m in 1:length(mom_list)){
    mds[[mom_list[m]]] <- unlist(matrix(dau_list[which(mom_list %in% mom_list[m])]))
  }#here it will automately join two similar mother into 1
  
  
  ##the analysis start from here
  result <- F_filter_tp_create_cut_red(mds, red.data, 20) #result[[1]] is the filtured.mds, result[[2]] is cut red.data
  f.mds <- result[[1]]
  f.red.data <- result[[2]]
  
  #Second, remove the cells that are fluctuating a lot (to do this, we only check the events 
  var.limit <- 12 #choose number of variance
  span <- 0.5 # choose number of span
  
  for (i in 1:length(f.mds)) {
    print(i)
    mom <- names(f.mds)[i]
    red1 <- f.red.data[,paste(f.mds[[i]][1],"_C4", sep="")]
    red2 <- f.red.data[,paste(f.mds[[i]][2],"_C4", sep="")]
    if (var(red1, na.rm=T) > var.limit | var(red2,na.rm=T) > var.limit) {
      result1 <- F_create_data_with_smoothed(red1, f.mds[[i]][1], span)
      result2 <- F_create_data_with_smoothed(red2, f.mds[[i]][2], span)
      data.plot <- rbind(result1, result2)
      pp <- F_making_plot(data.plot, 0,80,mom)
      print(pp)
      check <- readline(prompt="Enter ok/r/quit (r: remove, quit:to quit)")
      check <- as.character(check)
      if (check == "r"){
        f.mds[[i]] <- NA
      } 
      if (check == "quit"){
        break
      }
    }
  }
  save(f.mds,file=paste(OR_name,"total_pairs.RData", sep="")) #save this data! important
  
  #Third, check on the slope and difference for each daughter pairs: either they are different fate or just being different
  #a: different fates and different at the end
  slope.limit <- c(-0.3,0.3)
  change.limit <- 5
  span <- 0.5
  dt.mds.fate <- data.frame(matrix(nrow=0, ncol=12))
  for (i in 1:length(f.mds)) {
    if (length(f.mds[[i]])>1){
      red1 <- f.red.data[,paste(f.mds[[i]][1],"_C4", sep="")]
      red2 <- f.red.data[,paste(f.mds[[i]][2],"_C4", sep="")]
      smooth1 <- predict(loess(red1 ~ c(1:length(red1)),span=span))
      smooth2 <- predict(loess(red2 ~ c(1:length(red2)),span=span))
      change1 <- mean(tail(smooth1,10))-mean(head(smooth1,10))
      change2 <- mean(tail(smooth2,10))-mean(head(smooth2,10))
      
      diff <- mean(tail(smooth1,10))-mean(tail(smooth2,10))
      
      red1 <- red1[!is.na(red1)]/red1[!is.na(red1)][1]
      red2 <- red2[!is.na(red2)]/red2[!is.na(red2)][1]
      slope1 <- F_cal_slope(red1)
      slope2 <- F_cal_slope(red2)
      fate1 <- F_define_fate_2(change1,slope1,change.limit,slope.limit)
      fate2 <- F_define_fate_2(change2,slope2,change.limit,slope.limit)
      
      if (fate1 == fate2 & abs(diff) < 10) {
        comparison <- "similar"
      } else if (fate1 != fate2 & abs(diff) > 10) {
        comparison <- "different"
      } else {comparison <- "non-defined" }
      
      if (abs(diff) >10) {
        compare_diff <- "Yes"
      } else {compare_diff <- "No"}
      cell.data <- c(names(f.mds)[i],paste("Cell",names(f.mds)[i],"Cell",f.mds[[i]][1], sep="_"),slope1,change1,fate1,
                     paste("Cell",names(f.mds)[i],"Cell",f.mds[[i]][2], sep="_"),slope2,change2,fate2,abs(diff),comparison, compare_diff)
    }
    dt.mds.fate <- rbind(dt.mds.fate,cell.data)
  }
  colnames(dt.mds.fate) <- c("Mother","Dau_1","Slope_1","Change_1","Fate_1",
                             "Dau_2","Slope_2","Change_2","Fate_2", "Diff_btw_cell","Comparision_fate", "Comparision_diff")
  
  ##save final data and the export the final result
  write.csv(dt.mds.fate,paste(OR_name,"sister_analysis.csv", sep="_"))

  #to plot all the pairs used for quantification
  pdf(paste(OR_name,"plots.pdf", sep="_"))
  for (i in 1:length(f.mds)) {
    if (length(f.mds[[i]])>1) {
      mom <- names(f.mds)[i]
      red1 <- f.red.data[,paste(f.mds[[i]][1],"_C4", sep="")]
      red2 <- f.red.data[,paste(f.mds[[i]][2],"_C4", sep="")]
      result1 <- F_create_data_with_smoothed(red1, f.mds[[i]][1], span)
      result2 <- F_create_data_with_smoothed(red2, f.mds[[i]][2], span)
      data.plot <- rbind(result1, result2)
      pp <- F_making_plot(data.plot, 0,80,mom)
      ppp <- pp + annotate("text", x=min(which(!is.na(red1)))+5, y=70,col="black", size =5,
                           label=paste("Comparison:",
                                       dt.mds.fate$Comparision_fate[dt.mds.fate$Mother==mom],",Difference:",
                                       dt.mds.fate$Comparision_diff[dt.mds.fate$Mother==mom]))
      print(ppp)
    }
  }
  dev.off()
}


##----START the analysis from here
#Here we only check on 2 daughter cells, dont take into account mother cells



















