## Generates histograms like PlotViewer for each user-selected parameter




#########################################################
### Installing and loading required packages
#########################################################

if (!require("svDialogs")) {
  install.packages("svDialogs", dependencies = TRUE)
  library(svDialogs)
}

if (!require("flowCore")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("flowCore")
}

if (!require("tidyverse")) {
  install.packages("tidyverse", dependencies = TRUE)
  library(tidyverse)
}

if (!require("reshape2")) {
  install.packages("reshape2", dependencies = TRUE)
  library(reshape2)
}

if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}

if (!require("tcltk")) {
  install.packages("tcltk", dependencies = TRUE)
  library(tcltk)
}

#########################################################
### Script starts here
#########################################################

# Clear environment
rm(list = ls(all = TRUE))

# Data Import from file chosen by user

#library(svDialogs) # Moved to top 
# Get user input for file
testfile<-dlg_open()
# Convert to string value
testfile <- capture.output(testfile)[7]
{
  
  if ((testfile)=="character(0)")
    stop("File input cancelled")
  
  #Remove invalid characters from file input location
  testfile <- gsub("[\"]","",testfile)
  testfile<-substring (testfile,5)
  
  #Set file and directory
  filename <- basename (testfile)
  dir <- dirname (testfile)
  
  # Set working directory accoding to file chosen
  setwd(dir)
  
  #library(flowCore) #Moved to top
  
  # this read.FCS() function imports the flow data:
  raw_fcs<-read.FCS(filename, alter.names = TRUE)
  
  
  # Preparation work for arcsinh transform (columns is also used later for naming changes)
  # Create list of parameters
  columns<-colnames(raw_fcs)
  # Remove "Time" column to avoid it being transformed
  columns<-setdiff(columns,"Time")
  # Remove "Cell_Length" and Gaussians column to avoid it being transformed
  columns<-setdiff(columns,"Event_length")
  columns<-setdiff(columns,"Cell_length")
  columns<-setdiff(columns,"Center")
  columns<-setdiff(columns,"Offset")
  columns<-setdiff(columns,"Width")
  columns<-setdiff(columns,"Residual")
  ## Remove FSC and SSC
  removefscssc<-grep("FSC|SSC",columns,value=TRUE)
  columns<-columns[! columns %in% removefscssc]
  
  
  
  # Read data into a data frame
  FCSDATA <- as.data.frame(exprs(raw_fcs))
  
  ############ Optional Data Transform section
  
  #Remove comments from code lines to transform using logicle
  ## Automatically estimate the logicle transformation based on the data
  #lgcl <- estimateLogicle(raw_fcs, channels = c(columns))
  ## transform  parameters using the estimated logicle transformation
  #raw_fcs_trans <- transform(raw_fcs, lgcl)
  # Load into data frame
  #FCSDATA <- as.data.frame(exprs(raw_fcs_trans))
  
  ########### End of optional Data Transform section
  
  
  
  #Remove unnecessary parameter text
  names(FCSDATA)[-1] <- sub("Di", "", names(FCSDATA)[-1])
  names(FCSDATA)[-1] <- sub("Dd", "", names(FCSDATA)[-1])
  # Create list of channel / parameter descriptions 
  params<-parameters(raw_fcs)[["desc"]]
  # Replace parameters with descriptions, keeping things like Time, Event Length unchanged
  colnames(FCSDATA)[!is.na(params)] <- na.omit(params)
  
  # Determine whether data is CyTOF or Flow by presence of FSC
  # isflow will be 0 for a CyTOF or greater than 1 if flow
  isflow <-sum(grep("FSC",colnames(FCSDATA)))
  # Determine whether data is pre CyTOF 3 (Helios) by presence of "Cell_length", rather than "Event_length"
  isCyTOF2 <-sum(grep("Cell_length",colnames(FCSDATA)))
  
  ## Remove Time, Event_Length & Gaussian Parameters
  removecolumns <- c("Event_length", "Center", "Offset", "Width", "Residual", "Cell_length")
  FCSDATA <- FCSDATA[,!(names(FCSDATA) %in% removecolumns)]
  
  
  ## Remove FSC and SSC
  # library(tidyverse) # Moved to top
  FCSDATA <- FCSDATA %>% select(-contains("FSC"))
  FCSDATA <- FCSDATA %>% select(-contains("SSC"))
  
  
  
  
  #Calculate size of dataset
  DataSizeM <- (ncol(FCSDATA)*nrow(FCSDATA))/1000000
  #Subsample if dataset is large
  if (DataSizeM>3){
    #using random 10% of original rows
    #FCSDATA <- FCSDATA[sample(nrow(FCSDATA),nrow(FCSDATA)/10),]
    #OR
    #Subsample using a number of random rows, where the number is defined by numrows
    numrows <- 5000
    FCSDATA <- FCSDATA[sample(nrow(FCSDATA),numrows),]
  }
  
  
  
  
  
  # Ask user which parameters to plot
  markerlist<-tk_select.list(colnames(FCSDATA[-1]), multiple=TRUE,title="Select Markers to plot. Hit cancel to use all.") 
  
  # If user cancels dialog box, use all markers.
  if(length(markerlist)==0 ){
    markerlist <- colnames(FCSDATA[-1])
  }
  
  # Create list of positions of the user-selected markers
  # V and $ and gsub are used to ensure grep only matches exact / full names
  marker_cols<-NULL  
  for (m in markerlist){
    m <- paste ("^",m,"$")
    m <- gsub("\\s","",m)
    marker_cols<-c(marker_cols,grep(m,colnames(FCSDATA)))
    
  }
  
  # Remove markers that are not selected but keep first column (time)
  
  FCSDATA <- FCSDATA[,c(1,marker_cols)]
  
 
  
  
  # Melt the data into a continuous table, keeping Time for all values.
  # This allows plotting all parameters using facet_wrap in the next section
  # library(reshape2) # Moved to op
  fcsmelted <- melt(FCSDATA, id.var="Time", value.name = "intensity", variable.name="parameter")
  
  
  
  #use ggplot2 to draw dot plot
  # library(ggplot2) # Moved to top
  
  ## Add a tiny value to so we can plot on log scale
  fcsmelted$intensity <- fcsmelted$intensity + 0.9
  
  ## Generate log width bins
  ## Starting at the 0.9 minimum we set previously, ending at max intensity found in data.
  lseq <- function(from=1, to=100000, length.out=6) {
    exp(seq(log(from), log(to), length.out = length.out))
  }
  bins <- c(seq(0.9,99,1),lseq(100,max(fcsmelted$intensity),256))
  
  ## Generate x-axis labels
  xaxislabels<-c(0,10^seq(0,round(log10(max(fcsmelted$intensity)),0)))
 
  
  ## Plot histograms of intensity
  ggplot(fcsmelted, aes(x=intensity, alpha=0.9, fill=parameter)) +
    geom_histogram(bins=length(bins),breaks=bins)+
    scale_x_continuous(trans="log10",breaks=xaxislabels,labels=scales::scientific)+
    coord_cartesian(xlim=c(0.9, max(fcsmelted$intensity)))+
    # Repeat for all parameters...
    #facet_wrap("parameter")+
    # or with free scales
    facet_wrap("parameter",scales="free") +
    # Hide legend, make text smaller
    theme(legend.position = "none",axis.text.x = element_text(size=8),axis.text.y = element_text(size=8)) +
    ggtitle(filename)
    
  
    
    
  
} # End of file cancel loop

if ((testfile)=="character(0)"){
  stop("File selection cancelled")}




