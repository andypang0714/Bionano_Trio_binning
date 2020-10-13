readbnx2<-function(filename,mode="bz2") {
  bnx <- NULL
  try({
    
    con <- pipe(paste0("head -n 2000 ",filename," | grep -E '^0' "))
    bnx <- read.table(con,header=FALSE,sep="\t",stringsAsFactors=FALSE,nrows=2000)
    #colnames(bnx)<-c("LabelChannel", "MapID", "Length", "Avg_Intensity", "SNR")
    colnames(bnx)<-c("LabelChannel", "MoleculeID", "Length", "Avg_Intensity", "SNR", "Number_of_Labels", "OriginalMoleculeID", "ScanNumber", "ScenDirection", "ChipId", "Flowcell", "RunId", "GlobalScanNumber")
    class(bnx$MoleculeID) <- "character"
    #class(bnx$MapID) <- "character"
    
    con2 <- pipe(paste0("grep -E '^0' ",filename))
    
    bnx <- as.data.frame(scan(con2,what=bnx,comment.char="#",quiet=TRUE),stringsAsFactors=F)
    
    close(con2)
  })
  return(bnx)
}
