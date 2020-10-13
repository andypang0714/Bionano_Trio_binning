#utility functions to read standard files
# readCMap <- function(filename) {
#     cmap <- NULL
#     try( {
#
#         #message(filename)
#         cmap <- read.table(filename, comment.char="#", header=FALSE, stringsAsFactors=FALSE, nrows=50, fill = TRUE) #adjust nrows as needed
#         cmap <- as.data.frame(scan(filename, what=cmap, comment.char="#", fill = TRUE, quiet=TRUE))
#
#         if (ncol(cmap) == 9) {
#             colnames(cmap) <- c(  "CMapId", "ContigLength", "NumSites",
#                                   "SiteID", "LabelChannel", "Position",
#                                   "StdDev", "Coverage", "Occurrence"
#             )
#         } else if (ncol(cmap) == 11) {
#             colnames(cmap) <- c(  "CMapId", "ContigLength", "NumSites",
#                                   "SiteID", "LabelChannel", "Position",
#                                   "StdDev", "Coverage", "Occurrence",
#                                   "GmeanSNR", "lnSNRsd" #new columns
#             )
#         } else if (ncol(cmap) == 12) {
#         colnames(cmap) <- c(  "CMapId", "ContigLength", "NumSites",
#                               "SiteID", "LabelChannel", "Position",
#                               "StdDev", "Coverage", "Occurrence",
#                               "GmeanSNR", "lnSNRsd", #new columns
#                               "SNR" #new columns
#             )
#         }
#     } ) # try
#     return(cmap)
# } #readCMap
#
# readXMap <- function(filename) {
#     xmap <- NULL
#     try( {
#
#         #message(filename)
#         xmap <- read.table(filename, comment.char="#", header=FALSE, stringsAsFactors=FALSE, nrows=50)
#         xmap <- as.data.frame(scan(filename, what=xmap, comment.char="#", quiet=TRUE))
#
#         colnames(xmap) <- c( "XmapEntryID", "QryContigID", "RefcontigID",
#                              "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos",
#                              "Orientation", "Confidence", "HitEnum"
#                             )
#     },
#     silent = TRUE
#     ) #try
#     return(xmap)
# } #readXMap
#
readSMap <- function( filename ) {
    smap <- NULL;
    try ( {

        #message(filename)
        smap <- read.table(filename, comment.char="#", header=FALSE, stringsAsFactors=FALSE, nrows=1000000)

        if (ncol(smap) == 12) {
            colnames( smap ) <- c( "SmapEntryID", "QryContigID", "RefcontigID1", "RefcontigID2",
                                   "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos",
                                   "Confidence", #check for column for orientation
                                   "Type", "XmapID1", "XmapID2"
            )
        } else if (ncol(smap) == 13) {
            #01222015 EL
            colnames( smap ) <- c( "SmapEntryID", "QryContigID", "RefcontigID1", "RefcontigID2",
                                   "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos",
                                   "Confidence",
                                   "Type", "XmapID1", "XmapID2","LinkID"
            )

#             colnames( smap ) <- c( "SmapEntryID", "QryContigID", "RefcontigID1", "RefcontigID2",
#                                    "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos",
#                                    "Orientation", "Confidence", #check for column for orientation
#                                    "Type", "XmapID1", "XmapID2"
#             )
        } else if (ncol(smap)==17) {
            #03242015 EL
            colnames(smap) <- c("SmapEntryID","QryContigID","RefcontigID1","RefcontigID2",
                                "QryStartPos","QryEndPos","RefStartPos","RefEndPos",
                                "Confidence","Type","XmapID1","XmapID2",
                                "LinkID","QryStartIdx","QryEndIdx","RefStartIdx","RefEndIdx")
        } else if (ncol(smap)==18) {
          colnames(smap) <- c("SmapEntryID","QryContigID","RefcontigID1","RefcontigID2",
                              "QryStartPos","QryEndPos","RefStartPos","RefEndPos",
                              "Confidence","Type","XmapID1","XmapID2",
                              "LinkID","QryStartIdx","QryEndIdx","RefStartIdx","RefEndIdx",
                              "Zygosity")
        } else if (ncol(smap)==20) {
            colnames(smap) <- c("SmapEntryID","QryContigID","RefcontigID1","RefcontigID2",
                                "QryStartPos","QryEndPos","RefStartPos","RefEndPos",
                                "Confidence","Type","XmapID1","XmapID2",
                                "LinkID","QryStartIdx","QryEndIdx","RefStartIdx","RefEndIdx",
                                "Zygosity","Genotype","GenotypeGroup")
        } else if (ncol(smap)==21) {
          colnames(smap) <- c("SmapEntryID","QryContigID","RefcontigID1","RefcontigID2",
                              "QryStartPos","QryEndPos","RefStartPos","RefEndPos",
                              "Confidence","Type","XmapID1","XmapID2",
                              "LinkID","QryStartIdx","QryEndIdx","RefStartIdx",
                              "RefEndIdx","Zygosity","Genotype","GenotypeGroup",
                              "RawConfidence")
        } else if (ncol(smap)==24) {
          colnames(smap) <- c("SmapEntryID","QryContigID","RefcontigID1","RefcontigID2",
                              "QryStartPos","QryEndPos","RefStartPos","RefEndPos",
                              "Confidence","Type","XmapID1","XmapID2",
                              "LinkID","QryStartIdx","QryEndIdx","RefStartIdx",
                              "RefEndIdx","Zygosity","Genotype","GenotypeGroup",
                              "RawConfidence",
                              "RawConfidenceLeft","RawConfidenceRight","RawConfidenceCenter")
        }
    },
    silent = TRUE #added 03202014 to suppres read smap error there are no lines
    ) #try
    return(smap)
} #readSMap

#improved functions
#read a .map file
readmap <- function(filename) {
    map <- NULL
    try( {

        con1 <- pipe(paste0("cut -f1-18 ",filename))
        #map <- read.table(con1, skip=8, header=FALSE, stringsAsFactors=FALSE, fill=TRUE, col.names=1:18, nrows=50) #read in partial file
        map <- read.table(con1, skip=6, header=FALSE, stringsAsFactors=FALSE, fill=TRUE, col.names=1:18, nrows=50) #20150324 EL - format change

        colnames(map) <- c("MappedMoleculeId", "MoleculeId", "MoleculeIndex", "ContigId",
                           "Score", "Zscore", "Direction", "StartLocation",
                           "EndLocation", "StartMatchLocation", "EndMatchLocation", "DetectedLabelCount",
                           "TruePositiveLableCount", "FalsePositiveLableCount", "FalseNegativeLabelCount", "Log10Pvalue",
                           "LeftEndChim","RightEndChim")
        class(map$MoleculeId) <- "character" #change MoleculeId to "character" to prevent overflow

        con2 <- pipe(paste0("cut -f1-18 ",filename))
        #map <- as.data.frame(scan(con2, what=map, skip=8, quiet=TRUE),stringsAsFactors=F) #read entire file with scan
        map <- as.data.frame(scan(con2, what=map, skip=6, quiet=TRUE),stringsAsFactors=F) #20150324 EL - format change
        close(con2)

    } ) # try
    return(map)
} #readmap

#read a .map file (11132014 - RStudio Server version)
readmap1 <- function(filename) {
  map <- NULL
  try( {

    con1 <- pipe(paste0("head -n 50 ",filename," | cut -f1-18")) #cut has to finish what it's doing first
    #map <- read.table(con1, skip=8, header=FALSE, stringsAsFactors=FALSE, fill=TRUE, col.names=1:18)
    map <- read.table(con1, skip=6, header=FALSE, stringsAsFactors=FALSE, fill=TRUE, col.names=1:18, nrows=50) #20150324 EL - format change

    colnames(map) <- c("MappedMoleculeId", "MoleculeId", "MoleculeIndex", "ContigId",
                       "Score", "Zscore", "Direction", "StartLocation",
                       "EndLocation", "StartMatchLocation", "EndMatchLocation", "DetectedLabelCount",
                       "TruePositiveLableCount", "FalsePositiveLableCount", "FalseNegativeLabelCount", "Log10Pvalue",
                       "LeftEndChim","RightEndChim")
    class(map$MoleculeId) <- "character"

    con2 <- pipe(paste0("cut -f1-18 ",filename))
    #map <- as.data.frame(scan(con2, what=map, skip=8, quiet=TRUE),stringsAsFactors=F)
    map <- as.data.frame(scan(con2, what=map, skip=6, quiet=TRUE),stringsAsFactors=F) #20150324 EL - format change
    close(con2)

  } ) # try
  return(map)
} #readmap

#read an .xmap file
#112414 - modified to convert type for MoleculeId; assignment of colnames moved up
readxmap <- function(filename) {
    xmap <- NULL
    try( {

        xmap <- read.table(filename, comment.char="#", header=FALSE, stringsAsFactors=FALSE, nrows=50, sep="\t")

        if (ncol(xmap) == 10) {
          colnames(xmap) <- c("XmapEntryID", "QryContigID", "RefcontigID",
                              "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos",
                              "Orientation", "Confidence", "HitEnum")
        } else if (ncol(xmap) == 14) { #121514 - new columns
          colnames(xmap) <- c("XmapEntryID", "QryContigID", "RefcontigID",
                              "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos",
                              "Orientation", "Confidence", "HitEnum",
                              "QryLen","RefLen","LabelChannel","Alignment") #121514
        } else if (ncol(xmap) == 16) { #20160725 - new columns
          colnames(xmap) <- c("XmapEntryID", "QryContigID", "RefcontigID",
                              "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos",
                              "Orientation", "Confidence", "HitEnum",
                              "QryLen","RefLen","LabelChannel","Alignment",
                              "MapWt","MaxOutlierKb") #20160725
        }

        class(xmap$QryContigID) <- "character" #change QryContigID to "character" to prevent overflow

        xmap <- as.data.frame(scan(filename, what=xmap, comment.char="#", quiet=TRUE),stringsAsFactors=FALSE)

    } ) #try
    return(xmap)
} #readxmap

#read a .cmap file
readcmap <- function(filename) {
    cmap <- NULL
    try( {

        cmap <- read.table(filename, comment.char="#", header=FALSE, stringsAsFactors=FALSE, nrows=50) #adjust nrows as needed

        if (ncol(cmap) == 9) {
            colnames(cmap) <- c(  "CMapId", "ContigLength", "NumSites",
                                  "SiteID", "LabelChannel", "Position",
                                  "StdDev", "Coverage", "Occurrence"
            )
        } else if (ncol(cmap) == 11) {
            colnames(cmap) <- c(  "CMapId", "ContigLength", "NumSites",
                                  "SiteID", "LabelChannel", "Position",
                                  "StdDev", "Coverage", "Occurrence",
                                  "GmeanSNR", "lnSNRsd" #new columns
            )
        } else if (ncol(cmap) == 12) {
            colnames(cmap) <- c(  "CMapId", "ContigLength", "NumSites",
                                  "SiteID", "LabelChannel", "Position",
                                  "StdDev", "Coverage", "Occurrence",
                                  "GmeanSNR", "lnSNRsd", #new columns
                                  "SNR" #new columns
            )
        } else if (ncol(cmap) == 13) { #20150921
          colnames(cmap) <- c(  "CMapId", "ContigLength", "NumSites",
                                "SiteID", "LabelChannel", "Position",
                                "StdDev", "Coverage", "Occurrence", "ChimScore",
                                "GmeanSNR", "lnSNRsd", #new columns
                                "SNR" #new columns
          )
        }

        class(cmap$CMapId) <- "character"

        cmap <- as.data.frame(scan(filename, what=cmap, comment.char="#", quiet=TRUE),stringsAsFactors=F)

    } ) # try
    return(cmap)
} #readcmap

#20151103 EL
readcmap <- function(filename) {
  cmap <- NULL
  try({

    header0 <- readLines(filename,n=50)

    if (all(grepl("#",header0))) { #if all lines are comments, return NULL
      return(NULL)
    }

    header1 <- header0[regexpr("^#h", header0)>=0] #get column names

    header_splilt <- unlist(strsplit(header1,"\\s"))
    header <- tail(header_splilt,-1)

    cmap <- read.table(filename,comment.char="#",header=FALSE,stringsAsFactors=FALSE,nrows=50) #adjust nrows as needed

    ncols <- ncol(cmap)
    header <- head(header,ncols)

    colnames(cmap) <- header

    class(cmap$CMapId) <- "character"

    cmap <- as.data.frame(scan(filename,what=cmap,comment.char="#",quiet=TRUE),stringsAsFactors=F)

  }) # try

  return(cmap)

} #readcmap

#20160308 EL
readcmap1 <- function(filename) {
  con1 <- pipe(paste0("head -n 100 ",filename," | cut -f1-13")) #cut has to finish what it's doing first
  cmap <- read.table(con1,comment.char="#",header=FALSE,stringsAsFactors=FALSE,fill=TRUE,col.names=1:13,nrows=100) #20150324 EL - format change

  potential_colnames <- c("CMapId","ContigLength","NumSites",
                          "SiteID","LabelChannel","Position",
                          "StdDev","Coverage","Occurrence","ChimScore",
                          "GmeanSNR","lnSNRsd","SNR")

  ncols <- ncol(cmap)
  cmap_colnames <- potential_colnames[1:ncols]

  colnames(cmap) <-cmap_colnames

  class(cmap$CMapId) <- "character"

  con2 <- pipe(paste0("cut -f1-13 ",filename))
  cmap <- as.data.frame(scan(con2,what=cmap,comment.char="#",quiet=TRUE),stringsAsFactors=F)
  close(con2)

  return(cmap)
}

#20151103 EL
readcmap_w_headers <- function(filename,sep="") { #output a list with the header lines and the actual content
  header0 <- readLines(filename,n=100) #find header lines
  header <- header0[regexpr("^#",header0)>=0]

  header1 <- header0[regexpr("^#h", header0)>=0] #get column names
  header_splilt <- unlist(strsplit(header1,"\\s"))
  header_colnames <- tail(header_splilt,-1)

  cmap <- read.table(filename,comment.char="#",header=FALSE,stringsAsFactors=FALSE,nrows=50) #adjust nrows as needed

  ncols <- ncol(cmap)
  header_colnames <- head(header_colnames,ncols)

  colnames(cmap) <- header_colnames

  class(cmap$CMapId) <- "character"

  cmap <- as.data.frame(scan(filename,what=cmap,comment.char="#",quiet=TRUE),stringsAsFactors=F)

  outlist <- list("header"=header,
                  "cmap"=cmap)

  return(outlist)
}

#read a .bed file
readbed <- function(filename) {
    bed <- NULL
    try( {

        bed <- read.table(filename, comment.char="#", header=FALSE, stringsAsFactors=FALSE, fill=TRUE)

        if (ncol(bed)==5) {
            colnames(bed) <- c( "chr", "start", "end", "walkID", "orientation")
        } else if (ncol(bed)==3) {
            colnames(bed) <- c( "chr", "start", "end")
        }

    } ) #try
    return(bed)
} #readbed

#20150522 EL - read an .err file
readerr <- function(filename) {
  err <- NULL
  try({

    err <- read.table(filename,comment.char="#",header=TRUE,stringsAsFactors=FALSE,fill=FALSE)

    if (ncol(err)==19) {
      colnames(err) <- c("iteration","fp","fn","sf","sd","bpp","res","maps","log10LR","aligned",
                         "log10LR.aligned","bppSD","fp.rate","sr","se","label.density","resSD","mres","mresSD")

      use_columns <- c("log10LR","aligned","log10LR.aligned","label.density")
    } else {
      use_columns <- c("Log10LR..Maps.","GoodMaps","log10LR..GoodMaps.","LabelDensity..100kb.")
    }



    F <- names(err)%in%use_columns

    err1 <- tail(err,1) #last line
    err2 <- head(tail(err,2),1) #second to last line

    err1[F] <- err2[F]

    err <- err1
  }) #try
  return(err)
} #readerr
#err_file <- "/mnt/bionf_tmp/elam/20150519_ESHG/rubin_case.bnx_errEst.err"
#err <- readerr(err_file)

#20150803 - read a .align file (skipping the rows with aligned labels)
readalign <- function(filename) {
  align <- NULL
  try({

    con <- pipe(paste0("head -n 50 ",filename," | grep -E '>0' "))

    align <- read.table(con,comment.char="#",header=FALSE,stringsAsFactors=FALSE,nrows=50) #adjust nrows as needed
    #print(dim(align))

    if (ncol(align)==13) {
      #>0    AlignmentID    Mol0ID  Mol1ID  Score	CenterOffset	Overlap	Orientation	PvalueLog10	TrueOffset	TrueOverlapFraction	FileID0	FileID1
      colnames(align) <- c("lineID","AlignmentID","Mol0ID","Mol1ID",
                           "Score","CenterOffset","Overlap",
                           "Orientation","PvalueLog10","TrueOffset",
                           "TrueOverlapFraction","FileID0","FileID1")
    } else if (ncol(align)==14) {
      #>0    AlignmentID    Mol0ID  Mol1ID  Score  CenterOffset	Overlap	Orientation	PvalueLog10	TrueOffset	TrueOverlapFraction	FileID0	FileID1 Offset
      colnames(align) <- c("lineID","AlignmentID","Mol0ID","Mol1ID",
                           "Score","CenterOffset","Overlap",
                           "Orientation","PvalueLog10","TrueOffset",
                           "TrueOverlapFraction","FileID0","FileID1",
                           "Offset")
    }

    class(align$Mol0ID) <- "character"
    class(align$Mol1ID) <- "character"

    con2 <- pipe(paste0("grep -E '>0' ",filename))

    align <- as.data.frame(scan(con2,what=align,comment.char="#",quiet=TRUE),stringsAsFactors=F)
    align$lineID <- NULL #disard first column

    close(con2)
  },
  silent = TRUE) # try

  return(align)
} #readalign

#20150803 - read a .bnx file (skipping the rows with labels and QX values)
readbnx <- function(filename,mode="bz2") {
  bnx <- NULL
  try({

    con <- pipe(paste0("head -n 50 ",filename," | grep -E '^0' "))
    bnx <- read.table(con,header=FALSE,sep="\t",stringsAsFactors=FALSE,nrows=50)
    colnames(bnx)<-c("LabelChannel", "MapID", "Length", "StartLoc", "EndLoc")

    class(bnx$MapID) <- "character"

    con2 <- pipe(paste0("grep -E '^0' ",filename))

    bnx <- as.data.frame(scan(con2,what=bnx,comment.char="#",quiet=TRUE),stringsAsFactors=F)

    close(con2)
  })
  return(bnx)
}
