getqrycmapcoor<-function(labelstart, rcmap, xmap, refid){
  
  r<-rcmap
  x<-xmap
  
  refstart<-r[r$SiteID == labelstart& r$CMapId==refid,]
  refstartpos<-refstart$Position
  refstartmolcount<-refstart$Occurrence
  
  # refend<-r[r$SiteID == labelend,]
  # refendpos<-refend$Position
  # refendmolcount<-refend$Occurrence
  x<-x[x$RefcontigID == refid,]
  withalignments<-x[grepl(labelstart, x$Alignment),]
  #withalignmentsend<-x[grepl(labelend, x$Alignment),]
  
  correspondingmolcoor<-function(rownum, reflabelnum){
    print(rownum)
    permolalignment<-withalignments[rownum,]
    eachlabel<-strsplit(permolalignment$Alignment, split="[)(]")
    wantedqrylabels<-eachlabel
    #wantedqrylabels<-eachlabel[grepl(reflabelnum, eachlabel)]
    #print(length(wantedqrylabels))
    
    getperwantedqrylabel<-function(num, reflabelnum){
      #print(num)
      perwantedqrylabels<-wantedqrylabels[[num]]
      #print("pass1")
      the_align<-perwantedqrylabels[grep(reflabelnum, perwantedqrylabels)]
      the_align<-unlist(strsplit(the_align, split=". "), recursive=F)
      peralign<-strsplit(the_align, split=",")
      findcorrespond<-data.frame(matrix(unlist(peralign), length(peralign), byrow =T))
      good<-which(findcorrespond$X1 == reflabelnum)
      #peralign<-data.frame(matrix(unlist(peralign), nrow=length(peralign), byrow=T))
      if(length(good)>0){
        the_align<-the_align[good]
        qrylabel<-do.call(rbind, strsplit(the_align, split=","))[,2]
        #wantedqrylabel<-data.frame(strsplit(perwantedqrylabels, ","))
        
        print(paste0(length(qrylabel), " found labels"))
        #print("pass2")
        return(qrylabel)
      } else {
        return("NA")
      }
    }
    #print(length(wantedqrylabels))
    qrylabels<-lapply(seq(1, length(wantedqrylabels), 1), getperwantedqrylabel, reflabelnum)
    return(qrylabels)
  }
  
  qrycoor<-lapply(seq(1, nrow(withalignments), 1), correspondingmolcoor, labelstart)
  # qryends<-lapply(seq(1, nrow(withalignmentsend), 1), correspondingmolcoor, refend, withalignmentsend)
  #print(paste0(nrow(withalignments), " alignments"))
  #print(paste0(length(qrycoor)," qry coordinates"))
  
  
  result<-data.frame(withalignments$QryContigID, do.call(rbind, qrycoor))
  result$refPos<-refstartpos
  result$Qrylabel<-unlist(result$Qrylabel)
  colnames(result)<-c("QryContigID", "Qrylabel", "refpos")
  result$QryContigID<-as.numeric(as.character(result$QryContigID))
  result$Qrylabel<-unlist(result$Qrylabel)
  result<-result[!result$Qrylabel=="NA",]
  return(result)
}


