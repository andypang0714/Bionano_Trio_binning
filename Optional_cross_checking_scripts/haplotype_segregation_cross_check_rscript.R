#!/usr/bin/env Rscript
library("optparse")
source("/home/users/elam/rscripts/paste0.R")
source("/home/users/jlee/JLwrittenscripts/get_header.R")
source("~/JLwrittenscripts/get_corresponding_qrymap_coordinate.R")

#ftfpath<-"/home/users7/jlee/data/Rockefeller_Homo_sapiens_20190826/father/Father_selected_nonhap_ES_sdb_20190830_cloud/father_selected_nonhap_GM_to_father_ori_GM_xmapUnique0_IndelNoOverlap0//exp_refineFinal1_sv/merged_smaps/exp_refineFinal1_merged_filter_inversions.smap"
#ftfxmap<-"/home/users7/jlee/data/Rockefeller_Homo_sapiens_20190826/father/Father_selected_nonhap_ES_sdb_20190830_cloud/father_selected_nonhap_GM_to_father_ori_GM_xmapUnique0_IndelNoOverlap0//exp_refineFinal1_sv/merged_smaps/exp_refineFinal1_merged.xmap"
#ftfrmap<-"/home/users7/jlee/data/Rockefeller_Homo_sapiens_20190826/father/Father_selected_nonhap_ES_sdb_20190830_cloud/father_selected_nonhap_GM_to_father_ori_GM_xmapUnique0_IndelNoOverlap0//exp_refineFinal1_sv/merged_smaps/exp_refineFinal1_merged_r.cmap"
#ftfqmap<-"/home/users7/jlee/data/Rockefeller_Homo_sapiens_20190826/father/Father_selected_nonhap_ES_sdb_20190830_cloud/father_selected_nonhap_GM_to_father_ori_GM_xmapUnique0_IndelNoOverlap0//exp_refineFinal1_sv/merged_smaps/exp_refineFinal1_merged_q.cmap"
#mtfpath<-"/home/users7/jlee/data/Rockefeller_Homo_sapiens_20190826/mother/Mother_selected_nonhap_ES_sdb_20190830_cloud/mother_selected_nonhap_GM_to_father_ori_GM_xmapUnique0_IndelNoOverlap0//exp_refineFinal1_sv/merged_smaps/exp_refineFinal1_merged_filter_inversions.smap"
#mtfxmap<-"/home/users7/jlee/data/Rockefeller_Homo_sapiens_20190826/mother/Mother_selected_nonhap_ES_sdb_20190830_cloud/mother_selected_nonhap_GM_to_father_ori_GM_xmapUnique0_IndelNoOverlap0//exp_refineFinal1_sv/merged_smaps/exp_refineFinal1_merged.xmap"
#outputfile<-"/home/users7/jlee/data/Rockefeller_Homo_sapiens_20190826/father/Father_selected_nonhap_ES_sdb_20190830_cloud/father_selected_nonhap_GM_to_father_ori_GM_xmapUnique0_IndelNoOverlap0/cross_check_i1_result.txt"
#sizethres<-2000
#buffer<-10000
#types<-c("insertion", "deletion")
###zygosity, size filter at the mtf. indels only. Duplication may come later
###Not including indels that involve more than 1 match group e.g. those XmapID1 != XmapID2)
###What if mom has a large indels too that involves 2 matchgroups?

parser <- OptionParser()
parser <- add_option(parser, c("-f", "--ftfpath"), help="Parent1 allele smap for against Parent1 original")
parser <- add_option(parser, c("-x", "--ftfxmap"), help="Parent1 allele xmap for against Parent1 original")
parser <- add_option(parser, c("-r", "--ftfrmap"), help="Parent1 allele rcmap for against Parent1 original")
parser <- add_option(parser, c("-q", "--ftfqmap"), help="Parent1 allele qcmap for against Parent1 original")
parser <- add_option(parser, c("-m", "--mtfpath"), help="Parent2 allele smap for against Parent1 original")
parser <- add_option(parser, c("-y", "--mtfxmap"), help="Parent2 allele xmap for against Parent1 original")
parser <- add_option(parser, c("-o", "--outputfile"), help="cross check cuts file output")
parser <- add_option(parser, c("-s", "--sizethres"), help="sizethreshold", default= 2000)
parser <- add_option(parser, c("-b", "--buffer"), help="buffer for regions to be checked for SV", default=10000)
#parser <- add_option(parser, c("-t", "--types"), help="SV types to trigger cross checking", default=c("insertion", "deletion"))
args<-parse_args(parser)

#args<-parser$parse_args()

ftfpath<-args$ftfpath
ftfxmap<-args$ftfxmap
ftfrmap<-args$ftfrmap
ftfqmap<-args$ftfqmap
mtfpath<-args$mtfpath
mtfxmap<-args$mtfxmap
outputfile<-args$outputfile
sizethres<-args$sizethres
buffer<-args$buffer
#types<-args$types
types<-c("insertion", "deletion")

cross_checking<-function(ftfpath, ftfxmap, ftfrmap, ftfqmap, mtfpath, mtfxmap, outputfile, sizethres, buffer, types){
ftf<-readtable(ftfpath)
mtf<-readtable(mtfpath)
ftf$SVsize<-as.numeric(as.character(ftf$SVsize))
ftf$QryContigID<-as.numeric(as.character(ftf$QryContigID))
mtf$SVsize<-as.numeric(as.character(mtf$SVsize))
mtf$QryContigID<-as.numeric(as.character(mtf$QryContigID))
ftf$RefcontigID1<-as.numeric(as.character(ftf$RefcontigID1))
mtf$RefcontigID1<-as.numeric(as.character(mtf$RefcontigID1))
ftf$RefStartPos<-as.numeric(as.character(ftf$RefStartPos))
mtf$RefStartPos<-as.numeric(as.character(mtf$RefStartPos))
ftf$RefEndPos<-as.numeric(as.character(ftf$RefEndPos))
mtf$RefEndPos<-as.numeric(as.character(mtf$RefEndPos))
mtx<-readtable(mtfxmap)
ftx<-readtable(ftfxmap)
ftx$RefContigID<-as.numeric(as.character(ftx$RefContigID))
ftx$RefcontigID<-ftx$RefContigID
ftx$Alignment<-as.character(ftx$Alignment)
ftx$QryContigID<-as.numeric(as.character(ftx$QryContigID))
mtx$RefContigID<-as.numeric(as.character(mtx$RefContigID))
mtx$QryContigID<-as.numeric(as.character(mtx$QryContigID))
ftx$RefStartPos<-as.numeric(as.character(ftx$RefStartPos))
mtx$RefStartPos<-as.numeric(as.character(mtx$RefStartPos))
ftx$RefEndPos<-as.numeric(as.character(ftx$RefEndPos))
mtx$RefEndPos<-as.numeric(as.character(mtx$RefEndPos))
ftfq<-readtable(ftfqmap)
ftfq$SiteID<-as.numeric(as.character(ftfq$SiteID))
ftfq$`#h CMapId`<-as.numeric(as.character(ftfq$`#h CMapId`))
ftfr<-readtable(ftfrmap)
ftfr$`#h CMapId`<-as.numeric(as.character(ftfr$`#h CMapId`))
ftfr$CMapId<-ftfr$`#h CMapId`
ftfr$SiteID<-as.numeric(as.character(ftfr$SiteID))

ftfs<-ftf[ftf$SVsize>sizethres & !(ftf$Zygosity == "homozygous") & ftf$Type %in% types,]
ftfs<-ftfs[ftfs$Zygosity == "heterozygous" & (as.numeric(as.character(ftfs$XmapID1)) == as.numeric(as.character(ftfs$XmapID2))),]

chrtoberun<-ftfs$RefcontigID1[!duplicated(ftfs$RefcontigID1)]
#chr<-191
print(table(ftfs$RefcontigID1))

get_SV_alt_sample_at_region<-function(chr){
  ###Loop number of ref with SV
  print(paste0("chr", chr))
  ftfp<-ftfs[ftfs$RefcontigID1 == chr,]
  npchr<-nrow(ftfp)

  get_mother_SV_at_region<-function(npc){ 
  ##loop number of SV per anchor
    print(paste0("number", npc))
  ftfo<-ftfp[npc,]
  mtfp<-mtf[mtf$RefcontigID1 == chr & (mtf$RefStartPos>(ftfo$RefStartPos-buffer)) & (mtf$RefEndPos < (ftfo$RefEndPos+buffer)),]
  
  if(nrow(mtfp)==0){
   ###Check alignments at the region
    ftfxr<-ftx[ftx$RefContigID == chr & (ftx$RefStartPos <=ftfo$RefStartPos) & (ftx$RefEndPos >=ftfo$RefEndPos),]
    #ftxn<-nrow(ftfxr)###But large indels might involve more than 1 match group in the region. What to do with those? Check xmapID1 and xmapID2
    ftxn<-length(ftfxr[!duplicated(ftfxr$QryContigID),]$QryContigID)
    ftfpt<-ftfp[ftfp$RefcontigID1 == chr & ftfp$RefStartPos <=ftfo$RefStartPos & ftfp$RefEndPos >= ftfo$RefEndPos,]
  
    mtfxr<-mtx[mtx$RefContigID == chr & (mtx$RefStartPos <=ftfo$RefStartPos) & (mtx$RefEndPos >=ftfo$RefEndPos),]
    if(ftfo$Zygosity == "heterozygous" & ftxn==2 & nrow(mtfxr)<=2 & nrow(mtfxr) !=0 & nrow(ftfpt) ==1){
    ###Case where mother is a,a and father is a,b (with no call as it's same as anchor)
      print("Found SV at alternate sample")
      alteralign<-ftfxr[ftfxr$QryContigID != ftfo$QryContigID,]
      queriesstart<-getqrycmapcoor(as.numeric(as.character(ftfo$RefStartIdx)), ftfr, ftx, as.numeric(as.character(ftfo$RefcontigID1)))
      queriesend<-getqrycmapcoor(as.numeric(as.character(ftfo$RefEndIdx)), ftfr, ftx, as.numeric(as.character(ftfo$RefcontigID1)))
      queriesstartwanted<-queriesstart[queriesstart$QryContigID != ftfo$QryContigID[1],]
      queriesendwanted<-queriesend[queriesend$QryContigID != ftfo$QryContigID[1],]
      queriesstartwanted<-queriesstartwanted[queriesstartwanted$QryContigID %in% queriesendwanted$QryContigID,]
      queriesstartwanted<-queriesstartwanted[!duplicated(queriesstartwanted$QryContigID),]
      queriesendwanted<-queriesendwanted[queriesendwanted$QryContigID %in% queriesstartwanted$QryContigID,]
      queriesendwanted<-queriesendwanted[!duplicated(queriesendwanted$QryContigID),]
      querytotallabel<-ftfq[ftfq$`#h CMapId` == queriesstartwanted$QryContigID,]$NumSites[1]
      #QryStartPos<-ftfq[ftfq$`#h CMapId` == queriesstartwanted$QryContigID & ftfq$SiteID == queriesstartwanted$Qrylabel,]$Position
      #QryEndPos<-ftfq[ftfq$`#h CMapId` == queriesendwanted$QryContigID & ftfq$SiteID == queriesendwanted$Qrylabel,]$Position
      # result<-data.frame(ftfo$`#h SmapEntryID`, ftfxr$RefContigID[1], queriesstartwanted$refpos, queriesstartwanted$QryContigID, queriesstartwanted$Qrylabel, QryStartPos, queriesendwanted$QryContigID, queriesendwanted$Qrylabel, QryEndPos, "cut", querytotallabel, ftfo$Type, "Case1") 
      # colnames(result)<-c("SmapID", "RefContigID", "RefStartPos", "QryContigID", "QryStartId", "QryStartPos", "QryContigID", "QryEndID", "QryEndPos", "cutType", "querytotallabe", "Type", "Case")
      if(nrow(queriesstartwanted)==0 | nrow(queriesendwanted)==0){
        result<-data.frame(ftfo$`#h SmapEntryID`, ftfo$RefcontigID1, format(ftfo$RefStartPos,nsmall=1), ftfo$QryContigID, ftfo$QryStartIdx, ftfo$QryContigID, ftfo$QryEndIdx, "undefined", "-", ftfo$Type, "Case7")
        colnames(result)<-c("SmapID", "RefContigID", "RefStartPos", "QryContigID", "QryStartId", "QryContigID", "QryEndID", "cutType", "querytotallabe", "Type", "Case")
        return(result)
      }
      result<-data.frame(ftfo$`#h SmapEntryID`, ftfxr$RefContigID[1], queriesstartwanted$refpos, queriesstartwanted$QryContigID, queriesstartwanted$Qrylabel, queriesendwanted$QryContigID, queriesendwanted$Qrylabel, "cut", querytotallabel, ftfo$Type, "Case1") 
      colnames(result)<-c("SmapID", "RefContigID", "RefStartPos", "QryContigID", "QryStartId", "QryContigID", "QryEndID", "cutType", "querytotallabe", "Type", "Case")
      print("to Case1")
      return(result)
      } else if (ftfo$Zygosity == "heterozygous" & ftxn==2 & nrow(mtfxr)<=2 & nrow(mtfxr) !=0 & nrow(ftfpt) !=1){
    ###Case where mother is b, b and father is a,c (with no call in mother as it's same as anchor)
      print("to Case4")
      # result<-data.frame(ftfo$`#h SmapEntryID`, ftfo$RefcontigID1, format(ftfo$RefStartPos,nsmall=1), ftfo$QryContigID, ftfo$QryStartIdx, ftfo$QryStartPos, ftfo$QryContigID, ftfo$QryEndIdx, ftfo$QryEndPos, "alt_homo_het_diff", "-", as.character(ftfo$Type), "Case4")
      # colnames(result)<-c("SmapID", "RefContigID", "RefStartPos", "QryContigID", "QryStartId", "QryStartPos", "QryContigID", "QryEndID", "QryEndPos", "cutType", "querytotallabe", "Type", "Case")
        result<-data.frame(ftfo$`#h SmapEntryID`, ftfo$RefcontigID1, format(ftfo$RefStartPos,nsmall=1), ftfo$QryContigID, ftfo$QryStartPos, ftfo$QryContigID, ftfo$QryEndIdx, "alt_homo_het_diff", "-", as.character(ftfo$Type), "Case4")
        colnames(result)<-c("SmapID", "RefContigID", "RefStartPos", "QryContigID", "QryStartId", "QryContigID", "QryEndID", "cutType", "querytotallabe", "Type", "Case")
        return(result)
      } else {
        # result<-data.frame(ftfo$`#h SmapEntryID`, ftfo$RefcontigID1, format(ftfo$RefStartPos,nsmall=1), ftfo$QryContigID, ftfo$QryStartIdx, ftfo$QryStartPos, ftfo$QryContigID, ftfo$QryEndIdx, ftfo$QryEndPos, "undefined", "-", ftfo$Type, "Case8")
        # colnames(result)<-c("SmapID", "RefContigID", "RefStartPos", "QryContigID", "QryStartId", "QryStartPos", "QryContigID", "QryEndID", "QryEndPos", "cutType", "querytotallabe", "Type", "Case")
        result<-data.frame(ftfo$`#h SmapEntryID`, ftfo$RefcontigID1, format(ftfo$RefStartPos,nsmall=1), ftfo$QryContigID, ftfo$QryStartIdx, ftfo$QryContigID, ftfo$QryEndIdx, "undefined", "-", ftfo$Type, "Case8")
        colnames(result)<-c("SmapID", "RefContigID", "RefStartPos", "QryContigID", "QryStartId", "QryContigID", "QryEndID", "cutType", "querytotallabe", "Type", "Case")
        return(result)
      }
    
    } else if(nrow(mtfp)!=0){
    ###Case where mother is a,b and father is a,b
    ###Test if a-b in father is same as the a-b in mum
       print("No SV found at alternate sample")
       mtfps<-mtfp[mtfp$Type %in% types & (mtfp$RefStartPos <= ftfo$RefStartPos) & (mtfp$RefEndPos >= ftfo$RefEndPos),]
      #mtfps<-mtfp[mtfp$Type %in% types & (mtfp$RefStartPos <= ftfo$RefStartPos) & (mtfp$RefEndPos >= ftfo$RefEndPos),]
      if(as.numeric(as.character(length(mtfps$Zygosity[mtfps$Zygosity=="heterozygous"])))>0){
      ###size comparison
        ###Need fix for cases where father is a,b, mother is b,c, should be rare after trio binning tho as proband should only has 2 alleles and in cross checking.
        print("to Case2")
        sizecomp<-min(abs(mtfps$SVsize-ftfo$SVsize))
        if(sizecomp>500){   ##is this a problem if size comp is longer than 1 element?
            # result<-data.frame(ftfo$`#h SmapEntryID`, ftfo$RefcontigID1, format(ftfo$RefStartPos,nsmall=1), ftfo$QryContigID, ftfo$QryStartIdx, ftfo$QryStartPos, ftfo$QryContigID, ftfo$QryEndIdx, ftfo$QryEndPos, "both_het_diff", "-", as.character(ftfo$Type), "Case5")
            # colnames(result)<-c("SmapID", "RefContigID", "RefStartPos", "QryContigID", "QryStartId", "QryStartPos", "QryContigID", "QryEndID", "QryEndPos", "cutType", "querytotallabe", "Type", "Case")
            result<-data.frame(ftfo$`#h SmapEntryID`, ftfo$RefcontigID1, format(ftfo$RefStartPos,nsmall=1), ftfo$QryContigID, ftfo$QryStartIdx, ftfo$QryContigID, ftfo$QryEndIdx, "both_het_diff", "-", as.character(ftfo$Type), "Case5")
            colnames(result)<-c("SmapID", "RefContigID", "RefStartPos", "QryContigID", "QryStartId", "QryContigID", "QryEndID",  "cutType", "querytotallabe", "Type", "Case")
            return(result)
          } else {
            # result<-data.frame(ftfo$`#h SmapEntryID`, ftfo$RefcontigID1, format(ftfo$RefStartPos, nsmall=1), ftfo$QryContigID, ftfo$QryStartIdx, ftfo$QryStartPos, ftfo$QryContigID, ftfo$QryEndIdx, ftfo$QryEndPos, "both_het_same", "-", ftfo$Type, "Case2")
            # colnames(result)<-c("SmapID", "RefContigID", "RefStartPos", "QryContigID", "QryStartId", "QryStartPos", "QryContigID", "QryEndID", "QryEndPos", "cutType", "querytotallabe", "Type", "Case")
            result<-data.frame(ftfo$`#h SmapEntryID`, ftfo$RefcontigID1, format(ftfo$RefStartPos, nsmall=1), ftfo$QryContigID, ftfo$QryStartIdx, ftfo$QryContigID, ftfo$QryEndIdx, "both_het_same", "-", ftfo$Type, "Case2")
            colnames(result)<-c("SmapID", "RefContigID", "RefStartPos", "QryContigID", "QryStartId", "QryContigID", "QryEndID", "cutType", "querytotallabe", "Type", "Case")
            return(result)
          }
        return(result)
      } else if (as.numeric(as.character(length(mtfps$Zygosity[mtfps$Zygosity=="homozygous"])))>0){
        ###Case where mother is a,a and father is a,b (with homozgyous call as both a aligns to the b allele) 
        ###compare SV size between father and mother
        print("to Case3")
        sizecomp<-min(abs(mtfps$SVsize-ftfo$SVsize))
        if(sizecomp>500){
          # result<-data.frame(ftfo$`#h SmapEntryID`, ftfo$RefcontigID1, format(ftfo$RefStartPos,nsmall=1), ftfo$QryContigID, ftfo$QryStartIdx, ftfo$QryStartPos, ftfo$QryContigID, ftfo$QryEndIdx, ftfo$QryEndPos, "alt_homo_diff", "-", as.character(ftfo$Type), "Case6")
          # colnames(result)<-c("SmapID", "RefContigID", "RefStartPos", "QryContigID", "QryStartId", "QryStartPos", "QryContigID", "QryEndID", "QryEndPos", "cutType", "querytotallabe", "Type", "Case")
          result<-data.frame(ftfo$`#h SmapEntryID`, ftfo$RefcontigID1, format(ftfo$RefStartPos,nsmall=1), ftfo$QryContigID, ftfo$QryStartIdx, ftfo$QryContigID, ftfo$QryEndIdx, "alt_homo_diff", "-", as.character(ftfo$Type), "Case6")
          colnames(result)<-c("SmapID", "RefContigID", "RefStartPos", "QryContigID", "QryStartId", "QryContigID", "QryEndID", "cutType", "querytotallabe", "Type", "Case")
          return(result)
        } else if(sizecomp<500){
            # result<-data.frame(ftfo$`#h SmapEntryID`, ftfo$RefcontigID1, format(ftfo$RefStartPos, nsmall=1), ftfo$QryContigID, ftfo$QryStartIdx, ftfo$QryStartPos, ftfo$QryContigID, ftfo$QryEndIdx, ftfo$QryEndPos, "alt_homo_same_tobecut", "-", ftfo$Type, "Case3")
            # colnames(result)<-c("SmapID", "RefContigID", "RefStartPos", "QryContigID", "QryStartId", "QryStartPos", "QryContigID", "QryEndID", "QryEndPos", "cutType", "querytotallabe", "Type", "Case")
            # result$RefStartPos<-format(result$RefStartPos, nsmall=1)
          result<-data.frame(ftfo$`#h SmapEntryID`, ftfo$RefcontigID1, format(ftfo$RefStartPos, nsmall=1), ftfo$QryContigID, ftfo$QryStartIdx, ftfo$QryContigID, ftfo$QryEndIdx, "alt_homo_same_tobecut", "-", ftfo$Type, "Case3")
          colnames(result)<-c("SmapID", "RefContigID", "RefStartPos", "QryContigID", "QryStartId", "QryContigID", "QryEndID", "cutType", "querytotallabe", "Type", "Case")
          result$RefStartPos<-format(result$RefStartPos, nsmall=1)
          return(result)
        }
      } else {
          # result<-data.frame(ftfo$`#h SmapEntryID`, ftfo$RefcontigID1, format(ftfo$RefStartPos,nsmall=1), ftfo$QryContigID, ftfo$QryStartIdx, ftfo$QryStartPos, ftfo$QryContigID, ftfo$QryEndIdx, ftfo$QryEndPos, "undefined", "-", ftfo$Type, "Case7")
          # colnames(result)<-c("SmapID", "RefContigID", "RefStartPos", "QryContigID", "QryStartId", "QryStartPos", "QryContigID", "QryEndID", "QryEndPos", "cutType", "querytotallabe", "Type", "Case")
          result<-data.frame(ftfo$`#h SmapEntryID`, ftfo$RefcontigID1, format(ftfo$RefStartPos,nsmall=1), ftfo$QryContigID, ftfo$QryStartIdx, ftfo$QryContigID, ftfo$QryEndIdx, "undefined", "-", ftfo$Type, "Case7")
          colnames(result)<-c("SmapID", "RefContigID", "RefStartPos", "QryContigID", "QryStartId", "QryContigID", "QryEndID", "cutType", "querytotallabe", "Type", "Case")
          return(result)
      }
      return(result)
      print(result)
     }
  }
  perchrresult<-lapply(seq(1, npchr, 1), get_mother_SV_at_region)
  perchrresult_complete<-do.call(rbind,perchrresult)
}

allresults<-lapply(chrtoberun, get_SV_alt_sample_at_region) 
allresults_complete<-do.call(rbind, allresults)

print(paste0(nrow(ftfs), " heterozygous calls to be checked"))
print(paste0(nrow(allresults_complete), " regions checked"))
print(table(allresults_complete$cutType))
print(table(allresults_complete$Case))

write.table(allresults_complete, file=outputfile, row.names=F, col.names = T, sep="\t", quote=F)
print(paste0("written to ", outputfile))
check<-data.frame(table(allresults_complete$SmapID))
print(check[check$Freq>1,])   ###These are usually involving end of one of the query contigs

}

cross_checking(ftfpath, ftfxmap, ftfrmap, ftfqmap, mtfpath, mtfxmap, outputfile, sizethres, buffer, types)
 ##loop number of SV per anchor
 ###Loop number of ref with SV

###Cut at the query labelIDs but on the original query map instead of the q.cmap (Check # of labels to make sure they are the same between the input query and the q.cmap)



