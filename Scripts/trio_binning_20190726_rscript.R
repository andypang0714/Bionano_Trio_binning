#!/usr/bin/env Rscript
#install.packages("optparse")
library("optparse")

source("./Scripts/readV3bnx_2.R")
source("./Scripts/readmaps.R")

#a is bnx mapped to parent1, b is bnx mapped to parent2, x1 is alignment of all mol to parent1, x2 is alignment of all mol to parent2, all is autonoise1_rescaled.bnx in alignment to parent1, and all2 is autonoise1_rescaled.bnx in alignment to paretn2. 

parser <- OptionParser()
parser <- add_option(parser, c("-a", "--bnx1"), help="bnx mapped to parent1")
parser <- add_option(parser, c("-b", "--bnx2"), help="bnx mapped to parent2")
parser <- add_option(parser, c("-x", "--alignment1"), help="alignment of all mol to parent1")
parser <- add_option(parser, c("-y", "--alignment2"), help="alignment of all mol to parent2")
parser <- add_option(parser, c("-m", "--autonoise_bnx1"), help="autonoise1_rescaled.bnx in alignment to parent1")
parser <- add_option(parser, c("-n", "--autonoise_bnx2"), help="autonoise1_rescaled.bnx in alignment to parent2")
args<-parse_args(parser)

#args<-parser$parse_args()

a<-args$bnx1
b<-args$bnx2
x1<-args$alignment1
x2<-args$alignment2
all<-args$autonoise_bnx1
all2<-args$autonoise_bnx2


trio_binning<-function(a,b,x1,x2, all, all2){
  print("in loop")
  print(paste("bnx aligned to parent1:", a))
  print(paste("bnx aligned to parent2:", b))
  print(paste("alignment to parent1", x1))
  print(paste("alignment to parent2", x2))
  print(paste("proband autonoise_rescaled.bnx to parent1:", all))
  print(paste("proband autonoise_rescaled.bnx to parent2:", all2))
  prefix1<-strsplit(a, split="[.]")[[1]][1]
  prefix1<-strsplit(prefix1, split="_prebinning")[[1]][1]
  prefix2<-strsplit(b, split="[.]")[[1]][1]
  prefix2<-strsplit(prefix2, split="_prebinning")[[1]][1]

  a<-readbnx2(a)
  b<-readbnx2(b)
  x1<-readxmap(x1)
  x2<-readxmap(x2)
  all<-readbnx2(all)
  all2<-readbnx2(all2)

  print(paste0(nrow(all), " molecules (", sum(all$Length), " bp) input to align to parent1"))
  print(paste0(nrow(all2), " molecules (", sum(all2$Length), " bp) input to align to parent2"))

  a$MapID<-a$MoleculeID
  b$MapID<-b$MoleculeID
  all$MapID<-all$MoleculeID
  all2$MapID<-all2$MoleculeID
  au<-a[!a$MapID %in% b$MapID,]
  aul<-sum(au$Length)
  #print(paste0(nrow(a), " molecules aligned to parent1"))
  #print(paste0(nrow(au), " molecules (", aul, " bp) out of ", nrow(a), " molecules (", sum(a$Length), " bp) aligned to parent1 did not align to parent2"))
  bu<-b[!b$MapID %in% a$MapID,]
  bul<-sum(bu$Length)
  #print(paste0(nrow(b), " molecules aligned to parent2"))
  #print(paste0(nrow(bu), " molecules (", bul, " bp) out of ", nrow(b), " molecules (", sum(b$Length), " bp) aligned to parent2 did not align to parent1"))

  ###For molecules aligning to both sequences, see if the molecule alignment confidence is very different between the parents
  ai<-a[a$MapID %in% b$MapID,]
  ain<-a[!a$MapID %in% b$MapID,]
  print(paste0(nrow(ai), "molecules (", sum(ai$Length), " bp) aligned to parent1 also align to parent2"))
  print(paste0(nrow(ain), " molecules (", sum(ain$Length), " bp) aligned to parent1 does not align to parent2"))
  x1f<-x1[x1$QryContigID %in% ai$MapID,]
  bi<-b[b$MapID %in% a$MapID,]
  bin<-b[!b$MapID %in% a$MapID,]
  print(paste0(nrow(bi), " molecules (", sum(bi$Length), " bp) aligned to parent2 also align to parent1"))
  print(paste0(nrow(bin), " molecules (", sum(bin$Length), " bp) aligned to parent2 does not align to parent1"))
  x2f<-x2[x2$QryContigID %in% bi$MapID,]
  xaf<-merge(x1f, x2f, by="QryContigID")

  xaf$ratio1<-xaf$Confidence.x/xaf$Confidence.y
  xaf1<-xaf[xaf$ratio1>0.8,]
  xaf$ratio2<-xaf$Confidence.y/xaf$Confidence.x
  xaf2<-xaf[xaf$ratio2>0.8,]
  aia<-ai[ai$MapID %in% xaf1$QryContigID,]
  print(paste0(nrow(ai), " aligned to both, ", nrow(x1f), " in xmap", nrow(xaf), "in combined table"))
  print(paste0(nrow(bi), " aligned to both, ", nrow(x2f), " in xmap", nrow(xaf), "in combined table"))
  bia<-bi[bi$MapID %in% xaf2$QryContigID,]
  #print(paste0(nrow(aia), " molecules (", sum(aia$Length), " bp) aligned to parent1 also align to parent2 (after alignment score ratio filter of 0.8)"))
  #print(paste0(nrow(bia), " molecules (", sum(bia$Length), " bp) aligned to parent2 also align to parent1 (after alignment score ratio filter of 0.8)"))
  aa<-rbind(aia, ain)    ###These are molecules aligned to parent1 that does not align to parent2 or did not pass relative alignment score in parent2
  #print(paste0(nrow(aa), " molecules (", sum(aa$Length), " bp) aligned to parent1 and does not align to parent2 or did not pass relative alignment score in parent2"))
  bb<-rbind(bia, bin)    ###These are molecules aligned to parent2 that does not align to parent1 or did not pass relative alignment score in parent1
  #print(paste0(nrow(bb), " molecules (", sum(bb$Length), " bp) aligned to parent2 and does not align to parent1 or did not pass relative alignment score in parent1"))


  ###second way of filtering by alignment
  xaf$diff1<-xaf$Confidence.x-xaf$Confidence.y
  xaf$diff2<-xaf$Confidence.y-xaf$Confidence.x
  xaf1<-xaf[xaf$diff1>-2,]
  xaf2<-xaf[xaf$diff2>-2,]
  common<-xaf1[xaf1$QryContigID %in% xaf2$QryContigID,]

  cl<-sum(all[all$MapID %in% common$QryContigID,]$Length)
  lonep<-sum(all[all$MapID %in% xaf1$QryContigID,]$Length)
  ltwop<-sum(all2[all2$MapID %in% xaf2$QryContigID,]$Length)
  loneonly<-lonep-cl+aul
  ltwoonly<-ltwop-cl+bul
  print(paste0(nrow(xaf1), " molecules (", sum(all[all$MapID %in% xaf1$QryContigID,]$Length), " bp) aligned to parent1 (after alignment score difference filter of -2) - may align with low conf. or maynot align at all to parent2-prebpp"))
  print(paste0(nrow(xaf1), " molecules (", sum(all2[all2$MapID %in% xaf2$QryContigID,]$Length), " bp) aligned to parent2 (after alignment score difference filter of -2) - may align with low conf. or maynot align at all to parent1-prebpp"))
  print(paste0(nrow(common), " molecules (", sum(all[all$MapID %in% common$QryContigID,]$Length), " bp) align to both parents equally well and needs to be subsampled-prebpp"))
  print(paste0("Leaving ",   loneonly, " bp (", lonep-cl, "bp +",  aul, " bp) uniquely aligned to parent1"))
  print(paste0("Leaving ",   ltwoonly, " bp (", ltwop-cl," bp +", bul, " bp) uniquely aligned to parent2"))

  #print(paste0(nrow(xaf1), " molecules (", sum(xaf1$QryLen.x), " bp) aligned to parent1 (after alignment score difference filter of -2) - may align with low conf. or maynot align at all to parent2"))
  #print(paste0(nrow(xaf2), " molecules (", sum(xaf2$QryLen.y), " bp) aligned to parent2 (after alignment score difference filter of -2) - may align with low conf. or maynot align at all to parent1"))
  #print(paste0(nrow(common), " molecules (", sum(common$QryLen.x), " bp) align to both parents equally well and needs to be subsampled"))
  #print(paste0("Leaving ",   sum(xaf1$QryLen.x)-sum(common$QryLen.x) +aul, " bp (", sum(xaf1$QryLen.x)-sum(common$QryLen.x), "bp +",  aul, " bp) uniquely aligned to parent1"))
  #print(paste0("Leaving ",   sum(xaf2$QryLen.y) - sum(common$QryLen.y)+bul, " bp (", sum(xaf2$QryLen.y) - sum(common$QryLen.y)," bp +", bul, " bp) uniquely aligned to parent2"))

  s1<-sample(1:nrow(common), nrow(common)/2)
  tb1<-common[s1,]
  tb2<-common[!common$QryContigID %in% tb1$QryContigID,]
  xaf1n<-xaf1[!xaf1$QryContigID %in% common$QryContigID,]
  xaf2n<-xaf2[!xaf2$QryContigID %in% common$QryContigID,]

  #print(paste0(nrow(xaf1n), " molecules (", sum(xaf1n$QryLen.x), " bp) aligned to parent1 only"))
  #print(paste0(nrow(xaf2n), " molecules (", sum(xaf2n$QryLen.y), " bp) aligned to parent2 only"))

  xaf1g<-rbind(xaf1n, tb1)
  xaf2g<-rbind(xaf2n, tb2)
  aia<-ai[ai$MapID %in% xaf1g$QryContigID,]
  bia<-bi[bi$MapID %in% xaf2g$QryContigID,]
  #print(paste0(nrow(aia), " molecules (", sum(aia$Length), " bp) aligned to parent1 also align to parent2 (after alignment score difference filter of -2) and after downsampling to keep even coverage"))
  #print(paste0(nrow(bia), " molecules (", sum(bia$Length), " bp) aligned to parent2 also align to parent1 (after alignment score difference filter of -2) and after downsampling to keep even coverage"))
  aa<-rbind(aia, ain)
  #print(paste0(nrow(aa), " molecules (", sum(aa$Length), " bp) aligned to parent1 and does not align to parent2 or did not pass relative alignment score in parent2"))
  #print(paste0(nrow(aa), " molecules (", sum(aa$Length), " bp) aligned to parent1 and does or does not align to parent2 or did not pass relative alignment score in parent2, include those that align also to parent2 and have been downsampled"))

  bb<-rbind(bia, bin)
  #print(paste0(nrow(bb), " molecules (", sum(bb$Length), " bp) aligned to parent2 and does not align to parent1 or did not pass relative alignment score in parent1"))
  #print(paste0(nrow(bb), " molecules (", sum(bb$Length), " bp) aligned to parent2 and does or does not align to parent1 or did not pass relative alignment score in parent1, include those that align also to parent1 and have been downsampled"))


  noalign<-all[!(all$MapID %in% a$MapID | all$MapID %in% b$MapID),]
  noalign2<-all2[!(all2$MapID %in% a$MapID | all2$MapID %in% b$MapID),]
  print(paste0(nrow(noalign), " molecules (", sum(noalign$Length), " bp) do not align to parent1 also did not align to parent2, downsample them to half"))
  print(paste0(nrow(noalign2), " molecules (", sum(noalign2$Length), " bp) do not align to parent2 also did not align to parent1, downsample them to half"))

  noalignfor1<-noalign[!noalign$MapID %in% noalign2$MapID,]
  noalignfor2<-noalign2[!noalign2$MapID %in% noalign$MapID,]
  noalign_common<-noalign[noalign$MapID %in% noalign2$MapID,]
  nonaligned<-c(noalign$MapID , noalign2$MapID)
  nonalignedmol<-nonaligned[!duplicated(nonaligned)]
  print(paste0(nrow(noalign_common), " non-aligned molecules are in both bnx, ", nrow(noalignfor1), " non-aligned molecules only exist in paret1, and ", nrow(noalignfor2), " non-aligned molecules only exist in paretn2"))

  ###Only consider noalign_common
  noalignfor1split<-noalignfor1[sample(1:nrow(noalignfor1), nrow(noalignfor1)/2),]
  noalignfor1split_2<-noalignfor1[!noalignfor1$MoleculeID %in% noalignfor1split$MoleculeID,]
  noalignfor2split<-noalignfor2[sample(1:nrow(noalignfor2), nrow(noalignfor2)/2),]
  noalignfor2split_2<-noalignfor2[!noalignfor2$MoleculeID %in% noalignfor2split$MoleculeID,]
  noaligncommonsplit<-noalign_common[sample(1:nrow(noalign_common), nrow(noalign_common)/2),]
  noaligncommonsplit_2<-noalign_common[!noalign_common$MoleculeID %in% noaligncommonsplit$MoleculeID,]
  noalignsum1<-rbind(noalignfor1split, noalignfor2split, noaligncommonsplit)
  noalignsum2<-rbind(noalignfor1split_2, noalignfor2split_2, noaligncommonsplit_2)
  # noalign<-noalign[sample(1:nrow(noalign), nrow(noalign)/2),]
  # noalign2<-noalign2[sample(1:nrow(noalign2), nrow(noalign2)/2),]

  finala<-rbind(aa, noalignsum1)
  finalb<-rbind(bb, noalignsum2)

  unique_parent1_molID1<-xaf1[!xaf1$QryContigID %in% common$QryContigID,]$QryContigID
  unique_parent1_molID<-c(unique_parent1_molID1, au$MapID)
  unique_parent2_molID1<-xaf2[!xaf2$QryContigID %in% common$QryContigID,]$QryContigID
  unique_parent2_molID<-c(unique_parent2_molID1, bu$MapID)
  #nonaligned<-c(noalign$MapID , noalign2$MapID)

  if(nrow(all)>nrow(all2)){                
    #There are some molecules that only exist in one of the autonoise bnx due to different bpp between the sets. Usually one is inclusive of the other set, so we want to use the bigger (more inclusive) set for exporting.
    #Warning: There can be some molecules that did not pass 150 kbp in one of the alignments and falsely identified as sample 1 unique or vice versa. Nonetheless, these amount of molecules are minimal and does not affect the final assembly of the parental alleles.
    print("use parent1 autonoise bnx for exporting bnx ")
  } else if (nrow(all2)>nrow(all)){
    print ("use parent2 autonoise bnx for exporting bnx")
  }

  all_mol_ID<-c(common$QryContigID, unique_parent1_molID, unique_parent2_molID, nonaligned)
  all_mol_nonunique<-all_mol_ID[!duplicated(all_mol_ID)]
  check<-data.frame(table(all_mol_ID))
  dup<-check[check$Freq >1,]

  print("Summary:")
  print(paste0(length(all_mol_nonunique), " molecules_involved"))
  print(paste0(nrow(common), " molecules (", cl, " bp) align to both parents equally well and needs to be subsampled"))
  print(paste0("Leaving ",length(unique_parent1_molID), " molecules (", loneonly, " bp (", lonep-cl, "bp +",  aul, " bp)) uniquely aligned to parent1"))
  print(paste0("Leaving ", length(unique_parent2_molID), " molecules (",  ltwoonly, " bp (", ltwop-cl," bp +", bul, " bp)) uniquely aligned to parent2"))
  #print(paste0(nrow(noalign2)+nrow(noalign), " molecules (", sum(noalign2$Length)+sum(noalign$Length), " bp) do not align to either parent1 or parent2"))
  print(paste0(length(nonalignedmol), " molecules (", sum(noalignfor1$Length)+sum(noalignfor2$Length)+sum(noalign_common$Length), " bp) , ", nrow(noalignfor1)+nrow(noalignfor2)+nrow(noalign_common), " molecules do not align to either parent1 or parent2"))
  print(paste0(nrow(finala), " molecules (", sum(finala$Length), " bp) in bnx for parent1 assembly"))
  print(paste0(nrow(finalb), " molecules (", sum(finalb$Length), " bp) in bnx for parent2 assembly"))

  write.table(common$QryContigID, file=paste0(prefix1,"_shared_mol_ID.txt"), col.names=F, row.names=F, quote=F)
  write.table(unique_parent1_molID, file=paste0(prefix1,"_aligned_only_to_parent1_molID.txt"), col.names=F, row.names=F, quote=F, sep="\n")
  write.table(unique_parent2_molID, file=paste0(prefix2,"_aligned_only_to_parent2_molID.txt"), col.names=F, row.names=F, quote=F, sep="\n")
  write.table(nonalignedmol, file=paste0(prefix1,"_aligned_to_neither_molID.txt"), col.names=F, row.names=F, quote=F, sep="\n")
  write.table(all_mol_ID, file=paste0(prefix1,"_all_involved_molID.txt"), col.names=F, row.names=F, quote=F, sep="\n")

  write.table(finala$MapID, file=paste0(prefix1, "_postbinning_molID.txt"), col.names=F, row.names=F, quote=F)
  write.table(finalb$MapID, file=paste0(prefix2, "_postbinning_molID.txt"), col.names=F, row.names=F, quote=F)

}

trio_binning(a,b,x1,x2, all, all2)
