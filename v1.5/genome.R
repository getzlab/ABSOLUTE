## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

GetChrLens <- function(x=FALSE) {
  lens <- c(247249719, 242951149, 199501827, 191273063,
            180857866, 170899992, 158821424, 
            146274826, 140273252, 135374737, 134452384,
            132349534, 114142980, 106368585, 
            100338915, 88827254, 78774742, 76117153, 63811651,
            62435964, 46944323, 49691432, 154913754)
    
  if (x == FALSE) {
    lens <- lens[c(1:22)]
  }
    
  return(lens)
}

GetCentromerePos <- function(x=FALSE) {
  #data(ChrArmsDat, package = "ABSOLUTE")
  #load("/xchip/tcga/Tools/absolute/releases/v1.5/data/ChrArmsDat.RData")
  #this is loaded through RunAbsolute already
  chrarm_names <- paste(c(1:22), "q", sep = "")
  
  if (x) {
    chrarm_names <- c(chrarm_names, "Xq")
  }
  
  cent_pos <- chr.arms.dat[chrarm_names, "Start.bp"]
  
  return(cent_pos)
}
