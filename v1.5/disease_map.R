## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

DetermineGroup <- function(primary.disease) {

  return(NA) #fix for now

  data(diseaseMap)

  group <- try(get(primary.disease, disease_map))
  if (inherits(group, "try-error")) {
    ## It doesn't exist, just return the primary disease, it'll
    ## fall through
    return(primary.disease)
  } else {
    return(group)
  }
}
