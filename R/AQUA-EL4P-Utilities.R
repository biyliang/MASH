#################################################################
#
#   MASH/MBITES
#   Aquatic Ecology
#   EL4P Data Logging Utilities
#   R Version
#   Sean Wu
#   March 31, 2017
#
#################################################################


#################################################################
#
# Numerical tracking of aquatic populations:
# observe values of state varibles across each site in landscape
#
#################################################################

# trackEL4P_init: write output inside of "OUTPUT" folder in given directory.
# dir: directory to write in (usually set to MASH)
# fileName: name of file to write to
trackEL4P_init <- function(directory, fileName){
  if(!dir.exists(paste0(directory,"OUTPUT"))){
    dir.create(paste0(directory))
    dir.create(paste0(directory,"OUTPUT"))
  }
  if(file.exists(paste0(directory,"OUTPUT/",fileName))){
    stop("trackEL4P_init cannot init a file that already exists!")
  }
  conEL4P = file(description = paste0(directory,"OUTPUT/",fileName),open = "wt")
  writeLines(text =paste0(c("labels",paste0("ix",1:LANDSCAPE$nA)),collapse = ","),con = conEL4P,sep = "\n") # write header
  return(conEL4P)
}

# trackEL4P: write aquatic population state variables to .csv
# con: text connection to write to
trackEL4P <- function(con){
  el4p = sapply(LANDSCAPE$aquaSites,function(x){x$EL4P},simplify = "matrix")
  el4pNames = rownames(el4p)
  time = matrix(data = rep(tMax,length = LANDSCAPE$nA),ncol = LANDSCAPE$nA,nrow = 1)
  writeLines(text = paste0(c("time",time),collapse = ","),con = con,sep = "\n")
  for(i in 1:nrow(el4p)){
    writeLines(text = paste0(c(el4pNames[i],el4p[i,]),collapse = ","),con = con,sep = "\n")
  }
}

# importEL4P: import aquatic population counts
importEL4P <- function(directory, fileName){
  read.csv(file = paste0(directory,"OUTPUT/",fileName),header = TRUE)
}
