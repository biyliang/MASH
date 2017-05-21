##################################################################
##################################################################
##
##  M-BITES (Mosquito Bout-based and Individual-based Transmission Ecology Simulation)
##  Version 0.9
##  April 6, 2017
##
##  This version was designed and written by David L. Smith (aka.
##  Dave), Sean Wu, and Hector Sanchez.
##  Please send bug reports, comments, and suggestions to
##  <smitdave@gmail.com> or <slwu89@berkeley.edu>.
##
##  Robert C. Reiner, Jr. (aka Bobby) <bcreiner@uw.edu>, Hector
##  Sanchez Castellanos <sanchez.hmsc@gmail.com> Sean Wu
##  <slwu89@berkeley.edu>, and Amit Verma <amit.verma13@gmail.com>
##  helped with development, debugging and documentation of
##  version 1.0.
##
##  M-BITES (formerly DHM) was conceived of by David Smith, and it was inspired
##  by discussions with many people, including Bobby, Hector,
##  Sean, Amit, Arnaud Le Menach, Nick Ruktanonchai, Samson
##  Kiware, Gerry Killeen, Tom Scott, Ellis McKenzie, Steven W.
##  Lindsay, Willem Takken, Philip Eckhoff, Nicole Achee, Chris
##  Barker, Nakul Chitnis, Justin Cohen, Su Yun Kang, Audrey
##  Lenhart, John Marshall, Phil McCall, Catherine Moyes, Doug
##  Norris, Alex Perkins, Chris Stone, Edward Wenger, and Anne
##  Wilson.
##
##################################################################
##################################################################

#################################################################
#
#   MASH/MBITES
#   M-BITES Data Logging Utilities
#   R version
#   Sean Wu
#   April 4, 2017
#
#################################################################


#################################################################
# Coarse-grained tracking of adult counts via .csv output
#################################################################

# trackAdults_init: write output inside of "OUTPUT" folder in given directory.
# dir: directory to write in (usually set to MASH)
# fileName: name of file to write to
trackAdults_init <- function(directory, fileName){
  if(!dir.exists(paste0(directory,"OUTPUT"))){
    dir.create(paste0(directory))
    dir.create(paste0(directory,"OUTPUT"))
  }
  if(file.exists(paste0(directory,"OUTPUT/",fileName))){
    stop("trackAdults_init cannot init a file that already exists!")
  }
  conAdults = file(description = paste0(directory,"OUTPUT/",fileName),open = "wt")
  writeLines(text = paste0(c("time","Ff","Bf","Rf","Lf","Of","Mf","Sf","Ef","Mm","Sm","Rm"),collapse = ","),con = conAdults,sep = "\n") # write header
  return(conAdults)
}

# trackAdults: write aquatic population state variables to .csv
# con: text connection to write to
trackAdults <- function(con){

  ixM = which(sapply(MPopF$mosy,function(x){isAlive(x) & !is.null(x$id)}))
  ixF = which(sapply(MPopM$mosy,function(x){isAlive(x) & !is.null(x$id)}))

  Fstate = c("F"=0,"B"=0,"R"=0,"L"=0,"O"=0,"M"=0,"S"=0,"E"=0)
  Mstate = c("M"=0,"S"=0,"R"=0)

  Fcount = table(unlist(lapply(MPopF$mosy[ixF],function(x){x$state})))
  Fstate[names(Fcount)] = Fcount

  Mcount = table(unlist(lapply(MPopM$mosy[ixM],function(x){x$state})))
  Mstate[names(Mcount)] = Mcount

  txtOut = paste0(c(tMax,Fstate,Mstate),collapse = ",")
  writeLines(text = txtOut,con = con,sep = "\n")
}

# importAdults: import adult population counts
importAdults <- function(directory, fileName){
  read.csv(file = paste0(directory,"OUTPUT/",fileName),header = TRUE)
}


#################################################################
# Tracking of mosquito histories and bionomics via .json output
#################################################################

# clearOutput: erase all files in given OUTPUT directory (use with caution)
clearOutput <- function(directory){
  if(!dir.exists(paste0(directory,"OUTPUT"))){
    stop("directory does not exist; nothing to clear")
  }
  dirFiles = system(command = paste0("ls ",directory,"OUTPUT/"),intern = TRUE)
  xx = menu(c("Yes", "No"), title=paste0("There are ",length(dirFiles)," files in OUTPUT/.. are you sure you want to delete all?"))
  if(xx==1){
    file.remove(paste0(directory,"OUTPUT/",dirFiles))
  }
}

# trackBionomics: write bionomics to .json
# dir: directory to write in (usually set to MASH)
# fileName: name of file to write to
trackBionomics <- function(directory, fileName){
  con = file(description = paste0(directory,"OUTPUT/",fileName),open = "wt")
  deadM = sapply(MPopF$mosy,Negate(isAlive))
  deadID = sapply(MPopF$mosy[deadM],function(x){x$id})
  # bionomics
  listOut = lapply(MPopF$mosy[deadM],function(x){x$bionomics})
  names(listOut) = paste0("mosy",deadID)
  writeLines(text = jsonlite::toJSON(x = listOut,pretty = TRUE),con = con)
  close(con)
}

# trackHistory: write female histories to .json
# dir: directory to write in (usually set to MASH)
# fileName: name of file to write to
trackHistory <- function(directory, fileName){
  con = file(description = paste0(directory,"OUTPUT/",fileName),open = "wt")
  deadM = sapply(MPopF$mosy,Negate(isAlive))
  deadID = sapply(MPopF$mosy[deadM],function(x){x$id})
  # histories
  listOut = lapply(MPopF$mosy[deadM],function(x){x$history})
  names(listOut) = paste0("mosy",deadID)
  writeLines(text = jsonlite::toJSON(x = listOut,pretty = TRUE),con = con)
  close(con)
}

# trackHistory: write male histories to .json
# dir: directory to write in (usually set to MASH)
# fileName: name of file to write to
trackHistoryM <- function(directory, fileName){
  con = file(description = paste0(directory,"OUTPUT/",fileName),open = "wt")
  deadM = sapply(MPopM$mosy,Negate(isAlive))
  deadID = sapply(MPopM$mosy[deadM],function(x){x$id})
  # histories
  listOut = lapply(MPopM$mosy[deadM],function(x){x$history})
  names(listOut) = paste0("mosy",deadID)
  writeLines(text = jsonlite::toJSON(x = listOut,pretty = TRUE),con = con)
  close(con)
}


#################################################################
# Re-importation of .json files into R
#################################################################

# importBionomics: import female bionomics
importBionomics <- function(directory){
  dirFiles = system(command = paste0("ls ",directory,"OUTPUT/"),intern = TRUE)
  bionomics = grep("bionomics[[:digit:]]+.json",dirFiles)
  bionomicsOut = parallel::mclapply(X = dirFiles[bionomics],FUN = function(x){
    jsonlite::fromJSON(txt = paste0(directory,"OUTPUT/",x))
  },mc.cores = parallel::detectCores()-2)
  bionomicsOut = Reduce(f = c,x = bionomicsOut)
  # extract mosquito IDs and append to the history
  bionomicsIx = unname(sapply(names(bionomicsOut),function(x){sub(pattern = "mosy",replacement =  "",x = x)}))
  for(ix in 1:length(bionomicsOut)){
    bionomicsOut[[ix]]$id = bionomicsIx[ix]
  }
  return(bionomicsOut)
}

# importHistory: import female histories
importHistory <- function(directory){
  dirFiles = system(command = paste0("ls ",directory,"OUTPUT/"),intern = TRUE)
  historyF = grep("historyF[[:digit:]]+.json",dirFiles)
  if(length(historyF)==0){
    stop(paste0("no female histories found in",directory,"OUTPUT/.."))
  }
  historyOut = parallel::mclapply(X = dirFiles[historyF],FUN = function(x){
    jsonlite::fromJSON(txt = paste0(directory,"OUTPUT/",x))
  },mc.cores = parallel::detectCores()-2)
  historyOut = Reduce(f = c,x = historyOut)
  # extract mosquito IDs and append to the history
  historyIx = unname(sapply(names(historyOut),function(x){sub(pattern = "mosy",replacement =  "",x = x)}))
  for(ix in 1:length(historyOut)){
    historyOut[[ix]]$id = historyIx[ix]
  }
  return(historyOut)
}

# importHistoryM: import male histories
importHistoryM <- function(directory){
  dirFiles = system(command = paste0("ls ",directory,"OUTPUT/"),intern = TRUE)
  historyM = grep("historyM[[:digit:]]+.json",dirFiles)
  if(length(historyM)==0){
    stop(paste0("no male histories found in",directory,"OUTPUT/.."))
  }
  historyOut = parallel::mclapply(X = dirFiles[historyM],FUN = function(x){
    jsonlite::fromJSON(txt = paste0(directory,"OUTPUT/",x))
  },mc.cores = parallel::detectCores()-2)
  historyOut = Reduce(f = c,x = historyOut)
  # extract mosquito IDs and append to the history
  historyIx = unname(sapply(names(historyOut),function(x){sub(pattern = "mosy",replacement =  "",x = x)}))
  for(ix in 1:length(historyOut)){
    historyOut[[ix]]$id = historyIx[ix]
  }
  return(historyOut)
}
