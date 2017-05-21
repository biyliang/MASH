#################################################################
#
#   MASH/MBITES
#   Auxiliary Functions
#   R version
#   Sean Wu
#   January 24, 2017
#
#################################################################

replicateArgs <- function(n, f, ...){
  args = as.list(substitute(list(...)))[-1L]
  replicate(n = n,expr = do.call(f,args),simplify = FALSE)
}
