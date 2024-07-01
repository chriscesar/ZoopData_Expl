# scores.gllvm.R
## function provided by Bert van der Veen at:
# https://github.com/BertvanderVeen/GLLVM-workshop/blob/90a89dcf7d6cc0db2bae9a493dd1d535b28f8b57/3Wednesday/Tools.Rmd#L230
# to extract scores from a gllvm object

scores.gllvm <- function(x, display ="sites", choices = NULL){
  sol <- list()
  if(is.null(choices) || any(choices>(x$num.lv+x$num.lv.c+x$num.RR))){
    choices <- 1:(x$num.lv+x$num.lv.c+x$num.RR)
  }
  if(display%in%c("sites","both")){
    sol$sites <- gllvm::getLV(x)[,choices,drop=FALSE]
    if(display=="sites"){
      sol <- sol$sites
    }
  }
  if(display%in%c("species","both")){
    sol$species <- gllvm::getLoadings(x)[,choices,drop=FALSE]
    if(display=="species"){
      sol <- sol$species
    }
  }
  return(sol)
}