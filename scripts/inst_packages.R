
inst_packages <- function(pkgs){
  pkg <- sub(".*[/]","",unique(pkgs) )
  suppressMessages(source("http://bioconductor.org/biocLite.R",prompt.echo = F,echo=F))
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  inst.pkg <- pkg[(pkg %in% installed.packages()[, "Package"])]
  git.pkg <- pkgs[grep("[/]",pkgs)]
  new.pkg <- new.pkg[!(new.pkg %in% sub(".*[/]","",git.pkg))]

  if (length(git.pkg)) {
    if (!require(devtools)) { install.packages("devtools") }
    library(devtools)
    for (i in git.pkg){
      suppressWarnings(suppressMessages(install_github(i)))
  } }

  if (length(new.pkg)) {
    suppressWarnings(suppressMessages(install.packages(new.pkg, repos='http://cran.us.r-project.org', dependencies = TRUE, quiet = T)))
    if (!new.pkg %in% installed.packages()[, "Package"]) {
      Bio.pkg <- new.pkg[!(new.pkg %in% installed.packages()[, "Package"])]
      suppressWarnings(suppressMessages(biocLite(Bio.pkg,suppressUpdates = T,suppressAutoUpdate = T)))
  } }

  data.frame("Pkg_Status"=suppressMessages(sapply(pkg, require, character.only = TRUE)),
             "Pkg_Version"=sapply(pkg,function(x) paste0(x,"_",sessionInfo()$otherPkgs[[x]]$Version)),
             "New_installation_needed?"=ifelse((pkg %in% installed.packages()[, "Package"]),"No","Yes"),row.names = pkg)
  }


print_session_info <- function(){
  cat("\n\n\n\n")
  cat("\n##############################")
  cat("\n### SCRIPT RAN SUCESSFULLY ###")
  cat("\n##############################")
  cat("\n\n\n\n")
  
  cat("\n##########################")
  cat("\n### SYSTEM INFORMATION ###")
  cat("\n##########################")
  Sys.info()
  cat("\n\n\n\n")
  
  cat("\n###########################")
  cat("\n### SESSION INFORMATION ###")
  cat("\n###########################")
  sessionInfo()
}