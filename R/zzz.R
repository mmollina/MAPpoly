.onAttach <- function(libname, pkgname){
  packageStartupMessage(StartWelcomeMessage())
}

StartWelcomeMessage <- function(){
  paste(c("==============================\n",
          "MAPpoly Package",
          " [Version ", utils::packageDescription("mappoly")$Version,
           utils::packageDescription("mappoly")$Date, "]\n",
          "More information: https://github.com/mmollina/MAPpoly\n",
          "=============================="),
        sep="")
}
