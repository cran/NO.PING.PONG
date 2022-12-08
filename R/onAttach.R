
.onAttach<- function(libname, pkgname){
  packageStartupMessage("**************************************************************************************************\n",
                        pkgname," ",packageDescription("NO.PING.PONG")$Version,
                        "\n\nPlease contact Brian O'Connor at brian.oconnor@ubc.ca if you have questions or suggestions.\n",
                        "**************************************************************************************************", 
                        appendLF = TRUE)
}

NULL

