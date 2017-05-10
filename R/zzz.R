THISPKG <- 'CancerInSilico'
.onAttach <- function(libname, pkgname)
{
    version <- packageDescription('CancerInSilico', fields='Version')
	packageStartupMessage(paste("Welcome to CancerInSilico version ", version, "\n", sep = "" ) )
}
