

# Building manuals and namespace using roxygen2
library(devtools)
document()


# Remove .Rhistory
unlink(".Rhistory")
unlink("R/.Rhistory")

