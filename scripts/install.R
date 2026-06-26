if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", lib="~/my_R_libs", repos="https://cloud.r-project.org")

BiocManager::install("lefser",lib="~/my_R_libs")


install.packages("ggh4x", lib="~/my_R_libs", repos="https://cloud.r-project.org")


devtools::install_github("teunbrand/legendry", lib="~/my_R_libs")