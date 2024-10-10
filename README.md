This Shiny app provides an interactive platform for visualizing population genetic structure. Users can upload their own datasets or use example files to explore population relationships, genetic diversity, and ancestral components. The app dynamically reads PLINK binary files and integrates them with worldwide reference populations from the 1000 Genomes Project (1000G) and the Human Genome Diversity Project (HGDP), enabling basic population analyses such as Principal Component Analysis (PCA), ADMIXTURE, and FST analyses.
The app performs PCA for dimensionality reduction, producing scatterplots that reveal population clusters. ADMIXTURE analysis is employed to assess population ancestry components, with a geographic visualization feature that displays these components as pie charts on a world map—each slice representing an ancestral proportion, and chart sizes scaled according to population size. Finally, FST analysis enables users to quantify genetic differentiation between pairs of populations.
Designed for both researchers and students in population genetics, the app offers a user-friendly interface for exploring genetic structure through multiple methods. It facilitates interactive data exploration and produces publication-ready plots. Overall, the app integrates PCA, ADMIXTURE, and FST analyses, providing a comprehensive view of population dynamics and evolutionary relationships.

In order to run this app on desktop R you can install the following libraries:
run: install.packages(c("shiny", "qqman", "readr", "data.table", "ggplot2", "shinyjs", "dplyr", "tidyr", "rmarkdown", "knitr", "ggforce", "mapplots", "maps"))
 
The Reference .bed file for : BPGA-a-Shiny-app-to-perform-Basic-Population-Genetic-Analysis
        
can be downloaded from https://doi.org/10.6084/m9.figshare.27203142.v1. Save file at "REFERENCE" folder
