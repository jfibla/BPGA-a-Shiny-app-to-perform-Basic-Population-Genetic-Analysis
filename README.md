This Shiny app provides an interactive platform for visualizing population genetic structure. Users can upload their own datasets or use example files to explore population relationships, genetic diversity, and ancestral components. The app dynamically reads PLINK binary files and integrates them with worldwide reference populations from the 1000 Genomes Project (1000G) and the Human Genome Diversity Project (HGDP), enabling basic population analyses such as Principal Component Analysis (PCA), ADMIXTURE, and FST analyses.
The app performs PCA for dimensionality reduction, producing scatterplots that reveal population clusters. ADMIXTURE analysis is employed to assess population ancestry components, with a geographic visualization feature that displays these components as pie charts on a world map—each slice representing an ancestral proportion, and chart sizes scaled according to population size. Finally, FST analysis enables users to quantify genetic differentiation between pairs of populations.
Designed for both researchers and students in population genetics, the app offers a user-friendly interface for exploring genetic structure through multiple methods. It facilitates interactive data exploration and produces publication-ready plots. Overall, the app integrates PCA, ADMIXTURE, and FST analyses, providing a comprehensive view of population dynamics and evolutionary relationships.

References:
Alexander, D. H., Novembre, J. & Lange, K. Fast model-based estimation of ancestry in unrelated individuals. Genome Res 19, 1655–1664 (2009).
Auton, A. et al. A global reference for human genetic variation. Nature 526, 68–74 (2015).
Bergström, A. et al. Insights into human genetic variation and population history from 929 diverse genomes. Science 367, (2020).
Purcell, S. et al. PLINK: A Tool Set for Whole-Genome Association and Population-Based Linkage Analyses. The American Journal of Human Genetics 81, 559–575 (2007).
Yang, J., Lee, S. H., Goddard, M. E. & Visscher, P. M. GCTA: A Tool for Genome-wide Complex Trait Analysis. Am J Hum Genetics 88, 76–82 (2011).

Note: folder "www" app dependencies for iOS users. 
Other operating systems:
Download pages:
http://dalexander.github.io/admixture/download.html
https://www.cog-genomics.org/plink2/
https://yanglab.westlake.edu.cn/software/gcta/#Download

Available Docker container at https://doi.org/10.6084/m9.figshare.28009556.v1
