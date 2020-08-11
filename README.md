# Bee_euvir
### This repository contains the scripts used to analyse the Belgian bee virus sequencing data, and the results from the SRA screenings.

R scripts and seperate python scripts sit in the 'Scripts' folder.

BEAST XML files that were used to generate the phylogenetic trees, in the 'XML' folder
Notebooks are located in the 'Notebooks' folder and include:

  * Taxonomy parsing and network generation (minimzed nested block graph) for the Belgian bee viruses
  * The initial analysis of the mapping results against the SRA datasets
  * Taxonomy parsing and network generation (minimized nested block graph) for the SRA screening results
  * The code to generate the SRA abundance heatmap
  * The selection strategy & connected component analysis to get sets of sequences to include in phylogenies
  * The code used to determine virus sharing within Apidae and within hymenoptera samples

Furthermore, the fasta sequences for all viruses included in the phylogenies, as well as the complete contig set are available under the 'fasta' folder.

Some of the generated / included data sits under 'data' folder, the rest is hosted in Zenodo under (DOI 10.5281/zenodo.3979324):

The generated network CSV files are available under the 'networks' folder

The preferred order of running scripts are:

  * Rscripts
  * Tax.py
  * Taxonomy notebooks
 
