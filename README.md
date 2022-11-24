# New perspectives on an old grouping: the genomic and phenotypic differences of Oxalobacter formigenes
Scripts to reproduce the analyses in "New perspectives on an old grouping: the genomic and phenotypic differences of Oxalobacter formigenes"

## Contents
### 16S rRNA gene analysis
`01_custom_dada2_script.R`: DADA2 pipeline for reads resulting from the custom demultiplexing script

`02_filtering.R`: filters out low-abundance and rare ASVs + low-depth samples

`03_decontam.R`: script for the decontam pipeline

`04_adonis.R`: script for Adonis2 (PERMANOVA) analysis

`05_ancom_asv.R`: script to run ANCOM

`06_maaslin2_asv.R`: script to run MaAsLin2

`07_differential_abundance_fig.R`: script to generate differential abudance plot for Figure 5

`08_picrust2.R`: runs PICRUSt2 for metabolic inferencing

`09_barplots.R`: generates bar plots on Figure 5

`10_PCA.R`: generates PCA in Figure 5

### Genome analysis
`01_Figure1.R`: generates Figure 1

`02_Figure2.R`: generates Figure 2

`03_Figure3.R`: generates Figure 3

`RUCS_primer_filter.R`: scripts to filter primers from RUCS

`oxalobacter_genome_analysis.txt`: scripts for main genome analysis

## Data Availability
To use these scripts, please download the demultiplexed reads from the Sequence Read Archive (PRJNA836912) and/or the genomes from GenBank (PRJNA841018).

## Type strain availability
Type strains can be found at the following culutre collection centers.

<i>O. formigenes</i> OxB<sup>T</sup> --> <a href="https://www.atcc.org/products/35274">ATCC 35274</a>

<i>O. aliformigenes</i> Va3<sup>T</sup>

<i>O. paeniformigenes</i> OxGP1<sup>T</sup>

<i>O. paraformigenes</i> HOxBLS<sup>T</sup>

<i>O. vibrioformis</i> WoOx3<sup>T</sup> --> <a href="https://www.dsmz.de/collection/catalogue/details/culture/DSM-5502">DSM 5502</a>

## Citation
Chmiel JA, Carr C, Stuivenberg GA, Venema R, Chanyi RM, Al KF, Giguere D, Say H, Akouris PP, Domínguez Romero SA, Kwong A, Tai V, Koval SF, Razvi H, Bjazevic J, Burton JP (2022) New perspectives on an old grouping: The genomic and phenotypic variability of <i>Oxalobacter formigenes</i> and the implications for calcium oxalate stone prevention. Front. Microbiol. 13:1011102. doi: 10.3389/fmicb.2022.1011102
