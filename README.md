# kb_IMAG-viz

### Under active development as of July 2019

This script package uses genome annotation data as input to perform a taxa-specific assessment of genome, usually metagenome-assembled genome (MAG) quality.

The core steps are as follows:

1) Annotate your query genomes and reference genomes using the SAME tool for all genomes (e.g. RAST, Prokka, etc.). 
2) Munge the genome annotation data such that files are organized to have one genome per line, with the text string of each annotation separated by tab
3) Generate genome completeness and contamination estimates for all genomes using CheckM
4) Generate taxonomic identifications for all genomes using GTDB-tk
5) Combine annotation and genome quality information into a single large table
6) Identify a taxonomic level (e.g. Phylum) and divide the single large table into multiple tables where each table corresponds to a taxa at that level (e.g. p__Crenarchaea, p__Altiarchaeota, p__Euryarchaeota, etc.)
7. For each taxa, generate a presence/absence count table for the annotations represented
8. Perform dimensional reduction of annotation count tables
9. Generate plots of dimensional reduction results and color by taxonomy and shape by genome type (i.e. Isolate, SAG, MAG)

For the test data and instructions to run listed below, steps 1-4 have already been completed so effectively you are starting at step 5.


### Installation

Requirements (eventually the only requirement will be Docker, but for now one must install these Python and R packages manually):
* Python3
* Python packages: (pandas, numpy)
* R
* R packages (ggplot2, ggpubr)


### Running Instructions

1) Clone this repo

```git clone https://github.com/jungbluth/kb_iMAG-viz```

2) Set application location as a variable

```PATH_TO_KB_IMAG_VIZ="/Applications/ResearchSoftware/kb_iMAG-viz"```

3) Change permissions to executable

```chmod +x ${PATH_TO_KB_IMAG_VIZ}/kb_iMAG-viz-workflow.py```

4) Optional: if running on the test data, forgo the time-consuming count-table generation step by copying the count-tables to your local directory. If doing this, then in Step 5 set the --generate_count_tables flag to 'n'.

```cp ${PATH_TO_KB_IMAG_VIZ}/test/output/*count-data* ./```

5) Run application

```
/Applications/ResearchSoftware/kb_iMAG-viz/kb_iMAG-viz-workflow.py \
-i ${PATH_TO_KB_IMAG_VIZ}/test/query-genomes/TARA-MAGs_Delmont-Archaea-only-2017.RAST.txt \
--taxa_level Phylum \
--save_master_table Yes \
--path_to_kb_imagviz ${PATH_TO_KB_IMAG_VIZ} \
--generate_count_tables n \
--dimensional_reduction_method pca \
--plotting_method ggplot
```

6) Sweet, it worked! Grab a beer and review the newly-produced pdf files to learn something about your genomez. :)
