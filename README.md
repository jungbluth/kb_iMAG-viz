# kb_IMAG-viz

### Under active development as of July 2019 [use at your own risk!]

### Workflow Description

This script package uses genome annotation data as input to perform a taxa-specific assessment of genome quality.

The core steps are as follows:

1) Annotate your query genomes and reference genomes using the SAME tool for all genomes (e.g., RAST, Prokka, etc.).
2) Munge the genome annotation data such that files are organized to have one genome per line, with the text string of each annotation separated by tab.
3) Generate genome completeness and contamination estimates for all genomes using CheckM.
4) Generate taxonomic identifications for all genomes using GTDB-tk.
5) Combine annotation and genome quality information into a single large table.
6) Identify a taxonomic level (e.g., Phylum) and divide the single large table into multiple tables where each table corresponds to a taxa at that level (e.g., p__Crenarchaeaota, p__Altiarchaeota, p__Euryarchaeota, etc.).
7. For each taxa, generate a presence/absence count table for the annotations represented.
8. Perform dimensional reduction of annotation count tables.
9. Generate plots of dimensional reduction results and color by taxonomy and shape by genome type (i.e., isolate, SAG, MAG).

For the test data and instructions to run listed below, steps 1-4 were run previously, so effectively you are starting at step 5.


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


### Help Menu

>>> kb_iMAG-viz-workfloy.py -h
__________________________________________________________
|  _    _        _ __  __    _    ____            _      |
| | | _| |__    (_)  \/  |  / \  / ___|    __   _(_)____ |
| | |/ / '_ \   | | |\/| | / _ \| |  _ ____\ \ / / |_  / |
| |   <| |_) |  | | |  | |/ ___ \ |_| |_____\ V /| |/ /  |
| |_|\_\_.__/___|_|_|  |_/_/   \_\____|      \_/ |_/___| |
|          |_____|                                       |
|                                           version 0.1  |
|________________________________________________________|


usage: kb_iMAG-viz-workflow.py -i [query_annotation_table] --version

    kb_iMAG-viz is software to evaluate microbial genome quality against reference genomes.
    ----------------------------------------------------------------------------------------------------------------
    kb_iMAG-viz-workflow.py performs the following steps.
    The workflow goes as follows:
    STEP 1. Annotations data is read from a source (e.g. RAST), and converted into a tab-delimited text
    STEP 1a. Raw data is cleaned to remove extraneous characters
    STEP 2. Data are subsetted based on the taxonomic level of interest (e.g. domain, phylum, class, order, ...)
    STEP 3. Annotation data is converted to presence/absence data (for that taxonomic level)
    STEP 4: Dimensional reduction on annotation count data to identify major trends and outliers
    STEP 5: Plotting dimensional reduction results and associated genome quality information.

optional arguments:
  -h, --help            show this help message and exit
  -i QUERY_ANNOTATION_TABLE
                        Indicate a tab-delimited table where genome ID is the first column and subset columns correspond to gene names (default: /Applications/ResearchSoftware/kb_iMAG-viz/test/query-genomes/TARA-MAGs_Delmont-Archaea-only-2017.RAST.txt)
  --taxa_level TAXA_LEVEL
                        Indicate the taxonomic lineage. (default: Phylum)
  --save_master_table SAVE_MASTER_TABLE
                        Save the master merged table; warning it can be large as it contains an unflattened version of all the input data. (default: Yes)
  --generate_count_tables GENERATE_COUNT_TABLES
                        Regenerate annotation count tables (time-consuming step). (default: Yes)
  --dimensional_reduction_method DIMENSIONAL_REDUCTION_METHOD
                        Pick a dimensional reduction method. (default: pca)
  --plotting_method PLOTTING_METHOD
                        Pick a ploting method. (default: ggplot)
  --path_to_kb_imagviz PATH_TO_KB_IMAGVIZ
                        Indicate the path to kb_imagviz software
  --version             show program's version number and exit

