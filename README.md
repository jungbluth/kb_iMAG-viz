# kb_IMAG-viz

### Under active development as of July 2019

This script package uses genome annotation data as input to perform a taxa-specific assessment of genome, usually metagenome-assembled genome (MAG) quality.

Requirements (eventually the only requirement will be Docker, but for now one must install these Python and R packages manually):
* Python3
* Python packages: (pandas, numpy)
* R
* R packages (ggplot2, ggpubr)

To run:

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
