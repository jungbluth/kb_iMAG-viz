# kb_IMAG-viz

### Under active development as of July 2019

To run:

1) Clone this repo

```git clone https://github.com/jungbluth/kb_iMAG-viz```

2) Change permissions to executable

```chmod +x kb_iMAG-viz-workflow.py```

3) Set application location as a variable

```PATH_TO_KB_IMAG_VIZ="/Applications/ResearchSoftware/kb_iMAG-viz"```

4) Run application

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