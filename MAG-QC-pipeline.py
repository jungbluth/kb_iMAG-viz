#!/usr/bin/env python3

__author__ = "Sean P. Jungbluth"
__copyright__ = "Copyright 2019"
__license__ = "GPL 3.0"
__maintainer__ = "Sean P. Jungbluth"
__email__ = "jungbluth.sean@gmail.com"

import sys
import time
import pandas
import numpy
import os
import re
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt


# check to see if python3 installed and being used; should be automatic because of the shebang, but useful if someone calls explicitly with python 2
def check_for_py3():
    print("Running check_for_py3")
    time_start = time.time()
    try:
        if sys.version_info.major != 3:
            sys.stderr.write('\nError: python version 3 is required - you have python version %d.\n\n' % sys.version_info.major)
            sys.exit(-1)
    except Exception:
        sys.stderr.write("Unable to detect python version - assuming python 3 installed.\n\n")
    print('check_for_py3 done! Time elapsed: {} seconds'.format(time.time() - time_start) + '\n')

# From Miriam: function to return object ID based on object name
def get_obj_id(obj_name):
    """
    Example: get_obj_id("MAG-QC_Archaea.SAGs.Prokka")
    """
    from biokbase.workspace.client import Workspace
    import os
    ws = Workspace('https://kbase.us/services/ws')
    ws_name = os.environ['KB_WORKSPACE_ID']
    try:
        obj_id = ws.get_object_info3({'objects': [{'workspace': ws_name, 'name': obj_name}]})['paths'][0][0]
        return obj_id
    except:
        return False


# From Miriam: function to return object data based on object ID
def get_object_data(object_id):
    """
    Fetch data from the workspace. Example1: get_object_data(get_obj_id("MAG-QC_Archaea.SAGs.RAST")); Example2: get_object_data(u'43402/2132/1')
    """
    from biokbase.workspace.client import Workspace
    ws = Workspace('https://kbase.us/services/ws')
    return ws.get_objects([{"ref": object_id}])[0]

# todo: double check to see if indexed properly
def extract_annotations_from_genomeset(genomeset_name, export_filename):
    print("Running extract_annotations_from_genomeset")
    time_start = time.time()
    from biokbase.narrative.jobs.appmanager import AppManager
    ws = biokbase.narrative.clients.get("workspace")
    for genome in get_object_data(get_obj_id(genomeset_name)).values()[3]:
        functionList = ws.get_objects([{'ref': genome}])[0]['data']['cdss']
        df = []
        for function in range(len(functionList)):
            df.append(functionList[function]['functions'])
        with open(export_filename, 'a') as f:
            if get_object_data(genome).values()[0][1] != "2528311097___MAG-QC_Archaea__Isolates.RAST":  # this genome has something weird going on in the annotation; result is it breaks the JSON file structure, suspected ({,},') symbols, need to troubleshoot
                f.write(get_object_data(genome).values()[0][1] + '\t' + str(df) + '\n')
        f.close()
    print('extract_annotations_from_genomeset done! Time elapsed: {} seconds'.format(time.time() - time_start) + '\n')


# wont be needed when annotation data pulled directly from Genome Set Objects
def read_and_parse_rast_annotations(inputfile, outputfile):
    print("Running read_and_parse_rast_annotations")
    time_start = time.time()
    f1 = open(inputfile)
    if Path(outputfile).is_file():
        os.remove(outputfile)
    f2 = open(outputfile, "w+")
    f2 = open(outputfile, "a")
    ftrim = f1.readlines()[1:]
    for line in ftrim:
        rep = {"], [u": "|", "', u'": "|", '"': "", "[[u": ""}
        rep = dict((re.escape(k), v) for k, v in rep.items())
        pattern = re.compile("|".join(rep.keys()))
        line = pattern.sub(lambda m: rep[re.escape(m.group(0))], line)
        for j in range(len(line.split('\t')[1].split("|")) - 1):  # last entry needs to be removed
            if "''" != str(line.split('\t')[1].split("|")[j]): # write only non-empty lines
                if "hypothetical protein" not in line.split('\t')[1].split("|")[j]: # not working
                    if "unknown" not in line.split('\t')[1].split("|")[j]: # not working
                        if "conserved protein" not in line.split('\t')[1].split("|")[j]: # not working
                            if "Conserved protein" not in line.split('\t')[1].split("|")[j]: # not working
                                if "predicted protein" not in line.split('\t')[1].split("|")[j]: # not working
                                    if (str(line.split('\t')[1].split("|")[j])[0] == "'") and (str(line.split('\t')[1].split("|")[j])[-1] == "'"):
                                        f2.write(line.split('\t')[0].split('___')[0]+'\t'+line.split('\t')[1].split("|")[j][1:-1]+'\n')
                                    elif (str(line.split('\t')[1].split("|")[j])[0] == "'") and (str(line.split('\t')[1].split("|")[j])[-1] != "'"):
                                        f2.write(line.split('\t')[0].split('___')[0]+'\t'+line.split('\t')[1].split("|")[j][1:]+'\n')
                                    elif (str(line.split('\t')[1].split("|")[j])[0] != "'") and (str(line.split('\t')[1].split("|")[j])[-1] == "'"):
                                        f2.write(line.split('\t')[0].split('___')[0]+'\t'+line.split('\t')[1].split("|")[j][:-1]+'\n')
                                    elif (str(line.split('\t')[1].split("|")[j])[0] != '"') and (str(line.split('\t')[1].split("|")[j])[-1] == '"'):
                                        f2.write(line.split('\t')[0].split('___')[0]+'\t'+line.split('\t')[1].split("|")[j][1:-1]+'\n')
                                    else:
                                        f2.write(line.split('\t')[0].split('___')[0]+'\t'+line.split('\t')[1].split("|")[j]+'\n')
    f1.close()
    f2.close()
    print('read_and_parse_rast_annotations done! Time elapsed: {} seconds'.format(time.time() - time_start) + '\n')



def combine_external_checkm_and_taxonomy_info(file1, file2, file3, output1):
    print("Running combine_external_checkm_and_taxonomy_info")
    time_start = time.time()
    filenames = [file1, file2, file3]
    with open(output1, 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
    print('combine_external_checkm_and_taxonomy_info done! Time elapsed: {} seconds'.format(time.time() - time_start) + '\n')


def import_and_merge_tables(saveoutput):
    print("Running import_and_merge_tables")
    time_start = time.time()
    f1 = pandas.read_csv(file1, header=None, sep="\t")
    f2 = pandas.read_csv(file2, header=None, sep="\t")
    f3 = pandas.read_csv(file3, header=None, sep="\t")
    d1 = pandas.read_csv(data1, header=None, sep="\t")
    d2 = pandas.read_csv(data2, header=None, sep="\t")
    d3 = pandas.read_csv(data3, header=None, sep="\t")
    q1 = pandas.read_csv(query1, header=None, sep="\t")
    f123 = pandas.concat([f1, f2, f3], axis=0)
    d123q1 = pandas.concat([d1, d2, d3, q1], axis=0)
    global merge
    merge = f123.merge(d123q1, left_on=0, right_on=0)
    merge.columns = ['genomeID', 'genomeSet', 'V1', 'V2', 'V3', 'V4', 'GenomeType', 'V6', 'Completeness', 'Contamination', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'RAST_Annotation']
    if saveoutput == "Yes":
        merge.to_csv("_Master-table.tsv")
    print('import_and_merge_tables done! Time elapsed: {} seconds'.format(time.time() - time_start) + '\n')


def list_lineages_for_selected_level(taxa_level):
    if taxa_level == 'Domain':
        taxa_level_list = merge.Domain.unique()
        taxa_level_down1 = 'Phylum'
        taxa_level_down2 = 'Class'
    elif taxa_level == 'Phylum':
        taxa_level_list = merge.Phylum.unique()
        taxa_level_down1 = 'Class'
        taxa_level_down2 = 'Order'
    elif taxa_level == 'Class':
        taxa_level_list = merge.Class.unique()
        taxa_level_down1 = 'Order'
        taxa_level_down2 = 'Family'
    elif taxa_level == 'Order':
        taxa_level_list = merge.Order.unique()
        taxa_level_down1 = 'Family'
        taxa_level_down2 = 'Genus'
    elif taxa_level == 'Family':
        taxa_level_list = merge.Family.unique()
        taxa_level_down1 = 'Genus'
        taxa_level_down2 = 'Species'
    else:
        taxa_level_list = merge.Genus.unique()
        taxa_level_down1 = 'Species'
        taxa_level_down2 = ''
    return taxa_level_list, taxa_level_down1, taxa_level_down2


# merge example
# genomeID    genomeSet   V1  V2  V3  ... Order   Family  Genus   Species RAST_Annotation
# 1861359 GCA_002509575.1_ASM250957v1_genomic PRJNA348753-8000Genomes y   y   y   ... o__Methanomicrobiales   f__Methanocullaceae g__Methanoculleus   s__GCA_002508705.1  Acyl-CoA    synthetase  (NDP    forming)
# 1861360 GCA_002509575.1_ASM250957v1_genomic PRJNA348753-8000Genomes y   y   y   ... o__Methanomicrobiales   f__Methanocullaceae g__Methanoculleus   s__GCA_002508705.1  Iron(III)   dicitrate-binding   protein
# 1861361 GCA_002509575.1_ASM250957v1_genomic PRJNA348753-8000Genomes y   y   y   ... o__Methanomicrobiales   f__Methanocullaceae g__Methanoculleus   s__GCA_002508705.1  Iron(III)   dicitrate   transport   system  permease    ...
# 1861362 GCA_002509575.1_ASM250957v1_genomic PRJNA348753-8000Genomes y   y   y   ... o__Methanomicrobiales   f__Methanocullaceae g__Methanoculleus   s__GCA_002508705.1  Iron(III)   dicitrate   transport   ATP-binding protein
# 1861363 GCA_002509575.1_ASM250957v1_genomic PRJNA348753-8000Genomes y   y   y   ... o__Methanomicrobiales   f__Methanocullaceae g__Methanoculleus   s__GCA_002508705.1  Diadenosine 5'5'''-P1,P4-tetraphosphate pyroph...


def extract_lineages_from_table(taxa_level):
    print("Running extract_lineages_from_table")
    time_start = time.time()
    global merge
    global merge_reduced
    for lineage in range(len(list_lineages_for_selected_level(taxa_level)[0])):
        # print(lineage)
        print(list_lineages_for_selected_level(taxa_level)[0][lineage])
        # merge = merge.loc[merge[taxa_level] == str(list_lineages_for_selected_level(taxa_level)[lineage])]
    # print("merge " + str(type(merge)))
    # print(merge.head())
    merge = merge.loc[merge[taxa_level] == "p__Halobacterota"]
    # print(list_lineages_for_selected_level(taxa_level)[1])
    # down1list = merge[list_lineages_for_selected_level(taxa_level)[1]].to_string(index=False)
    merge_reduced = merge.drop(columns="RAST_Annotation").drop_duplicates() # per genomeID, extra information. table is the length of the object used for dimensional reduction
    # print(merge_reduced)
    # list_taxa_down_one_level = merge_reduced[list_lineages_for_selected_level(taxa_level)[1]]
    # list_genomeid = merge_reduced["genomeID"]
    # print("down1list length"+str(len(down1list)))
    # print(merge.head)
    print('extract_lineages_from_table done! Time elapsed: {} seconds'.format(time.time() - time_start) + '\n')


def count_annotation_data_for_level(merge, outputfile1):
    print("Running count_annotation_data_for_level")
    time_start = time.time()
    if Path(outputfile1).is_file():
        os.remove(outputfile1)
    f1 = open(outputfile1, "w+")
    f1 = open(outputfile1, "a")
    genomelist = merge.genomeID.unique()
    annotationlist = merge.RAST_Annotation.unique()
    f1.write("Annotation" + "\t")  # write header line
    for genome in range(len(genomelist)):
        f1.write(str(genomelist[genome]))
        if genome != (len(genomelist) - 1):
            f1.write("\t")
        else:
            f1.write("\n")
    # for annotation in range(len(annotationlist)): # takes long to run (up to hours)
    for annotation in range(0, 20):
        f1.write(str(annotationlist[annotation]) + "\t")
        temp_merge = merge.loc[merge['RAST_Annotation'] == annotationlist[annotation]]
        dat = temp_merge['genomeID'].value_counts().rename_axis('genomeID').reset_index(name='counts')
        for genome in range(len(genomelist)):
            dat = temp_merge['genomeID'].rename_axis('V1').reset_index(name='genomeID')
            f1.write(str(dat['genomeID'].str.contains(str(genomelist[genome])).sum()))
            if genome != (len(genomelist) - 1):
                f1.write("\t")
            else:
                f1.write("\n")
    print('count_annotation_data_for_level done! Time elapsed: {} seconds'.format(time.time() - time_start) + '\n')


def import_count_and_combine_with_genome_metadata():
    global dat
    global dat_combined
    with open(inputfile1) as f:
        ncols = len(f.readline().split('\t'))
    dat = numpy.loadtxt(inputfile1, delimiter="\t", skiprows=1, usecols=range(1, ncols))
    dat = numpy.transpose(dat)  # transpose data
    print("***" + str(len(merge_reduced)))
    print("***" + str(len(dat)))


def run_tsne_dimensional_reduction(inputfile1):
    print("Running run_tsne_dimensional_reduction")
    time_start = time.time()
    from sklearn.manifold import TSNE
    tsne = TSNE(n_components=2, perplexity=40.0, early_exaggeration=12.0, learning_rate=200.0, n_iter=300, n_iter_without_progress=300, min_grad_norm=1e-07, metric="euclidean", init="random", verbose=1, random_state=None, method="barnes_hut", angle=0.5)
    tsne_results = tsne.fit_transform(dat)  # fit into an embedded space
    global dat1
    dat1 = pandas.DataFrame(data=tsne_results, columns=['tsne-2d-one', 'tsne-2d-two'])
    dat1_combined = pandas.DataFrame(numpy.concatenate((dat1, merge_reduced),axis = 1)) # combine count table with genome metadata
    dat1_combined.columns = ['tsne-2d-one', 'tsne-2d-two', 'genomeID', 'genomeSet', 'V1', 'V2', 'V3', 'V4', 'GenomeType', 'V6', 'Completeness', 'Contamination', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    print('run_tsne_dimensional_reduction done! Time elapsed: {} seconds'.format(time.time() - time_start) + '\n')


def run_principal_component_analysis(inputfile1):
    print("Running run_principal_component_analysis")
    time_start = time.time()
    from sklearn.decomposition import PCA
    pca = PCA(n_components=3)
    pca_result = pca.fit_transform(dat)
    global dat1_combined
    dat1 = pandas.DataFrame(data=pca_result, columns=['PC1', 'PC2', 'PC3'])
    dat1_combined = pandas.DataFrame(numpy.concatenate((dat1, merge_reduced),axis = 1)) # combine count table with genome metadata
    dat1_combined.columns = ['PC1', 'PC2', 'PC3', 'genomeID', 'genomeSet', 'V1', 'V2', 'V3', 'V4', 'GenomeType', 'V6', 'Completeness', 'Contamination', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    print(dat1_combined)
    print(dat1_combined.shape)
    print('run_principal_component_analysis done! Time elapsed: {} seconds'.format(time.time() - time_start) + '\n')


def plot_dimensional_reduction_results_seaborn(mode):
    print("Running plot_dimensional_reduction_results_seaborn")
    time_start = time.time()
    import seaborn
    # list_genomeid = merge_reduced["genomeID"]
    if mode == "tsne":
        sns_plot = seaborn.scatterplot(x="tsne-2d-one", y="tsne-2d-two", hue="Class", style=None, size=None, data=dat1_combined, palette=seaborn.color_palette("hls", len(dat1_combined["Class"].unique())), hue_order=None, hue_norm=None, sizes=None, size_order=None, size_norm=None, markers=dat1_combined["Order"], style_order=None, x_bins=None, y_bins=None, units=None, estimator=None, ci=95, n_boot=1000, alpha=0.3, x_jitter=None, y_jitter=None, legend='full', ax=None)
    else:
        # print(type(merge_reduced["Class"]))
        # print(type(merge_reduced))
        df = pandas.DataFrame(merge_reduced)
        # print(type(df))
        # print(df)
        #sns_plot = seaborn.scatterplot(x="PC1", y="PC2", hue=df["Class"], style=None, size=None, data=dat1, palette=seaborn.color_palette("hls", len(df["Class"].unique())), hue_order=None, hue_norm=None, sizes=None, size_order=None, size_norm=None, markers=True, style_order=None, x_bins=None, y_bins=None, units=None, estimator=None, ci=95, n_boot=1000, alpha=0.3, x_jitter=None, y_jitter=None, legend='full', ax=None)
        sns_plot = seaborn.scatterplot(x="PC1", y="PC2", hue="Class", style=None, size=None, data=dat1_combined, palette=seaborn.color_palette("hls", len(dat1_combined["Class"].unique())), hue_order=None, hue_norm=None, sizes=None, size_order=None, size_norm=None, markers=dat1_combined["Order"], style_order=None, x_bins=None, y_bins=None, units=None, estimator=None, ci=95, n_boot=1000, alpha=0.3, x_jitter=None, y_jitter=None, legend='full', ax=None)
    fig = sns_plot.get_figure()
    fig.savefig("output.png")
    # plot_ex.show()
    # plot_ex.savefig('test.png')
    # plot_ex.savefig('test.pdf')
    print('plot_dimensional_reduction_results_seaborn done! Time elapsed: {} seconds'.format(time.time() - time_start) + '\n')


def plot_dimensional_reduction_results_rggplot(mode):
    print()
    import rpy2
    from rpy2 import robjects
    from rpy2.robjects.packages import importr

    # The R 'print' function
    ggplot2 = importr('ggplot2')
    # grdevices = importr('grDevices')
    # base = importr('base')
    # datasets = importr('datasets')

    grid.activate()


def test_text():
    print(merge.tail(n=5))
    print(merge.shape)
    print(merge.shape[0])
    for col in merge.columns:
        print(col)
    # print(type(merge['RAST_Annotation'].value_counts()))
    # print(f123.columns)
    # print(d123q1.columns)
    # print(d1.tail(n=5))
    # print(d2.tail(n=5))
    # print(d3.tail(n=5))
    # print("f123\n"+str(f123.shape))
    # print("d123q1\n"+str(d123q1.shape))
    # print(merge.tail(n=5))
    # print(merge.shape)
    # print(merge.Domain.unique())


# run main MAG-QC-pipeline.py script
if __name__ == "__main__":

    inputfile = os.path.realpath(sys.argv[1])
    outputfile = os.path.realpath(sys.argv[2])
    outputfile1 = os.path.realpath(sys.argv[3])
    inputfile1 = os.path.realpath(sys.argv[3])


    # files needed for function combine_external_checkm_and_taxonomy_info AND import_and_merge_tables
    file1 = os.path.realpath("Delmont_genomeQC-data.tsv")
    file2 = os.path.realpath("IMG_genomeQC-data.tsv")
    file3 = os.path.realpath("Other_genomeQC-data.tsv")
    output1 = os.path.realpath("_All_genomeQC-data.tsv")

    # files needed for function import_and_merge_tables
    data1 = os.path.realpath("MAG-QC_Archaea.Isolates_clean.RAST.txt")
    data2 = os.path.realpath("MAG-QC_Archaea.MAGs_clean.RAST.txt")
    data3 = os.path.realpath("MAG-QC_Archaea.SAGs_clean.RAST.txt")
    query1 = os.path.realpath("TARA-MAGs_Delmont-Archaea-only-2017_clean.RAST.txt")


    # declare variables
    saveoutput = "Yes"
    taxa_level = "Phylum"


    check_for_py3()

    # read_and_parse_rast_annotations(inputfile, outputfile)

    # combine_external_checkm_and_taxonomy_info(file1, file2, file3, output1)

    # import_and_merge_tables(saveoutput)

    # extract_lineages_from_table(taxa_level)

    # # count_annotation_data_for_level(merge, outputfile1)

    # import_count_and_combine_with_genome_metadata()

    # run_tsne_dimensional_reduction(inputfile1)

    # run_principal_component_analysis(inputfile1)

    # plot_dimensional_reduction_results_seaborn(mode="pca")

    plot_dimensional_reduction_results_rggplot()

    # test_text()
