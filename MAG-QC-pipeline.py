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


#todo
#convert matrix to presence/absence

# check to see if python3 installed and being used; should be automatic because of the shebang, but useful if someone calls explicitly with python 2
def check_for_py3():
    print("\n"+"Running check_for_py3")
    time_start = time.time()
    try:
        if sys.version_info.major != 3:
            sys.stderr.write('\nError: python version 3 is required - you have python version %d.\n\n' % sys.version_info.major)
            sys.exit(-1)
    except Exception:
        sys.stderr.write("Unable to detect python version - assuming python 3 installed.\n\n")
    print('check_for_py3 done! Time elapsed: ' + '{}'.format(time.time() - time_start)[:7] + ' seconds\n')

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
    print("\n"+"Running extract_annotations_from_genomeset")
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
    print('extract_annotations_from_genomeset done! Time elapsed: ' + '{}'.format(time.time() - time_start)[:7] + ' seconds\n')


# wont be needed when annotation data pulled directly from Genome Set Objects
def read_and_parse_rast_annotations(inputfile, outputfile):
    print("\n"+"Running read_and_parse_rast_annotations")
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
    print('read_and_parse_rast_annotations done! Time elapsed: ' + '{}'.format(time.time() - time_start)[:7] + ' seconds\n')



def combine_external_checkm_and_taxonomy_info(file1, file2, file3):
    print("\n"+"Running combine_external_checkm_and_taxonomy_info")
    time_start = time.time()
    filenames = [file1, file2, file3]
    with open(output1, 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
    print("\nTotal number of genomes in metadata table: " + str(len(open(output1).readlines())))
    print('combine_external_checkm_and_taxonomy_info done! Time elapsed: ' + '{}'.format(time.time() - time_start)[:7] + ' seconds\n')


def import_and_merge_tables(saveoutput):
    print("\n"+"Running import_and_merge_tables")
    time_start = time.time()
    f1 = pandas.read_csv(file1, header=None, sep="\t", index_col=False, dtype=str)
    global f2
    f2 = pandas.read_csv(file2, header=None, sep="\t", index_col=False, dtype=str)
    f3 = pandas.read_csv(file3, header=None, sep="\t", index_col=False, dtype=str)
    d1 = pandas.read_csv(data1, header=None, sep="\t", index_col=False, dtype=str)
    d2 = pandas.read_csv(data2, header=None, sep="\t", index_col=False, dtype=str)
    d3 = pandas.read_csv(data3, header=None, sep="\t", index_col=False, dtype=str)
    q1 = pandas.read_csv(query1, header=None, sep="\t", index_col=False, dtype=str)
    f123 = pandas.concat([f1, f2, f3], axis=0)
    d123q1 = pandas.concat([d1, d2, d3, q1], axis=0)
    global merge
    merge = f123.merge(d123q1, how='outer', left_on=0, right_on=0)
    merge.columns = ['genomeID', 'genomeSet', 'V1', 'V2', 'V3', 'V4', 'GenomeType', 'V6', 'Completeness', 'Contamination', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'RAST_Annotation']
    if saveoutput == "Yes":
        merge.to_csv("MAG-QC-output_all-merged-data.tsv")
    print('import_and_merge_tables done! Time elapsed: ' + '{}'.format(time.time() - time_start)[:7] + ' seconds\n')


def extract_lineages_for_selected_level(taxa_level):
    print("\n"+"Running extract_lineages_for_selected_level")
    time_start = time.time()
    if taxa_level == 'Domain':
        taxa_level_list = merge.Domain.unique()
        taxa_level_down1 = 'Phylum'
        taxa_level_down2 = 'Class'
        #query_unique = f2.Domain.unique()
    elif taxa_level == 'Phylum':
        taxa_level_list = merge.Phylum.unique()
        taxa_level_down1 = 'Class'
        taxa_level_down2 = 'Order'
        #query_unique = f2.Phylum.unique()
    elif taxa_level == 'Class':
        taxa_level_list = merge.Class.unique()
        taxa_level_down1 = 'Order'
        taxa_level_down2 = 'Family'
        #query_unique = f2.Class.unique()
    elif taxa_level == 'Order':
        taxa_level_list = merge.Order.unique()
        taxa_level_down1 = 'Family'
        taxa_level_down2 = 'Genus'
        #query_unique = f2.Order.unique()
    elif taxa_level == 'Family':
        taxa_level_list = merge.Family.unique()
        taxa_level_down1 = 'Genus'
        taxa_level_down2 = 'Species'
        #query_unique = f2.Family.unique()
    else:
        taxa_level_list = merge.Genus.unique()
        taxa_level_down1 = 'Species'
        taxa_level_down2 = ''
        #query_unique = f2.Genus.unique()
    return taxa_level_list, taxa_level_down1, taxa_level_down2
    print('extract_lineages_for_selected_level done! Time elapsed: ' + '{}'.format(time.time() - time_start)[:7] + ' seconds\n')


# merge example
# genomeID    genomeSet   V1  V2  V3  ... Order   Family  Genus   Species RAST_Annotation
# 1861359 GCA_002509575.1_ASM250957v1_genomic PRJNA348753-8000Genomes y   y   y   ... o__Methanomicrobiales   f__Methanocullaceae g__Methanoculleus   s__GCA_002508705.1  Acyl-CoA    synthetase  (NDP    forming)
# 1861360 GCA_002509575.1_ASM250957v1_genomic PRJNA348753-8000Genomes y   y   y   ... o__Methanomicrobiales   f__Methanocullaceae g__Methanoculleus   s__GCA_002508705.1  Iron(III)   dicitrate-binding   protein
# 1861361 GCA_002509575.1_ASM250957v1_genomic PRJNA348753-8000Genomes y   y   y   ... o__Methanomicrobiales   f__Methanocullaceae g__Methanoculleus   s__GCA_002508705.1  Iron(III)   dicitrate   transport   system  permease    ...
# 1861362 GCA_002509575.1_ASM250957v1_genomic PRJNA348753-8000Genomes y   y   y   ... o__Methanomicrobiales   f__Methanocullaceae g__Methanoculleus   s__GCA_002508705.1  Iron(III)   dicitrate   transport   ATP-binding protein
# 1861363 GCA_002509575.1_ASM250957v1_genomic PRJNA348753-8000Genomes y   y   y   ... o__Methanomicrobiales   f__Methanocullaceae g__Methanoculleus   s__GCA_002508705.1  Diadenosine 5'5'''-P1,P4-tetraphosphate pyroph...


def subset_data_by_lineage(lineage):
    print("\n\t"+"Running subset_data_by_lineage")
    time_start = time.time()
    global merge_reduced
    global merge_reduced_trim
    merge_reduced = merge.loc[merge[taxa_level] == lineage]
    merge_reduced_trim = merge_reduced.drop(columns="RAST_Annotation").drop_duplicates() # per genomeID, extra information.
    print('\tsubset_data_by_lineage done! Time elapsed: ' + '{}'.format(time.time() - time_start)[:7] + ' seconds\n')


def count_annotation_data_for_level(merge_reduced, outputfile1, lineage):
    print("\n\t"+"Running count_annotation_data_for_level")
    time_start = time.time()
    if Path(outputfile1+"_"+lineage+".tsv").is_file():
        os.remove(outputfile1+"_"+lineage+".tsv")

    f1 = open(outputfile1+"_"+lineage+".tsv", "w+")
    f1 = open(outputfile1+"_"+lineage+".tsv", "a")
    genomelist = merge_reduced.genomeID.unique()
    annotationlist = merge_reduced.RAST_Annotation.unique()
    f1.write("Annotation" + "\t")  # write header line
    for genome in range(len(genomelist)):
        f1.write(str(genomelist[genome]))
        if genome != (len(genomelist) - 1):
            f1.write("\t")
        else:
            f1.write("\n")
    for annotation in range(len(annotationlist)): # takes long to run (up to hours)
    #for annotation in range(0, 20):
        f1.write(str(annotationlist[annotation]) + "\t")
        temp_merge = merge_reduced.loc[merge_reduced['RAST_Annotation'] == annotationlist[annotation]]
        #dat = temp_merge['genomeID'].value_counts().rename_axis('genomeID').reset_index(name='counts')
        for genome in range(len(genomelist)):
            dat = temp_merge['genomeID'].value_counts().rename_axis('genomeID').reset_index(name='counts')
            f1.write(str(dat['genomeID'].astype(str).str.contains(str(genomelist[genome])).sum()))
            if genome != (len(genomelist) - 1):
                f1.write("\t")
            else:
                f1.write("\n")
    print('\tcount_annotation_data_for_level done! Time elapsed: ' + '{}'.format(time.time() - time_start)[:7] + ' seconds\n')


def import_count_and_combine_with_genome_metadata(lineage):
    global dat
    global dat_combined
    global genome_number
    with open(inputfile1+"_"+lineage+".tsv") as f:
        ncols = len(f.readline().split('\t'))
    dat = numpy.loadtxt(inputfile1+"_"+lineage+".tsv", delimiter="\t", skiprows=1, usecols=range(1, ncols))
    dat = numpy.transpose(dat)  # transpose data
    genome_number = len(merge_reduced_trim)
    print(dat)
    print(merge_reduced_trim)
    print("\tNumber of genomes: " + str(len(merge_reduced_trim)))


def run_tsne_dimensional_reduction():
    print("\n\t"+"Running run_tsne_dimensional_reduction")
    time_start = time.time()
    from sklearn.manifold import TSNE
    tsne = TSNE(n_components=2, perplexity=40.0, early_exaggeration=12.0, learning_rate=200.0, n_iter=300, n_iter_without_progress=300, min_grad_norm=1e-07, metric="euclidean", init="random", verbose=1, random_state=None, method="barnes_hut", angle=0.5)
    tsne_results = tsne.fit_transform(dat)  # run tsne, fit into an embedded space
    global dat1
    dat1 = pandas.DataFrame(data=tsne_results, columns=['tsne-2d-one', 'tsne-2d-two'])
    dat1_combined = pandas.DataFrame(numpy.concatenate((dat1, merge_reduced_trim),axis = 1)) # combine count table with genome metadata
    dat1_combined.columns = ['tsne-2d-one', 'tsne-2d-two', 'genomeID', 'genomeSet', 'V1', 'V2', 'V3', 'V4', 'GenomeType', 'V6', 'Completeness', 'Contamination', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    print('\trun_tsne_dimensional_reduction done! Time elapsed: ' + '{}'.format(time.time() - time_start)[:7] + ' seconds\n')


def run_principal_component_analysis():
    print("\n\t"+"Running run_principal_component_analysis")
    time_start = time.time()
    from sklearn.decomposition import PCA
    pca = PCA(n_components=3)
    pca_result = pca.fit_transform(dat)
    global dat1_combined
    dat1 = pandas.DataFrame(data=pca_result, columns=['PC1', 'PC2', 'PC3'])
    dat1_combined = pandas.DataFrame(numpy.concatenate((dat1, merge_reduced_trim),axis = 1)) # combine count table with genome metadata
    dat1_combined.columns = ['PC1', 'PC2', 'PC3', 'genomeID', 'genomeSet', 'V1', 'V2', 'V3', 'V4', 'GenomeType', 'V6', 'Completeness', 'Contamination', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    #print(dat1_combined)
    #print(dat1_combined.shape)
    dat1_combined.to_csv("MAG-QC-output_PCA-input-data_"+str(lineage)+".tsv", sep='\t')
    print('\trun_principal_component_analysis done! Time elapsed: ' + '{}'.format(time.time() - time_start)[:7] + ' seconds\n')


def plot_dimensional_reduction_results_seaborn(lineage, mode):
    print("\n\t"+"Running plot_dimensional_reduction_results_seaborn")
    time_start = time.time()
    import seaborn
    # list_genomeid = merge_reduced["genomeID"]
    if mode == "tsne":
        sns_plot = seaborn.scatterplot(x="tsne-2d-one", y="tsne-2d-two", hue="Class", style=None, size=None, data=dat1_combined, palette=seaborn.color_palette("hls", len(dat1_combined["Class"].unique())), hue_order=None, hue_norm=None, sizes=None, size_order=None, size_norm=None, markers=dat1_combined["Class"], style_order=None, x_bins=None, y_bins=None, units=None, estimator=None, ci=95, n_boot=1000, alpha=0.3, x_jitter=None, y_jitter=None, legend='full', ax=None)
    else:
        # print(type(merge_reduced["Class"]))
        df = pandas.DataFrame(merge_reduced_trim)
        #sns_plot = seaborn.scatterplot(x="PC1", y="PC2", hue=df["Class"], style=None, size=None, data=dat1, palette=seaborn.color_palette("hls", len(df["Class"].unique())), hue_order=None, hue_norm=None, sizes=None, size_order=None, size_norm=None, markers=True, style_order=None, x_bins=None, y_bins=None, units=None, estimator=None, ci=95, n_boot=1000, alpha=0.3, x_jitter=None, y_jitter=None, legend='full', ax=None)
        sns_plot = seaborn.scatterplot(x="PC1", y="PC2", hue="Class", style=None, size=None, data=dat1_combined, palette=seaborn.color_palette("hls", len(dat1_combined["Class"].unique())), hue_order=None, hue_norm=None, sizes=None, size_order=None, size_norm=None, markers=dat1_combined["Class"], style_order=None, x_bins=None, y_bins=None, units=None, estimator=None, ci=95, n_boot=1000, alpha=0.3, x_jitter=None, y_jitter=None, legend='full', ax=None)
    fig = sns_plot.get_figure()
    fig.savefig("MAG-QC-output_PCA_"+str(lineage)+".png")
    print('\tplot_dimensional_reduction_results_seaborn done! Time elapsed: ' + '{}'.format(time.time() - time_start)[:7] + ' seconds\n')


def plot_dimensional_reduction_results_rggplot(lineage, mode):
    print("\n\t"+"Running plot_dimensional_reduction_results_rggplot")
    import shutil
    import os
    shutil.copy("MAG-QC-output_PCA-input-data_"+str(lineage)+".tsv", "R-input-data.tsv")
    time_start = time.time()
    import subprocess
    command = "/Applications/ResearchSoftware/kb_MAG-QC/make_ggplot_scatterplot.R"
    process=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process.wait()
    out, err = process.communicate()
    if process.returncode != 0: sys.exit("*Error generating figure with ggplot2")
    shutil.move("R-output-plot.pdf","MAG-QC-output_PCA_"+str(lineage)+".pdf")
    os.remove("R-input-data.tsv")
    print('\tplot_dimensional_reduction_results_rggplot done! Time elapsed: ' + '{}'.format(time.time() - time_start)[:7] + ' seconds\n')


def test_text():
    print(merge.tail(n=5))
    print(merge.shape)
    print(merge.shape[0])
    for col in merge.columns:
        print(col)
    # print("f123\n"+str(f123.shape))
    # print("d123q1\n"+str(d123q1.shape))
    # print(merge.tail(n=5))
    # print(merge.shape)
    # print(merge.Domain.unique())


# run main MAG-QC-pipeline.py script
if __name__ == "__main__":

    inputfile = os.path.realpath(sys.argv[1])
    outputfile = "MAG-QC-output_query-flattened-annotation-data.tsv"
    outputfile1 = "MAG-QC-output_annotation-count-data"
    inputfile1 = "MAG-QC-output_annotation-count-data"


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
    #lineage = "p__Euryarchaeota"


    check_for_py3()

    read_and_parse_rast_annotations(inputfile, outputfile)

    combine_external_checkm_and_taxonomy_info(file1, file2, file3)

    import_and_merge_tables(saveoutput)

    taxalist = extract_lineages_for_selected_level(taxa_level)[0]  # get list of lineages to iterate over


    for lineage_number in range(len(taxalist)):
        lineage = taxalist[lineage_number]
        #if (str(lineage) == "p__Nanoarchaeota") or (str(lineage) == "p__Micrarchaeota") or (str(lineage) == "p__Euryarchaeota") or (str(lineage) == "p__Asgardarchaeota"):
        #if (str(lineage) == "p__Asgardarchaeota"):

        print("Starting with lineage: "+str(lineage))

        subset_data_by_lineage(lineage)

        count_annotation_data_for_level(merge_reduced, outputfile1, lineage)

        import_count_and_combine_with_genome_metadata(lineage)

        if genome_number > 3:

            run_tsne_dimensional_reduction()

            run_principal_component_analysis()

            plot_dimensional_reduction_results_seaborn(lineage, mode="pca")

            plot_dimensional_reduction_results_rggplot(lineage, mode="pca")
        else:
            print("\nWarning: not enough data (genomes) to run a meaningful dimensional reduction. Select a different group.")


        print("Finished with lineage: "+str(lineage))


    # test_text()
