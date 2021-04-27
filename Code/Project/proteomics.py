import os
import random
import math
import streamlit
import gensim
import altair_saver
import pandas as pd
import numpy as np
import altair as alt
import datetime as dt
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from gensim.models import Word2Vec
from sklearn.cluster import KMeans, OPTICS
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn import preprocessing, model_selection, metrics
from scipy.spatial.distance import jaccard, pdist, squareform


# Defining paths for each and every omic

path_root_data = os.path.join(
    "..", "..", "Data", "Extracted", "First source", "Databases"
)

path_all_fasta = os.path.join(path_root_data, "fasta_files", "AllBins")
path_genomics_78 = os.path.join(path_root_data, "fasta_files",
                                "rmags_filtered")
path_genomics_kegg = os.path.join(path_root_data, "Annotations", "KEGG")
path_normalised_metabolomics = os.path.join(
    path_root_data, "Metabolomics", "Normalised_Tables"
)
path_proteomics_78 = os.path.join(path_root_data, "Proteomics", "set_of_78")
path_physico_chemical = os.path.join(path_root_data, "PhysicoChemical")
path_second_source = os.path.join("..", "..", "Data", "Extracted",
                                  "Second source")

path_model_save_root = os.path.join("..", "Saved_models")
path_figures_save_root = os.path.join("..", "Output_figures")

num_of_mags = len([i for i in os.listdir(path_genomics_78) if
                   i.endswith("fa")])
num_of_proteomics = len([i for i in os.listdir(path_proteomics_78) if 
                         i.endswith("faa")])

SEED = 42
END = num_of_mags
ALL_DAYS = 51
MAX_ROWS = 15000
EPOCHS = 10
NUM_OF_WORKERS = 8
START_DATE = dt.datetime.strptime("2011-03-21", "%Y-%m-%d")
random.seed(SEED)
np.random.seed(SEED)
alt.data_transformers.enable(
    "default", max_rows=MAX_ROWS
)  # Important if you want to visualize datasets with >5000 samples


# Functions below are shared among different omics Function that saves charts
# from list_of_charts with names from list_of_names
def save_charts(list_of_chart, list_of_names):

    for chart, name in zip(list_of_chart, list_of_names):
        print(chart, name)
        # altair_saver.save(chart, os.path.join (path_figures_save_root, name))
        chart.save(os.path.join(path_figures_save_root, name))


# This function creates new dataframe with column that represent season
# according to date It also concatenates important types with metabolite names
def season_data(data, temporal_column):
    new_df = data
    new_df["season"] = new_df[temporal_column].dt.month % 12 // 3 + 1

    # important_types = [metabolite_column] + important_types
    # new_df['new_name'] = df[important_types].agg('\n'.join, axis=1)

    return new_df


def create_temporal_column(list_of_days, start_date, end):

    list_of_dates = []

    # This is specific to the metaomics data set I am using Creating list of
    # dates for every rMAG
    for i in list_of_days[:end]:

        tmp_datetime = start_date + dt.timedelta(weeks=int(i[1:3]))

        if tmp_datetime not in list_of_dates:
            list_of_dates.append(tmp_datetime)

        else:
            tmp_datetime = tmp_datetime.replace(day=tmp_datetime.day + 1)
            list_of_dates.append(tmp_datetime)

    return list_of_dates


# This functions are used for PROTEOMIC


def import_proteomics(end=25, path_proteomics=path_proteomics_78):

    print("Importing proteomics data")

    # There are 78 FASTA files I have to traverse every FASTA file, and in each
    # file every protein sequence

    fasta_files = [i for i in os.listdir(path_proteomics) if (i[-3:] == "faa")]
    tmp_all = []

    # This was done so that I could work with first 100 FASTA files only.
    # Otherwise, I should just remove: i, and enumerate
    for i, fasta_file_name in enumerate(fasta_files):

        if i == end:
            break

        else:
            with open(
                os.path.join(path_proteomics, fasta_file_name), "r"
            ) as input_file:

                one_mag_list = []
                for fasta_string in SeqIO.parse(input_file, "fasta"):

                    # Analyzing protein (peptide) and creating list of values
                    # for one MAG
                    sequence = str(fasta_string.seq)

                    if "*" in sequence:
                        continue

                    else:

                        sequence_analysis = ProteinAnalysis(sequence)

                        tmp_list = [
                            sequence_analysis.molecular_weight(),
                            sequence_analysis.gravy(),
                            sequence_analysis.aromaticity(),
                            sequence_analysis.instability_index(),
                            sequence_analysis.isoelectric_point(),
                        ]

                        tmp_sec_str = sequence_analysis.secondary_structure_fraction()
                        tmp_list += [tmp_sec_str[0], tmp_sec_str[1], tmp_sec_str[2]]
                        tmp_list.append(
                            sequence.count("K")
                            + sequence.count("R")
                            - sequence.count("D")
                            - sequence.count("E")
                        )  # Electricity

                        amino_acid_perc = sequence_analysis.get_amino_acids_percent()

                        tmp_list.append(sum([amino_acid_perc[aa] for aa in
                                             "AGILPV"]))
                        tmp_list.append(sum([amino_acid_perc[aa] for aa in
                                             "STNQ"]))
                        tmp_list.append(sum([amino_acid_perc[aa] for aa in
                                             "QNHSTYCMW"]))
                        tmp_list.append(sum([amino_acid_perc[aa] for aa in
                                             "AGILPVF"]))
                        tmp_list.append(sum([amino_acid_perc[aa] for aa in
                                             "HKR"]))
                        tmp_list.append(sum([amino_acid_perc[aa] for aa in
                                             "CM"]))
                        tmp_list.append(sum([amino_acid_perc[aa] for aa in
                                             "DE"]))
                        tmp_list.append(sum([amino_acid_perc[aa] for aa in
                                             "NQ"]))
                        tmp_list.append(sum([amino_acid_perc[aa] for aa in
                                             "ST"]))

                        # Now I put all these values in one_mag_list as a numpy
                        # arrays
                        one_mag_list.append(np.asarray(tmp_list))

                # Now I put one mag values, aggregated by mean, into the all
                # mag list
                tmp_all.append(np.asarray(one_mag_list).mean(axis=0))

    COLUMN_LIST = [
        "Molecular weight",
        "Gravy",
        "Aromaticity",
        "Instability index",
        "Isoelectric point",
        "Secondary structure fraction 0",
        "Secondary structure fraction 1",
        "Secondary structure fraction 2",
        "Electricity",
        "Fraction aliphatic",
        "Fraction uncharged polar",
        "Fraction polar",
        "Fraction hydrophobic",
        "Fraction positive",
        "Fraction sulfur",
        "Fraction negative",
        "Fraction amide",
        "Fraction alcohol",
    ]
    all_mag_df = pd.DataFrame(tmp_all, columns=COLUMN_LIST)

    print("Finished importing")

    return all_mag_df


def visualize_proteomics(data):

    # Adding another column that replaces temporal data for now
    if "Index_tmp" not in data.columns:
        data.insert(0, "Index_tmp", data.index.values)

    # Create repeated chart
    chart = (alt.Chart(data).mark_area().encode(
             alt.X("Index_tmp", type="quantitative"),
             alt.Y(alt.repeat("row"), type="quantitative"),
             ).properties(width=1200).repeat(row=data.columns.values))
    # .resolve_scale(size = 'independent')#.interactive()

    return chart
