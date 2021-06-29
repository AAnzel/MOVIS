import os
import math
import random
import visualize
import pandas as pd
import numpy as np
import altair as alt
import datetime as dt
import streamlit as st
from sklearn.cluster import KMeans, OPTICS
from sklearn import preprocessing
from sklearn import metrics
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from gensim.models import Word2Vec
from scipy.spatial.distance import jaccard, pdist, squareform


# Defining paths for each and every omic

path_root_data = os.path.join(
    "..", "..", "Data", "Extracted", "First source", "Databases"
)

path_cached_save_root = 'cached'

path_all_fasta = os.path.join(path_root_data, "fasta_files", "AllBins")
path_genomics_78 = os.path.join(path_root_data, "fasta_files",
                                "rmags_filtered")
path_genomics_kegg = os.path.join(path_root_data, "Annotations", "KEGG")
path_genomics_bins = os.path.join(path_root_data, "Annotations", "Bins")

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
                   (i.endswith("fa") and i.startswith("D"))])
num_of_proteomics = len([i for i in os.listdir(path_proteomics_78) if
                         (i.endswith("faa") and i.startswith("D"))])

SEED = 42
END = num_of_mags
ALL_DAYS = 51
MAX_ROWS = 15000
MAX_RANGE = 100
EPOCHS = 10
NUM_OF_WORKERS = 8
START_DATE = dt.datetime.strptime("2011-03-21", "%Y-%m-%d")
EX_1 = 1
EX_2 = 2
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


# Everything below is used for proteomics data set exclusively
def import_proteomics(end=25, path_proteomics=path_proteomics_78):

    print("Importing proteomics data")

    # There are 78 FASTA files I have to traverse every FASTA file, and in each
    # file every protein sequence

    fasta_files = os.listdir(path_proteomics)
    fasta_files.sort()
    tmp_all = []

    # This was done so that I could work with first 100 FASTA files only.
    # Otherwise, I should just remove: i, and enumerate
    for i, fasta_file_name in enumerate(fasta_files):

        if i == end:
            break

        else:
            with open(os.path.join(path_proteomics, fasta_file_name), "r")\
                     as input_file:

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

                        tmp_sec_str =\
                            sequence_analysis.secondary_structure_fraction()
                        tmp_list += [tmp_sec_str[0], tmp_sec_str[1],
                                     tmp_sec_str[2]]
                        tmp_list.append(
                            sequence.count("K")
                            + sequence.count("R")
                            - sequence.count("D")
                            - sequence.count("E")
                        )  # Electricity

                        amino_acid_perc =\
                            sequence_analysis.get_amino_acids_percent()

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


# Everything below is used for genomics data set exclusively
# Function that splits each genome into k-mers thus creating even longer
# sentence (MAG) It returns tokenized genome i.e. [kmer, kmer,...]
def split_genome(genome, k=5):
    new_genome = []
    n = len(genome)

    if n - k <= 0:
        return genome
    else:
        for i in range(n - k):
            new_genome.append(genome[i: i + k])

        return new_genome


def vectorize_one_mag(one_mag, w2v_model):

    # We have to generate vectors for each word in one MAG and then create
    # vector representation of that MAG by averaging vectors of its words
    zero_vector = np.zeros(w2v_model.vector_size)
    word_vectors = []
    one_mag_vector = []

    for sentence in one_mag:
        for word in sentence:
            if word in w2v_model.wv:
                try:
                    word_vectors.append(w2v_model.wv[word])
                except KeyError:
                    print("Key Error")
                    continue

    if word_vectors:
        word_vectors = np.asarray(word_vectors)
        one_mag_vector = word_vectors.mean(axis=0)

    else:
        one_mag_vector = zero_vector

    return one_mag_vector


# Function that vectorizes a MAG (document) with a pretrained word2vec model.
# It returns vector representation of a given MAG Vectorization is done by
# averaging word (k-mer) vectors for the whole document (MAG)
def vectorize_mags(w2v_model, path_fasta=path_genomics_78, end=25):

    print("Vectorizing MAGs")

    fasta_files = os.listdir(path_fasta)
    list_of_mag_vectors = []

    # This was done so that I could work with first 'end' FASTA files only.
    # Otherwise, I should just remove: i, and enumerate
    for i, fasta_file_name in enumerate(fasta_files):

        if i == end:
            break

        else:
            with open(os.path.join(path_fasta, fasta_file_name), "r")\
                    as input_file:

                one_mag = []
                for fasta_string in SeqIO.parse(input_file, "fasta"):

                    # Get kmers of a genome and create a sentence (list of
                    # words)
                    temp_kmers = split_genome(str(fasta_string.seq))

                    # Create a document (list of sentences)
                    one_mag.append(temp_kmers)

                # Vectorizing MAGs one by one
                list_of_mag_vectors.append(vectorize_one_mag(one_mag,
                                                             w2v_model))

    print("Finished vectorizing")

    return list_of_mag_vectors


# If one wants to import MAGs to train word2vec model, one should use only end
# argument, so that first 'end' MAGs are used for training
def import_mags_and_build_model(end=25, path_fasta=path_genomics_78):

    print("Importing MAGs and building model")

    # There are 1364 MAGs enclosed in FASTA files of the first dataset I have
    # to traverse every FASTA file, and in each file every sequence

    fasta_files = os.listdir(path_fasta)
    fasta_files.sort()
    fasta_ids = []

    # This was done so that I could work with first 100 FASTA files only.
    # Otherwise, I should just remove: i, and enumerate
    for i, fasta_file_name in enumerate(fasta_files):

        if i == end:
            break

        else:
            with open(os.path.join(path_fasta, fasta_file_name), "r")\
                    as input_file:

                one_mag = []
                one_mag_ids = []
                for fasta_string in SeqIO.parse(input_file, "fasta"):

                    # Get kmers of a genome and create a sentence (list of
                    # words)
                    temp_kmers = split_genome(str(fasta_string.seq))

                    # Create a document (list of sentences)
                    one_mag.append(temp_kmers)
                    # Save FASTA ids for every MAG
                    one_mag_ids.append(str(fasta_string.id))

                # Save list of ids for one MAG in global list
                fasta_ids.append(one_mag_ids)

                # If we do not have a model, we build one
                if i == 0:
                    print("Building w2v model")
                    # We build our model on the first MAG
                    w2v_model = Word2Vec(
                        sentences=one_mag, vector_size=100,
                        workers=NUM_OF_WORKERS, seed=SEED)

                # Else we just expand its vocabulary
                else:
                    # Now we expand our vocabulary
                    w2v_model.build_vocab(one_mag, update=True)

    print("Finished building")

    return w2v_model, fasta_files


def train_model(w2v_model, epochs=EPOCHS, path_fasta=path_genomics_78, end=25):

    print("Starting model training")

    # There are 1364 MAGs enclosed in FASTA files I have to traverse every
    # FASTA file, and in each file every sequence
    fasta_files = os.listdir(path_fasta)
    fasta_files.sort()

    # This was done so that I could work with first 100 FASTA files only.
    # Otherwise, I should just remove: i, and enumerate
    for i, fasta_file_name in enumerate(fasta_files):

        if i == end:
            break

        else:
            with open(os.path.join(path_fasta, fasta_file_name), "r")\
                    as input_file:

                one_mag = []
                for fasta_string in SeqIO.parse(input_file, "fasta"):

                    # Get kmers of a genome and create a sentence (list of
                    # words)
                    temp_kmers = split_genome(str(fasta_string.seq))

                    # Create a document (list of sentences)
                    one_mag.append(temp_kmers)

                w2v_model.train(
                    one_mag, total_examples=w2v_model.corpus_count,
                    epochs=epochs)

    print("Model training finished")

    return w2v_model


def example_1_import_kegg_and_create_df(end=51, path_fasta=path_genomics_78,
                                        path_all_keggs=path_genomics_kegg):

    print("Importing KEGG data")

    # There are 51 files for each day, in which there are KEGG IDs for each
    # genome collected that day I have to traverse every KEGG file, and create
    # DataFrame for each and every one

    kegg_files = [i for i in os.listdir(path_all_keggs) if
                  (i.endswith("besthits") and i.startswith("D"))]
    kegg_files.sort()
    rmags_78_names = [os.path.splitext(i)[0] for i in os.listdir(path_fasta) if
                      (i.endswith("fa") and i.startswith("D"))]
    rmags_78_names.sort()

    kegg_data_list = []

    # This was done so that I could work with first 100 files only. Otherwise,
    # I should just remove: i, and enumerate
    for i, kegg_file_name in enumerate(kegg_files):

        if i == end:
            break

        else:
            # Now I create a DataFrame out of it and save it in the list of
            # DataFrames
            tmp_df = pd.read_csv(os.path.join(path_genomics_kegg,
                                              kegg_file_name), delimiter="\t")
            tmp_filter = (tmp_df["Gene"].apply(lambda x: str(x).split("_")[0]
                                               + "_" + str(x).split("_")[1])
                                        .isin(rmags_78_names))

            tmp_df = tmp_df[tmp_filter]

            tmp_df["Gene"] = tmp_df["Gene"].apply(
                lambda x: str(x).split("_")[0] + "_" + str(x).split("_")[1]
            )
            tmp_df["ID"] = tmp_df["ID"].apply(lambda x: str(x).split(":")[1])
            tmp_df.drop(["maxScore", "hitNumber"], axis=1, inplace=True)
            tmp_df.reset_index(drop=True, inplace=True)

            kegg_data_list.append(tmp_df)

    print("Finished importing")
    return example_1_create_kegg_matrix(kegg_data_list, path_fasta)


def example_1_create_kegg_matrix(list_data, path_fasta=path_genomics_78):

    print("Creating KEGG matrix")

    rmags_78_names = [
        os.path.splitext(i)[0]
        for i in os.listdir(path_fasta)
        if (i.endswith("fa") and i.startswith("D"))
    ]
    rmags_78_names.sort()

    result_matrix_df = pd.DataFrame(columns=rmags_78_names)

    for i in list_data:
        tmp_df = i.value_counts().reset_index()

        for i, row in tmp_df.iterrows():
            result_matrix_df.at[row["ID"], row["Gene"]] = row[0]

    result_matrix_df.fillna(0, inplace=True)
    result_matrix_df = result_matrix_df.transpose()

    print("Finished creating")
    return result_matrix_df.sort_index()


def get_number_of_clusters(data):

    if 'DateTime' in data.columns.to_list():
        data = data.drop('DateTime', axis=1)

    mag_scaler = preprocessing.StandardScaler()
    scaled_data = mag_scaler.fit_transform(data)

    num_of_entries = data.shape[0]  # Getting number of rows
    k_range_end = int(math.sqrt(num_of_entries))  # Usually it is sqrt(#)
    k_range = range(1, k_range_end)

    k_mean_models = [KMeans(n_clusters=i, random_state=SEED) for i in k_range]
    k_scores = [
        k_mean_model.fit(scaled_data).score(scaled_data)
        for k_mean_model in k_mean_models
    ]
    k_model_data = pd.DataFrame({"k_range": k_range, "k_scores": k_scores})

    return k_model_data


def cluster_data(data, num_of_clusters, min_samples, model_name):

    if 'DateTime' in data.columns.to_list():
        data = data.drop('DateTime', axis=1)

    if model_name == 'K-Means':
        model = KMeans(n_clusters=num_of_clusters, random_state=SEED)

    elif model_name == 'OPTICS':
        model = OPTICS(min_samples=min_samples, n_jobs=NUM_OF_WORKERS)

    else:
        pass

    predicted_values = model.fit_predict(data)

    return predicted_values


def evaluate_clustering(data, predicted):
    return metrics.silhouette_score(data, predicted)


def create_pairwise_jaccard(data):

    tmp_data = data.clip(0, 1)
    result = squareform(pdist(tmp_data.astype(bool), jaccard))

    return pd.DataFrame(result, index=data.index, columns=data.index)


def example_1_create_annotated_data_set():

    rmags_names = [os.path.splitext(i)[0] for i in os.listdir(path_genomics_78)
                   if (i.endswith("fa") and i.startswith("D"))]
    rmags_names.sort()
    relevant_files = [i for i in os.listdir(path_genomics_bins)
                      if i.startswith(tuple(rmags_names))]
    relevant_files.sort()

    # Creating nested dictionary where each rmag has a dict of products and
    # their number of occurence for that rmag
    final_dict = {}
    rmags_names.sort()
    for i in rmags_names:
        final_dict[i] = {}
    # Traversing every annotation file, line by line, and saving only 'product'
    # column
    for annotation_file in relevant_files:
        with open(os.path.join(path_genomics_bins, annotation_file), 'r')\
                as input_file:

            rmag_name = os.path.splitext(annotation_file)[0]
            for line in input_file:
                product = line.split('product=')[-1].rstrip()
                if product not in final_dict[rmag_name]:
                    final_dict[rmag_name][product] = 0
                else:
                    final_dict[rmag_name][product] += 1

    result_df = pd.DataFrame.from_dict(final_dict).fillna(0).transpose()

    # I also save only first 10 products in regard to the number of occurence
    # This is done for easier visualization
    sorted_columns = result_df.sum(axis=0).sort_values(ascending=False)
    sorted_columns = sorted_columns.index.tolist()
    result_df = result_df[sorted_columns]
    top_10_result_df = result_df.iloc[:, :10]
    top_10_result_df['Other'] = result_df.iloc[:, 10:].sum(axis=1)

    return result_df, top_10_result_df


def cache_dataframe(dataframe, num_of_example, name):
    if num_of_example == EX_1:
        dataframe.to_pickle(os.path.join(path_cached_save_root,
                                         "example_1", "data_frames", name +
                                         "_dataframe.pkl"))
    elif num_of_example == EX_2:
        dataframe.to_pickle(os.path.join(path_cached_save_root,
                                         "example_2", "data_frames", name +
                                         "_dataframe.pkl"))
    else:
        dataframe.to_pickle(name)

    return None


@st.cache(allow_output_mutation=True)
def get_cached_dataframe(num_of_example, name):
    if num_of_example == EX_1:
        return pd.read_pickle(os.path.join(path_cached_save_root,
                                           "example_1", "data_frames", name +
                                           "_dataframe.pkl")).convert_dtypes()
    elif num_of_example == EX_2:
        return pd.read_pickle(os.path.join(path_cached_save_root,
                                           "example_2", "data_frames", name +
                                           "_dataframe.pkl")).convert_dtypes()
    else:
        return pd.read_pickle(name).convert_dtypes()


def fix_dataframe_columns(dataframe):

    old_columns = dataframe.columns.to_list()
    old_columns = [str(i) for i in old_columns]
    new_columns_map = {}
    bad_symbols = ['[', ']', '.']

    for column in old_columns:
        if any(char in column for char in bad_symbols):
            new_column = column
            for i in bad_symbols:
                new_column = new_column.replace(i, '_')

        else:
            new_column = column

        new_columns_map[column] = new_column

    return dataframe.rename(columns=new_columns_map)


def create_temporal_column(list_of_days, start_date, end, day_or_week):

    # In this case we just have to extract file names
    if day_or_week == 'TIMESTAMP':
        return [dt.datetime.strptime(i.split('.')[0], "%Y-%m-%d")
                for i in list_of_days]

    else:
        list_of_dates = []
        list_of_days.sort()

        # This is specific to the metaomics data set I am using Creating list
        # of dates for every rMAG. IT IS SAMPLED WEEKLY ! ! ! ! !! !
        for i in list_of_days[:end]:

            # Taking the number after D or W in D03.fa, without .fa so 03
            tmp_number = int(i.split('.')[0][1:])
            if day_or_week == 'W':
                tmp_datetime = start_date + dt.timedelta(weeks=tmp_number)
            else:
                tmp_datetime = start_date + dt.timedelta(days=tmp_number)

            while tmp_datetime in list_of_dates:
                if day_or_week == 'W':
                    tmp_datetime = tmp_datetime + dt.timedelta(days=1)
                else:
                    tmp_datetime = tmp_datetime + dt.timedelta(hours=1)

            list_of_dates.append(tmp_datetime)

        return list_of_dates


# ---
# # GENOMIC ANALYSIS
# ---
# Important only if we have filename D01.fa as in example 1
# This function fixes those types of file names
def example_1_fix_archive_file_names(start_date, unpack_archive_path):

    # I will first remove every FASTA file that doesn't start with 'D'
    files_list_old = os.listdir(unpack_archive_path)
    files_list_new = [i for i in files_list_old if i.startswith('D')]
    files_list_remove = [i for i in files_list_old if i not in files_list_new]

    for i in files_list_remove:
        os.remove(os.path.join(unpack_archive_path, i))

    # Sorting files and fixing the name. Originally they start with D but
    # represent weeks, so I will replace D with W
    files_list_new.sort()
    files_list_pass = [i.replace('D', 'W') for i in files_list_new]
    files_list_pass = [i.split('_')[0] + '.fa' for i in files_list_pass]

    list_of_dates = create_temporal_column(
        files_list_pass, start_date, len(files_list_pass),
        files_list_pass[0][0])

    list_of_new_names = [i.strftime('%Y-%m-%d') for i in list_of_dates]

    for i in range(len(files_list_pass)):
        os.rename(
            os.path.join(unpack_archive_path, files_list_new[i]),
            os.path.join(unpack_archive_path, list_of_new_names[i] + '.fa'))

    return None


def example_1_calc_genomics():

    kegg_matrix_df = example_1_import_kegg_and_create_df(
        end=ALL_DAYS, path_fasta=path_genomics_78,
        path_all_keggs=path_genomics_kegg)

    fasta_names = [i for i in os.listdir(path_genomics_78) if
                   (i.endswith("fa") and i.startswith("D"))]
    fasta_names.sort()
    list_of_dates = create_temporal_column(fasta_names, START_DATE, END, 'W')
    temporal_kegg_matrix_df = kegg_matrix_df.copy()
    temporal_kegg_matrix_df.insert(0, 'DateTime', list_of_dates)

    cache_dataframe(temporal_kegg_matrix_df, EX_1, 'genomics_kegg_temporal')

    ##########

    # KEGG examination but with pairwise Jaccard distance matrix(as
    #   seen in paper)
    # kegg_pairwise = create_pairwise_jaccard(kegg_matrix_df)
    # kegg_mds_chart = visualize.visualize_with_mds(
    #    kegg_pairwise, START_DATE, END, path_genomics_78)

    #######

    final_model, fasta_names =\
        import_mags_and_build_model(
            end=END, path_fasta=path_genomics_78)

    # Train model. It tooks ~10 minutes for END = 25 amount of MAGs
    final_model = train_model(final_model, epochs=EPOCHS, end=END)

    final_model.save(
        os.path.join(path_model_save_root, "genomics_model_78.saved"))

    # Now I should vectorize documents with this model. For further use, I
    # could save this model's weights, and use it to vectorize all mags. That
    # would take a lot, but every MAG will have its vector representation
    # > This could be done by importing one MAG at a time, then tokenizing
    # > it(like before), then getting vector representations of that MAG's
    # > sentences(genomes) and then finding the vector representation of the
    # > whole MAG(document). If I do that for one MAG at a time, There is no
    # > need to worry about memory
    #################
    list_of_mag_vectors = vectorize_mags(
        final_model, path_fasta=path_genomics_78, end=END)
    mags_df = pd.DataFrame(list_of_mag_vectors)
    cache_dataframe(mags_df, EX_1, 'genomics_mags')

    mags_df = get_cached_dataframe(EX_1, 'genomics_mags')

    list_of_dates = create_temporal_column(fasta_names, START_DATE, END, 'W')
    temporal_mags_df = mags_df.copy()
    temporal_mags_df.insert(0, 'DateTime', list_of_dates)
    cache_dataframe(temporal_mags_df, EX_1, 'genomics_mags_temporal')

    #################
    annotated_mags_df, top_10_annotated_mags_df =\
        example_1_create_annotated_data_set()

    list_of_dates = create_temporal_column(fasta_names, START_DATE, END, 'W')
    temporal_annotated_mags_df = annotated_mags_df.copy()
    temporal_annotated_mags_df.insert(0, 'DateTime', list_of_dates)
    temporal_top_10_annotated_mags_df = top_10_annotated_mags_df.copy()
    temporal_top_10_annotated_mags_df.insert(0, 'DateTime', list_of_dates)

    cache_dataframe(temporal_annotated_mags_df, EX_1,
                    'genomics_mags_annotated_temporal')
    cache_dataframe(temporal_top_10_annotated_mags_df, EX_1,
                    'genomics_mags_top_10_annotated_temporal')

    return None


def example_1_cluster_elbow():
    # There is already a function cluster_data
    # This should only be a wraper

    return None


# ---
# # METABOLOMIC ANALYSIS
# ---
def example_1_calc_metabolomics():
    # ## Importing Metabolomic data
    metabolomics_file_name = os.path.join(
        path_normalised_metabolomics,
        os.listdir(path_normalised_metabolomics)[0])
    metabolomics_df = pd.read_csv(metabolomics_file_name, delimiter="\t")

    # ## Data preprocessing

    metabolomics_df["date"] = pd.to_datetime(metabolomics_df["date"])
    metabolomics_df.insert(0, "date", metabolomics_df.pop("date"))
    metabolomics_df.sort_values("date", inplace=True, ignore_index=True)

    cache_dataframe(metabolomics_df, EX_1, "metabolomics")

    # Changing metabolite name if it is unknown
    metabolomics_df.loc[
        metabolomics_df["known_type"].eq("unknown"), "Metabolite"
    ] = np.nan

    print("Dataset uniqueness:")
    print("\t1. Timestamps:", len(metabolomics_df["date"].unique()))
    print("\t2. Metabolites:", len(metabolomics_df["Metabolite"].unique()))
    print("\t3. Types:", len(metabolomics_df["type"].unique()))
    print("\t4. Known types:", len(metabolomics_df["known_type"].unique()))
    print("\t5. Ns:", len(metabolomics_df["N"].unique()))
    print("\t6. Type 2s:", len(metabolomics_df["type2"].unique()))
    print("\t7. Measurements:", len(metabolomics_df["measurement"].unique()))

    # Saving the name column and removing unnecessairy columns metabolite_names
    # = metabolomics_df['Metabolite'] metabolomics_df.drop(labels =
    # ['Metabolite', 'tp', 'KEGG.Compound.ID', 'Chebi.Name',
    # 'Chebi.Name_combined'], axis = 1, inplace = True)
    metabolomics_df.drop(
        labels=["tp", "KEGG.Compound.ID", "Chebi.Name", "Chebi.Name_combined"],
        axis=1,
        inplace=True,
    )

    # Dummy eencoding categorical data
    scaled_metabolomics_df = pd.get_dummies(metabolomics_df,
                                            columns=["type", "known_type", "N",
                                                     "type2", "measurement"])

    # Standardizing data
    metabolomics_scaler = preprocessing.StandardScaler()
    scaled_metabolomics_df[
        ["means", "medians", "sds", "se", "ci"]
    ] = metabolomics_scaler.fit_transform(
        metabolomics_df[["means", "medians", "sds", "se", "ci"]]
    )

    metabolomics_df.dropna(inplace=True)
    metabolomics_df.reset_index(drop=True, inplace=True)

    return None


# ---
# # PROTEOMIC ANALYSIS
# ---
def example_1_calc_proteomics():
    # ## Importing Proteomic data

    # I could create something similar to Fig. 5 of the original paper, where I
    # would calculate mean of different proteomic feature values for each rMAG
    # calculated by days So I would have a table: date | feature 1 | feature 2
    # | ... Where each feature is mean of all values for one day of each MAG in
    # that rMAG

    proteomics_data = import_proteomics(end=num_of_proteomics)

    # I have to add temporality to this data set, according to file names
    fasta_files = [i for i in os.listdir(path_proteomics_78)
                   if (i.endswith("faa"))]
    fasta_files.sort()
    list_of_dates = create_temporal_column(fasta_files, START_DATE, END, 'W')
    proteomics_data.insert(0, 'DateTime', list_of_dates)

    # I will save this dataframe to show to the end-user
    cache_dataframe(proteomics_data, EX_1, "proteomics")

    return None


# ---
# # PHYSICO-CHEMICAL ANALYSIS
# ---
def example_1_calc_phy_che():
    # ## Importing Physico-chemical data
    phy_che_file_name = os.path.join(
        path_physico_chemical,
        [
            i
            for i in os.listdir(path_physico_chemical)
            if (i.endswith((".tsv", ".csv")))
        ][1],
    )
    phy_che_df = pd.read_csv(phy_che_file_name, decimal=",")

    # ## Data preprocessing
    phy_che_df.drop(index=0, axis=1, inplace=True)
    phy_che_df["Date"] = pd.to_datetime(phy_che_df["Date"])
    phy_che_df["Time"] = pd.to_timedelta(phy_che_df["Time"], unit="h")

    filtered_phy_che_df = phy_che_df[
                                     (phy_che_df["Date"] >= "2011-03-21") &
                                     (phy_che_df["Date"] <= "2012-05-03")]
    tmp_column = pd.Series(filtered_phy_che_df["Date"] +
                           filtered_phy_che_df["Time"])

    filtered_phy_che_df.drop(["Date", "Time"], axis=1, inplace=True)
    filtered_phy_che_df.reset_index(inplace=True, drop=True)
    filtered_phy_che_df = filtered_phy_che_df.apply(
        lambda x: pd.to_numeric(x.astype(str).str.replace(",", "."))
    )  # , errors='coerce'))
    filtered_phy_che_df.insert(0, "DateTime", tmp_column.values)

    # I will save this dataframe to show to the end-user
    cache_dataframe(filtered_phy_che_df, EX_1, "phy_che")

    return None


########
def show_data_set(df):
    with st.spinner('Showing the data set and related info'):
        st.markdown('First 100 entries')
        st.dataframe(df.head(100))
        tmp_df = df.describe()
        if len(tmp_df.columns.to_list()) > 1:
            st.dataframe(df.describe())

    return None


# This function is used when working with proteomics data i.e. FASTA files
# It is used to show calculated features of those FASTA files
def show_calculated_data_set(df, text_info):
    with st.spinner('Calculating features and showing the data set'):
        if len(df.columns.to_list()) > 50 or len(df.columns.to_list()) == 1:
            st.markdown('First 50 entries and first 8 features (columns). '
                        + '**' + text_info + '**')
            st.dataframe(df.iloc[:50, :8])
            # st.dataframe(df.describe())
        else:
            st.markdown('First 100 entries ' + '**' + text_info + '**')
            st.dataframe(df.head(100))
            st.dataframe(df.describe())

    return None


def show_clustering_info(df, key_suffix):

    clustering_methods = st.multiselect(
        'Choose clustering method:', ['K-Means', 'OPTICS'],
        key='choose_clus' + key_suffix
        )

    elbow_vis_col, num_input_col = st.beta_columns([3, 1])

    tmp_df = get_number_of_clusters(df)
    elbow_vis_col.altair_chart(visualize.elbow_rule(tmp_df),
                               use_container_width=True)

    help_text = '''Choose the number according to the elbow rule. The number of
                   clusters should be the number on the x-axis of the Elbow
                   chart where is the "elbow".'''

    if all(i in clustering_methods for i in ['K-Means', 'OPTICS']):
        cluster_number = num_input_col.slider(
            'Select a number of clusters for K-Means using the elbow rule:',
            min_value=1, max_value=15, value=1, step=1, format='%d',
            key='slider_cluster_Kmeans_' + key_suffix, help=help_text)
        cluster_samples = num_input_col.slider(
            'Select a minimum number of samples for OPTICS to be considered as\
            a core point:', min_value=1, max_value=15, value=1, step=1,
            format='%d', key='slider_cluster_Optics_' + key_suffix,
            help=help_text)

    elif 'K-Means' in clustering_methods:
        cluster_number = num_input_col.slider(
            'Select a number of clusters for K-Means using the elbow rule:',
            min_value=1, max_value=15, value=1, step=1, format='%d',
            key='slider_cluster_Kmeans_' + key_suffix, help=help_text)
        cluster_samples = 0

    elif 'OPTICS' in clustering_methods:
        cluster_number = 0
        cluster_samples = num_input_col.slider(
            'Select a minimum number of samples for OPTICS to be considered as\
            a core point:', min_value=1, max_value=15, value=1, step=1,
            format='%d', key='slider_cluster_Optics_' + key_suffix,
            help=help_text)

    else:
        pass

    # We create new columns that hold labels for each chosen method
    # it holds pairs (name of method, labels)
    labels_list = []
    for i in clustering_methods:
        labels_list.append((i, cluster_data(
            df, cluster_number, cluster_samples, i)))

    return labels_list


def show_folder_structure(uploaded_folder_path):

    space = '    '
    tee = '├── '
    max_lines = 7

    uploaded_folder_name = os.path.basename(
        os.path.normpath(uploaded_folder_path))

    folder_structure_text = uploaded_folder_name + '\n' + space
    files_list = os.listdir(uploaded_folder_path)
    files_list.sort()

    for i in range(0, min(len(files_list), max_lines)):
        folder_structure_text += tee + files_list[i] + '\n' + space
    folder_structure_text += '...'

    with st.spinner('Showing folder structure'):
        st.code(folder_structure_text)

    return files_list[len(files_list)-1][0]


# This functions will hold different fixes for uploaded data sets
def fix_data_set(df):

    # FIX 1:
    # We want to change ',' to '.' for all columns exept datetime eventhough
    # this is important only for float columns
    columns_to_fix = df.select_dtypes(
        exclude=[np.datetime64, 'datetime', 'datetime64',
                 np.timedelta64, 'timedelta', 'timedelta64',
                 'category', 'datetimetz']).columns.to_list()

    for column in columns_to_fix:
        df[column] = df[column].apply(
            lambda x: str(x).replace(",", "."))

    # FIX 2:

    return df


def fix_archive_file_names(start_date, unpack_archive_path):

    # I will first remove every FASTA file that doesn't start with 'D'
    files_list_old = os.listdir(unpack_archive_path)
    files_list_new = [i for i in files_list_old if i.startswith('D')]
    files_list_remove = [i for i in files_list_old if i not in files_list_new]

    for i in files_list_remove:
        os.remove(os.path.join(unpack_archive_path, i))

    files_list_new.sort()

    list_of_dates = create_temporal_column(
        files_list_new, start_date, len(files_list_new),
        files_list_new[0][0])

    list_of_new_names = [i.strftime('%Y-%m-%d') for i in list_of_dates]

    for i in range(len(files_list_new)):
        os.rename(
            os.path.join(unpack_archive_path, files_list_new[i]),
            os.path.join(unpack_archive_path, list_of_new_names[i] + '.fa'))

    return None


def find_temporal_feature(df):
    feature_list = df.columns.astype(str).to_list()
    temporal_feature = None
    datetime_strings = ['date', 'time']
    temporal_columns = []

    try:
        for i in feature_list:
            # This means that i is datetime column
            if any(datetime in i.lower() for datetime in datetime_strings):
                temporal_columns.append(i)

        if len(temporal_columns) == 0:
            raise ValueError

        # This means that we have only one temporal column with date or date
        # and time combined
        elif len(temporal_columns) == 1:
            df[temporal_columns[0]] = pd.to_datetime(df[temporal_columns[0]])
            temporal_feature = temporal_columns[0]

            st.success('Detected 1 temporal column named **' + temporal_feature
                       + '**.')

        # This means that we have 2 temporal columns, one for date and the
        # other one for time
        elif len(temporal_columns) == 2:
            for i in temporal_columns:
                if i.lower() == 'date':
                    tmp_date = pd.to_datetime(df[i])
                    df[i] = tmp_date
                else:
                    tmp_time = pd.to_timedelta(pd.to_numeric(df[i]), unit='h')
                    df[i] = tmp_time

            tmp_combined = tmp_date + tmp_time
            df.insert(0, 'DateTime', tmp_combined.values,
                      allow_duplicates=True)
            temporal_feature = 'DateTime'
            df.drop(temporal_columns, axis=1, inplace=True)

            st.success('Detected 2 temporal columns. Successfully merged them\
                        into one called **DateTime**.')

        # We are not providing any functionality if there are >=3 columns
        else:
            raise ValueError

        feature_list = df.columns.to_list()
        feature_list.remove(temporal_feature)
        df[feature_list] = df[feature_list].apply(
            pd.to_numeric, errors='ignore')
        df = df.convert_dtypes()

        return temporal_feature, feature_list

    except ValueError:
        st.error('Datetime column not detected.')
        st.stop()
        return None, None


def visualize_data_set(df, temporal_feature, feature_list, key_suffix):

    chosen_charts = []

    if key_suffix.startswith('Cluster'):
        # TODO: Implement t-SNE reduction
        visualizations = st.multiselect('Choose your visualization',
                                        ['PCA visualization',
                                         'MDS visualization'],
                                        key='vis_data_' + key_suffix)

    else:
        visualizations = st.multiselect('Choose your visualization',
                                        ['Feature through time',
                                         'Two features scatter-plot',
                                         'Scatter-plot matrix',
                                         'Multiple features parallel chart',
                                         'Heatmap',
                                         'Top 10 count through time'],
                                        key='vis_data_' + key_suffix)

    for i in visualizations:
        # I have to check which clustering method was used and visualize it
        if i == 'PCA visualization':
            chosen_charts.append(
                (visualize.visualize_clusters(df, temporal_feature,
                                              feature_list, 'PCA'),
                 i + '_' + key_suffix + '_PCA'))

        elif i == 'MDS visualization':
            chosen_charts.append(
                (visualize.visualize_clusters(df, temporal_feature,
                                              feature_list, 'MDS'),
                 i + '_' + key_suffix + '_MDS'))

        elif i == 'Feature through time' and temporal_feature is not None:
            selected_feature = st.selectbox('Select feature to visualize',
                                            feature_list)
            selected_color = None

            if str(df[selected_feature].dtype) not in ['string']:  # ,'Int64']:
                selected_color = st.color_picker(
                    'Select line color', value='#000000')

            chosen_charts.append(
                (visualize.time_feature(df, selected_feature, temporal_feature,
                                        selected_color), i + '_' + key_suffix))

        elif i == 'Two features scatter-plot':
            feature_1 = None
            feature_2 = None

            feature_1 = st.selectbox('Select 1. feature', feature_list)

            if feature_1 in feature_list:
                feature_list.remove(feature_1)

            feature_2 = st.selectbox('Select 2. feature', feature_list)

            if feature_1 is None or feature_2 is None:
                st.stop()

            chosen_charts.append(
                (visualize.two_features(df, feature_1, feature_2),
                 i + '_' + key_suffix))

        elif i == 'Scatter-plot matrix':
            target_feature = st.selectbox(
                'Select target feature for color', feature_list)

            list_of_features = st.multiselect('Choose at least 2 features',
                                              feature_list)

            if len(list_of_features) < 2:
                st.stop()

            list_of_features.append(target_feature)

            chosen_charts.append(
                    (visualize.scatter_matrix(df, list_of_features,
                                              target_feature),
                     i + '_' + key_suffix))

        elif i == 'Multiple features parallel chart':
            if temporal_feature is not None:
                target_feature = st.selectbox(
                    'Select target feature for color',
                    feature_list + [temporal_feature],
                    index=len(feature_list))
            else:
                target_feature = st.selectbox(
                    'Select target feature for color',
                    feature_list,
                    index=len(feature_list))

            list_of_features = st.multiselect('Choose at least 2 features',
                                              feature_list)

            if len(list_of_features) < 2:
                st.stop()

            list_of_features.append(target_feature)

            chosen_charts.append(
                    (visualize.parallel_coordinates(df, list_of_features,
                                                    target_feature),
                     i + '_' + key_suffix))
            # st.altair_chart(
            #    visualize.parallel_coordinates(df, list_of_features,
            #                                target_feature),
            #    use_container_width=True)
            #

        elif i == 'Heatmap':
            chosen_charts.append((visualize.heatmap(df), i + '_' + key_suffix))

        elif i == 'Top 10 count through time' and temporal_feature is not None:
            chosen_charts.append((visualize.top_10_time(df, feature_list,
                                                        temporal_feature),
                                  i + '_' + key_suffix))

        else:
            pass

    return chosen_charts
