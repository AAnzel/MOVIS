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

    # important_types = [metabolite_column]
    #                   + important_types new_df['new_name']
    # = df[important_types].agg('\n'.join, axis=1)

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


# This functions are used for GENOMICS

# Function that splits each genome into k-mers thus creating even longer
# sentence (MAG) It returns tokenized genome i.e. [kmer, kmer,...]
def split_genome(genome, k=5):
    new_genome = []
    n = len(genome)

    if n - k <= 0:
        return genome
    else:
        for i in range(n - k):
            new_genome.append(genome[i : i + k])

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

    fasta_files = [i for i in os.listdir(path_fasta) if
                   (i.endswith("fa") and i.startswith("D"))]
    list_of_mag_vectors = []

    # This was done so that I could work with first 'end' FASTA files only.
    # Otherwise, I should just remove: i, and enumerate
    for i, fasta_file_name in enumerate(fasta_files):

        if i == end:
            break

        else:
            with open(os.path.join(path_fasta, fasta_file_name), "r") as input_file:

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

    fasta_files = [i for i in os.listdir(path_fasta) if
                   (i.endswith("fa") and i.startswith("D"))]
    fasta_ids = []

    # This was done so that I could work with first 100 FASTA files only.
    # Otherwise, I should just remove: i, and enumerate
    for i, fasta_file_name in enumerate(fasta_files):

        if i == end:
            break

        else:
            with open(os.path.join(path_fasta, fasta_file_name), "r") as input_file:

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
                        sentences=one_mag, size=100, workers=NUM_OF_WORKERS,
                        seed=SEED)

                # Else we just expand its vocabulary
                else:
                    # Now we expand our vocabulary
                    w2v_model.build_vocab(one_mag, update=True)

    print("Finished building")

    return w2v_model, fasta_files, fasta_ids


def train_model(w2v_model, epochs, path_fasta=path_genomics_78, end=25):

    print("Starting model training")

    # There are 1364 MAGs enclosed in FASTA files I have to traverse every
    # FASTA file, and in each file every sequence

    fasta_files = [i for i in os.listdir(path_fasta) if
                   (i.endswith("fa") and i.startswith("D"))]

    # This was done so that I could work with first 100 FASTA files only.
    # Otherwise, I should just remove: i, and enumerate
    for i, fasta_file_name in enumerate(fasta_files):

        if i == end:
            break

        else:
            with open(os.path.join(path_fasta, fasta_file_name), "r") as input_file:

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


def visualize_with_pca(data, labels, centers):

    pca_model = PCA(n_components=2, random_state=SEED)
    data_transformed = pca_model.fit_transform(data)

    data_transformed = pd.DataFrame(data_transformed)
    data_transformed.columns = ["PC_1", "PC_2"]
    data_transformed["Labels"] = labels

    chart_data = (alt.Chart(data_transformed).mark_circle(opacity=1).encode(
            alt.X("PC_1:Q"),
            alt.Y("PC_2:Q"),
            alt.Color("Labels:N", legend=alt.Legend())
        ))

    # This means we are visualising centroids from k_means (there are less
    # centroids that data points)
    if labels.shape[0] != centers.shape[0]:

        centers_transformed = pca_model.fit_transform(centers)
        centers_transformed = pd.DataFrame(centers_transformed)
        centers_transformed.columns = ["PC_1", "PC_2"]

        chart_centers = (alt.Chart(centers_transformed)
                         .mark_point(shape="diamond", color="black", size=50,
                                     opacity=0.7).encode(
                                                         alt.X("PC_1:Q"),
                                                         alt.Y("PC_2:Q"),
            ))

        return chart_data + chart_centers

    # For DBSCAN there are no centroids
    else:
        return chart_data


def visualize_temporal_mags(data, list_of_days, start_date, end):

    list_of_dates = create_temporal_column(list_of_days, start_date, end)

    pca_model = PCA(n_components=2, random_state=SEED)
    data_transformed = pca_model.fit_transform(data)

    data_transformed = np.hstack(
        ((np.asarray(list_of_dates))[:, np.newaxis], data_transformed)
    )
    data_transformed = pd.DataFrame(
        data_transformed, columns=["DateTime", "PCA_1", "PCA_2"]
    )

    data_transformed = season_data(data_transformed, "DateTime")

    chart_data = (
        alt.Chart(data_transformed).mark_circle(opacity=1).encode(
            alt.X("PCA_1:Q"),
            alt.Y("PCA_2:Q"),
            alt.Color("season:N",
                      scale=alt.Scale(range=["blue", "green", "orange",
                                             "brown"])),
        ).properties(width=1200)
    )

    return chart_data


def import_kegg_and_create_df(end=51, path_fasta=path_genomics_78,
                              path_all_keggs=path_genomics_kegg):

    print("Importing KEGG data")

    # There are 51 files for each day, in which there are KEGG IDs for each
    # genome collected that day I have to traverse every KEGG file, and create
    # DataFrame for each and every one

    kegg_files = [i for i in os.listdir(path_all_keggs) if
                  (i.endswith("besthits") and i.startswith("D"))]

    rmags_78_names = [os.path.splitext(i)[0] for i in os.listdir(path_fasta) if
                      (i.endswith("fa") and i.startswith("D"))]

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
    return create_kegg_matrix(kegg_data_list, path_fasta)


def create_kegg_matrix(list_data, path_fasta=path_genomics_78):

    print("Creating KEGG matrix")

    rmags_78_names = [
        os.path.splitext(i)[0]
        for i in os.listdir(path_fasta)
        if (i.endswith("fa") and i.startswith("D"))
    ]
    result_matrix_df = pd.DataFrame(columns=rmags_78_names)

    for i in list_data:
        tmp_df = i.value_counts().reset_index()

        for i, row in tmp_df.iterrows():
            result_matrix_df.at[row["ID"], row["Gene"]] = row[0]

    result_matrix_df.fillna(0, inplace=True)

    print("Finished creating")
    return result_matrix_df.T


def create_pairwise_jaccard(data):

    tmp_data = data.clip(0, 1)
    result = squareform(pdist(tmp_data.astype(bool), jaccard))

    return pd.DataFrame(result, index=data.index, columns=data.index)


def visualize_with_mds(data, start_date, end, path_fasta=path_genomics_78):

    mds_model = MDS(
        n_components=2,
        random_state=SEED,
        dissimilarity="precomputed",
        n_jobs=NUM_OF_WORKERS,
    )
    mds_pos = mds_model.fit_transform(data)

    list_of_days = [i for i in os.listdir(path_fasta) if
                    (i.endswith("fa") and i.startswith("D"))]

    temporal_column = create_temporal_column(list_of_days, start_date, end)

    data_transformed = pd.DataFrame(mds_pos)
    data_transformed.columns = ["MDS_1", "MDS_2"]
    data_transformed = np.hstack(
        ((np.asarray(temporal_column))[:, np.newaxis], data_transformed)
    )
    data_transformed = pd.DataFrame(
        data_transformed, columns=["DateTime", "MDS_1", "MDS_2"]
    )

    data_transformed = season_data(data_transformed, "DateTime")

    chart_data = (alt.Chart(data_transformed).mark_circle(opacity=1).encode(
            alt.X("MDS_1:Q"),
            alt.Y("MDS_2:Q"),
            alt.Color("season:N", scale=alt.Scale(range=["blue", "green",
                                                         "orange", "brown"])),
        ))

    return chart_data
