import os
import random
import base64
import json
import altair_saver
import pandas as pd
import numpy as np
import altair as alt
import datetime as dt
import streamlit as st
# from sklearn.decomposition import PCA
# from sklearn.manifold import MDS
from sklearn import preprocessing
# from sklearn import metrics

import calculate

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
path_cached_save_root = os.path.join("cached")

num_of_mags = len([i for i in os.listdir(path_genomics_78) if
                  i.endswith("fa")])
num_of_proteomics = len(
    [i for i in os.listdir(path_proteomics_78) if i.endswith("faa")]
)
SEED = 42
END = num_of_mags
ALL_DAYS = 51
MAX_ROWS = 15000
EPOCHS = 10
NUM_OF_WORKERS = 8
START_DATE = dt.datetime.strptime("2011-03-21", "%Y-%m-%d")
EX_1 = 1
EX_2 = 2
random.seed(SEED)
np.random.seed(SEED)
alt.data_transformers.enable("default", max_rows=MAX_ROWS)
# Important if you want to visualize datasets with >5000 samples


# Functions below are shared among different omics Function that saves charts
# from list_of_charts with names from list_of_names
def save_charts(list_of_chart, list_of_names):

    for chart, name in zip(list_of_chart, list_of_names):
        altair_saver.save(chart, os.path.join(path_figures_save_root, name))
        # chart.save(os.path.join(path_figures_save_root, name))

    return None


def cache_dataframe(dataframe, num_of_example, name):
    if num_of_example == 1:
        dataframe.to_pickle(os.path.join(path_cached_save_root,
                                         "example_1", "data_frames", name +
                                         "_dataframe.pkl"))
    else:
        dataframe.to_pickle(os.path.join(path_cached_save_root,
                                         "example_2", "data_frames", name +
                                         "_dataframe.pkl"))
    return None


@st.cache
def get_cached_dataframe(num_of_example, name):
    if num_of_example == EX_1:
        return pd.read_pickle(os.path.join(path_cached_save_root,
                                           "example_1", "data_frames", name +
                                           "_dataframe.pkl")).convert_dtypes()
    else:
        return pd.read_pickle(os.path.join(path_cached_save_root,
                                           "example_2", "data_frames", name +
                                           "_dataframe.pkl")).convert_dtypes()


def fix_dataframe_columns(dataframe):

    old_columns = dataframe.columns.to_list()
    old_columns = [str(i) for i in old_columns]
    new_columns_map = {}
    bad_symbols = ['[', ']']

    for column in old_columns:
        if any(char in column for char in bad_symbols):
            new_column = column
            for i in bad_symbols:
                new_column = new_column.replace(i, '(')

        else:
            new_column = column

        new_columns_map[column] = new_column

    return dataframe.rename(columns=new_columns_map)


def get_hashed_name(list_of_strings):
    return (base64.urlsafe_b64encode(
        ('_'.join(list_of_strings)).encode("utf-8"))).decode("utf-8")


def cache_chart(chart, num_of_example, name):
    if num_of_example == 1:
        chart.save(os.path.join(path_cached_save_root, "example_1", "charts",
                                name + "_chart.json"))
    else:
        chart.save(os.path.join(path_cached_save_root, "example_2", "charts",
                                name + "_chart.json"))
    return None


def get_cached_chart(num_of_example, name):
    if num_of_example == 1:
        tmp_path = os.path.join(path_cached_save_root, "example_1", "charts",
                                name + "_chart.json")

    else:
        tmp_path = os.path.join(path_cached_save_root, "example_2", "charts",
                                name + "_chart.json")

    # If we cached the chart, we return it. If not, we return None as a signal
    if os.path.exists(tmp_path):
        with open(tmp_path) as f:
            tmp_json_data = json.load(f)
            
        print(tmp_json_data)
        return alt.Chart.from_dict(tmp_json_data)
    else:
        return None


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
    list_of_days.sort()

    # This is specific to the metaomics data set I am using Creating list of
    # dates for every rMAG. IT IS SAMPLED WEEKLY ! ! ! ! !! !
    for i in list_of_days[:end]:

        tmp_datetime = start_date + dt.timedelta(weeks=int(i[1:3]))

        if tmp_datetime not in list_of_dates:
            list_of_dates.append(tmp_datetime)

        else:
            tmp_datetime = tmp_datetime.replace(day=tmp_datetime.day + 1)
            list_of_dates.append(tmp_datetime)

    return list_of_dates


# ---
# # GENOMIC ANALYSIS
# ---
def example_1_calc_genomics():
    # **Important**: I should review the way I look at MAGs. The names of all
    # fasta files beggining with 'D_##' represent the days those MAGs were
    # obtained. Therefore, I should look at this also as timeseries data. Also,
    # maybe I should only consider 78 MAGs, and not all ~1300. After some
    # consideration, I conclude that I should definetly use only 78 MAGs,
    # because that way I wouldn't be tied to meta-omics data only. I also
    # thinked about what should I visualize in that case. One idea is that I
    # should also encode those MAGs with word2wec, and then make a 3D chart
    # where one dimension is time, and other two dimensions would be PCA
    # dimensions of those MAGs. I could also use this function to visualize
    # proteomics data if I want.

    # Another important thing is that I should actually include FASTA headers
    # and possibly use only them. That way, I could make figures like in a
    # relevant paper where MAGs are groupped according to their taxonomy etc. I
    # should look more into this.
    # ## MAG examination
    # ### KEGG examination

    #######################################################################
    # CREATE DIFFERENT FUNCTIONS FOR KEGG WORK AND FASTA WORK
    # CACHE KEGG MATRIX, FASTA MODEL AND MAGS_DF
    # PUT VISUALIZATION PARTS IN visualize.py
    #######################################################################
    '''
    kegg_matrix_df = calculate.import_kegg_and_create_df(
        end=ALL_DAYS, path_fasta=path_genomics_78,
        path_all_keggs=path_genomics_kegg)

    fasta_names = [i for i in os.listdir(path_genomics_78) if
                   (i.endswith("fa") and i.startswith("D"))]
    list_of_dates = create_temporal_column(fasta_names, START_DATE, END)
    temporal_kegg_matrix_df = kegg_matrix_df.copy()
    temporal_kegg_matrix_df.insert(0, 'DateTime', list_of_dates)

    cache_dataframe(temporal_kegg_matrix_df, EX_1, 'genomics_kegg_temporal')
    '''
    # ### KEGG examination but with pairwise Jaccard distance matrix(as
    # seen in paper)
    # kegg_pairwise = calculate.create_pairwise_jaccard(kegg_matrix)
    # kegg_mds_chart = visualize.visualize_with_mds(kegg_pairwise, START_DATE,
    #                                             END, path_genomics_78)
    # ---
    # # VAZNO:
    # Sledece sto treba da se uradi je da se nadje transcriptomic data set i da
    # se obradi i on u potpunosti. Nakon toga, treba da se sve podeli po
    # skriptama i da se odluci o dizajnu. Posle ostaje jos da se napravi front
    # end.
    #
    # ---

    # FOR CLUSTERING I SHOULD CREATE A DATAFRAME WITH MAGs INDEXES AND THEIR
    # VECTOR REPRESENTATIONS
    '''
    final_model, fasta_names, fasta_ids =\
        calculate.import_mags_and_build_model(end=END,
                                              path_fasta=path_genomics_78)

    # Train model. It tooks ~10 minutes for END = 25 amount of MAGs
    final_model = calculate.train_model(final_model, epochs=EPOCHS, end=END)

    final_model.wv.save_word2vec_format(
        os.path.join(path_model_save_root, "genomics_model_78.bin"),
        binary=True)
    '''
    # Now I should vectorize documents with this model. For further use, I
    # could save this model's weights, and use it to vectorize all mags. That
    # would take a lot, but every MAG will have its vector representation
    # > This could be done by importing one MAG at a time, then tokenizing
    # > it(like before), then getting vector representations of that MAG's
    # > sentences(genomes) and then finding the vector representation of the
    # > whole MAG(document). If I do that for one MAG at a time, There is no
    # > need to worry about memory
    #
    '''
    list_of_mag_vectors = calculate.vectorize_mags(final_model,
                                                   path_fasta=path_genomics_78,
                                                   end=END)
    mags_df = pd.DataFrame(list_of_mag_vectors)
    cache_dataframe(mags_df, EX_1, 'genomics_mags')

    mags_df = get_cached_dataframe(EX_1, 'genomics_mags')

    fasta_names = [i for i in os.listdir(path_genomics_78) if
                   (i.endswith("fa") and i.startswith("D"))]
    list_of_dates = create_temporal_column(fasta_names, START_DATE, END)
    temporal_mags_df = mags_df.copy()
    temporal_mags_df.insert(0, 'DateTime', list_of_dates)
    cache_dataframe(temporal_mags_df, EX_1, 'genomics_mags_temporal')

    # ## Data preprocessing
    mag_scaler = preprocessing.StandardScaler()
    scaled_mags_df = mag_scaler.fit_transform(mags_df)

    # PCA for visualizing MAGs
    pca_model = PCA(n_components=2, random_state=SEED)
    temporal_mags_df = pd.DataFrame(
        pca_model.fit_transform(scaled_mags_df), columns=['PCA_1', 'PCA_2'])

    temporal_mags_df.insert(0, 'DateTime', list_of_dates)

    cache_dataframe(temporal_mags_df, EX_1, 'genomics_mags_temporal_PCA')
    '''
    # THERE IS AN ERROR HERE, CHECK
    # MDS for visualizing KEGG
    # kegg_matrix_df = get_cached_dataframe(EX_1, 'genomics_kegg')

    # mds_model = MDS(n_components=2, random_state=SEED,
    #                dissimilarity="precomputed", n_jobs=NUM_OF_WORKERS)
    # kegg_matrix_transformed_df = pd.DataFrame(
    #    mds_model.fit_transform(kegg_matrix_df), columns=['MDS_1', 'MDS_2'])

    # kegg_matrix_transformed_df.insert(0, 'DateTime', list_of_dates)

    # cache_dataframe(kegg_matrix_transformed_df, EX_1,
    #                'genomics_kegg_temporal_MDS')
    '''
    annotated_mags_df, top_10_annotated_mags_df =\
            calculate.create_annotated_data_set()

    fasta_names = [i for i in os.listdir(path_genomics_78) if
                   (i.endswith("fa") and i.startswith("D"))]
    list_of_dates = create_temporal_column(fasta_names, START_DATE, END)
    temporal_annotated_mags_df = annotated_mags_df.copy()
    temporal_annotated_mags_df.insert(0, 'DateTime', list_of_dates)
    temporal_top_10_annotated_mags_df = top_10_annotated_mags_df.copy()
    temporal_top_10_annotated_mags_df.insert(0, 'DateTime', list_of_dates)

    cache_dataframe(temporal_annotated_mags_df, EX_1,
                    'genomics_mags_annotated_temporal')
    cache_dataframe(temporal_top_10_annotated_mags_df, EX_1,
                    'genomics_mags_top_10_annotated_temporal')
    '''
    return None


def example_1_cluster_genomics():
    # There is already a function calculate.cluster_data
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

    # ## Time series examination
    # metabolites_chart = visualize.visualize_metabolites(
    #   metabolomics_df, "date", "Metabolite", ["type", "type2", "measurement",
    #                                            "N"])

    # save_charts([metabolites_chart], ['metabolomics_metabolites_chart.png'])

    # ## Clustering
    # Deep learning temporal clustering Should I even do this? Previous
    # visualizations are descriptive enough. It would be a lot of work for not
    # much benefit
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

    proteomics_data = calculate.import_proteomics(end=num_of_proteomics)

    # I have to add temporality to this data set, according to file names
    fasta_files = [i for i in os.listdir(path_proteomics_78)
                   if (i.endswith("faa"))]
    list_of_dates = create_temporal_column(fasta_files, START_DATE, END)
    proteomics_data.insert(0, 'DateTime', list_of_dates)

    # I will save this dataframe to show to the end-user
    cache_dataframe(proteomics_data, EX_1, "proteomics")

    # chart_proteomics = visualize.visualize_proteomics(proteomics_data)
    # save_charts([chart_proteomics], ['proteomics_chart_proteomics.png'])
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

    # Visualize temperature, air_temperature, conductivity, inflow_pH, nitrate,
    # oxygen, pH
    # chart_phy_che = visualize.visualize_phy_che(
    #   filtered_phy_che_df, "DateTime", filtered_phy_che_df.columns.values[4:]
    # )
    # chart_phy_che_corr = visualize.visualize_phy_che_heatmap(
    #    filtered_phy_che_df)

    # save_charts([chart_phy_che_corr, chart_phy_che],
    # ['physico_chemical_chart_psy_che_corr.png',
    # 'physico_chemical_chart_psy_che.png'])
    return None
