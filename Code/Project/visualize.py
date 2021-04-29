import os
import random
import numpy as np
import pandas as pd
import altair as alt
import datetime as dt
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn.manifold import MDS


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


def time_feature(df, selected_column, temporal_column):

    if str(df[selected_column].dtype) == 'string':
        chart = alt.Chart(df).mark_bar().encode(
            alt.X(temporal_column, type='temporal',
                  scale=alt.Scale(nice=True)),
            alt.Y(selected_column, type='nominal'))
    else:
        chart = alt.Chart(df).mark_line().encode(
            alt.X(temporal_column, type='temporal',
                  scale=alt.Scale(nice=True)),
            alt.Y(selected_column, type='quantitative'))

    return chart.interactive()


def two_features(df, feature_1, feature_2):

    if (str(df[feature_1].dtype) == 'string' and
            str(df[feature_2].dtype) != 'string'):
        chart = alt.Chart(df).mark_bar().encode(
            alt.X(feature_1, type='nominal', scale=alt.Scale(nice=True)),
            alt.Y(feature_2, type='quantitative'))

    elif (str(df[feature_1].dtype) != 'string' and
            str(df[feature_2].dtype) == 'string'):
        chart = alt.Chart(df).mark_bar().encode(
            alt.X(feature_2, type='nominal'),
            alt.Y(feature_1, type='quantitative', scale=alt.Scale(nice=True)))

    elif (str(df[feature_1].dtype) == 'string' and
            str(df[feature_2].dtype) == 'string'):
        chart = alt.Chart(df).mark_point().encode(
            alt.X(feature_1, type='nominal'),
            alt.Y(feature_2, type='nominal'))
    else:
        chart = alt.Chart(df).mark_point().encode(
            alt.X(feature_1, type='quantitative'),
            alt.Y(feature_2, type='quantitative'))

    return chart.interactive()


def parallel_coordinates(df, list_of_features, target):
    '''
    new_df = df[list_of_features].reset_index().melt(id_vars=['index', target])

    chart = alt.Chart(new_df).mark_line().encode(
        alt.X('variable:N'),
        alt.Y('value:Q'),
        alt.Color(target, type='nominal'),
        alt.Detail('index:N'),
        opacity=alt.value(0.4)
    )
    '''
    chart = px.parallel_coordinates(
        df, color=range(0, len(df[target])), dimensions=list_of_features,
        color_continuous_scale=px.colors.diverging.Tealrose,
        color_continuous_midpoint=2)

    return chart


def scatter_matrix(df, list_of_features, target_feature):

    list_of_features.remove(target_feature)

    chart = alt.Chart(df).mark_circle().encode(
        alt.X(alt.repeat("column"), type='quantitative'),
        alt.Y(alt.repeat("row"), type='quantitative'),
        color=alt.Color(target_feature, type='quantitative')
    ).properties(
        width=150,
        height=150
    ).repeat(
        row=list_of_features,
        column=list_of_features
    )

    return chart


# Everything below is used for metabolomics data set exclusively
def visualize_metabolites(data, temporal_column, metabolite_column,
                          type_columns):

    data_seasoned = season_data(data, temporal_column)

    # Extract columns with float values
    float_columns = []

    for i in data_seasoned.columns:
        if (data_seasoned[i].dtypes == "float64"
                or data_seasoned[i].dtypes == "float32"):
            float_columns.append(i)

    # Create repeated chart with varying size encodings
    chart = (alt.Chart(data_seasoned).mark_point(opacity=1).encode(
            alt.X(temporal_column, type="temporal",
                  scale=alt.Scale(nice=True)),
            alt.Y(metabolite_column, type="nominal"),
            alt.Size(alt.repeat("row"), type="quantitative"),
            alt.Color("season:N",
                      scale=alt.Scale(range=["blue", "green", "orange",
                                             "brown"])),
            alt.Tooltip(type_columns, type="nominal"),
        ).properties(width=1200).repeat(row=float_columns)
        .resolve_scale(color="independent", size="independent")
    )
    # .interactive()

    return chart


# Everything below is used for phy_che data set exclusively
def visualize_phy_che(data, temporal_column, list_of_columns):

    # Create repeated chart
    chart = (
        alt.Chart(data).mark_line().encode(
            alt.X(temporal_column, type="temporal"),  # , timeUnit = 'month'),
            alt.Y(alt.repeat("row"), type="quantitative"),
        ).properties(width=1200).repeat(row=list_of_columns)
    )
    # .resolve_scale(color = 'independent')#.interactive()

    return chart


def visualize_phy_che_heatmap(data):

    new_data = data.drop("DateTime", axis=1)
    corr = new_data.corr().reset_index().melt("index")
    corr.columns = ["var_1", "var_2", "correlation"]

    # Create correlation chart
    chart = alt.Chart(corr).mark_rect().encode(
            alt.X("var_1", title=None, axis=alt.Axis(labelAngle=-45)),
            alt.Y("var_2", title=None),
            alt.Color(
                "correlation",
                legend=None,
                scale=alt.Scale(scheme="redblue", reverse=True),
            ),
        ).properties(width=alt.Step(30), height=alt.Step(30))

    chart += chart.mark_text(size=8).encode(
        alt.Text("correlation", format=".2f"),
        color=alt.condition("abs(datum.correlation) > 0.5",
                            alt.value("white"), alt.value("black"))
    )

    # This returns only lower triangle
    return chart.transform_filter("datum.var_1 < datum.var_2").interactive()


# Everything below is used for proteomics data set exclusively
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


# Everything below is used for genomics data set exclusively
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
