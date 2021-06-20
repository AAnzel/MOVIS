import os
import random
import numpy as np
import pandas as pd
import altair as alt
import datetime as dt
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn import preprocessing
# import plotly.express as px


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
    new_data = data
    new_data["season"] = new_data[temporal_column].dt.month % 12 // 3 + 1

    # important_types = [metabolite_column]
    #                   + important_types new_data['new_name']
    # = data[important_types].agg('\n'.join, axis=1)

    return new_data


def time_feature(data, selected_column, temporal_column, selected_color):

    selected_column_type = str(data[selected_column].dtype)
    if selected_column_type == 'string':  # or selected_column_type == 'Int64':
        chart = alt.Chart(data).mark_bar().encode(
            alt.X(temporal_column, type='temporal',
                  scale=alt.Scale(nice=True)),
            alt.Y('count(' + selected_column + ')', type='nominal'),
            alt.Color(selected_column, type='nominal'),
            alt.Tooltip('count(' + selected_column + ')', type='nominal'))

    else:
        chart = alt.Chart(data).mark_line(stroke=selected_color).encode(
            alt.X(temporal_column, type='temporal',
                  scale=alt.Scale(nice=True)),
            alt.Y(selected_column, type='quantitative'))

    return chart.interactive()


def two_features(data, feature_1, feature_2):

    if (str(data[feature_1].dtype) == 'string' and
            str(data[feature_2].dtype) != 'string'):
        chart = alt.Chart(data).mark_bar().encode(
            alt.X(feature_1, type='nominal', scale=alt.Scale(nice=True)),
            alt.Y(feature_2, type='quantitative'))

    elif (str(data[feature_1].dtype) != 'string' and
            str(data[feature_2].dtype) == 'string'):
        chart = alt.Chart(data).mark_bar().encode(
            alt.X(feature_2, type='nominal'),
            alt.Y(feature_1, type='quantitative', scale=alt.Scale(nice=True)))

    elif (str(data[feature_1].dtype) == 'string' and
            str(data[feature_2].dtype) == 'string'):
        chart = alt.Chart(data).mark_point().encode(
            alt.X(feature_1, type='nominal'),
            alt.Y(feature_2, type='nominal'))
    else:
        chart = alt.Chart(data).mark_point().encode(
            alt.X(feature_1, type='quantitative'),
            alt.Y(feature_2, type='quantitative'))

    return chart.interactive()


def parallel_coordinates(data, list_of_features, target_feature):

    # TODO: Implement normalization before creating a chart
    # See: https://github.com/AAnzel/DataVis_Supplementary_Material/

    new_data = data[list_of_features].reset_index().melt(
        id_vars=['index', target_feature])

    chart = alt.Chart(new_data).mark_line().encode(
        alt.X('variable:N'),
        alt.Y('value:Q'),
        alt.Color(target_feature, type='quantitative'),
        alt.Detail('index:N'),
        opacity=alt.value(0.4)
    )
    '''
    chart = px.parallel_coordinates(
        data, color=range(0, len(data[target_feature])),
        dimensions=list_of_features,
        color_continuous_scale=px.colors.sequential.Inferno,
        color_continuous_midpoint=2)
    '''

    return chart.interactive()


def scatter_matrix(data, list_of_features, target_feature):

    list_of_features.remove(target_feature)

    chart = alt.Chart(data).mark_circle().encode(
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


def heatmap(data):

    new_data = data.select_dtypes(include=np.number)
    corr = new_data.corr().reset_index().melt("index")
    corr.columns = ["var_1", "var_2", "correlation"]

    # Create correlation chart
    chart = alt.Chart(corr).mark_rect().encode(
            alt.X("var_1", title=None, axis=alt.Axis(labelAngle=-45)),
            alt.Y("var_2", title=None),
            alt.Color(
                "correlation",
                legend=alt.Legend(tickCount=5),
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


def top_10_time(data, list_of_features, temporal_column):

    # I want to create a stacked bar chart where on x axis I will have time
    # and on y axis I will have stacked precentages of a whole
    # Example: https://altair-viz.github.io/gallery/bar_rounded.html
    new_data = data.reset_index().melt(id_vars=['index', temporal_column])

    brush = alt.selection(type='interval')

    chart = alt.Chart(new_data).mark_bar().encode(
        alt.X(temporal_column, type='temporal', scale=alt.Scale(domain=brush)),
        alt.Y('value:Q'),  # , stack='normalize'),
        alt.Color('variable:N', scale=alt.Scale(scheme='category10')),
        # legend=alt.Legend(orient='top', direction='vertical')),
        tooltip=['value']
    )

    interval_chart = alt.Chart(new_data).mark_line().encode(
        alt.X(temporal_column, type='temporal'),
        alt.Y('sum(value):Q')
    ).add_selection(brush).properties(height=60)

    # IMPORTANT: There is a streamlit bug that prevents vconcatenated chart
    # to fill the full width of the screen area
    return alt.vconcat(interval_chart, chart)


def elbow_rule(data):

    chart = alt.Chart(
        data, title='Elbow rule chart', description='Description'
    ).mark_line().encode(
        alt.X("k_range:Q"), alt.Y("k_scores:Q"))

    return chart


# Everything below is used for metabolomics data set exclusively
def visualize_metabolites(data, temporal_column, metabolite_column,
                          type_columns):

    data_seasoned = season_data(data, temporal_column)

    # Extract columns with float values
    # NOT TESTED, MIGHT NOT WORK AS EXPECTED
    float_columns = data.select_dtypes(include=[np.float32,
                                                np.float64]).columns.tolist()

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


# Everything below is used for proteomics data set exclusively
def visualize_proteomics(data):

    # Adding another column that replaces temporal data for now
    if "Index_tmp" not in data.columns:
        data.insert(0, "Index_tmp", data.index.tolist())

    # Create repeated chart
    chart = (alt.Chart(data).mark_area().encode(
             alt.X("Index_tmp", type="quantitative"),
             alt.Y(alt.repeat("row"), type="quantitative"),
             ).properties(width=1200).repeat(row=data.columns.tolist()))
    # .resolve_scale(size = 'independent')#.interactive()

    return chart


# Everything below is used for genomics data set exclusively
def visualize_clusters(data, temporal_feature, labels_feature, method):

    temporal_series = data[temporal_feature]
    tmp_data = data.drop(temporal_feature, axis=1)

    labels_series = tmp_data[labels_feature]
    tmp_data = tmp_data.drop(labels_feature, axis=1)

    scaler = preprocessing.StandardScaler()
    scaled_data = scaler.fit_transform(tmp_data)

    if method == 'PCA':
        pca_model = PCA(n_components=2, random_state=SEED)
        tmp_data = pd.DataFrame(
            pca_model.fit_transform(scaled_data), columns=['PCA_1', 'PCA_2'])

    elif method == 'MDS':
        mds_model = MDS(n_components=2, random_state=SEED,
                        dissimilarity="euclidean", n_jobs=NUM_OF_WORKERS)
        tmp_data = pd.DataFrame(
            mds_model.fit_transform(scaled_data), columns=['MDS_1', 'MDS_2'])

    else:
        pass

    tmp_data.insert(0, temporal_feature, temporal_series)
    tmp_data.insert(1, labels_feature, labels_series)
    # tmp_data['New_DateTime'] =\
    #     tmp_data[temporal_feature].apply(lambda x: x.value)

    # time_start = tmp_data['New_DateTime'].min()
    # time_end = tmp_data['New_DateTime'].max()
    # slider = alt.binding_range(min=time_start, max=time_end, step=100)
    # select_time = alt.selection_single(
    #     fields=['New_DateTime'], bind=slider)

    chart = alt.Chart(tmp_data).mark_circle(opacity=1).encode(
            alt.X(str(tmp_data.columns[2]), type='quantitative'),
            alt.Y(str(tmp_data.columns[3]), type='quantitative'),
            alt.Color(str(tmp_data.columns[1]), type='nominal'),
            alt.Tooltip(str(tmp_data.columns[0]), type='temporal')
        )
    # .add_selection(select_time).transform_filter(select_time)

    return chart.interactive()


def visualize_seasonal_clusters(data, temporal_column):

    data_transformed = season_data(data, temporal_column)

    chart = alt.Chart(data_transformed).mark_circle(opacity=1).encode(
            alt.X(str(data.columns[1]), type='quantitative'),
            alt.X(str(data.columns[2]), type='quantitative'),
            alt.Color("season:N",
                      scale=alt.Scale(range=["blue", "green", "orange",
                                             "brown"])),
        )

    return chart
