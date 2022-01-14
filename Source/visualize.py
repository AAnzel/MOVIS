import random
import numpy as np
import pandas as pd
import altair as alt
from sklearn.decomposition import PCA
from sklearn.manifold import MDS, TSNE
from sklearn import preprocessing


__author__ = 'Aleksandar Anžel'
__copyright__ = ''
__credits__ = ['Aleksandar Anžel', 'Georges Hattab']
__license__ = 'GNU General Public License v3.0'
__version__ = '1.0'
__maintainer__ = 'Aleksandar Anžel'
__email__ = 'aleksandar.anzel@uni-marburg.de'
__status__ = 'Dev'


SEED = 42
MAX_ROWS = 15000
MAX_COLUMNS = 100
EPOCHS = 10
NUM_OF_WORKERS = 8
random.seed(SEED)
np.random.seed(SEED)
alt.data_transformers.enable(
    "default", max_rows=MAX_ROWS
)  # Important if you want to visualize datasets with >5000 samples


# This function creates new dataframe with column that represent season
# according to date It also concatenates important types with metabolite names
def season_data(data, temporal_column):
    new_data = data
    new_data["season"] = new_data[temporal_column].dt.month % 12 // 3 + 1

    # important_types = [metabolite_column]
    #                   + important_types new_data['new_name']
    # = data[important_types].agg('\n'.join, axis=1)

    return new_data


# This function is used to check whether the data set is too big to visualize
# If that is the case, it is resized to a preditermined value
def shrink_data(data, num_of_columns=MAX_COLUMNS):
    shrink_signal = False

    if len(data.columns) > MAX_COLUMNS:
        data.drop(data.columns.tolist()[MAX_COLUMNS:], axis=1, inplace=True)
        shrink_signal = True

    # if num_of_columns == 10:
    #     # Returning the TOP 10 features based on counts
    #     # TODO: There is a bug here. We have to include only numerical
    #     #  columns, without temporal and other column types
    #     sorted_columns = data.sum(axis=0).sort_values(ascending=False)
    #     sorted_columns = sorted_columns.index.tolist()
    #     data = data[sorted_columns]

    #     other_series = data.iloc[:, 10:].sum(axis=1)

    #     data.drop(sorted_columns[10:], axis=1, inplace=True)
    #     data['Other'] = other_series

    if num_of_columns != MAX_COLUMNS:
        data.drop(data.columns.tolist()[num_of_columns:], axis=1, inplace=True)
        shrink_signal = True
    else:
        pass

    return data, shrink_signal


# TODO: Implement time sampling (yearly, monthly, daily)
def time_feature(data, selected_column, temporal_column, target_feature):

    selected_column_type = str(data[selected_column].dtype)

    if selected_column_type == 'string':  # or selected_column_type == 'Int64':
        chart = alt.Chart(
            data, title=selected_column + ' through time').mark_bar().encode(
                alt.X(temporal_column, type='temporal',
                      scale=alt.Scale(nice=True)),
                alt.Y('count(' + selected_column + ')', type='nominal'),
                alt.Color(selected_column, type='nominal'),
                alt.Tooltip('count(' + selected_column + ')', type='nominal'))

    else:
        if target_feature is not None:
            chart = alt.Chart(
                data,
                title=selected_column + ' through time').mark_line().encode(
                        alt.X(temporal_column, type='temporal',
                              scale=alt.Scale(nice=True)),
                        alt.Y(selected_column, type='quantitative'),
                        alt.Color(target_feature, type='nominal'),
                        alt.Tooltip([temporal_column, selected_column,
                                     target_feature]))
        else:
            # TODO: https://altair-viz.github.io/gallery/select_mark_area.html
            # brush = alt.selection_interval(encodings=['x'])
            # color=alt.condition(brush, alt.value('#4c78a8'),
            #                                alt.value('lightgray')),
            #            opacity=alt.condition(brush, alt.value(1.0),
            #                                  alt.value(0.2))
            #            ).add_selection(brush)

            # Creating a selection that chooses the nearest point and selects
            # it based on x-value
            nearest = alt.selection(
                type='single', nearest=True, on='mouseover',
                fields=[temporal_column], empty='none')

            # The line part of the chart
            line = alt.Chart().mark_line(interpolate='basis').encode(
                alt.X(temporal_column, type='temporal',
                      scale=alt.Scale(nice=True)),
                alt.Y(selected_column, type='quantitative'),
            )

            # Transparent selectors across the chart. This is what tells us
            # the x-value of the cursor
            selectors = alt.Chart().mark_point().encode(
                alt.X(temporal_column, type='temporal'),
                opacity=alt.value(0),
            ).add_selection(
                nearest
            )

            # Draw points on the line, and highlight based on selection
            points = line.mark_point().encode(
                opacity=alt.condition(nearest, alt.value(1), alt.value(0))
            )

            # Draw text labels near the points, and highlight based on
            # selection
            text_1 = line.mark_text(
                align='left', dx=10, dy=15, stroke='#222222',
                fill='#ffffff', filled=False, size=12).encode(
                    text=alt.condition(nearest, selected_column + ':Q',
                                       alt.value(' '))
            )

            text_2 = line.mark_text(
                align='left', dx=10, dy=35, stroke='#222222',
                fill='#ffffff', filled=False, size=12).encode(
                    text=alt.condition(nearest, temporal_column + ':T',
                                       alt.value(' '))
            )

            # Draw a rule at the location of the selection
            rules = alt.Chart().mark_rule(color='gray').encode(
                alt.X(temporal_column, type='temporal'),
            ).transform_filter(
                nearest
            )

            chart = alt.layer(
                line, selectors, points, rules, text_1, text_2, data=data,
                title=selected_column + ' through time')

    return chart


def two_features(data, feature_1, feature_2, temporal_feature):

    title_text = feature_1 + ' in function of ' + feature_2
    brush = alt.selection_interval()

    if (str(data[feature_1].dtype) == 'string' and
            str(data[feature_2].dtype) != 'string'):
        chart = alt.Chart(data, title=title_text).mark_bar().encode(
            alt.X(feature_1, type='nominal', scale=alt.Scale(nice=True)),
            alt.Y(feature_2, type='quantitative'),
            alt.Tooltip([feature_1, feature_2, temporal_feature]),
            color=alt.condition(
                brush, alt.value('#4c78a8'), alt.value('lightgray')),
            opacity=alt.condition(brush, alt.value(1.0), alt.value(0.2))
        ).add_selection(brush)

    elif (str(data[feature_1].dtype) != 'string' and
            str(data[feature_2].dtype) == 'string'):
        chart = alt.Chart(data, title=title_text).mark_bar().encode(
            alt.X(feature_2, type='nominal'),
            alt.Y(feature_1, type='quantitative', scale=alt.Scale(nice=True)),
            alt.Tooltip([feature_1, feature_2, temporal_feature]),
            color=alt.condition(
                brush, alt.value('#4c78a8'), alt.value('lightgray')),
            opacity=alt.condition(brush, alt.value(1.0), alt.value(0.2))
        ).add_selection(brush)

    elif (str(data[feature_1].dtype) == 'string' and
            str(data[feature_2].dtype) == 'string'):
        chart = alt.Chart(data, title=title_text).mark_circle().encode(
            alt.X(feature_1, type='nominal'),
            alt.Y(feature_2, type='nominal'),
            alt.Tooltip([feature_1, feature_2, temporal_feature]),
            color=alt.condition(
                brush, alt.value('#4c78a8'), alt.value('lightgray')),
            opacity=alt.condition(brush, alt.value(1.0), alt.value(0.2))
        ).add_selection(brush)

    else:
        chart = alt.Chart(data, title=title_text).mark_circle().encode(
            alt.X(feature_1, type='quantitative'),
            alt.Y(feature_2, type='quantitative'),
            alt.Tooltip([feature_1, feature_2, temporal_feature]),
            color=alt.condition(
                brush, alt.value('#4c78a8'), alt.value('lightgray')),
            opacity=alt.condition(brush, alt.value(1.0), alt.value(0.2))
        ).add_selection(brush)

    return chart


def parallel_coordinates(data, list_of_features, target_feature):

    # TODO: Implement normalization before creating a chart
    # Use https://altair-viz.github.io/gallery/normed_parallel_coordinates.html
    selected_column_type = str(data[target_feature].dtype)

    color_type = 'nominal' if selected_column_type == 'string'\
        else 'quantitative'

    new_data = data[list_of_features].reset_index().melt(
        id_vars=['index', target_feature])

    # chart = alt.Chart(new_data).transform_window(
    #     key='count()'
    # ).transform_fold(
    #     list_of_features
    # ).transform_joinaggregate(
    #     min='min(value)',
    #     max='max(value)',
    #     groupby=['variable']
    # ).transform_calculate(
    #     minmax_value=(alt.datum.value - alt.datum.min)/(alt.datum.max -
    #                                                     alt.datum.min),
    #     mid=(alt.datum.min + alt.datum.max)/2
    # ).mark_line().encode(
    #     alt.X('variable:N', axis=alt.Axis(labelAngle=0)),
    #     alt.Y('minmax_value:Q'),
    #     alt.Color(target_feature, type=color_type),
    #     alt.Detail('key:N'),
    #     opacity=alt.value(0.4)
    # )

    chart = alt.Chart(
        new_data,
        title='Parallel coordinates chart of selected features').mark_line(
        ).encode(
        alt.X('variable:N', axis=alt.Axis(labelAngle=0)),
        alt.Y('value:Q'),
        alt.Color(target_feature, type=color_type),
        alt.Detail('index:N'),
        alt.Tooltip(['value', target_feature]),
        opacity=alt.value(0.4)
    )

    return chart.interactive()


def scatter(data, selected_feature, temporal_feature):

    selected_column_type = str(data[selected_feature].dtype)

    selected_column_type = 'nominal' if selected_column_type == 'string'\
        else 'quantitative'

    chart = alt.Chart(data).mark_circle().encode(
        alt.X(temporal_feature, type='temporal'),
        alt.Y(selected_feature, type=selected_column_type)
    ).interactive()

    return chart


def scatter_matrix(data, list_of_features, target_feature, temporal_feature):

    list_of_features.remove(target_feature)
    selected_column_type = str(data[target_feature].dtype)

    color_type = 'N' if selected_column_type == 'string' else 'Q'

    brush = alt.selection_interval()

    chart = alt.Chart(data).mark_circle().encode(
        alt.X(alt.repeat("column"), type='quantitative'),
        alt.Y(alt.repeat("row"), type='quantitative'),
        alt.Tooltip(list_of_features + [target_feature, temporal_feature]),
        color=alt.condition(
            brush, target_feature + ':' + color_type, alt.value('lightgray')),
        opacity=alt.condition(brush, alt.value(1.0), alt.value(0.2))
    ).properties(
        width=150,
        height=150
    ).repeat(
        row=list_of_features,
        column=list_of_features,
        title='Scatter matrix chart of selected features'
    ).add_selection(brush)

    return chart


def correlation_heatmap(data):

    new_data = data.copy()
    new_data = data.select_dtypes(include=np.number)
    new_data, shrink_signal = shrink_data(data, num_of_columns=MAX_COLUMNS)
    corr = new_data.corr().reset_index().melt("index")
    corr.columns = ["var_1", "var_2", "Correlation"]

    if shrink_signal:
        title_text = 'Heatmap chart of numerical features (first '\
                     + str(MAX_COLUMNS) + ' features only)'
    else:
        title_text = 'Heatmap chart of numerical features'

    # Create correlation chart
    chart = alt.Chart(
        corr, title=title_text).mark_rect().encode(
            alt.X("var_1", title=None, axis=alt.Axis(labelAngle=-45)),
            alt.Y("var_2", title=None),
            alt.Color("Correlation", legend=alt.Legend(tickCount=5),
                      scale=alt.Scale(scheme="redblue", reverse=True)),
            alt.Tooltip(['var_1', 'var_2', 'Correlation'])
        )

    # chart += chart.mark_text(size=8).encode(
    #     alt.Text("correlation", format=".2f"),
    #     color=alt.condition("abs(datum.correlation) > 0.5",
    #                         alt.value("white"), alt.value("black"))
    # )

    # This returns only lower triangle
    return chart.transform_filter("datum.var_1 < datum.var_2").interactive()


def time_heatmap(data, feature_1, feature_2, color_feature, temporal_feature):
    # TODO: Check if timestamps span over months, days or minutes
    # And use yearmonth, monthdate etc. acordingly
    # https://altair-viz.github.io/user_guide/transform/timeunit.html
    if color_feature == 'Sequential Single-Hue':
        color_scheme = 'greys'
    else:
        color_scheme = 'redblue'

    if feature_2 == temporal_feature:
        chart = alt.Chart(
            data,
            title=feature_1 + ' time heatmap').mark_rect().encode(
                alt.X('date(' + temporal_feature + ')', type='ordinal'),
                alt.Y('yearmonth(' + temporal_feature + ')', type='ordinal'),
                alt.Color(feature_1, type='quantitative',
                          scale=alt.Scale(scheme=color_scheme)),
                alt.Tooltip([temporal_feature, feature_1]))

    else:
        feature_2_type = str(data[feature_2].dtype)
        feature_2_type = 'nominal' if feature_2_type == 'string'\
            else 'quantitative'

        chart = alt.Chart(
            data,
            title=feature_1 + ' time heatmap').mark_rect().encode(
                alt.X('yearmonthdate(' + temporal_feature + ')',
                      type='ordinal', axis=alt.Axis(labelAngle=-45)),
                alt.Y(feature_2, type=feature_2_type),
                alt.Color(feature_1, type='quantitative',
                          scale=alt.Scale(scheme=color_scheme, reverse=True)),
                alt.Tooltip([temporal_feature, feature_1, feature_2]))

    return chart.interactive()


def whisker(summary_data, temporal_column):
    # Summary data frame contains the following columns:
    # 'Q1','Q3', 'IQR', 'LowerLimit', 'UpperLimit', 'Mean'

    bar_chart = alt.Chart(summary_data).mark_bar(size=10).encode(
        alt.X(temporal_column, type='temporal'),
        alt.Y('Q1:Q', title=None),
        alt.Y2('Q3:Q', title=None),
        alt.Tooltip(['Q1:Q', 'Q3:Q', 'LowerLimit:Q', 'UpperLimit:Q', 'Mean:Q'])
    )

    whiskers_chart = alt.Chart(summary_data).mark_rule().encode(
        alt.X(temporal_column, type='temporal'),
        alt.Y('LowerLimit:Q', scale=alt.Scale(zero=False), title=None),
        alt.Y2('UpperLimit:Q', title=None),
        alt.Tooltip(['Q1:Q', 'Q3:Q', 'LowerLimit:Q', 'UpperLimit:Q', 'Mean:Q'])
    )

    mean_chart = alt.Chart(summary_data).mark_tick(
        color='black', height=80, opacity=1).encode(
        alt.X(temporal_column, type='temporal'),
        alt.Y('Mean:Q', title=None),
        alt.Tooltip(['Q1:Q', 'Q3:Q', 'LowerLimit:Q', 'UpperLimit:Q', 'Mean:Q'])
    )

    final_chart = (whiskers_chart + bar_chart + mean_chart).configure_scale(
        bandPaddingInner=0.2)

    return final_chart.interactive()


def top_10_time(data, list_of_features, temporal_column):

    # I want to create a stacked bar chart where on x axis I will have time
    # and on y axis I will have stacked precentages of a whole
    # Example: https://altair-viz.github.io/gallery/bar_rounded.html
    new_data = data.copy()
    new_data, shrink_signal = shrink_data(new_data, 10)
    new_data = new_data.reset_index().melt(id_vars=['index', temporal_column])

    brush = alt.selection(type='interval')

    chart = alt.Chart(
        new_data, title='Top 10 share through time').mark_bar(
        ).encode(
        alt.X(temporal_column, type='temporal', scale=alt.Scale(domain=brush)),
        alt.Y('value:Q'),  # , stack='normalize'),
        alt.Color('variable:N', scale=alt.Scale(scheme='category10')),
        # legend=alt.Legend(orient='top', direction='vertical')),
        tooltip=['value', temporal_column]
    )

    interval_chart = alt.Chart(new_data).mark_line().encode(
        alt.X(temporal_column, type='temporal'),
        alt.Y('sum(value):Q')
    ).add_selection(brush).properties(height=60)

    # IMPORTANT: There is a streamlit bug that prevents vconcatenated chart
    # to fill the full width of the screen area
    return alt.vconcat(interval_chart, chart)


def elbow_rule(data):

    chart = alt.Chart(data, title='Elbow rule chart').mark_line().encode(
            alt.X("k_range:Q"),
            alt.Y("k_scores:Q"),
            alt.Tooltip(['k_range', 'k_scores'])
        )

    return chart


# Everything below is used for genomics data set exclusively
def visualize_clusters(data, temporal_feature, feature_list, method):
    tmp_data = data.copy()
    temporal_series = tmp_data[temporal_feature]
    tmp_data = tmp_data.drop(temporal_feature, axis=1)

    labels_feature = feature_list[0]
    labels_series = tmp_data[labels_feature]
    tmp_data = tmp_data.drop(labels_feature, axis=1)

    scaler = preprocessing.StandardScaler()
    scaled_data = scaler.fit_transform(tmp_data)

    if method == 'PCA':
        tmp_title = '2 dimensional PCA scatter plot'
        pca_model = PCA(n_components=2, random_state=SEED)
        tmp_result = pca_model.fit_transform(scaled_data)

        axis_name_1 = 'PCA_1_' + str(np.round(
            100 * pca_model.explained_variance_ratio_[0],
            2)).replace('.', '_') + '%'
        axis_name_2 = 'PCA_2_' + str(np.round(
            100 * pca_model.explained_variance_ratio_[1],
            2)).replace('.', '_') + '%'

        tmp_data = pd.DataFrame(
            tmp_result,
            columns=[axis_name_1, axis_name_2])

    elif method == 'MDS':
        tmp_title = '2 dimensional MDS scatter plot'
        mds_model = MDS(n_components=2, random_state=SEED,
                        dissimilarity="euclidean", n_jobs=NUM_OF_WORKERS)
        tmp_data = pd.DataFrame(
            mds_model.fit_transform(scaled_data), columns=['MDS_1', 'MDS_2'])

    elif method == 't-SNE':
        tmp_title = '2 dimensional t-SNE scatter plot'
        t_sne_model = TSNE(n_components=2, random_state=SEED,
                           metric="euclidean", n_jobs=NUM_OF_WORKERS)
        tmp_data = pd.DataFrame(
            t_sne_model.fit_transform(scaled_data),
            columns=['t-SNE_1', 't-SNE_2'])

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

    chart = alt.Chart(tmp_data, title=tmp_title).mark_point(opacity=1).encode(
            alt.X(str(tmp_data.columns[2]), type='quantitative'),
            alt.Y(str(tmp_data.columns[3]), type='quantitative'),
            alt.Shape(str(tmp_data.columns[1]), type='nominal'),
            alt.Color(str(tmp_data.columns[0]), type='temporal',
                      scale=alt.Scale(scheme='greys')),
            alt.Tooltip([temporal_feature, labels_feature]),
        ).configure_legend(columns=2)
    # .add_selection(select_time).transform_filter(select_time)

    return chart


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
