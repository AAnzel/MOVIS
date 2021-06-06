import streamlit as st
import numpy as np
import pandas as pd
import omics_run

import visualize


def show_data_set(df):
    
    with st.spinner('Showing the data set and related info'):
        st.markdown('First 100 entries')
        st.dataframe(df.head(100))
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


def find_temporal_feature(df):
    feature_list = list(df.columns.values)

    temporal_feature = None
    temporal_feature = str(
        df.select_dtypes(include=[np.datetime64]).columns[0])

    if temporal_feature is None:
        st.text('Datetime column not detected.')
        # TODO: Implement choosing time column and modifying dataframe

    feature_list.remove(temporal_feature)

    return temporal_feature, feature_list


def choose_columns(df, key_prefix):

    new_columns = st.multiselect('Choose columns to visualize:',
                                 list(df.columns.values),
                                 key='choose_col' + key_prefix)

    return new_columns


def visualize_data_set(df, temporal_feature, feature_list, key_prefix):

    visualizations = st.multiselect('Choose your visualization',
                                    ['Feature through time',
                                     'Two features scatter-plot',
                                     'Scatter-plot matrix',
                                     'Multiple features parallel chart',
                                     'Heatmap',
                                     'Top 10 count through time'],
                                    key='vis_data_' + key_prefix)

    for i in visualizations:
        if i == 'Feature through time':
            selected_feature = st.selectbox('Select feature to visualize',
                                            feature_list)
            selected_color = st.color_picker('Select line color')

            with st.spinner('Visualizing...'):
                st.altair_chart(
                    visualize.time_feature(df, selected_feature,
                                           temporal_feature, selected_color),
                    use_container_width=True)

        elif i == 'Two features scatter-plot':
            feature_1 = st.selectbox('Select 1. feature', feature_list)

            if feature_1 in feature_list:
                feature_list.remove(feature_1)

            feature_2 = st.selectbox('Select 2. feature', feature_list)

            with st.spinner('Visualizing...'):
                st.altair_chart(
                    visualize.two_features(df, feature_1, feature_2),
                    use_container_width=True)

        elif i == 'Scatter-plot matrix':
            target_feature = st.selectbox('Select target feature for color',
                                          feature_list)

            list_of_features = st.multiselect('Choose features', feature_list)

            if len(list_of_features) < 2:
                st.stop()

            list_of_features.append(target_feature)

            with st.spinner('Visualizing...'):
                st.altair_chart(
                    visualize.scatter_matrix(df, list_of_features,
                                             target_feature),
                    use_container_width=True)

        elif i == 'Multiple features parallel chart':
            target_feature = st.selectbox(
                    'Select target feature for color',
                    feature_list + [temporal_feature],
                    index=len(feature_list))

            list_of_features = st.multiselect('Choose features', feature_list)

            if len(list_of_features) < 2:
                st.stop()

            list_of_features.append(target_feature)

            with st.spinner('Visualizing...'):
                '''
                st.altair_chart(
                    visualize.parallel_coordinates(df, list_of_features,
                                                   target_feature),
                    use_container_width=True)
                '''
                st.plotly_chart(
                    visualize.parallel_coordinates(df, list_of_features,
                                                   target_feature),
                    use_container_width=True)

        elif i == 'Heatmap':
            with st.spinner('Visualizing...'):
                st.altair_chart(
                    visualize.heatmap(df), use_container_width=True)

        elif i == 'Top 10 count through time':
            with st.spinner('Visualizing...'):
                st.altair_chart(
                    visualize.top_10_time(df, feature_list,
                                          temporal_feature),
                    use_container_width=True)

        else:
            pass

    return None


def create_main_example_1_genomics():
    st.header('Genomics')
    with st.spinner('Showing folder structure'):
        st.code('''rmags_filtered/
            ├── D03_O1.31.2.fa
            ├── D04_G2.13.fa
            ├── D04_G2.5.fa
            ├── D04_L6.fa
            ├── D04_O2.19.fa
            ├── D05_G3.14.1.fa
            ├── D05_G3.4.fa
            ├── D05_L3.17.fa
            ├── D08_G1.16.fa
            ├── D08_O6.fa
            ...
                ''')

    # Here I should implement multiple select where I provide user with
    # different choices for what kind of chart/computation the user wants
    data_set_list = ['W2V embedded MAGs', 'KEGG matrix',
                     'Product-annotated MAGs']
    choose_data_set = st.multiselect('Which data set do you want to see:',
                                     data_set_list)

    for i in choose_data_set:
        if i == 'W2V embedded MAGs':
            with st.spinner('Embedding MAGs into vectors...'):
                df_1 = omics_run.get_cached_dataframe(omics_run.EX_1,
                                                      'genomics_mags')
            show_calculated_data_set(df_1, 'Embedded MAGs')

        elif i == 'KEGG matrix':
            with st.spinner('Creating KO dataset...'):
                df_2 = omics_run.get_cached_dataframe(omics_run.EX_1,
                                                      'genomics_kegg_temporal')
            show_calculated_data_set(df_2, 'KO matrix')

        else:
            with st.spinner('Creating annotation dataset...'):
                df_3 = omics_run.get_cached_dataframe(
                    omics_run.EX_1, 'genomics_mags_annotated_temporal')
                df_4 = omics_run.get_cached_dataframe(
                    omics_run.EX_1, 'genomics_mags_top_10_annotated_temporal')
            show_calculated_data_set(df_3, 'Product annotations')

            temporal_feature, feature_list = find_temporal_feature(df_4)
            visualize_data_set(df_4, temporal_feature, feature_list,
                               'genomics')

    #################################################
    # TODO:
    # I SHOULD ESCAPE ALL BRACKETS IN COLUMN NAMES BECAUSE ALTAIR CANNOT
    # PROCESS THEM AND GIVES BLANK CHART
    # I SHOULD DO THIS FOR ALL USER UPLOADED DATA SETS AS WELL
    ###################################################

    # I should put cluster charts here, however I have to run it first
    # because I have rendered images and not altair charts
    # st.altair_chart()

    return None


def create_main_example_1_metabolomics():
    st.header('Metabolomics')

    # Here I show the head() of the data set and some summary() and info()
    with st.spinner('Getting data set...'):
        df = omics_run.get_cached_dataframe(omics_run.EX_1, 'metabolomics')
    show_data_set(df)

    temporal_feature, feature_list = find_temporal_feature(df)

    visualize_data_set(df, temporal_feature, feature_list, 'metabolomics')

    # Here I should implement multiple select where I provide user with
    # different choices for what kind of chart/computation the user wants

    return None


def create_main_example_1_proteomics():
    st.header('Proteomics')

    with st.spinner('Showing folder structure'):
        st.code('''set_of_78/
            ├── D03_O1.31.2.faa
            ├── D04_G2.13.faa
            ├── D04_G2.5.faa
            ├── D04_L6.faa
            ├── D04_O2.19.faa
            ├── D05_G3.14.1.faa
            ├── D05_G3.4.faa
            ├── D05_L3.17.faa
            ├── D08_G1.16.faa
            ├── D08_O6.faa
            ...
                ''')

    # Here I show the head() of the data set and some summary() and info()
    with st.spinner('Creating additional data set...'):
        df = omics_run.get_cached_dataframe(omics_run.EX_1, 'proteomics')
    show_calculated_data_set(df, 'Protein properties')

    temporal_feature, feature_list = find_temporal_feature(df)

    visualize_data_set(df, temporal_feature, feature_list, 'proteomics')

    # Here I should implement multiple select where I provide user with
    # different choices for what kind of chart/computation the user wants

    return None


def create_main_example_1_phy_che():
    st.header('Physico-chemical')

    # Here I show the head() of the data set and some summary() and info()
    df = omics_run.get_cached_dataframe(omics_run.EX_1, 'phy_che')
    show_data_set(df)

    temporal_feature, feature_list = find_temporal_feature(df)

    visualize_data_set(df, temporal_feature, feature_list, 'phy_che')

    # Here I should implement multiple select where I provide user with
    # different choices for what kind of chart/computation the user wants

    return None


def create_main_example_1():

    col_1, col_2 = st.beta_columns(2)

    col_1.info('''
            This data set comes from the following paper:

            **Herold, M., Martínez Arbas, S., Narayanasamy, S. et al.\
            Integration of time-series meta-omics data reveals how microbial\
            ecosystems respond to disturbance. Nat Commun 11, 5281(2020).\
            https://doi.org/10.1038/s41467-020-19006-2**. Analyzed samples\
            were collected from a biological wastewater treatment plant in\
            Schifflange, Luxembourg. A precise location is shown on the map\
            located on the right.

            It contains **genomics**, **metabolomics**, **proteomics**, and\
            **physico-chemical** data. The code used to parse the data can be\
            found here: [GitLab]
            (https://git-r3lab.uni.lu/malte.herold/laots_niche_ecology_analysis)
            ''')

    col_2.map(pd.DataFrame({'lat': [49.513414], 'lon': [6.017925]}),
              zoom=13, use_container_width=True)

    example_1_omics_list = ['Genomics', 'Metabolomics', 'Proteomics',
                            'Physico-chemical']
    choose_omics = st.multiselect('Which omic do you want to see:',
                                  example_1_omics_list)

    num_of_columns = len(choose_omics)

    with st.beta_expander('Show/hide data sets and related info', expanded=True):
        if num_of_columns >= 2:
            column_list = st.beta_columns(num_of_columns)
            curr_pos = 0

            for i in choose_omics:
                if i == 'Genomics':
                    with column_list[curr_pos]:
                        curr_pos += 1
                        create_main_example_1_genomics()
                elif i == 'Metabolomics':
                    with column_list[curr_pos]:
                        curr_pos += 1
                        create_main_example_1_metabolomics()
                elif i == 'Proteomics':
                    with column_list[curr_pos]:
                        curr_pos += 1
                        create_main_example_1_proteomics()
                else:
                    with column_list[curr_pos]:
                        curr_pos += 1
                        create_main_example_1_phy_che()

        else:
            for i in choose_omics:
                if i == 'Genomics':
                    create_main_example_1_genomics()
                elif i == 'Metabolomics':
                    create_main_example_1_metabolomics()
                elif i == 'Proteomics':
                    create_main_example_1_proteomics()
                else:
                    create_main_example_1_phy_che()
    
    return None
