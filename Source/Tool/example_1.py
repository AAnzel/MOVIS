import streamlit as st
import pandas as pd
import common

import visualize


def get_data_set(omic_name):

    spinner_text_map =\
        {'genomics_mags': 'Embedding MAGs into vectors...',
         'genomics_mags_temporal': 'Embedding MAGs into vectors...',
         'genomics_mags_temporal_PCA': 'Embedding MAGs into vectors...',
         'genomics_kegg_temporal': 'Creating KO dataset...',
         'genomics_mags_annotated_temporal': 'Creating annotation dataset...',
         'genomics_mags_top_10_annotated_temporal': 'Creating annotation\
                                                     dataset...',
         'proteomics': 'Creating additional data set...',
         'metabolomics': 'Creating additional data set...',
         'phy_che': 'Creating additional data set...'}

    with st.spinner(spinner_text_map[omic_name]):
        df = common.get_cached_dataframe(common.EX_1, omic_name)

    # Fixing dataframe columns
    df = common.fix_dataframe_columns(df)

    return df


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


def show_clustering_info(df, key_suffix):

    clustering_methods = st.multiselect(
        'Choose clustering method:', ['K-Means', 'OPTICS'],
        key='choose_clus' + key_suffix
        )

    elbow_vis_col, num_input_col = st.beta_columns([3, 1])

    tmp_df = common.get_number_of_clusters(df)
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
        labels_list.append((i, common.cluster_data(
            df, cluster_number, cluster_samples, i)))

    return labels_list


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

    chosen_charts = []
    for i in choose_data_set:
        if i == 'W2V embedded MAGs':
            df_1 = get_data_set('genomics_mags_temporal')
            show_calculated_data_set(df_1, 'Embedded MAGs')
            labels_list = show_clustering_info(df_1, 'Genomics_1')

            # Traversing pairs in list
            for i in labels_list:
                temporal_feature, feature_list = common.find_temporal_feature(
                    df_1)
                feature_list = i[0]
                df_1[i[0]] = i[1]
                chosen_charts += common.visualize_data_set(
                        df_1, temporal_feature, feature_list,
                        'Genomics_1_' + i[0])

        elif i == 'KEGG matrix':
            df_2 = get_data_set('genomics_kegg_temporal')
            show_calculated_data_set(df_2, 'KO matrix')
            labels_list = show_clustering_info(df_2, 'Genomics_2')

            # Traversing pairs in list
            for i in labels_list:
                temporal_feature, feature_list = common.find_temporal_feature(
                    df_2)
                feature_list = i[0]
                df_2[i[0]] = i[1]
                chosen_charts += common.visualize_data_set(
                        df_2, temporal_feature, feature_list,
                        'Genomics_2_' + i[0])

        else:
            tmp_df_3 = get_data_set('genomics_mags_annotated_temporal')
            df_3 = get_data_set('genomics_mags_top_10_annotated_temporal')
            show_calculated_data_set(tmp_df_3, 'Product annotations')

            temporal_feature, feature_list = common.find_temporal_feature(df_3)
            chosen_charts += common.visualize_data_set(
                df_3, temporal_feature, feature_list, 'Genomics_3')

    # I should put cluster charts here, however I have to run it first
    # because I have rendered images and not altair charts
    # st.altair_chart()

    return chosen_charts


def create_main_example_1_metabolomics():
    st.header('Metabolomics')

    # Here I show the head() of the data set and some summary() and info()
    df = get_data_set('metabolomics')
    show_data_set(df)

    temporal_feature, feature_list = common.find_temporal_feature(df)

    chosen_charts = common.visualize_data_set(
        df, temporal_feature, feature_list, 'Metabolomics')

    # Here I should implement multiple select where I provide user with
    # different choices for what kind of chart/computation the user wants

    return chosen_charts


def create_main_example_1_proteomics():
    st.header('Proteomics')

    # TODO: Implement the same data set creation as with the genomics data
    # Create a data set with w2v and then implement all data set visulizations
    # as with the genomics data
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
    df = get_data_set('proteomics')
    show_calculated_data_set(df, 'Protein properties')

    temporal_feature, feature_list = common.find_temporal_feature(df)

    chosen_charts = common.visualize_data_set(
        df, temporal_feature, feature_list, 'Proteomics')

    # Here I should implement multiple select where I provide user with
    # different choices for what kind of chart/computation the user wants

    return chosen_charts


def create_main_example_1_phy_che():
    st.header('Physico-chemical')

    # Here I show the head() of the data set and some summary() and info()
    df = get_data_set('phy_che')
    show_data_set(df)

    temporal_feature, feature_list = common.find_temporal_feature(df)

    chosen_charts = common.visualize_data_set(
        df, temporal_feature, feature_list, 'Physico-Chemical')

    # Here I should implement multiple select where I provide user with
    # different choices for what kind of chart/computation the user wants

    return chosen_charts


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

    charts = []  # An empty list to hold all pairs (visualizations, key)

    with st.beta_expander('Show/hide data sets and related info',
                          expanded=True):
        if num_of_columns >= 2:
            column_list = st.beta_columns(num_of_columns)
            curr_pos = 0

            for i in choose_omics:
                if i == 'Genomics':
                    with column_list[curr_pos]:
                        curr_pos += 1
                        charts += create_main_example_1_genomics()
                elif i == 'Metabolomics':
                    with column_list[curr_pos]:
                        curr_pos += 1
                        charts += create_main_example_1_metabolomics()
                elif i == 'Proteomics':
                    with column_list[curr_pos]:
                        curr_pos += 1
                        charts += create_main_example_1_proteomics()
                else:
                    with column_list[curr_pos]:
                        curr_pos += 1
                        charts += create_main_example_1_phy_che()

        else:
            for i in choose_omics:
                if i == 'Genomics':
                    charts += create_main_example_1_genomics()
                elif i == 'Metabolomics':
                    charts += create_main_example_1_metabolomics()
                elif i == 'Proteomics':
                    charts += create_main_example_1_proteomics()
                else:
                    charts += create_main_example_1_phy_che()

    with st.beta_expander('Show/hide visualizations', expanded=True):
        for i in charts:
            type_of_chart = type(i[0])

            with st.spinner('Visualizing...'):
                if 'altair' in str(type_of_chart):
                    st.altair_chart(i[0], use_container_width=True)
                else:
                    pass

    return None
