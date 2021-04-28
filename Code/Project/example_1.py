import os
import io
import random
import math
import gensim
import altair_saver
import streamlit as st
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

import omics_run
import metabolomics

URANDOM_LENGTH = 5


def show_data_set(df):
    with st.beta_expander('Show the data set and related info', expanded=True):
        st.markdown('First 100 entries')
        st.dataframe(df.head(100))
        st.dataframe(df.describe())

    return None


def choose_columns(df):

    new_columns = st.multiselect('Choose columns to visualize:',
                                 list(df.columns.values),
                                 key=os.urandom(URANDOM_LENGTH))
    return new_columns


def create_main_example_1_genomics():
    st.header('Genomics')
    with st.beta_expander('Show folder structure', expanded=True):
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
    feature_list = list(df.columns.values)

    temporal_feature = None
    temporal_feature = str(
        df.select_dtypes(include=[np.datetime64]).columns[0])

    if temporal_feature is None:
        st.text('Datetime column not detected.')
        # TODO: Implement choosing time column and modifying dataframe

    feature_list.remove(temporal_feature)


    visualizations = st.multiselect('Choose your visualization',
                                    ['Feature through time',
                                     'Two features scatter-plot',
                                     'Scatter-plot matrix',
                                     'Multiple features parallel chart'])

    for i in visualizations:
        if i == 'Feature through time':
            selected_feature = st.selectbox('Select feature to visualize',
                                            feature_list)

            with st.spinner('Visualizing...'):
                st.altair_chart(
                    metabolomics.visualize_time_feature(df, selected_feature,
                                                        temporal_feature),
                    use_container_width=True)

        elif i == 'Two features scatter-plot':
            col_1, col_2 = st.beta_columns(2)
            feature_1 = col_1.selectbox('Select 1. feature',
                                        feature_list)

            if feature_1 in feature_list:
                feature_list.remove(feature_1)

            feature_2 = col_2.selectbox('Select 2. feature',
                                        feature_list)

            with st.spinner('Visualizing...'):
                st.altair_chart(
                    metabolomics.visualize_two_features(df, feature_1,
                                                        feature_2),
                    use_container_width=True)

        elif i == 'Scatter-plot matrix':
            pass

        elif i == 'Multiple features parallel chart':
            target_feature = st.selectbox('Select target feature for color',
                                          feature_list)

            list_of_features = st.multiselect('Choose features', feature_list)

            if len(list_of_features) < 2:
                st.stop()

            list_of_features.append(target_feature)

            with st.spinner('Visualizing...'):
                '''
                st.altair_chart(
                    metabolomics.visualize_parallel(df, list_of_features,
                                                    target_feature),
                    use_container_width=True)
                '''
                st.plotly_chart(
                    metabolomics.visualize_parallel(df, list_of_features,
                                                    target_feature),
                    use_container_width=True)
                    
        else:
            pass

    # Here I should implement multiple select where I provide user with
    # different choices for what kind of chart/computation the user wants

    return None


def create_main_example_1_proteomics():
    st.header('Proteomics')

    # Here I show the head() of the data set and some summary() and info()

    # Here I should implement multiple select where I provide user with
    # different choices for what kind of chart/computation the user wants

    return None


def create_main_example_1_phy_che():
    st.header('Physico-chemical')

    # Here I show the head() of the data set and some summary() and info()
    df = omics_run.get_cached_dataframe(omics_run.EX_1, 'phy_che')
    show_data_set(df)
    
    with st.beta_expander('Show a correlation matrix', expanded=True):
        st.altair_chart(omics_run.phy_che.visualize_phy_che_heatmap(df),
                        use_container_width=True)

    # Here I should implement multiple select where I provide user with
    # different choices for what kind of chart/computation the user wants

    return None


def create_main_example_1():

    st.info('''
            This data set comes from the following paper:

            **Herold, M., Martínez Arbas, S., Narayanasamy, S. et al.\
            Integration of time-series meta-omics data reveals how microbial\
            ecosystems respond to disturbance. Nat Commun 11, 5281(2020).\
            https://doi.org/10.1038/s41467-020-19006-2**

            It contains **genomics**, **metabolomics**, **proteomics**, and\
            **physico-chemical** data. The code used to parse the data can be\
            found here: [GitHub](put_link)
            ''')

    example_1_omics_dict = {'Genomics': 0, 'Metabolomics': 0, 'Proteomics': 0,
                            'Physico-chemical': 0}
    choose_omics = st.multiselect('Which omic do you want to see:',
                                  list(example_1_omics_dict.keys()))

    num_of_columns = 0
    for i in choose_omics:
        example_1_omics_dict[i] = 1
        num_of_columns += 1

    if num_of_columns >= 2:
        column_list = st.beta_columns(num_of_columns)
        curr_pos = 0

        for i in list(example_1_omics_dict.keys()):
            if example_1_omics_dict[i] != 0:
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
        for i in list(example_1_omics_dict.keys()):
            if example_1_omics_dict[i] != 0:
                if i == 'Genomics':
                    create_main_example_1_genomics()
                elif i == 'Metabolomics':
                    create_main_example_1_metabolomics()
                elif i == 'Proteomics':
                    create_main_example_1_proteomics()
                else:
                    create_main_example_1_phy_che()
    return None
