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

def create_main_shared():
    
    # TODO: Be careful when placing values in beta_columns
    # These are hardcoded values for column width
    tmp_col_1, tmp_col_2, tmp_col_3 = st.beta_columns([1.5, 2, 1])
    tmp_col_2.title('Tool title - Multi-omics time series')
    st.markdown(' ')
    st.markdown(' ')
    st.markdown('---')

    return None

def show_data_set(df):
    with st.beta_expander('Show the data set and related info', expanded = True):
        st.markdown('First 100 entries')
        st.dataframe(df.head(100))
        st.dataframe(df.describe())

    return None

def create_main_example_2():
    create_main_shared()

    return None