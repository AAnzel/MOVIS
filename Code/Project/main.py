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
import example_1
import example_2
import upload

def create_sidebar_and_main():
    st.sidebar.markdown('Above should be a logo. Below should be a github repo\
                        logo that is clickable, and on the right should be\
                        the link to the paper.')
    
    # We need 5 columns so that we have nicely aligned images
    col_1, col_2, col_3, col_4, col_5 = st.sidebar.beta_columns([1, 2, 1, 2, 1])
    col_2.image(os.path.join('images', 'GitHub-Mark-120px-plus.png'), width = 52)
    col_2.markdown('[GitHub](https://github.com/AAnzel/Multi-omics_platform)')
    col_4.markdown('Paper doi with journal logo')
    
    st.sidebar.markdown('---')
    st.sidebar.markdown('Here put some info about the app.')
    st.sidebar.markdown('---')

    st.sidebar.markdown('**Choose your data set:**')
    choice_data_set = st.sidebar.radio('', ('Example 1', 'Example 2', 'Upload'),
                                       index = 0)

    if choice_data_set == 'Example 1':
        example_1.create_main_example_1()
    elif choice_data_set == 'Example 2':
        example_2.create_main_example_2()
    else:
        upload.create_main_upload()
    
    return None

def main():
    
    st.set_page_config(layout = 'wide')

    # I need to run this just once, so I create cache
    #omics_run.example_1_calc_phy_che()
    #omics_run.example_1_calc_metabolomics()
    
    create_sidebar_and_main()

    return None


if __name__ == '__main__':
    main()