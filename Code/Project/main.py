import os
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


def main_and_sidebar():
    st.title('Tool title - ML cruncher')
    st.info ('''
            * Describe what kind of data we allow.
            * Possibly add example data set in a table
            * By using this tool and uploading your files you agree that you
            are accepting the licence agreement (put link).
            ''')
    
    st.subheader('Dataset')

    imported_file = st.file_uploader ('Upload your data set here. Maximum size is 200MB.',
    type = ['csv', 'xlsx', 'xls'])
    delimiter = st.selectbox ('Select the delimiter in your data set', ['Comma (,)', 'Semicolon (;)', 'Tab (\\t)', 'Excel File'])




def main():
    main_and_sidebar()




if __name__ == '__main__':
    main()