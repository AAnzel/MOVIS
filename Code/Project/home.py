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

def create_home():
    st.markdown('''
                Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do\
                eiusmod tempor incididunt ut labore et dolore magna aliqua.\
                Mattis molestie a iaculis at erat. Magnis dis parturient montes\
                nascetur ridiculus mus mauris. Enim ut tellus elementum\
                sagittis vitae et leo duis ut. Nulla posuere sollicitudin\
                aliquam ultrices. Quisque non tellus orci ac. Ut venenatis\
                tellus in metus vulputate eu. Eget nunc scelerisque viverra\
                mauris in. Imperdiet sed euismod nisi porta lorem mollis\
                aliquam ut. Morbi tristique senectus et netus. Egestas sed sed\
                risus pretium quam vulputate dignissim suspendisse. Odio eu\
                feugiat pretium nibh. Eu tincidunt tortor aliquam nulla\
                facilisi cras fermentum odio. Netus et malesuada fames ac. Nibh\
                praesent tristique magna sit amet purus gravida quis. Urna
                condimentum mattis pellentesque id nibh tortor id aliquet\
                lectus. Id donec ultrices tincidunt arcu non. Pellentesque sit\
                amet porttitor eget dolor morbi non arcu.
                
                Amet mauris commodo quis imperdiet massa. Id volutpat lacus\
                laoreet non curabitur gravida arcu. Laoreet suspendisse\
                interdum consectetur libero id faucibus nisl tincidunt.\
                Congue nisi vitae suscipit tellus mauris a. Sed velit\
                dignissim sodales ut. Pretium fusce id velit ut tortor pretium\
                viverra suspendisse. Morbi quis commodo odio aenean sed\
                adipiscing diam. Purus in massa tempor nec feugiat nisl.\
                Nunc non blandit massa enim nec dui nunc mattis. Dictum at\
                tempor commodo ullamcorper a lacus vestibulum sed.
                ''')
    
    return None
