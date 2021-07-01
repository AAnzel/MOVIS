import streamlit as st


def show_data_set(df):
    with st.beta_expander('Show the data set and related info', expanded=True):
        st.markdown('First 100 entries')
        st.dataframe(df.head(100))
        st.dataframe(df.describe())

    return None


def create_main_example_2():

    return None
