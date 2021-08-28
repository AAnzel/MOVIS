import streamlit as st


def create_home():
    st.markdown('''
                Welcome to **MOVIS** - a time-series multi-omics visualization
                tool. MOVIS allows you to explore multiple time-series omics
                data sets at once, so that you can easily pin-point anomalies
                and/or patterns in your temporal data.

                Currently, MOVIS supports:
                *genomics*, *proteomics*, *metabolomics*, *transcriptomics*,
                and *physico-chemical* data sets.



                If you have any problems with MOVIS, please make an issue
                ticket so that we can fix it for you. Click on
                [this link](https://github.com/AAnzel/MOVIS/issues/new/choose)
                and then create a new *Bug report*.

                If you want to propose a new feature to MOVIS, click on
                [this link](https://github.com/AAnzel/MOVIS/issues/new/choose)
                and then create a new *Feature request*.

                If you have technical experience and want collaborate and
                improve MOVIS, check our
                [GitHub repo](https://github.com/AAnzel/MOVIS).
                ''')

    return None
