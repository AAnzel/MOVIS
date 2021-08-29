import streamlit as st


def create_home():
    st.info('''Welcome to **MOVIS** - a time-series multi-omics visualization
               tool. MOVIS allows you to explore multiple time-series omics
               data sets at once so that you can easily pinpoint anomalies
               and/or patterns in your temporal data.''')

    st.markdown('''
                Currently, MOVIS supports *genomics*, *proteomics*,
                *metabolomics*, *transcriptomics*, and *physico-chemical*
                time-series data sets.

                You can see the examples of what kind of data sets MOVIS
                supports on *Example 1* and *Example 2* pages on the sidebar.
                Example 1 consists of genomics, proteomics, metabolomics, and
                physico-chemical time-series data. Example 2 contains only
                transcriptomics time-series data.

                ---

                If you do not know how to use MOVIS, or you just want to see
                how MOVIS works, check out our documentation (link on the
                sidebar). If you want to start exploring your data, head out to
                the *Upload* page.
                ''')

    st.markdown('''
                   ---

                   ### Changelog
                   | Date | Info |
                   | --- | --- |
                   | 03.09.2021. | MOVIS v1.0.0 released    🎉 🎈 🍾|
                ''')

    return None
