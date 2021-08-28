import streamlit as st


def create_home():
    st.info('''Welcome to **MOVIS** - a time-series multi-omics visualization
               tool. MOVIS allows you to explore multiple time-series omics
               data sets at once so that you can easily pinpoint anomalies
               and/or patterns in your temporal data.''')

    st.markdown('''
                Currently, MOVIS supports *genomics*, *proteomics*,
                *metabolomics*, *transcriptomics*, and *physico-chemical*
                data sets.

                You can see the examples of what kind of data sets MOVIS
                supports on *Example 1* and *Example 2* pages on the left.
                Example 1 consists of genomics, proteomics, metabolomics, and
                physico-chemical data. Example 2 contains only transcriptomics
                data.

                ### Data set formats

                MOVIS supports various data set formats for each omic and we
                are constantly working on expanding the support to even more
                formats in the future. Current support is shown below:
                ''')

    col_1, col_2 = st.beta_columns(2)
    col_1.markdown('''
                    - Genomics
                        - Archived FASTA files
                        - Archived BIN files
                        - Tabular file
                    - Proteomics
                        - Archived FASTA files
                        - Tabular file
                    ''')

    col_2.markdown('''
                    - Metabolomics
                        - Tabular file
                    - Transcriptomics
                        - One or many tabular files
                    - Physico-chemical
                        - Tabular file
                    ''')
    ################################
    # TODO: Fix uploading CALCULATED PROTEOMICS DATA SET
    ################################
    st.warning('''For archived data sets, you must follow the naming convection
                  for archived files. Each file should be named as shown on the
                  *Upload* page.
                  **Archived** data sets are supported in the
                  following formats: TARBZ2, TBZ2, TARGZ, TGZ, TAR, TARXZ, TXZ,
                  ZIP. **Tabular** files are supported in CSV or TSV format.
                  ''')

    st.markdown('''
                   ---

                   ### Changelog
                   | Date | Info |
                   | --- | --- |
                   | 03.09.2021. | MOVIS v1.0.0 released    🎉 🎈 🍾|
                ''')

    return None
