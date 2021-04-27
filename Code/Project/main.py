import os
import streamlit as st

import example_1
import example_2
import upload
import home


def create_sidebar_and_main():

    # TODO: Be careful when placing values in beta_columns
    # These are hardcoded values for column width
    tmp_col_1, tmp_col_2, tmp_col_3 = st.beta_columns([1.5, 2, 1])
    tmp_col_2.title('Tool title - Multi-omics time series')
    st.markdown(' ')
    st.markdown(' ')
    st.markdown('---')

    st.sidebar.markdown('Above should be a logo. Below should be a github repo\
                        logo that is clickable, and on the right should be\
                        the link to the paper.')

    # We need 5 columns so that we have nicely aligned images
    col_1, col_2, col_3, col_4, col_5 = st.sidebar.beta_columns([1, 2, 1,
                                                                 2, 1])
    col_2.image(os.path.join('images', 'GitHub-Mark-120px-plus.png'), width=52)
    col_2.markdown('[GitHub](https://github.com/AAnzel/Multi-omics_platform)')
    col_4.markdown('Paper doi with journal logo')

    st.sidebar.markdown('---')
    st.sidebar.markdown('Here put some info about the app.')
    st.sidebar.markdown('---')

    st.sidebar.markdown('**Navigation:**')
    choice_data_set = st.sidebar.radio('', ('Home', 'Example 1', 'Example 2',
                                            'Upload'), index=0)

    if choice_data_set == 'Home':
        home.create_home()
    elif choice_data_set == 'Example 1':
        example_1.create_main_example_1()
    elif choice_data_set == 'Example 2':
        example_2.create_main_example_2()
    else:
        upload.create_main_upload()

    return None


def main():

    st.set_page_config(page_title='Tool title', layout='wide',
                       initial_sidebar_state='auto')
    # , page_icon='path_to_image')

    # I need to run this just once, so I create cache
    # omics_run.example_1_calc_phy_che()
    # omics_run.example_1_calc_metabolomics()

    create_sidebar_and_main()

    return None


if __name__ == '__main__':
    main()
