import os
import streamlit as st

import example_1
import example_2
import example_3
import upload
import home


# TODO: Implement transcriptomics in the upload page


def create_sidebar_and_main():

    tmp_col_1, tmp_col_2, tmp_col_3, tmp_col_4, tmp_col_5 = st.beta_columns(5)
    tmp_col_3.title('MOVIS')
    st.markdown(' ')
    st.markdown(' ')
    st.markdown('---')

    st.sidebar.image(os.path.join('images', 'movis_logo_transparent.png'),
                     use_column_width=True)

    # TODO: Uncomment this after publication
    # We need 5 columns so that we have nicely aligned images
    # col_1, col_2, col_3, col_4, col_5 = st.sidebar.beta_columns([1, 2, 1,
    #                                                             2, 1])
    # col_2.image(os.path.join('images', 'GitHub-Mark-120px-plus.png'),
    #             width=52)
    # col_2.markdown('[GitHub](https://github.com/AAnzel/MOVIS)')
    # col_4.markdown('Paper doi with journal logo')

    st.sidebar.markdown('---')
    col_1, col_2 = st.sidebar.beta_columns([1, 3])
    col_1.image(os.path.join('images', 'GitHub-Mark-120px-plus.png'), width=52)
    col_1.markdown('[GitHub](https://github.com/AAnzel/MOVIS)')
    col_2.markdown('**MOVIS** - **M**ulti-**O**mics **VIS**ualization tool\
                        for time series data sets')
    st.sidebar.markdown('---')

    st.sidebar.markdown('**Navigation:**')

    # TODO: Add 'Example 3' below when the data is ready
    choice_data_set = st.sidebar.radio(
        '', ('Home', 'Example 1', 'Example 2', 'Upload'), index=0)

    # Deleting old user-uploaded cached data
    upload.remove_cached_data()

    if choice_data_set == 'Home':
        home.create_home()
    elif choice_data_set == 'Example 1':
        example_1.create_main_example_1()
    elif choice_data_set == 'Example 2':
        example_2.create_main_example_2()
    elif choice_data_set == 'Example 3':
        example_3.create_main_example_3()
    else:
        upload.create_main_upload()

    return None


def main():

    st.set_page_config(page_title='MOVIS', layout='wide',
                       initial_sidebar_state='auto',
                       page_icon=os.path.join('images',
                                              'movis_logo_transparent.png'))

    create_sidebar_and_main()

    return None


if __name__ == '__main__':
    main()
