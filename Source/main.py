import os
import streamlit as st

import example_1
import example_2
import use_case
import upload
import home


__author__ = 'Aleksandar Anžel'
__copyright__ = ''
__credits__ = ['Aleksandar Anžel', 'Georges Hattab']
__license__ = 'GNU General Public License v3.0'
__version__ = '1.0'
__maintainer__ = 'Aleksandar Anžel'
__email__ = 'aleksandar.anzel@uni-marburg.de'
__status__ = 'Dev'


def create_sidebar_and_main():

    tmp_col_1, tmp_col_2, tmp_col_3, tmp_col_4, tmp_col_5 = st.columns(5)
    tmp_col_3.title('MOVIS')
    st.markdown(' ')
    st.markdown(' ')
    st.markdown('---')

    st.sidebar.image(os.path.join('images', 'movis_logo_transparent.png'),
                     use_column_width=True)

    # TODO: Uncomment this after publication
    # We need 5 columns so that we have nicely aligned images
    # col_1, col_2, col_3, col_4, col_5 = st.sidebar.columns([1, 2, 1,
    #                                                             2, 1])
    # col_2.image(os.path.join('images', 'GitHub-Mark-120px-plus.png'),
    #             width=52)
    # col_2.markdown('[GitHub](https://github.com/AAnzel/MOVIS)')
    # col_4.markdown('Paper doi with journal logo')

    st.sidebar.markdown('---')
    col_1, col_2 = st.sidebar.columns([3, 1])
    col_2.image(os.path.join('images', 'GitHub-Mark-120px-plus.png'), width=54)
    col_2.markdown('[GitHub](https://github.com/AAnzel/MOVIS)')
    col_1.markdown('''
                      **MOVIS**

                      **M**ulti-**O**mics **VIS**ualization tool for
                      time-series data sets''')

    st.sidebar.markdown('''
                            ---

                            **Navigation:**''')

    # TODO: Add 'Example 3' below when the data is ready
    choice_data_set = st.sidebar.radio(
        '', ('Home', 'Example 1', 'Example 2', 'Use case', 'Upload'), index=0)

    st.sidebar.markdown('''
        ---
        **Bug report**: [report here](https://github.com/AAnzel/MOVIS/issues/new?assignees=AAnzel&labels=bug&template=bug_report.md&title=)

        **Feature request**: [propose here](https://github.com/AAnzel/MOVIS/issues/new?assignees=AAnzel&labels=enhancement&template=feature_request.md&title=)

        **Documentation**: [see here](https://github.com/AAnzel/MOVIS/wiki)

        ---
        ''') # noqa

    # Deleting old user-uploaded cached data
    upload.remove_cached_data()

    if choice_data_set == 'Home':
        home.create_home()
    elif choice_data_set == 'Example 1':
        example_1.create_main_example_1()
    elif choice_data_set == 'Example 2':
        example_2.create_main_example_2()
    elif choice_data_set == 'Use case':
        use_case.create_main_use_case()
    elif choice_data_set == 'Upload':
        upload.create_main_upload()
    else:
        pass

    return None


def main():

    st.set_page_config(page_title='MOVIS', layout='wide',
                       initial_sidebar_state='expanded',
                       page_icon=os.path.join('images',
                                              'movis_logo_transparent.png'))

    # TODO: unsafe_allow_html might be removed in the future and replaced
    # with a proper solution. Beware
    st.markdown('''
                    <style>
                    #MainMenu {visibility: hidden;}
                    footer {visibility: hidden;}
                    details {
                                display: None;
                            }
                    </style>

                ''', unsafe_allow_html=True)

    create_sidebar_and_main()

    return None


if __name__ == '__main__':
    main()
