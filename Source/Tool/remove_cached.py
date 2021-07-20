import os
import shutil


# TODO: Change this for production
path_uploaded = 'uploaded_try'
path_uploaded_genomics = os.path.join(path_uploaded, 'genomics')
path_uploaded_proteomics = os.path.join(path_uploaded, 'proteomics')
path_uploaded_transcriptomics = os.path.join(path_uploaded, 'transcriptomics')
path_uploaded_metabolomics = os.path.join(path_uploaded, 'metabolomics')
path_uploaded_phy_che = os.path.join(path_uploaded, 'phy_che')

path_uploaded_dict = {
        'Genomics': path_uploaded_genomics,
        'Proteomics': path_uploaded_proteomics,
        'Transcriptomics': path_uploaded_transcriptomics,
        'Metabolomics': path_uploaded_metabolomics,
        'Physico-chemical': path_uploaded_phy_che
    }


def remove_cached_data(path_omic):

    for file_name in os.listdir(path_omic):

        file_path = os.path.join(path_omic, file_name)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except OSError as e:
            print('Failed to delete ' + file_path + '. Reason: ' + e)

    return None


def main():

    # Traversing every cache directory and deleting everything from it
    for omic in path_uploaded_dict:
        remove_cached_data(path_uploaded_dict(omic))
        print('Successfully removed everything from ' + omic + ' data set')

    return None


if __name__ == '__main__':
    main()
