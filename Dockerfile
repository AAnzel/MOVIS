# syntax=docker/dockerfile:1

FROM continuumio/miniconda3

WORKDIR .

RUN /opt/conda/bin/conda install --yes \
        numpy==1.21.2 \
        pandas==1.3.5 \
        scikit-learn==1.0.2 \
        scipy==1.7.3 \
        protobuf==3.19.1 \
        python-levenshtein==0.12.2 \
        nomkl

RUN /opt/conda/bin/conda config --add channels conda-forge

RUN /opt/conda/bin/conda install --yes \
        altair==4.1.0 \
        biopython==1.78 \
        gensim==4.0.1 \
        altair_saver==0.5.0

RUN /opt/conda/bin/conda clean --all --yes

RUN find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete

RUN pip install streamlit==1.2.0

COPY . .

WORKDIR ./Source/

EXPOSE 8501

CMD streamlit run main.py --browser.gatherUsageStats False
