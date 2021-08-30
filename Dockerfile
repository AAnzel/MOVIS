# syntax=docker/dockerfile:1

FROM continuumio/miniconda3

WORKDIR .

RUN /opt/conda/bin/conda install --yes \
        numpy==1.20.3 \
        pandas==1.3.2 \
        scikit-learn==0.24.2 \
        scipy==1.6.2 \
        protobuf==3.14.0 \
        nomkl

RUN /opt/conda/bin/conda config --add channels conda-forge

RUN /opt/conda/bin/conda install --yes \
        streamlit==0.81.1 \
        altair==4.1.0 \
        biopython==1.78 \
        gensim==4.0.1

RUN /opt/conda/bin/conda clean --all --yes

RUN find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete

COPY . .

WORKDIR ./Source/

EXPOSE 8501

CMD streamlit run main.py --browser.gatherUsageStats False
