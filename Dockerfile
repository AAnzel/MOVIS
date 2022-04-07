# syntax=docker/dockerfile:1

FROM continuumio/miniconda3
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

RUN apt update
RUN apt install -y fontconfig
RUN export PATH=~/anaconda3/bin:$PATH

WORKDIR .

RUN export FONTCONFIG_FILE=/etc/fonts/fonts.conf
RUN export FONTCONFIG_PATH=/etc/fonts/
RUN export FILECONFIG_PATH=/etc/fonts

RUN conda install --yes \
        numpy==1.21.2 \
        pandas==1.3.5 \
        scikit-learn==1.0.2 \
        scipy==1.7.3 \
        protobuf==3.19.1 \
        python-levenshtein==0.12.2

RUN conda config --add channels conda-forge

RUN conda install --yes \
        altair==4.1.0 \
        biopython==1.78 \
        gensim==4.0.1 \
	vega-cli==5.17.0 \
	vega-lite-cli==4.17.0 \
        altair_saver==0.5.0 \
	click==7.1.2

RUN conda clean --all --yes

RUN pip install streamlit==1.5.1

COPY . .

WORKDIR ./Source/

EXPOSE 8501

CMD streamlit run main.py --browser.gatherUsageStats False
