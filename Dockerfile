# Dockerfile for metilene3
# https://github.com/zzhu1372/metilene3

FROM condaforge/miniforge3:24.3.0-0

LABEL org.opencontainers.image.source=https://github.com/zzhu1372/metilene3
LABEL org.opencontainers.image.description="metilene3: multi-condition DMR detection with auto-classification"
LABEL org.opencontainers.image.licenses=GPL-2.0

ARG METILENE3_COMMIT=v3.1.1

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        git \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN mamba create -y -n metilene3 -c bioconda -c conda-forge \
        python==3.10.0 \
        pandas \
        numpy \
        matplotlib \
        pandarallel \
        scikit-learn \
        seaborn \
        biopython \
        gseapy \
        r-base \
        bioconductor-ChIPseeker \
        bioconductor-org.Hs.eg.db \
        bioconductor-txdb.hsapiens.ucsc.hg19.knowngene \
        bioconductor-txdb.hsapiens.ucsc.hg38.knowngene \
    && mamba clean -afy

RUN git clone https://github.com/zzhu1372/metilene3.git /opt/metilene3 \
    && cd /opt/metilene3 \
    && git checkout ${METILENE3_COMMIT} \
    && make

WORKDIR /data

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "metilene3", "python", "/opt/metilene3/metilene3.py"]
CMD ["--help"]
