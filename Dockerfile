FROM docker.1ms.run/library/python:3.9-slim

LABEL maintainer="siyangming"
LABEL software="eccsplorer_slim"
LABEL software.version="1.0.1"
LABEL description="Minimal ECCsplorer analysis scripts (Python+R only, no bioinformatics tools)"
LABEL org.opencontainers.image.source="https://github.com/crimBubble/ECCsplorer"

# Install R + bedtools (candidate_extract/coverage_profile) + build deps → compile → cleanup
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base r-cran-ggplot2 r-cran-ggrepel r-cran-gridextra r-cran-dplyr \
    bedtools procps build-essential gfortran \
    && pip install --no-cache-dir -i https://pypi.tuna.tsinghua.edu.cn/simple \
        numpy scipy biopython \
    && pip cache purge \
    && apt-get remove -y build-essential gfortran \
    && apt-get autoremove -y && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /root/.cache

# Copy scripts
COPY bin/ /opt/eccsplorer_slim/bin/
RUN chmod +x /opt/eccsplorer_slim/bin/*.py /opt/eccsplorer_slim/bin/*.R
ENV PATH="/opt/eccsplorer_slim/bin:${PATH}"

WORKDIR /data
CMD ["python", "-c", "import numpy, scipy; from Bio import SeqIO; print('ECCsplorer_slim ready')"]
