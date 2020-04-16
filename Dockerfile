FROM ubuntu:18.04

# For opencontainers label definitions, see:
#    https://github.com/opencontainers/image-spec/blob/master/annotations.md
LABEL org.opencontainers.image.title="HyP3 Geocode"
LABEL org.opencontainers.image.description="HyP3 plugin to create geocoded products from Sentinel-1 granules"
LABEL org.opencontainers.image.vendor="Alaska Satellite Facility"
LABEL org.opencontainers.image.authors="ASF APD/Tools Team <uaf-asf-apd@alaska.edu>"
LABEL org.opencontainers.image.licenses="BSD-3-Clause"
LABEL org.opencontainers.image.url="https://github.com/asfadmin/hyp3-geocode"
LABEL org.opencontainers.image.source="https://github.com/asfadmin/hyp3-geocode"
# LABEL org.opencontainers.image.documentation=""

# Dynamic lables to define at build time via `docker build --label`
# LABEL org.opencontainers.image.created=""
# LABEL org.opencontainers.image.version=""
# LABEL org.opencontainers.image.revision=""

ARG DEBIAN_FRONTEND=noninteractive
ENV PYTHONDONTWRITEBYTECODE=true

RUN apt-get update && apt-get upgrade -y && apt-get install -y software-properties-common && \
    add-apt-repository -y ppa:ubuntugis/ubuntugis-unstable && apt-get update && \
    apt-get install -y unzip vim wget curl gdal-bin libgdal-dev gimp \
    gnuplot gnuplot-qt libblas-dev libblas3 libfftw3-dev \
    libgtk2.0-bin libgtk2.0-common libgtk2.0-dev libhdf5-dev \
    liblapack-dev liblapack3 python3-dev python3-pip python3-h5py python3-matplotlib python3-scipy && \
    apt-get clean

RUN pip3 install --upgrade pip && \
    python3 -m pip install --upgrade numpy scipy statsmodels scikit-image

COPY GAMMA_SOFTWARE-20170707 /opt/gamma/

ARG S3_PYPI_HOST

RUN python3 -m pip install --no-cache-dir hyp3_geocode \
    --trusted-host "${S3_PYPI_HOST}" \
    --extra-index-url "http://${S3_PYPI_HOST}"

ARG CONDA_GID=1000
ARG CONDA_UID=1000

RUN groupadd -g "${CONDA_GID}" --system conda && \
    useradd -l -u "${CONDA_UID}" -g "${CONDA_GID}" --system -d /home/conda -m  -s /bin/bash conda && \
    chown -R conda:conda /opt && \
    echo "export GAMMA_HOME=/opt/gamma" >> /home/conda/.bashrc && \
    echo "export MSP_HOME=$GAMMA_HOME/MSP" >> /home/conda/.bashrc && \
    echo "export ISP_HOME=$GAMMA_HOME/ISP" >> /home/conda/.bashrc && \
    echo "export DIFF_HOME=$GAMMA_HOME/DIFF" >> /home/conda/.bashrc && \
    echo "export DISP_HOME=$GAMMA_HOME/DISP" >> /home/conda/.bashrc && \
    echo "export LAT_HOME=$GAMMA_HOME/LAT" >> /home/conda/.bashrc && \
    echo "export IPTA_HOME=$GAMMA_HOME/IPTA" >> /home/conda/.bashrc && \
    echo "export GEO_HOME=$GAMMA_HOME/GEO" >> /home/conda/.bashrc && \
    echo "export PATH=$PATH:$MSP_HOME/bin:$ISP_HOME/bin:$DIFF_HOME/bin:$LAT_HOME/bin:$DISP_HOME/bin:$IPTA_HOME/bin:$GEO_HOME/bin" >> /home/conda/.bashrc && \
    echo "export PATH=$PATH:$MSP_HOME/scripts:$ISP_HOME/scripts:$DIFF_HOME/scripts:$LAT_HOME/scripts:$IPTA_HOME/scripts"  >> /home/conda/.bashrc && \
    echo "export GAMMA_RASTER=BMP" >> /home/conda/.bashrc

USER ${CONDA_UID}
SHELL ["/bin/bash", "-l", "-c"]
ENV PYTHONDONTWRITEBYTECODE=true
WORKDIR /home/conda/

ENTRYPOINT ["/usr/local/bin/hyp3_geocode"]
CMD ["-v"]
