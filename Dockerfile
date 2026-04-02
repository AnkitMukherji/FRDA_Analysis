FROM bioconductor/bioconductor_docker:RELEASE_3_22

# ---- system dependencies ----
RUN apt-get update && apt-get install -y \
    git \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# ---- renv library ----
ENV RENV_PATHS_CACHE=/renv/cache
ENV RENV_PATHS_LIBRARY=/renv/library

RUN mkdir -p /renv/cache /renv/library \
    && chown -R rstudio:rstudio /renv

# ---- project directory ----
WORKDIR /workspace/frda_analysis

COPY renv.lock renv.lock
COPY renv/activate.R renv/activate.R
COPY .Rprofile .Rprofile

# ---- GitHub access (optional) ----
ARG GITHUB_PAT
ENV GITHUB_PAT=$GITHUB_PAT

# ---- install renv ----
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org'); \
          remotes::install_version('renv', version='1.1.7', repos='https://cloud.r-project.org')"

# ---- restore environment ----
RUN R -e "renv::restore(prompt = FALSE)"

# ---- remove token from final image ----
ENV GITHUB_PAT=""