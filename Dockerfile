FROM cyversevice/shiny-geospatial

LABEL maintainer="Daniel Leon Perinan <dleoper@upo.es>"

RUN apt-get update && apt-get install -y --no-install-recommends \
    sudo \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages('shinydashboard')"
RUN R -e "install.packages('shinycssloaders')"
RUN R -e "install.packages('plotly')"
RUN R -e "install.packages('spdep')"
RUN R -e "install.packages('forecast')"
RUN R -e "install.packages('future')"
RUN R -e "install.packages('ipc')"
RUN R -e "install.packages('segclust2d')"
RUN R -e "install.packages('dplR')"
RUN R -e "install.packages('tsmp')"
RUN R -e "install.packages('mclust')"
RUN R -e "install.packages('tsmp')"
RUN R -e "install.packages('lattice')"
RUN R -e "install.packages('TTR')"
RUN R -e "install.packages('egg')"
RUN R -e "install.packages('readxl')"
RUN R -e "install.packages('RCurl')"
RUN R -e "install.packages('shinythemes')"
RUN R -e "install.packages('broom')"
RUN R -e "install.packages('ggpubr')"
RUN R -e "install.packages('htmlwidgets')"
RUN R -e "install.packages('rapportools')"
RUN R -e "install.packages('segmenTier')"
RUN R -e "devtools::install_github('romainfrancois/tie')"
RUN R -e "install.packages('VLTimeCausality')"
RUN R -e "install.packages('network')"
RUN R -e "install.packages('sna')"
RUN R -e "install.packages('GGally')"
RUN R -e "if (!requireNamespace ('BiocManager', quietly = TRUE)) { install.packages('BiocManager') }"
RUN R -e "BiocManager::install('graph'); BiocManager::install('RBGL')"
RUN sudo apt-get update -y; sudo apt-get install -y libglpk-dev
RUN R -e "install.packages('igraph'); install.packages('pcalg')"
RUN R -e "install.packages('networkD3')"
RUN R -e "install.packages('PerformanceAnalytics')"

RUN addgroup --system app \
    && adduser --system --ingroup app app

WORKDIR /home/chromo

COPY app app
COPY examples examples

RUN chown app:app -R /home/chromo

USER app

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/home/chromo/app', host = '0.0.0.0', port = 3838)"]
