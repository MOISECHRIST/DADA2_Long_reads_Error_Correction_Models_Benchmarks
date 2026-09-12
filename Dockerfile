FROM quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0

LABEL maintainer="Moise Meka <moise.meka@students.unibe.ch>"
LABEL description="R environment for DADA2 pipeline (QC/trim + learning error + dada())"

USER root

ENV DEBIAN_FRONTEND=noninteractive

RUN R -e 'install.packages("remotes", repos = "https://cloud.r-project.org")'
RUN R -e ' \
    remotes::install_version("tidyr", version = "1.3.2", \
                              repos = "https://cloud.r-project.org", upgrade = "never"); \
    if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr install failed"); \
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr install failed (dependance de tidyr)") \
    '
    
RUN R -e ' \
    remotes::install_version("wavethresh", version = "4.7.3", \
                              repos = "https://cloud.r-project.org", upgrade = "never"); \
    if (!requireNamespace("wavethresh", quietly = TRUE)) stop("wavethresh install failed") \
    '

RUN R -e ' \
    stopifnot( \
      as.character(packageVersion("dada2")) == "1.38.0", \
      as.character(packageVersion("tidyr")) == "1.3.2", \
      packageVersion("dplyr") >= "1.0.0" \
    ); \
    library(dada2); \
    library(ggplot2); \
    library(tidyr); \
    cat("OK -- versions verifiees :\n"); \
    print(packageVersion("dada2")); \
    print(packageVersion("ggplot2")); \
    print(packageVersion("tidyr")); \
    print(packageVersion("dplyr")) \
    '
ENV R_ENVIRON_USER=/dev/null
ENV R_PROFILE_USER=/dev/null
ENV R_LIBS_USER=""
ENV TZ=UTC
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

WORKDIR /workdir

CMD ["R"]