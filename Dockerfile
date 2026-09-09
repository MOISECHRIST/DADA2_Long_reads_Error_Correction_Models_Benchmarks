FROM quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0

LABEL maintainer="Moise Meka <moise.meka@students.unibe.ch>"
LABEL description="R environment for DADA2 pipeline (QC/trim + learning error + dada()) -- dada2 1.38.0 (biocontainer bioconda, figé), tidyr 1.3.2 ajouté par-dessus"

USER root

RUN R -e 'install.packages("remotes", repos = "https://cloud.r-project.org")'
RUN R -e ' \
    remotes::install_version("tidyr", version = "1.3.2", \
                              repos = "https://cloud.r-project.org", upgrade = "never"); \
    if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr install failed"); \
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr install failed (dependance de tidyr)") \
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

WORKDIR /workdir

CMD ["R"]