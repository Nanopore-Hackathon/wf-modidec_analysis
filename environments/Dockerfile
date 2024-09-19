FROM docker.io/mambaorg/micromamba:1.5-jammy
USER root

COPY environment.yml /opt/conda.yml

RUN micromamba install -y -n base -f /opt/conda.yml && micromamba clean -afy

ENV PATH "${MAMBA_ROOT_PREFIX}/bin:${PATH}"

RUN micromamba install -n base conda-forge::gcc

RUN micromamba run -n base pip install ont-remora

USER $MAMBA_USER






