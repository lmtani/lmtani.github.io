FROM mambaorg/micromamba:0.15.3

COPY --chown=micromamba:micromamba environment.yml /tmp/env.yaml

RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes