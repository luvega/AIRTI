FROM mambaorg/micromamba:2.9.0@sha256:b62ed0c54940e3c801642d72ba7d2462f06356c378eba92d603d39f2ce5e4a0d

LABEL org.opencontainers.image.title="AIRTI unified reverse target fishing"
LABEL org.opencontainers.image.version="0.1.0"
LABEL org.opencontainers.image.source="https://github.com/luvega/AIRTI"
LABEL org.opencontainers.image.description="QuickVina2, fpocket, Meeko, Boltz-2, AmberTools, ParmEd and CUDA GROMACS in one audited runtime"

ENV PATH=/opt/conda/bin:${PATH} \
    NVIDIA_VISIBLE_DEVICES=all \
    NVIDIA_DRIVER_CAPABILITIES=compute,utility \
    CONDA_OVERRIDE_CUDA=12.8 \
    BOLTZ_CACHE=/models/boltz \
    AIRTI_ARTIFACT_ROOT=/data/airti-target-fishing \
    AIRTI_CACHE_ROOT=/models/airti-target-fishing \
    NUMBA_CACHE_DIR=/tmp/airti-cache/numba \
    MPLCONFIGDIR=/tmp/airti-cache/matplotlib \
    XDG_CACHE_HOME=/tmp/airti-cache/xdg \
    TRITON_CACHE_DIR=/tmp/airti-cache/triton \
    PYTHONUNBUFFERED=1

USER root
RUN mkdir -p /tmp/airti-cache/numba \
        /tmp/airti-cache/matplotlib \
        /tmp/airti-cache/xdg \
        /tmp/airti-cache/triton \
    && chmod 1777 /tmp/airti-cache \
        /tmp/airti-cache/numba \
        /tmp/airti-cache/matplotlib \
        /tmp/airti-cache/xdg \
        /tmp/airti-cache/triton
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER containers/environment.lock.yml /tmp/environment.lock.yml
RUN micromamba install --yes --name base --file /tmp/environment.lock.yml \
    && micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER containers/boltz.constraints.txt /tmp/boltz.constraints.txt
RUN /opt/conda/bin/python -m pip install --no-cache-dir \
      torch==2.7.1 \
      --index-url https://download.pytorch.org/whl/cu128 \
    && /opt/conda/bin/python -m pip install --no-cache-dir \
      --constraint /tmp/boltz.constraints.txt \
      'boltz[cuda]==2.2.1'

WORKDIR /opt/airti
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml README.md ./
COPY --chown=$MAMBA_USER:$MAMBA_USER src ./src
COPY --chown=$MAMBA_USER:$MAMBA_USER configs ./configs
COPY --chown=$MAMBA_USER:$MAMBA_USER templates ./templates
RUN /opt/conda/bin/python -m pip install --no-cache-dir --no-deps . \
    && /opt/conda/bin/python -m pip check

WORKDIR /work
ENTRYPOINT []
CMD ["airti-tf", "--help"]
