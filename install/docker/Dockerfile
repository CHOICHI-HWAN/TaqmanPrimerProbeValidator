# Image for `setup.sh --backend docker`: the pipeline environment plus this installer. The panels, their data and results
# live in the host folder mounted at /work, so they survive rebuilding or removing the container.
FROM mambaorg/micromamba:1.5.10
USER root
RUN apt-get update && apt-get install -y --no-install-recommends git ca-certificates curl && rm -rf /var/lib/apt/lists/*
COPY install/environment.yml install/environment.portable.yml /tmp/
# The exact versions are built for linux-64; on arm64 (Apple Silicon) fall back to the unpinned file.
ARG TARGETARCH
RUN if [ "$TARGETARCH" = "amd64" ] || [ -z "$TARGETARCH" ]; then f=/tmp/environment.yml; else f=/tmp/environment.portable.yml; fi \
    && micromamba create -y -n chichipes_ui -f "$f" && micromamba clean -a -y && chmod -R a+rX /opt/conda
COPY . /opt/chichipes
ENV PATH=/opt/conda/envs/chichipes_ui/bin:$PATH CHICHIPES_HOST=0.0.0.0 HOME=/work/.home USER=chichipes LOGNAME=chichipes
WORKDIR /work
ENTRYPOINT ["bash", "/opt/chichipes/install/docker/entrypoint.sh"]
