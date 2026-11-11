# Container for the Nipah edition: only the ENVIRONMENT (snakemake, mafft, pixy ...) is baked in; the project folder is mounted at /app
# (docker-compose.yml), so your config/, data/ and output/ stay on your computer. mfeprimer is not included (see README).
# Built for linux/amd64 on purpose: it is also what Docker Desktop runs on Windows and (via emulation) on Apple-silicon Macs,
# and mfeprimer's official binary is linux-amd64.
FROM --platform=linux/amd64 mambaorg/micromamba:1.5.10
USER root
RUN apt-get update && apt-get install -y --no-install-recommends git ca-certificates curl && rm -rf /var/lib/apt/lists/*
COPY environment.yml /tmp/environment.yml
RUN micromamba create -y -n primer_validation -f /tmp/environment.yml && micromamba clean -a -y
ENV PATH=/opt/conda/envs/primer_validation/bin:$PATH \
    MFEPRIMER_BIN=/usr/local/bin/mfeprimer \
    PYTHONUNBUFFERED=1 \
    HOME=/tmp
WORKDIR /app
EXPOSE 8767
# the container listens on all its interfaces; docker-compose.yml publishes the port on 127.0.0.1 only (this computer)
CMD ["python3", "webapp/server.py", "--host", "0.0.0.0", "--port", "8767"]
