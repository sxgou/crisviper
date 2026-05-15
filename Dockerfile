# CrisViper — CARLIN sequence analysis pipeline
# Multi-stage Docker build for minimal image size

FROM python:3.10-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY pyproject.toml ./
COPY crisviper/ crisviper/

RUN pip install --no-cache-dir --user .

FROM python:3.10-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    procps && \
    rm -rf /var/lib/apt/lists/*

COPY --from=builder /root/.local /root/.local

WORKDIR /data
ENV PATH=/root/.local/bin:$PATH

ENTRYPOINT ["crisviper"]
CMD ["--help"]
