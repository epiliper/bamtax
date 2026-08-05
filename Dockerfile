# syntax=docker/dockerfile:1

FROM ubuntu:22.04 AS builder

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        build-essential \
        ca-certificates \
        clang \
        cmake \
        curl \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libssl-dev \
        pkg-config \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN curl --proto '=https' --tlsv1.2 --fail --silent --show-error \
        https://sh.rustup.rs \
        | sh -s -- -y --profile minimal --default-toolchain stable

ENV PATH="/root/.cargo/bin:${PATH}"

WORKDIR /build/source
ADD https://github.com/epiliper/bamtax.git#release .
RUN --mount=type=cache,target=/usr/local/cargo/registry \
    cargo build --locked --release

FROM ubuntu:22.04

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        ca-certificates \
        libbz2-1.0 \
        libcurl4 \
        liblzma5 \
        libssl3 \
        procps \
        zlib1g \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/source/target/release/bamtax /usr/local/bin/bamtax

WORKDIR /work

CMD ["bamtax", "--help"]
