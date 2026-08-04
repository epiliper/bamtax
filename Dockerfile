# syntax=docker/dockerfile:1

FROM lukemathwalker/cargo-chef:latest-rust-1-bookworm AS source

WORKDIR /build/source
ADD https://github.com/epiliper/bamtax.git#release .

FROM lukemathwalker/cargo-chef:latest-rust-1-bookworm AS planner

WORKDIR /build/source
COPY --from=source /build/source .
RUN cargo chef prepare --recipe-path recipe.json

FROM lukemathwalker/cargo-chef:latest-rust-1-bookworm AS builder

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        clang \
        cmake \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libssl-dev \
        pkg-config \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build/source
COPY --from=planner /build/source/recipe.json recipe.json
RUN --mount=type=cache,target=/usr/local/cargo/registry \
    cargo chef cook --locked --release --recipe-path recipe.json

COPY --from=source /build/source .
RUN --mount=type=cache,target=/usr/local/cargo/registry \
    cargo build --locked --release

FROM debian:bookworm-slim

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
