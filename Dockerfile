FROM rust:1.85-slim AS builder

WORKDIR /build
COPY Cargo.toml Cargo.lock* ./
COPY src/ src/

RUN cargo build --release

FROM debian:bookworm-slim

LABEL org.opencontainers.image.source="https://github.com/ayobi/emits"
LABEL org.opencontainers.image.description="EMITS: EM abundance estimation for fungal ITS from long-read sequencing"
LABEL org.opencontainers.image.licenses="MIT"

RUN apt-get update && apt-get install -y --no-install-recommends \
    minimap2 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/target/release/emits /usr/local/bin/emits

ENTRYPOINT ["emits"]
CMD ["--help"]
