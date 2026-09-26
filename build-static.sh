#!/usr/bin/env bash
# Build the static, reproducible Linux binary.
#
#   ./build-static.sh            -> target/x86_64-unknown-linux-musl/release/expansionscout
#
# Same source, same Cargo.lock and the compiler pinned in rust-toolchain.toml
# give the same bytes, wherever the checkout is: the absolute paths the
# compiler would otherwise embed are remapped to fixed prefixes, and the final
# link uses the rust-lld that ships inside that toolchain rather than the
# system's cc and ld, whose version would otherwise change the output from
# one distribution to the next. (A system cc is still needed: it links the
# dependencies' build scripts and proc-macros, which run during the build and
# are not part of the binary.) The binary is statically linked against musl
# (whose startup objects also come with the toolchain), so it runs on any
# x86-64 Linux with nothing installed, and it carries the STRchive release it
# was built with.
#
# Check a binary against a published one with `sha256sum`.

set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cargo_home="${CARGO_HOME:-$HOME/.cargo}"

cd "$here"
RUSTFLAGS="-C linker=rust-lld -C linker-flavor=ld.lld --remap-path-prefix=$here=/build --remap-path-prefix=$cargo_home=/cargo" \
    cargo build --release --locked --target x86_64-unknown-linux-musl

bin="$here/target/x86_64-unknown-linux-musl/release/expansionscout"
"$bin" --version
sha256sum "$bin"
