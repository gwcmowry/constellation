# Prebuilt Linux releases

The initial supported binary is `x86_64-unknown-linux-musl`: a static Linux
x86_64 executable built for the baseline x86-64 CPU target. It has no dynamic
loader or shared-library dependencies. It includes gzip support; optional zstd
streams still require a separate `zstd` executable. macOS, Windows and ARM
packages have not been validated or enabled.

## What users receive

- `constellation-VERSION-x86_64-unknown-linux-musl.tar.gz` and its `.sha256` file.
- The executable, a `~/.local/bin` installer, and the tiny runnable example.
- Project and dependency license/copyright notices, with source copies where
  included by the packaging policy.
- `BUILD-INFO.json`: source commit, dirty flag, source-file hashes, toolchain,
  build flags, binary SHA-256 and smoke-test result.
- `MANIFEST.sha256` for all extracted files, checked by the installer.

Use the archive's README for installation. No source checkout or compiler is
needed by the user. The published package names omit the `-dev-` suffix used for
local snapshots. Checksums detect changed downloads; they are not signatures.

## Local packaging

Maintainers need Python 3.12+, GCC/binutils, `musl-gcc` (the `musl-tools` package
on Debian/Ubuntu), and Rust 1.93.1 with its musl target and documentation/license
files. The package script uses only Python's standard library.

```bash
rustup toolchain install 1.93.1 --profile minimal \
  --component rustfmt --component rust-docs --target x86_64-unknown-linux-musl
CC_x86_64_unknown_linux_musl=musl-gcc \
CARGO_TARGET_X86_64_UNKNOWN_LINUX_MUSL_LINKER=musl-gcc \
  python3 scripts/package_release.py --allow-dirty
```

Artifacts are written to ignored `dist/`. Without `--tag`, packages are
development snapshots, including on a clean tree. `--allow-dirty` records the
uncommitted state explicitly; it cannot be combined with a tagged release.
Normal `cargo build --release` keeps debug symbols for profiling; the separate
`dist` profile strips the distributable binary.

The script builds with locked dependencies and generic CPU flags, checks for
dynamic dependencies, runs both conservative and USA counting on the example,
creates the archive, extracts it, verifies its checksums, installs it into a
temporary prefix containing spaces, and repeats the example with the installed
binary. It never runs the 100M/500M benchmarks or writes to the user's usual
installation directory.

## Publishing a version

1. Review and commit the intended code, documentation, fixtures and release
   files. Bump package versions and `Cargo.lock` when preparing a new version.
2. Push the committed snapshot and an annotated `vVERSION` tag whose version
   matches `crates/constellation-cli/Cargo.toml`, for example `v0.1.0`.
3. The [release workflow](../.github/workflows/release.yml) runs formatting and
   the workspace tests on the musl target, then builds and smoke-tests the
   package. A separate job with release-write permission publishes the tested
   archive and checksum to that existing tag. A mismatched tag or dirty source
   tree fails; an existing release is not overwritten.

Manual **Run workflow** builds a development artifact for inspection without
creating a GitHub release. Tag pushes publish only after successful checks.
This setup does not itself commit, push, tag or publish the current working tree.

## Dependency notices

The inventory follows normal/build dependency edges from the CLI; it may include
build-time crates and feature-unified workspace dependencies as well as linked
code. Notices bundled in the original crate are copied. Missing upstream notices
are retained under `release/third-party-notices` with source URLs, immutable
revisions and checksums. The recorded `fxhash` lookup has no files; no notice
from that lookup is used.

Where a crate and its upstream revision provide only a compatible MIT license
declaration, the package includes its unmodified source/author metadata and the
declared standard MIT text without inventing a copyright holder. MPL-covered
crate source is also included. Rust distribution notices and the musl package
copyright file accompany the static runtime. Changing dependencies may require
updating these recorded notices; unresolved cases fail packaging.

Originals: [GitHub Rust workflows](https://docs.github.com/en/actions/tutorials/build-and-test-code/rust),
[GitHub releases API](https://docs.github.com/en/rest/releases/releases), and
the per-file upstream URLs in
[notice provenance](../release/third-party-notices/provenance.json).
