[![Rust](https://github.com/stuart-lab/f2m/actions/workflows/rust.yml/badge.svg)](https://github.com/stuart-lab/f2m/actions/workflows/rust.yml)

# f2m: fragments to matrix

Create region x cell matrix from a fragment file.

```
f2m -f <fragments.tsv.gz> -b <peaks.bed> -c <cells.txt> -o <output>
```

## Installation

Clone the git repo:

```
git clone git@github.com:stuart-lab/f2m.git
cd f2m
```

Compile:

```
cargo build --release
```
