[![Rust](https://github.com/stuart-lab/f2m/actions/workflows/rust.yml/badge.svg)](https://github.com/stuart-lab/f2m/actions/workflows/rust.yml)

# f2m: fragments to matrix

Create region x cell matrix from a fragment file.

```
f2m -f <fragments.tsv.gz> -b <peaks.bed> -c <cells.txt> -o <matrix.mtx>
```

We also provide a simple tool to select cell barcodes from the fragment file 
according to their total count:

```
cellselect -f <fragments.tsv.gz> -o <barcode_counts.tsv> -t <threshold> > barcodes.txt
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
