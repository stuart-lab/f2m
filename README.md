[![Rust](https://github.com/stuart-lab/f2m/actions/workflows/rust.yml/badge.svg)](https://github.com/stuart-lab/f2m/actions/workflows/rust.yml)

# f2m: fragments to matrix

Create region x cell matrix from a fragment file.

```
f2m -f <fragments.tsv.gz> -b <peaks.bed> -c <cells.txt> -o <matrix.mtx>
```

Processes at around 1 second per million lines in the fragment file.