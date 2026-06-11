# bitags

Fast and flexible suite of operations for barcoded NGS reads.

> **Warning**: This package is in early alpha. 
The API may change without notice and it is not yet recommended for production use.

## About
Bitags is a package designed to simplify working with barcoded NGS reads.
It relies on an adapted version of the bitap algorithm to perform fuzzy
matching of tags composing a barcode (usually from a combinatorial barcoding method).
Moreover, it provides additional utility functions for barcoded reads.

Currently supported features:

- **Barcoding** — identify and assign barcodes to reads using fuzzy bitap matching
- **Trimming** — trim reads based on a barcoding schema
- **Classification** — classify reads according to a barcoding schema
- **Manipulation** — embed barcodes, convert paired to unpaired reads
- **I/O** — read and write FASTQ files via lazy Polars DataFrames
- **Visualization** — visualize read structure in terminal, html or png.

Bitags can be used as a command line tool or through the Python API.

## Installation

```bash
pip install bitags
```

Or from source by cloning and using maturin to build.  

## Usage

*Proper documentation will be available soon, for now please rely on docstrings.*

## Cite

*At the moment, there is no paper on Bitags. If you find this package useful for
your work, please cite it by providing the link to the [GitHub page](https://github.com/StefanoCretti/Bitags).*

## License

BSD 3-Clause License. See [LICENSE](LICENSE) for details.
