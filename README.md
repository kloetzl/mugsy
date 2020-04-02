Mugsy - multiple whole genome alignment tool

This is a fork of the original sourceforge project [1]. It features bug fixes and performance improvements. However, at the moment compiling is still not ideal. Use with care.

## Installation

Good luck.

    make


## Usage

Mugsy requires two additional renvironment variables.

    PERL5LIB=$PWD MUGSY_INSTALL=$PWD ./mugsy genome1.fa genome2.fa genome3.fa

The output is written to `/tmp/tmp.maf`.


## License

This code is a chimera of at least the following projects:

- seqan
- mummer
- multiz
- mugsy  Artistic License 2.0


## Citation

Angiuoli SV, Salzberg SL. Mugsy: Fast multiple alignment of closely related whole genomes. Bioinformatics. 2010 Dec 9.


1: https://sourceforge.net/projects/mugsy/
