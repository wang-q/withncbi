# Name

withncbi - [egaz](https://github.com/wang-q/egaz) and [alignDB](https://github.com/wang-q/alignDB)
work with external (NCBI/EBI) data.

[TOC levels=1-3]: #

# Table of Contents
- [Name](#name)
- [Purpose](#purpose)
- [Directory organization](#directory-organization)
- [Conventions](#conventions)
  - [fasta](#fasta)
  - [fastq](#fastq)
- [Concepts](#concepts)
  - [Runlists in YAML](#runlists-in-yaml)
  - [Blocked fasta files](#blocked-fasta-files)
  - [Ranges and links of ranges](#ranges-and-links-of-ranges)
- [Author](#author)
- [Copyright and license](#copyright-and-license)


# Purpose

Fetch sequences, generate reports and build alignments according to various NCBI databases.

For more detailed, check `README.md` in each sub-directories.

# Directory organization

* [`db/`](db/): turn NCBI genome reports and assembly reports into a query-able MySQL database.

* [`ensembl/`](ensembl/): Ensembl related scripts.

* [`misc/`](misc/): miscellaneous projects.

* [`pop/`](pop/): build alignments on an whole Eukaryotes genus.

* [`taxon/`](taxon/): process (small) genomes according to NCBI Taxonomy.

* [`util/`](util/): miscellaneous utilities.

# Conventions

## fasta

* `.fa` - genomic sequences
* `.fas` - blocked fasta files
* `.fasta` - normal/miscellaneous fasta files

## fastq

Use `.fq` over `.fastq`

# Concepts

## Runlists in YAML

[App::RL](https://github.com/wang-q/App-RL)

[jrunlist](https://github.com/egateam/jrunlist)

## Blocked fasta files

Examples in [`example.fas`](https://github.com/wang-q/App-Fasops/blob/master/t/example.fas)

```text
>S288c.I(+):13267-13287|species=S288c
TCGTCAGTTGGTTGACCATTA
>YJM789.gi_151941327(-):5668-5688|species=YJM789
TCGTCAGTTGGTTGACCATTA
>RM11.gi_61385832(-):5590-5610|species=RM11
TCGTCAGTTGGTTGACCATTA
>Spar.gi_29362400(+):2477-2497|species=Spar
TCATCAGTTGGCAAACCGTTA

```

![blocked-fasta-files](doc/blocked-fasta-files.png)

[App::Fasops](https://github.com/wang-q/App-Fasops)

## Ranges and links of ranges

[App::Rangeops](https://github.com/wang-q/App-Rangeops)

[jrange](https://github.com/egateam/jrange)

# Author

Qiang Wang &lt;wang-q@outlook.com&gt;

# Copyright and license

This software is copyright (c) 2015 by Qiang Wang.

This is free software; you can redistribute it and/or modify it under the same terms as the Perl 5
programming language system itself.
