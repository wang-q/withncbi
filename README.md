# withncbi

withncbi - batch operations according to various NCBI/EBI databases.

## Purpose

Fetch sequences, generate reports and build alignments according to
various NCBI databases.

For more detailed, check `README.md` in each sub-directories.

## Directory organization

* [`db/`](db/): turn NCBI genome reports and assembly reports into a
  query-able MySQL database.

* [`ensembl/`](ensembl/): Ensembl related scripts.

* [`misc/`](misc/): miscellaneous projects.

* [`pop/`](pop/): build alignments on an whole Eukaryotes genus.

* [`taxon/`](taxon/): process (small) genomes according to NCBI Taxonomy.

* [`util/`](util/): miscellaneous utilities.

## Concepts

### IntSpans

[AlignDB::IntSpan](https://github.com/wang-q/AlignDB-IntSpan) 

[jintspan](https://github.com/egateam/jintspan) 

### Positions

[`S288c.txt`](https://github.com/wang-q/App-RL/blob/master/t/S288c.txt) 

```text
I:1-100
I(+):90-150
S288c.I(-):190-200
II:21294-22075
II:23537-24097
```

### Runlists in YAML

[App::RL](https://github.com/wang-q/App-RL) 

[jrunlist](https://github.com/egateam/jrunlist) 

### Blocked fasta files

[App::Fasops](https://github.com/wang-q/App-Fasops) 

### Ranges and links of ranges

[App::Rangeops](https://github.com/wang-q/App-Rangeops) 

[jrange](https://github.com/egateam/jrange) 

## Author

Qiang Wang &lt;wang-q@outlook.com&gt;

## Copyright and license

This software is copyright (c) 2015 by Qiang Wang.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.
