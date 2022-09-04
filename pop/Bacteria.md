# All genomes of *Archaea*, species by species

## Install `nwr` and initiate local databases

```shell
brew install wang-q/tap/nwr # 0.5.4 or above
brew install sqlite         # 3.34 or above

nwr download
nwr txdb

nwr ardb
nwr ardb --genbank

```

## Strain info

[Bacteria](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2)

### NCBI statistics

```shell
curl -L "https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=statistics&uncultured=hide&unspecified=hide" |
    pup 'table[bgcolor=#CCCCFF] tr td text{}' |
    grep '\S' |
    paste -d $'\t' - - - - - - |
    head -n 9 |
    mlr --itsv --omd cat


```

| Ranks:        | higher taxa |   genus | species | lower taxa |   total |
|---------------|------------:|--------:|--------:|-----------:|--------:|
| Archaea       |         536 |     240 |     830 |          0 |   1,606 |
| Bacteria      |       5,446 |   4,802 |  23,909 |        925 |  35,082 |
| Eukaryota     |      65,206 |  96,518 | 500,097 |     35,417 | 697,238 |
| Fungi         |       5,686 |   7,213 |  53,319 |      1,518 |  67,736 |
| Metazoa       |      47,431 |  68,706 | 262,013 |     17,646 | 395,796 |
| Viridiplantae |       8,205 |  16,756 | 170,584 |     15,888 | 211,433 |
| Viruses       |       1,934 |   2,525 |   7,158 |         65 |  11,682 |
| All taxa      |      73,152 | 104,087 | 531,980 |     36,407 | 745,626 |

### List all ranks

```shell
mkdir -p ~/data/Bacteria
cd ~/data/Bacteria

nwr member Bacteria |
    grep -v " sp." |
    tsv-summarize -H -g 3 --count |
    mlr --itsv --omd cat

```

| rank             | count  |
|------------------|--------|
| superkingdom     | 1      |
| phylum           | 169    |
| class            | 121    |
| order            | 278    |
| family           | 711    |
| no rank          | 6277   |
| species          | 108809 |
| genus            | 4808   |
| clade            | 132    |
| strain           | 40793  |
| varietas         | 24     |
| isolate          | 456    |
| subspecies       | 652    |
| subclass         | 5      |
| forma            | 4      |
| species group    | 92     |
| species subgroup | 28     |
| suborder         | 7      |
| biotype          | 7      |
| serotype         | 252    |
| serogroup        | 138    |
| subphylum        | 1      |
| subgenus         | 1      |
| tribe            | 2      |
| pathogroup       | 5      |
| subfamily        | 1      |

