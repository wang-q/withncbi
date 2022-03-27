# *Pseudomonas* HGT

- [*Pseudomonas* HGT](#pseudomonas-hgt)
    * [Software](#software)
    * [Strain info](#strain-info)
        + [List all ranks](#list-all-ranks)
        + [Species with assemblies](#species-with-assemblies)
        + [Outgroups](#outgroups)
    * [Download all assemblies](#download-all-assemblies)
    * [BioSample](#biosample)
    * [Count and group strains](#count-and-group-strains)
    * [NCBI taxonomy](#ncbi-taxonomy)
    * [Raw phylogenetic tree by MinHash](#raw-phylogenetic-tree-by-minhash)
        + [Tweak the mash tree](#tweak-the-mash-tree)
    * [Collect proteins](#collect-proteins)
        + [`all.pro.fa`](#allprofa)
        + [`all.replace.fa`](#allreplacefa)
        + [`all.info.tsv`](#allinfotsv)
    * [Phylogenetics with 40 single-copy genes](#phylogenetics-with-40-single-copy-genes)
        + [Find corresponding proteins by `hmmsearch`](#find-corresponding-proteins-by-hmmsearch)
        + [Create a valid marker gene list](#create-a-valid-marker-gene-list)
        + [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree)
        + [Tweak the concat tree](#tweak-the-concat-tree)
    * [Phylogenetics with bac120](#phylogenetics-with-bac120)
        + [Find corresponding proteins by `hmmsearch`](#find-corresponding-proteins-by-hmmsearch-1)
        + [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree-1)
        + [Tweak the concat tree](#tweak-the-concat-tree-1)
    * [Protein domains and families](#protein-domains-and-families)
        + [Proteins in Pseudomonas strains](#proteins-in-pseudomonas-strains)
        + [Scrap PFAM domains](#scrap-pfam-domains)
        + [Scan every domain](#scan-every-domain)
        + [InterProScan](#interproscan)
    * [Search for gene families with uneven members](#search-for-gene-families-with-uneven-members)
        + [Within species](#within-species)
        + [Among species](#among-species)
        + [IPR005999 - Glycerol kinase](#ipr005999---glycerol-kinase)
        + [IPR004800 - Phosphosugar isomerase, KdsD/KpsF-type](#ipr004800---phosphosugar-isomerase-kdsdkpsf-type)
        + [IPR005593 - Xylulose 5-phosphate/Fructose 6-phosphate phosphoketolase](#ipr005593---xylulose-5-phosphatefructose-6-phosphate-phosphoketolase)
        + [IPR035461 - GmhA/DiaA](#ipr035461---gmhadiaa)
    * [InterProScan on all proteins of typical strains](#interproscan-on-all-proteins-of-typical-strains)
        + [IPR007416 - YggL 50S ribosome-binding protein](#ipr007416---yggl-50s-ribosome-binding-protein)
    * [Collect CDS](#collect-cds)
        + [`all.cds.fa`](#allcdsfa)
        + [`YggL.cds.fa`](#ygglcdsfa)

The genus Pseudomonas includes the conditionally pathogenic bacteria Pseudomonas aeruginosa, plant
pathogens, plant beneficial bacteria, and soil bacteria. Microorganisms of Pseudomonas are extremely
rich in metabolic diversity, and it is thought that this diversity also allows them to survive in a
very wide range of ecological niches.

The metabolism of carbohydrates is a fundamental biochemical process that ensures a continuous
supply of energy to living cells.

The biodegradation of RNA is also an important part of metabolism, which is accomplished by the
degradosomes. However, the diversity of degradosomes in different environments has not been fully
investigated.

According to a recent [paper](https://doi.org/10.1128/mSystems.00543-20), there are some order-level
changes in Gammaproteobacteria. We include both old and new orders.

* Old ones: Cellvibrionales, Oceanospirillales, Pseudomonadales, and Alteromonadales
* New ones: Moraxellales, Kangiellales, and Pseudomonadales

## Software

* Install `nwr` and create a local taxonomy and assembly database.

```shell
brew install wang-q/tap/nwr # 0.5.5 or above
brew install sqlite         # 3.34 or above

nwr download
nwr txdb

nwr ardb
nwr ardb --genbank

```

* Other packages

```shell
brew install hmmer
brew install brewsci/bio/muscle
brew install brewsci/bio/fasttree
brew install brewsci/bio/newick-utils
brew install brewsci/bio/trimal

brew install datamash
brew install miller
brew install wang-q/tap/tsv-utils

brew install librsvg
brew install jq
brew install pup

# KaKs
brew install brewsci/bio/clustal-w
brew install brewsci/bio/paml

cpanm Bio::Tools::Run::Alignment::Clustalw
cpanm https://github.com/wang-q/Bio-Tools-Phylo-PAML.git

```

## Strain info

* [Pseudomonas](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=286)
* [Acinetobacter](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=469)

### List all ranks

```shell
mkdir -p ~/data/Pseudomonas
cd ~/data/Pseudomonas

nwr member Pseudomonas |
    grep -v " sp." |
    tsv-summarize -H -g 3 --count |
    mlr --itsv --omd cat

nwr member Acinetobacter |
    grep -v " sp." |
    tsv-summarize -H -g 3 --count |
    mlr --itsv --omd cat

nwr member Pseudomonas Acinetobacter -r "species group" -r "species subgroup" |
    tsv-select -f 1-3 |
    keep-header -- tsv-sort -k3,3 -k2,2 |
    sed 's/Pseudomonas /P. /g' |
    sed 's/Acinetobacter /A. /g' |
    mlr --itsv --omd cat

```

| rank             | count |
|------------------|-------|
| genus            | 1     |
| species          | 405   |
| strain           | 747   |
| subspecies       | 12    |
| no rank          | 120   |
| species group    | 6     |
| species subgroup | 5     |
| isolate          | 1     |

| rank             | count |
|------------------|-------|
| genus            | 1     |
| species group    | 2     |
| species subgroup | 3     |
| species          | 113   |
| strain           | 1110  |
| no rank          | 2     |
| subspecies       | 1     |
| isolate          | 2     |

| #tax_id | sci_name                                 | rank             |
|---------|------------------------------------------|------------------|
| 909768  | A. calcoaceticus/baumannii complex       | species group    |
| 2839056 | A. Taxon 24                              | species group    |
| 136841  | P. aeruginosa group                      | species group    |
| 136842  | P. chlororaphis group                    | species group    |
| 136843  | P. fluorescens group                     | species group    |
| 136845  | P. putida group                          | species group    |
| 136846  | P. stutzeri group                        | species group    |
| 136849  | P. syringae group                        | species group    |
| 2839060 | A. Taxon 24C                             | species subgroup |
| 2839057 | A. Taxon 24D                             | species subgroup |
| 2839061 | A. Taxon 24E                             | species subgroup |
| 627141  | P. nitroreducens/multiresinivorans group | species subgroup |
| 1232139 | P. oleovorans/pseudoalcaligenes group    | species subgroup |
| 578833  | P. stutzeri subgroup                     | species subgroup |
| 251695  | P. syringae group genomosp. 1            | species subgroup |
| 251698  | P. syringae group genomosp. 2            | species subgroup |

### Species with assemblies

Check also the order Pseudomonadales.

* Old ones: Cellvibrionales, Oceanospirillales, Pseudomonadales, and Alteromonadales
* New ones: Moraxellales, Kangiellales, and Pseudomonadales

```shell
cd ~/data/Pseudomonas

SPECIES=$(
    nwr member -r species \
        Cellvibrionales Oceanospirillales Alteromonadales \
        Moraxellales Kangiellales Pseudomonadales |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        grep -v " sp." |
        grep -v -E "\bbacterium\b" |
        grep -v -E "\bsymbiont\b" |
        sed '1d' |
        cut -f 1 |
        sort |
        uniq
)

for S in $SPECIES; do
    RS=$(
        echo "
            SELECT
                COUNT(*)
            FROM ar
            WHERE 1=1
                AND species_id = $S
            " |
            sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
    )

    CHR=$(
        echo "
            SELECT
                COUNT(*)
            FROM ar
            WHERE 1=1
                AND species_id = $S
                AND assembly_level IN ('Complete Genome', 'Chromosome')
            " |
            sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
    )

    if [[ ${RS} -gt 0 ]]; then
        echo -e "$S\t$RS\t$CHR"
    fi
done |
    nwr append stdin |
    tsv-select -f 1,4,2-3 |
    tsv-sort -k3,3nr -k4,4nr -k2,2 |
    (echo -e 'species_id\tspecies\tRS\tCHR' && cat) \
    > species.count.tsv

cat species.count.tsv |
    tsv-filter -H --ge CHR:5 |
    tsv-filter -H --invert --str-in-fld species:Pseudomonas --lt RS:30 |
    tsv-filter -H --invert --str-in-fld species:Acinetobacter --lt RS:30 |
    sed 's/Pseudomonas /P. /g' |
    sed 's/Acinetobacter /A. /g' |
    mlr --itsv --omd cat

```

| species_id | species                         | RS   | CHR |
|------------|---------------------------------|------|-----|
| 287        | P. aeruginosa                   | 6417 | 477 |
| 470        | A. baumannii                    | 5955 | 335 |
| 33069      | P. viridiflava                  | 1536 | 7   |
| 317        | P. syringae                     | 540  | 41  |
| 48296      | A. pittii                       | 322  | 32  |
| 294        | P. fluorescens                  | 257  | 39  |
| 480        | Moraxella catarrhalis           | 210  | 16  |
| 303        | P. putida                       | 190  | 50  |
| 106654     | A. nosocomialis                 | 160  | 11  |
| 316        | P. stutzeri                     | 133  | 30  |
| 29438      | P. savastanoi                   | 116  | 5   |
| 38313      | Shewanella algae                | 108  | 22  |
| 587753     | P. chlororaphis                 | 99   | 60  |
| 756892     | A. indicus                      | 91   | 19  |
| 47877      | P. amygdali                     | 86   | 8   |
| 380021     | P. protegens                    | 73   | 23  |
| 40215      | A. junii                        | 65   | 9   |
| 1530123    | A. seifertii                    | 59   | 25  |
| 29430      | A. haemolyticus                 | 55   | 14  |
| 76759      | P. monteilii                    | 46   | 9   |
| 40214      | A. johnsonii                    | 43   | 19  |
| 296        | P. fragi                        | 38   | 6   |
| 28090      | A. lwoffii                      | 31   | 11  |
| 43657      | Pseudoalteromonas luteoviolacea | 25   | 5   |
| 28108      | Alteromonas macleodii           | 24   | 9   |
| 34062      | Moraxella osloensis             | 21   | 10  |
| 43662      | Pseudoalteromonas piscicida     | 19   | 6   |
| 314275     | Alteromonas mediterranea        | 16   | 16  |
| 24         | Shewanella putrefaciens         | 15   | 9   |
| 62322      | Shewanella baltica              | 14   | 11  |
| 386891     | Moraxella bovoculi              | 9    | 7   |
| 1697053    | Thiopseudomonas alkaliphila     | 7    | 7   |

### Outgroups

Use these model organisms as outgroups.

```shell
cd ~/data/Pseudomonas

GENUS=$(
    nwr member Bacteria -r genus |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        sed '1d' |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
.headers ON

    SELECT
        *
    FROM ar
    WHERE 1=1
        AND genus_id IN ($GENUS)
        AND refseq_category IN ('reference genome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    > reference.tsv

cat reference.tsv |
    sed '1s/^/#/' |
    nwr append stdin -r phylum -r class |
    tsv-select -H -f 1,2,phylum,class |
    parallel --col-sep "\t" -j 1 '
        if [[ "{3}" != "Proteobacteria" ]]; then
            printf "%s\t%s\t%s\n" {1} {2} {3}
        else
            printf "%s\t%s\t%s\n" {1} {2} {4}
        fi
    ' |
    mlr --itsv --omd cat

```

| #tax_id | organism_name                                                    | phylum                |
|---------|------------------------------------------------------------------|-----------------------|
| 565050  | Caulobacter vibrioides NA1000                                    | Alphaproteobacteria   |
| 192222  | Campylobacter jejuni subsp. jejuni NCTC 11168 = ATCC 700819      | Epsilonproteobacteria |
| 208964  | Pseudomonas aeruginosa PAO1                                      | Gammaproteobacteria   |
| 871585  | Acinetobacter pittii PHEA-2                                      | Gammaproteobacteria   |
| 511145  | Escherichia coli str. K-12 substr. MG1655                        | Gammaproteobacteria   |
| 386585  | Escherichia coli O157:H7 str. Sakai                              | Gammaproteobacteria   |
| 1125630 | Klebsiella pneumoniae subsp. pneumoniae HS11286                  | Gammaproteobacteria   |
| 99287   | Salmonella enterica subsp. enterica serovar Typhimurium str. LT2 | Gammaproteobacteria   |
| 198214  | Shigella flexneri 2a str. 301                                    | Gammaproteobacteria   |
| 227377  | Coxiella burnetii RSA 493                                        | Gammaproteobacteria   |
| 272561  | Chlamydia trachomatis D/UW-3/CX                                  | Chlamydiae            |
| 93061   | Staphylococcus aureus subsp. aureus NCTC 8325                    | Firmicutes            |
| 224308  | Bacillus subtilis subsp. subtilis str. 168                       | Firmicutes            |
| 169963  | Listeria monocytogenes EGD-e                                     | Firmicutes            |
| 83332   | Mycobacterium tuberculosis H37Rv                                 | Actinobacteria        |

## Download all assemblies

```shell
cd ~/data/Pseudomonas

# Pseudomonas aeruginosa PAO1 is in the reference list
cat reference.tsv |
    tsv-select -H -f organism_name,species,genus,ftp_path,assembly_level \
    > raw.tsv

# Species with 2 or more genomes were retained.
SPECIES=$(
    cat species.count.tsv |
        tsv-filter -H --ge CHR:2 |
        sed '1d' |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

# Pseudomonadales
echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND assembly_level IN ('Complete Genome', 'Chromosome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    tsv-filter --invert --str-eq 2:"Pseudomonas aeruginosa" --str-eq 5:"Chromosome" |
    tsv-filter --invert --str-eq 2:"Acinetobacter baumannii" --str-eq 5:"Chromosome" \
    >> raw.tsv

# Also includes representative strains of Gammaproteobacteria.
# families not in our orders
FAMILY=$(
    nwr member Gammaproteobacteria -r family |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        sed '1d' |
        cut -f 1 |
        nwr append stdin -r order |
        tsv-filter --str-ne 2:"Cellvibrionales" --str-ne 2:"Oceanospirillales" --str-ne 2:"Alteromonadales" |
        tsv-filter --str-ne 2:"Moraxellales" --str-ne 2:"Kangiellales" --str-ne 2:"Pseudomonadales" |
        tsv-filter --str-ne 2:"NA" |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND family_id IN ($FAMILY)
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND refseq_category IN ('representative genome')
        AND assembly_level IN ('Complete Genome', 'Chromosome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    grep -v -i "symbiont " |
    tsv-filter --str-not-in-fld 1:"[" \
    >> raw.tsv

# Create abbr.
cat raw.tsv |
    grep -v '^#' |
    tsv-uniq |
    perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "1,2,3" -s '\t' -m 3 --shortsub |
    (echo -e '#name\tftp_path\torganism\tassembly_level' && cat ) |
    perl -nl -a -F"," -e '
        BEGIN{my %seen};
        /^#/ and print and next;
        /^organism_name/i and next;
        $seen{$F[3]}++; # ftp_path
        $seen{$F[3]} > 1 and next;
        $seen{$F[5]}++; # abbr_name
        $seen{$F[5]} > 1 and next;
        printf qq{%s\t%s\t%s\t%s\n}, $F[5], $F[3], $F[1], $F[4];
        ' |
    keep-header -- sort -k3,3 -k1,1 \
    > Pseudomonas.assembly.tsv

# find potential duplicate strains or assemblies
cat Pseudomonas.assembly.tsv |
    tsv-uniq -f 1 --repeated

# Edit .tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Pseudomonas.assembly.tsv
# cp Pseudomonas.assembly.tsv ~/Scripts/withncbi/pop

# Comment out unneeded strains

# Cleaning
rm raw*.*sv

```

```shell
cd ~/data/Pseudomonas

perl ~/Scripts/withncbi/taxon/assembly_prep.pl \
    -f ~/Scripts/withncbi/pop/Pseudomonas.assembly.tsv \
    -o ASSEMBLY

# Remove dirs not in the list
find ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
    tr "/" "\t" |
    cut -f 2 |
    tsv-join --exclude -k 1 -f ASSEMBLY/rsync.tsv -d 1 |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        echo Remove {}
        rm -fr ASSEMBLY/{}
    '

# Run
proxychains4 bash ASSEMBLY/Pseudomonas.assembly.rsync.sh

bash ASSEMBLY/Pseudomonas.assembly.collect.sh

# md5
cat ASSEMBLY/rsync.tsv |
    tsv-select -f 1 |
    parallel -j 4 --keep-order '
        echo "==> {}"
        cd ASSEMBLY/{}
        md5sum --check md5checksums.txt
    ' |
    grep -v ": OK"

```

## BioSample

ENA's BioSample missed many strains, so NCBI's was used.

```shell
cd ~/data/Pseudomonas

mkdir -p biosample

ulimit -n `ulimit -Hn`

cat ASSEMBLY/Pseudomonas.assembly.collect.csv |
    tsv-select -H -d, -f BioSample |
    grep "^SAM" |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        if [ ! -s biosample/{}.txt ]; then
            >&2 echo {}
            curl -fsSL "https://www.ncbi.nlm.nih.gov/biosample/?term={}&report=full&format=text" -o biosample/{}.txt
#            curl -fsSL "https://www.ebi.ac.uk/biosamples/samples/{}" -o biosample/{}.json
        fi
    '

find biosample -name "SAM*.txt" | wc -l
# 1956

find biosample -name "SAM*.txt" |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        cat {} |
            perl -nl -e '\''
                print $1 if m{\s+\/([\w_ ]+)=};
            '\''
    ' |
    tsv-uniq --at-least 50 | # ignore rare attributes
    grep -v "^INSDC" |
    grep -v "^ENA" \
    > attributes.lst

cat attributes.lst |
    (echo -e "BioSample" && cat) |
    tr '\n' '\t' |
    sed 's/\t$/\n/' \
    > Pseudomonas.biosample.tsv

find biosample -name "SAM*.txt" |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        >&2 echo {/.}
        cat {} |
            perl -nl -MPath::Tiny -e '\''
                BEGIN {
                    our @keys = grep {/\S/} path(q{attributes.lst})->lines({chomp => 1});
                    our %stat = ();
                }

                m(\s+\/([\w_ ]+)=\"(.+)\") or next;
                my $k = $1;
                my $v = $2;
                if ( $v =~ m(\bNA|missing|Not applicable|not collected|not available|not provided|N\/A|not known|unknown\b)i ) {
                    $stat{$k} = q();
                } else {
                    $stat{$k} = $v;
                }

                END {
                    my @c;
                    for my $key ( @keys ) {
                        if (exists $stat{$key}) {
                            push @c, $stat{$key};
                        }
                        else {
                            push @c, q();
                        }
                    }
                    print join(qq{\t}, q{{/.}}, @c);
                }
            '\''
    ' \
    >> Pseudomonas.biosample.tsv

```

## Count and group strains

* Check N50 of assemblies

* Some strains were anomalously labeled and identified by the `mash` tree.
    * Pseudom_flu_GCF_900636635_1
    * Pseudom_chl_GCF_001023535_1
    * Pseudom_syr_GCF_004006335_1
    * Pseudom_puti_GCF_003228315_1 and Pseudom_puti_GCF_020172705_1

```shell
cd ~/data/Pseudomonas

for dir in $(find ASSEMBLY -maxdepth 1 -mindepth 1 -type d | sort); do
    1>&2 echo "==> ${dir}"
    name=$(basename ${dir})

    find ${dir} -type f -name "*_genomic.fna.gz" |
        grep -v "_from_" | # exclude CDS and rna
        xargs cat |
        faops n50 -C -S stdin |
        (echo -e "name\t${name}" && cat) |
        datamash transpose
done |
    tsv-uniq |
    tee ASSEMBLY/n50.tsv

cat ASSEMBLY/n50.tsv |
    tsv-filter \
        -H --or \
        --le 4:100 \
        --ge 2:100000 |
    tsv-filter -H --ge 3:1000000 |
    tr "\t" "," \
    > ASSEMBLY/n50.pass.csv

wc -l ASSEMBLY/n50*
#  1958 ASSEMBLY/n50.pass.csv
#  1959 ASSEMBLY/n50.tsv

tsv-join \
    ASSEMBLY/Pseudomonas.assembly.collect.csv \
    --delimiter "," -H --key-fields 1 \
    --filter-file ASSEMBLY/n50.pass.csv |
    tsv-filter --delimiter "," -H --str-not-in-fld 1:GCF_900636635 |
    tsv-filter --delimiter "," -H --str-not-in-fld 1:GCF_001023535 |
    tsv-filter --delimiter "," -H --str-not-in-fld 1:GCF_004006335 |
    tsv-filter --delimiter "," -H --str-not-in-fld 1:GCF_003228315 |
    tsv-filter --delimiter "," -H --str-not-in-fld 1:GCF_020172705 \
    > ASSEMBLY/Pseudomonas.assembly.pass.csv

wc -l ASSEMBLY/Pseudomonas.assembly*csv
#  1959 ASSEMBLY/Pseudomonas.assembly.collect.csv
#  1953 ASSEMBLY/Pseudomonas.assembly.pass.csv

```

* Order

```shell
cd ~/data/Pseudomonas

# Group by order
cat ASSEMBLY/Pseudomonas.assembly.pass.csv |
    sed -e '1d' |
    tsv-select -d, -f 3 |
    tsv-uniq |
    nwr append stdin -r order |
    tsv-select -f 2 |
    tsv-uniq \
    > order.lst

cat order.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        n_species=$(cat ASSEMBLY/Pseudomonas.assembly.pass.csv |
            sed "1d" |
            tsv-select -d, -f 3 |
            nwr append stdin -r order -r species |
            grep {} |
            tsv-select -f 1,3 |
            tsv-uniq |
            wc -l)

        n_strains=$(cat ASSEMBLY/Pseudomonas.assembly.pass.csv |
            sed "1d" |
            tsv-select -d, -f 3 |
            nwr append stdin -r order |
            grep {} |
            wc -l)

        printf "%s\t%d\t%d\n" {} ${n_species} ${n_strains}
    ' |
    nwr append stdin --id |
    tsv-select -f 5,4,2,3 |
    tsv-sort -k2,2 |
    (echo -e '#tax_id\tfamily\t#species\t#strains' && cat) |
    mlr --itsv --omd cat

```

| #tax_id | family                | #species | #strains |
|---------|-----------------------|----------|----------|
| 1692040 | Acidiferrobacterales  | 3        | 3        |
| 135624  | Aeromonadales         | 18       | 18       |
| 135622  | Alteromonadales       | 52       | 118      |
| 1385    | Bacillales            | 3        | 3        |
| 213849  | Campylobacterales     | 1        | 1        |
| 135615  | Cardiobacteriales     | 2        | 2        |
| 204458  | Caulobacterales       | 1        | 1        |
| 1706369 | Cellvibrionales       | 2        | 4        |
| 51291   | Chlamydiales          | 1        | 1        |
| 135613  | Chromatiales          | 32       | 32       |
| 85007   | Corynebacteriales     | 1        | 1        |
| 91347   | Enterobacterales      | 179      | 179      |
| 1934945 | Immundisolibacterales | 1        | 1        |
| 118969  | Legionellales         | 18       | 18       |
| 135618  | Methylococcales       | 14       | 14       |
| 2887326 | Moraxellales          | 50       | 523      |
| 1775403 | Nevskiales            | 1        | 1        |
| 135619  | Oceanospirillales     | 9        | 18       |
| 1240482 | Orbales               | 3        | 3        |
| 135625  | Pasteurellales        | 34       | 34       |
| 72274   | Pseudomonadales       | 183      | 846      |
| 72273   | Thiotrichales         | 32       | 32       |
| 135623  | Vibrionales           | 46       | 46       |
| 135614  | Xanthomonadales       | 53       | 53       |

* strains

```shell
cd ~/data/Pseudomonas

# list strains
mkdir -p taxon

rm taxon/* strains.lst *.tmp
cat ASSEMBLY/Pseudomonas.assembly.pass.csv |
    sed -e '1d' |
    tr "," "\t" |
    tsv-select -f 1,2,3 |
    nwr append stdin -c 3 -r species -r genus -r family -r order |
    parallel --col-sep "\t" --no-run-if-empty --linebuffer -k -j 4 '
        echo {1} >> strains.lst

        echo {4} >> species.tmp
        echo {5} >> genus.tmp
        echo {6} >> family.tmp

        echo {7} >> order.tmp
        echo {1} >> taxon/{7}

        printf "%s\t%s\t%d\t%s\t%s\t%s\t%s\n" {1} {2} {3} {4} {5} {6} {7}
    ' \
    > strains.taxon.tsv

cat species.tmp | tsv-uniq > species.lst
cat genus.tmp | tsv-uniq > genus.lst
cat family.tmp | tsv-uniq > family.lst
cat order.tmp | tsv-uniq > order.lst

# Omit strains without protein annotations
for STRAIN in $(cat strains.lst); do
    if ! compgen -G "ASSEMBLY/${STRAIN}/*_protein.faa.gz" > /dev/null; then
        echo ${STRAIN}
    fi
    if ! compgen -G "ASSEMBLY/${STRAIN}/*_cds_from_genomic.fna.gz" > /dev/null; then
        echo ${STRAIN}
    fi
done |
    tsv-uniq \
    > omit.lst
# All OK

rm *.tmp

```

## NCBI taxonomy

Done by `bp_taxonomy2tree.pl` from BioPerl.

```shell
mkdir -p ~/data/Pseudomonas/tree
cd ~/data/Pseudomonas/tree

bp_taxonomy2tree.pl -e \
    $(
        cat ../genus.lst |
            tr " " "_" |
            parallel echo '-s {}'
    ) \
    > ncbi.nwk

nw_display -s -b 'visibility:hidden' -w 600 -v 30 ncbi.nwk |
    rsvg-convert -o Pseudomonas.ncbi.png

```

## Raw phylogenetic tree by MinHash

```shell
mkdir -p ~/data/Pseudomonas/mash
cd ~/data/Pseudomonas/mash

for strain in $(cat ../strains.lst ); do
    2>&1 echo "==> ${strain}"

    if [[ -e ${strain}.msh ]]; then
        continue
    fi

    find ../ASSEMBLY/${strain} -name "*_genomic.fna.gz" |
        grep -v "_from_" |
        xargs cat |
        mash sketch -k 21 -s 100000 -p 8 - -I "${strain}" -o ${strain}
done

mash triangle -E -p 8 -l <(
    cat ../strains.lst | parallel echo "{}.msh"
    ) \
    > dist.tsv

# fill matrix with lower triangle
tsv-select -f 1-3 dist.tsv |
    (tsv-select -f 2,1,3 dist.tsv && cat) |
    (
        cut -f 1 dist.tsv |
            tsv-uniq |
            parallel -j 1 --keep-order 'echo -e "{}\t{}\t0"' &&
        cat
    ) \
    > dist_full.tsv

cat dist_full.tsv |
    Rscript -e '
        library(readr);
        library(tidyr);
        library(ape);
        pair_dist <- read_tsv(file("stdin"), col_names=F);
        tmp <- pair_dist %>%
            pivot_wider( names_from = X2, values_from = X3, values_fill = list(X3 = 1.0) )
        tmp <- as.matrix(tmp)
        mat <- tmp[,-1]
        rownames(mat) <- tmp[,1]

        dist_mat <- as.dist(mat)
        clusters <- hclust(dist_mat, method = "ward.D2")
        tree <- as.phylo(clusters)
        write.tree(phy=tree, file="tree.nwk")

        group <- cutree(clusters, h=0.4) # k=5
        groups <- as.data.frame(group)
        groups$ids <- rownames(groups)
        rownames(groups) <- NULL
        groups <- groups[order(groups$group), ]
        write_tsv(groups, "groups.tsv")
    '

```

### Tweak the mash tree

```shell
mkdir -p ~/data/Pseudomonas/tree
cd ~/data/Pseudomonas/tree

nw_reroot ../mash/tree.nwk Bac_subti_subtilis_168 Sta_aure_aureus_NCTC_8325 |
    nw_order -c n - \
    > mash.reroot.newick

# rank::col
ARRAY=(
    'order::7'
    'family::6'
    'genus::5'
    'species::4'
)

rm mash.condensed.map
CUR_TREE=mash.reroot.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/Scripts/withncbi/taxon/condense_tree.sh ${CUR_TREE} ../strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick mash.${GROUP_NAME}.newick
    cat condense.map >> mash.condensed.map

    CUR_TREE=mash.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 600 -v 30 mash.species.newick |
    rsvg-convert -o Pseudomonas.mash.png

```

## Collect proteins

### `all.pro.fa`

```shell script
cd ~/data/Pseudomonas

mkdir -p PROTEINS

find ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
    sort |
    grep 'ASSEMBLY/' |
    wc -l
# 1958

find ASSEMBLY -type f -name "*_protein.faa.gz" |
    wc -l
# 1958

cat strains.lst |
    wc -l
# 1952

for STRAIN in $(cat strains.lst); do
    gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz
done \
    > PROTEINS/all.pro.fa

cat PROTEINS/all.pro.fa |
    perl -nl -e '
        BEGIN { our %seen; our $h; }

        if (/^>/) {
            $h = (split(" ", $_))[0];
            $seen{$h}++;
            $_ = $h;
        }
        print if $seen{$h} == 1;
    ' \
    > PROTEINS/all.uniq.fa

# counting proteins
cat PROTEINS/all.pro.fa |
    grep "^>" |
    wc -l
#8754303

cat PROTEINS/all.pro.fa |
    grep "^>" |
    tsv-uniq |
    wc -l
#3995205

# annotations may be different
cat PROTEINS/all.uniq.fa |
    grep "^>" |
    wc -l
#3944568

# ribonuclease
cat PROTEINS/all.pro.fa |
    grep "ribonuclease" |
    grep -v "deoxyribonuclease" |
    perl -nl -e 's/^>\w+\.\d+\s+//g; print' |
    perl -nl -e 's/\s+\[.+?\]$//g; print' |
    perl -nl -e 's/MULTISPECIES: //g; print' |
    sort |
    uniq -c |
    sort -nr

```

### `all.replace.fa`

```shell
cd ~/data/Pseudomonas

rm PROTEINS/all.strain.tsv
for STRAIN in $(cat strains.lst); do
    gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
        grep "^>" |
        cut -d" " -f 1 |
        sed "s/^>//" |
        STRAIN=${STRAIN} perl -nl -e '
            $n = $_;
            $s = $n;
            $s =~ s/\.\d+//;
            printf qq{%s\t%s_%s\t%s\n}, $n, $ENV{STRAIN}, $s, $ENV{STRAIN};
        ' \
    > PROTEINS/${STRAIN}.replace.tsv

    cut -f 2,3 PROTEINS/${STRAIN}.replace.tsv >> PROTEINS/all.strain.tsv

    faops replace -s ASSEMBLY/${STRAIN}/*_protein.faa.gz <(cut -f 1,2 PROTEINS/${STRAIN}.replace.tsv) stdout

    rm PROTEINS/${STRAIN}.replace.tsv
done \
    > PROTEINS/all.replace.fa

cat PROTEINS/all.replace.fa |
    grep "^>" |
    wc -l
#8754303

(echo -e "#name\tstrain" && cat PROTEINS/all.strain.tsv)  \
    > temp &&
    mv temp PROTEINS/all.strain.tsv

faops size PROTEINS/all.replace.fa > PROTEINS/all.replace.sizes

(echo -e "#name\tsize" && cat PROTEINS/all.replace.sizes) > PROTEINS/all.size.tsv

rm PROTEINS/all.replace.sizes

```

### `all.info.tsv`

```shell
cd ~/data/Pseudomonas

for STRAIN in $(cat strains.lst); do
    gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
        grep "^>" |
        sed "s/^>//" |
        perl -nl -e '/\[.+\[/ and s/\[/\(/; print' |
        perl -nl -e '/\].+\]/ and s/\]/\)/; print' |
        perl -nl -e 's/\s+\[.+?\]$//g; print' |
        perl -nl -e 's/MULTISPECIES: //g; print' |
        STRAIN=${STRAIN} perl -nl -e '
            /^(\w+)\.\d+\s+(.+)$/ or next;
            printf qq{%s_%s\t%s\n}, $ENV{STRAIN}, $1, $2;
        '
done \
    > PROTEINS/all.annotation.tsv

cat PROTEINS/all.annotation.tsv |
    wc -l
#8754303

(echo -e "#name\tannotation" && cat PROTEINS/all.annotation.tsv) \
    > temp &&
    mv temp PROTEINS/all.annotation.tsv

# check differences
cat PROTEINS/all.size.tsv |
    grep -F -f <(cut -f 1 PROTEINS/all.annotation.tsv) -v

tsv-join \
    PROTEINS/all.strain.tsv \
    --data-fields 1 \
    -f PROTEINS/all.size.tsv \
    --key-fields 1 \
    --append-fields 2 \
    > PROTEINS/all.strain_size.tsv

tsv-join \
    PROTEINS/all.strain_size.tsv \
    --data-fields 1 \
    -f PROTEINS/all.annotation.tsv \
    --key-fields 1 \
    --append-fields 2 \
    > PROTEINS/all.info.tsv

cat PROTEINS/all.info.tsv |
    wc -l
#8754304

```

## Phylogenetics with 40 single-copy genes

### Find corresponding proteins by `hmmsearch`

* Download HMM models as described in [`hmm/README.md`](../hmm/README.md)

* The `E_VALUE` was manually adjusted to 1e-20 to reach a balance between sensitivity and
  speciality.

```shell
E_VALUE=1e-20

cd ~/data/Pseudomonas

## example
#gzip -dcf ASSEMBLY/Ac_axa_ATCC_25176/*_protein.faa.gz |
#    hmmsearch -E 1e-20 --domE 1e-20 --noali --notextw ~/data/HMM/scg40/bacteria_and_archaea_dir/BA00001.hmm - |
#    grep '>>' |
#    perl -nl -e '/>>\s+(\S+)/ and print $1'

# Find all genes
for marker in BA000{01..40}; do
    >&2 echo "==> marker [${marker}]"

    mkdir -p PROTEINS/${marker}

    for ORDER in $(cat order.lst); do
        >&2 echo "==> ORDER [${ORDER}]"

        cat taxon/${ORDER} |
            parallel --no-run-if-empty --linebuffer -k -j 8 "
                gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                    hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw ~/data/HMM/scg40/bacteria_and_archaea_dir/${marker}.hmm - |
                    grep '>>' |
                    perl -nl -e ' m{>>\s+(\S+)} and printf qq{%s\t%s\n}, \$1, {}; '
            " \
            > PROTEINS/${marker}/${ORDER}.replace.tsv
    done

    echo
done

```

### Create a valid marker gene list

* `hmmsearch` may identify more than one copy for some marker genes.

    * BA00004: translation initiation factor EF-2
    * BA00005: translation initiation factor IF-2
    * BA00008: signal recognition particle protein

* Acinetobacter
    * BA00028

Compare proteins and strains.

```shell script
cd ~/data/Pseudomonas

for marker in BA000{01..03} BA000{06..07} BA000{09..27} BA000{29..40}; do
    echo ${marker}
done > marker.lst

for marker in $(cat marker.lst); do
    echo "==> marker [${marker}]"

    for ORDER in $(cat order.lst); do
        cat PROTEINS/${marker}/${ORDER}.replace.tsv |
            cut -f 2 |
            diff - taxon/${ORDER}
    done

    echo
done

```

### Align and concat marker genes to create species tree

* Strains within a species share a large proportion of identical protein sequences.

* Use `trimal -automated1` to remove gaps.

```shell script
cd ~/data/Pseudomonas

# Extract sequences
cat marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        >&2 echo "==> marker [{}]"

        for ORDER in $(cat order.lst); do
            cat PROTEINS/{}/${ORDER}.replace.tsv
        done \
            > PROTEINS/{}/{}.replace.tsv

        faops some PROTEINS/all.uniq.fa <(
            cat PROTEINS/{}/{}.replace.tsv |
                cut -f 1 |
                tsv-uniq
            ) stdout \
            > PROTEINS/{}/{}.pro.fa
    '

# Align each markers with `muscle`
cat marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        >&2 echo "==> marker [{}]"

        muscle -quiet -in PROTEINS/{}/{}.pro.fa -out PROTEINS/{}/{}.aln.fa
    '

for marker in $(cat marker.lst); do
    >&2 echo "==> marker [${marker}]"

    # 1 name to many names
    cat PROTEINS/${marker}/${marker}.replace.tsv |
        parallel --no-run-if-empty --linebuffer -k -j 4 "
            faops replace -s PROTEINS/${marker}/${marker}.aln.fa <(echo {}) stdout
        " \
        > PROTEINS/${marker}/${marker}.replace.fa
done

# Concat marker genes
for marker in $(cat marker.lst); do
    # sequences in one line
    faops filter -l 0 PROTEINS/${marker}/${marker}.replace.fa stdout

    # empty line for .fas
    echo
done \
    > PROTEINS/scg40.aln.fas

fasops concat PROTEINS/scg40.aln.fas strains.lst -o PROTEINS/scg40.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in PROTEINS/scg40.aln.fa -out PROTEINS/scg40.trim.fa -automated1

# FastTree produces NJ trees to simulate ML ones
FastTree PROTEINS/scg40.trim.fa > PROTEINS/scg40.trim.newick

```

### Tweak the concat tree

```shell script
cd ~/data/Pseudomonas/tree

nw_reroot ../PROTEINS/scg40.trim.newick Bac_subti_subtilis_168 Sta_aure_aureus_NCTC_8325 |
    nw_order -c n - \
    > scg40.reroot.newick

# rank::col
ARRAY=(
#    'order::7'
#    'family::6'
    'genus::5'
    'species::4'
)

rm scg40.condensed.map
CUR_TREE=scg40.reroot.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/Scripts/withncbi/taxon/condense_tree.sh ${CUR_TREE} ../strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick scg40.${GROUP_NAME}.newick
    cat condense.map >> scg40.condensed.map

    CUR_TREE=scg40.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 600 -v 30 scg40.species.newick |
    rsvg-convert -o Pseudomonas.scg40.png

```

## Phylogenetics with bac120

### Find corresponding proteins by `hmmsearch`

```shell
E_VALUE=1e-20

cd ~/data/Pseudomonas

# Find all genes
for marker in $(cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1); do
    >&2 echo "==> marker [${marker}]"

    mkdir -p PROTEINS/${marker}

    for ORDER in $(cat order.lst); do
        >&2 echo "==> ORDER [${ORDER}]"

        cat taxon/${ORDER} |
            parallel --no-run-if-empty --linebuffer -k -j 8 "
                gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                    hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw ~/data/HMM/bac120/HMM/${marker}.HMM - |
                    grep '>>' |
                    perl -nl -e ' m{>>\s+(\S+)} and printf qq{%s\t%s\n}, \$1, {}; '
            " \
            > PROTEINS/${marker}/${ORDER}.replace.tsv
    done

    echo
done

```

### Align and concat marker genes to create species tree

```shell script
cd ~/data/Pseudomonas

cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        cat PROTEINS/{}/*.replace.tsv | wc -l
    ' |
    tsv-summarize --median 1
# 1951

cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo {}
        cat PROTEINS/{}/*.replace.tsv | wc -l
    ' |
    paste - - |
    tsv-filter --invert --ge 2:1500 --le 2:2500 |
    cut -f 1 \
    > PROTEINS/bac120.omit.lst

# Extract sequences
# Multiple copies slow down the alignment process
cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    grep -v -Fx -f PROTEINS/bac120.omit.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        >&2 echo "==> marker [{}]"

        for ORDER in $(cat order.lst); do
            cat PROTEINS/{}/${ORDER}.replace.tsv
        done \
            > PROTEINS/{}/{}.replace.tsv

        faops some PROTEINS/all.uniq.fa <(
            cat PROTEINS/{}/{}.replace.tsv |
                cut -f 1 |
                tsv-uniq
            ) stdout \
            > PROTEINS/{}/{}.pro.fa
    '

# Align each markers with muscle
cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        >&2 echo "==> marker [{}]"
        if [ ! -s PROTEINS/{}/{}.pro.fa ]; then
            exit
        fi
        if [ -s PROTEINS/{}/{}.aln.fa ]; then
            exit
        fi

        muscle -quiet -in PROTEINS/{}/{}.pro.fa -out PROTEINS/{}/{}.aln.fa
    '

for marker in $(cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1); do
    >&2 echo "==> marker [${marker}]"
    if [ ! -s PROTEINS/${marker}/${marker}.pro.fa ]; then
        continue
    fi

    # 1 name to many names
    cat PROTEINS/${marker}/${marker}.replace.tsv |
        parallel --no-run-if-empty --linebuffer -k -j 4 "
            faops replace -s PROTEINS/${marker}/${marker}.aln.fa <(echo {}) stdout
        " \
        > PROTEINS/${marker}/${marker}.replace.fa
done

# Concat marker genes
for marker in $(cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1); do
    if [ ! -s PROTEINS/${marker}/${marker}.pro.fa ]; then
        continue
    fi

    # sequences in one line
    faops filter -l 0 PROTEINS/${marker}/${marker}.replace.fa stdout

    # empty line for .fas
    echo
done \
    > PROTEINS/bac120.aln.fas

fasops concat PROTEINS/bac120.aln.fas strains.lst -o PROTEINS/bac120.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in PROTEINS/bac120.aln.fa -out PROTEINS/bac120.trim.fa -automated1

faops size PROTEINS/bac120.*.fa |
    tsv-uniq -f 2 |
    cut -f 2
#52756
#25526

# To make it faster
FastTree -fastest -noml PROTEINS/bac120.trim.fa > PROTEINS/bac120.trim.newick

```

### Tweak the concat tree

```shell script
cd ~/data/Pseudomonas/tree

nw_reroot ../PROTEINS/bac120.trim.newick Bac_subti_subtilis_168 Sta_aure_aureus_NCTC_8325 |
    nw_order -c n - \
    > bac120.reroot.newick

rm bac120.condensed.map

# rank::col
ARRAY=(
#    'order::7'
#    'family::6'
    'genus::5'
    'species::4'
)

rm bac120.condensed.map
CUR_TREE=bac120.reroot.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/Scripts/withncbi/taxon/condense_tree.sh ${CUR_TREE} ../strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick bac120.${GROUP_NAME}.newick
    cat condense.map >> bac120.condensed.map

    CUR_TREE=bac120.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 600 -v 30 bac120.species.newick |
    rsvg-convert -o Pseudomonas.bac120.png

```

## Protein domains and families

### Proteins in Pseudomonas strains

* [`GO:0005975` carbohydrate metabolic process](https://www.ebi.ac.uk/QuickGO/GTerm?id=GO:0005975)
    * http://www.cazy.org/Glycoside-Hydrolases.html

* [`GO:0016837` carbon-oxygen lyase activity, acting on polysaccharides](https://www.ebi.ac.uk/QuickGO/GTerm?id=GO:0016837)
    * http://www.cazy.org/Polysaccharide-Lyases.html

* Search `https://www.pseudomonas.com/goterms` for `GO:0005975`
    * 107 - Pseudomonas aeruginosa PAO1
    * 7360 - Pseudomonas putida KT2440
    * 479 - Pseudomonas chlororaphis subsp. aureofaciens 30-84
    * 116 - Pseudomonas fluorescens SBW25
    * 113 - Pseudomonas protegens Pf-5
    * 123 - Pseudomonas stutzeri A1501
    * 112 - Pseudomonas syringae pv. syringae B728a
    * 114 - Pseudomonas savastanoi pv. phaseolicola 1448A
    * 117 - Pseudomonas entomophila L48
    * 109 - Pseudomonas aeruginosa UCBPP-PA14

```shell
cd ~/data/Pseudomonas/
mkdir -p ~/data/Pseudomonas/DOMAINS/ref_strain

for ID in 107 7360 479 116 113 123 112 114 117 109 ; do
    URL=$(printf 'https://www.pseudomonas.com/goterms/list?accession=GO:0005975&strain_id=%d&format=TAB' $ID)
    curl -L ${URL}
done  |
    sed '1d' |
    tsv-select -f 1 \
    > DOMAINS/ref_strain/locus.lst

gzip -dcf \
    ASSEMBLY/Pseudom_aeru_PAO1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_puti_KT2440_GCF_000007565_2/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_chl_aureofaciens_30_84_GCF_000281915_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_fluo_SBW25_GCF_000009225_2/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_prot_Pf_5_GCF_000012265_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_stu_A1501_GCF_000013785_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_syr_pv_syringae_B728a_GCF_000012245_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_sav_pv_phaseolicola_1448A_GCF_000012205_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_entomophi_L48_GCF_000026105_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_aeru_UCBPP_PA14_GCF_000014625_1/*_genomic.gff.gz |
    grep -v "^#" |
    grep -F -w -f DOMAINS/ref_strain/locus.lst |
    tsv-filter --str-eq 3:gene |
    perl -nl -e 'print $1 if /\bID=(.+?);/i' \
    > DOMAINS/ref_strain/gene.lst

gzip -dcf \
    ASSEMBLY/Pseudom_aeru_PAO1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_puti_KT2440_GCF_000007565_2/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_chl_aureofaciens_30_84_GCF_000281915_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_fluo_SBW25_GCF_000009225_2/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_prot_Pf_5_GCF_000012265_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_stu_A1501_GCF_000013785_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_syr_pv_syringae_B728a_GCF_000012245_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_sav_pv_phaseolicola_1448A_GCF_000012205_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_entomophi_L48_GCF_000026105_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_aeru_UCBPP_PA14_GCF_000014625_1/*_genomic.gff.gz |
    grep -v "^#" |
    grep -F -w -f DOMAINS/ref_strain/gene.lst |
    tsv-filter --str-eq 3:CDS |
    perl -nl -e 'print $1 if /\bName=(.+?);/i' \
    > DOMAINS/ref_strain/pro.lst

wc -l DOMAINS/ref_strain/*.lst
#  497 DOMAINS/ref_strain/gene.lst
#  513 DOMAINS/ref_strain/locus.lst
#  497 DOMAINS/ref_strain/pro.lst

faops some PROTEINS/all.uniq.fa DOMAINS/ref_strain/pro.lst stdout \
    > DOMAINS/ref_strain/pro.fa

faops size DOMAINS/ref_strain/pro.fa | wc -l
# 495

# prepare an HMM database for faster hmmscan searches
#hmmpress ~/data/HMM/PFAM/Pfam-A.hmm

E_VALUE=1e-5

hmmscan --cpu 4 -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw \
    ~/data/HMM/PFAM/Pfam-A.hmm DOMAINS/ref_strain/pro.fa |
    grep '>>' |
    perl -nl -e '/>>\s+(\S+)/ and print $1' |
    tee DOMAINS/ref_strain/domain.lst

# query ID against .hmm.dat to find AC
cat ~/data/HMM/PFAM/Pfam-A.hmm.dat |
    grep -F -w -f DOMAINS/ref_strain/domain.lst -A 2 |
    grep -E ' (ID|AC) ' |
    perl -nl -e 'print substr($_, 10) ' |
    paste -d $'\t' - - |
    perl -nlp -e 's/\.\d+$//g' |
    tsv-select -f 2,1 \
    > DOMAINS/ref_strain/pfam_domain.tsv

wc -l < DOMAINS/ref_strain/pfam_domain.tsv
# 107


```

### Scrap PFAM domains

* Perform keyword search in `pfam` and save the result page as `html only`.
    * [GO:0005975](https://pfam.xfam.org/search/keyword?query=GO%3A0005975)
    * [Glyco_hyd](http://pfam.xfam.org/search/keyword?query=Glyco_hyd)

```shell
cd ~/data/Pseudomonas/

cp DOMAINS/ref_strain/pfam_domain.tsv raw.tsv

cat GO_0005975.htm |
    pup 'table.resultTable tr td text{}' |
    grep '\S' |
    paste -d $'\t' - - - - - |
    tsv-select -f 2-4 \
    >> raw.tsv

cat Glyco_hyd.htm |
    pup 'table.resultTable tr td text{}' |
    grep '\S' |
    paste -d $'\t' - - - - - |
    tsv-select -f 2-4 \
    >> raw.tsv

#* [carbohydrate metabolic](http://pfam.xfam.org/search/keyword?query=carbohydrate+metabolic)
#cat carbohydrate_metabolic.htm |
#    pup 'table.resultTable tr td text{}' |
#    grep '\S' |
#    paste -d $'\t' - - - - - - - - |
#    tsv-select -f 2-4 \
#    >> raw.tsv

#* [Glycosyl_hydrolase](https://pfam.xfam.org/search/keyword?query=Glycosyl+hydrolase)
#cat Glycosyl_hydrolase.htm |
#    pup 'table.resultTable tr td text{}' |
#    grep '\S' |
#    paste -d $'\t' - - - - - - - - - |
#    tsv-filter --ne 5:10000000 | # Text fields of Pfam entries are not empty
#    tsv-select -f 2-4 \
#    >> raw.tsv

wc -l < raw.tsv
#290

cat raw.tsv |
    tsv-filter --str-not-in-fld 2:"DUF" |
    tsv-uniq -f 1 |
    tsv-sort -k2,2 \
    > pfam_domain.tsv

wc -l < pfam_domain.tsv
#216

mkdir -p ~/data/Pseudomonas/DOMAINS/HMM

cat pfam_domain.tsv |
    parallel --col-sep "\t" --no-run-if-empty --linebuffer -k -j 4 '
        >&2 echo {2}
        curl -L http://pfam.xfam.org/family/{1}/hmm > DOMAINS/HMM/{2}.hmm
    '

find DOMAINS/HMM -type f -name "*.hmm" |
    wc -l
#216

```

### Scan every domain

* The `E_VALUE` was adjusted to 1e-5 to capture all possible sequences.

```shell script
E_VALUE=1e-5

cd ~/data/Pseudomonas/

for domain in $(cat pfam_domain.tsv | cut -f 2 | sort); do
    >&2 echo "==> domain [${domain}]"

    if [ -e DOMAINS/${domain}.replace.tsv ]; then
        continue;
    fi

    for ORDER in $(cat order.lst); do
        >&2 echo "==> ORDER [${ORDER}]"

        cat taxon/${ORDER} |
            parallel --no-run-if-empty --linebuffer -k -j 8 "
                gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                    hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw DOMAINS/HMM/${domain}.hmm - |
                    grep '>>' |
                    perl -nl -e '
                        m{>>\s+(\S+)} or next;
                        \$n = \$1;
                        \$s = \$n;
                        \$s =~ s/\.\d+//;
                        printf qq{%s\t%s_%s\n}, \$n, {}, \$s;
                    '
            "
    done \
        > DOMAINS/${domain}.replace.tsv

    >&2 echo
done

for domain in $(cat pfam_domain.tsv | cut -f 2 | sort); do
    wc -l DOMAINS/${domain}.replace.tsv
done |
    datamash reverse -W |
    tsv-filter --ge 2:2000 |
    tsv-sort -k2,2nr |
    (echo -e "Domain\tCount" && cat) |
    mlr --itsv --omd cat

```

| Domain                              | Count |
|-------------------------------------|-------|
| DOMAINS/Epimerase.replace.tsv       | 40284 |
| DOMAINS/Hydrolase.replace.tsv       | 35532 |
| DOMAINS/HAD.replace.tsv             | 22725 |
| DOMAINS/CBS.replace.tsv             | 21493 |
| DOMAINS/GDP_Man_Dehyd.replace.tsv   | 21275 |
| DOMAINS/HAD_2.replace.tsv           | 20959 |
| DOMAINS/RmlD_sub_bind.replace.tsv   | 19014 |
| DOMAINS/Glycos_transf_2.replace.tsv | 18461 |
| DOMAINS/F420_oxidored.replace.tsv   | 15900 |
| DOMAINS/Glyco_tranf_2_3.replace.tsv | 15734 |
| DOMAINS/3Beta_HSD.replace.tsv       | 15218 |
| DOMAINS/Hydrolase_3.replace.tsv     | 14970 |
| DOMAINS/Hydrolase_like.replace.tsv  | 14372 |
| DOMAINS/PfkB.replace.tsv            | 11621 |
| DOMAINS/SIS.replace.tsv             | 11589 |
| DOMAINS/AAA_33.replace.tsv          | 10988 |
| DOMAINS/Glyco_trans_4_4.replace.tsv | 10216 |
| DOMAINS/Alpha-amylase.replace.tsv   | 8387  |
| DOMAINS/AP_endonuc_2.replace.tsv    | 7594  |
| DOMAINS/ApbA.replace.tsv            | 7173  |
| DOMAINS/Glyco_trans_2_3.replace.tsv | 7143  |
| DOMAINS/Phos_pyr_kin.replace.tsv    | 6989  |
| DOMAINS/CTP_transf_like.replace.tsv | 6561  |
| DOMAINS/Polysacc_deac_1.replace.tsv | 6150  |
| DOMAINS/FGGY_C.replace.tsv          | 6079  |
| DOMAINS/PGM_PMM_I.replace.tsv       | 5353  |
| DOMAINS/FGGY_N.replace.tsv          | 5288  |
| DOMAINS/PGM_PMM_II.replace.tsv      | 5240  |
| DOMAINS/PGM_PMM_III.replace.tsv     | 5234  |
| DOMAINS/PGM_PMM_IV.replace.tsv      | 4831  |
| DOMAINS/Glyco_transf_21.replace.tsv | 4791  |
| DOMAINS/NAD_Gly3P_dh_N.replace.tsv  | 4291  |
| DOMAINS/SKI.replace.tsv             | 3889  |
| DOMAINS/QRPTase_C.replace.tsv       | 3824  |
| DOMAINS/Glyco_hydro_3.replace.tsv   | 3698  |
| DOMAINS/Aldose_epim.replace.tsv     | 3536  |
| DOMAINS/CBM_48.replace.tsv          | 3407  |
| DOMAINS/DctQ.replace.tsv            | 3160  |
| DOMAINS/LamB_YcsF.replace.tsv       | 2784  |
| DOMAINS/Glyco_transf_28.replace.tsv | 2446  |
| DOMAINS/Glyco_tran_28_C.replace.tsv | 2409  |
| DOMAINS/TAL_FSA.replace.tsv         | 2340  |
| DOMAINS/F_bP_aldolase.replace.tsv   | 2239  |
| DOMAINS/Ribul_P_3_epim.replace.tsv  | 2124  |
| DOMAINS/Chitin_synth_2.replace.tsv  | 2029  |

Check each domain, such as `https://pfam.xfam.org/family/Epimerase`. Some domains are not directly
related to carbohydrate metabolism.

* ATP, AMP, GDP, CTP, and NADP associated
    * AAA_33
    * GDP_Man_Dehyd
    * CBS*
    * CTP_transf_like
* Binding to other small molecules
    * HAD* 卤素
    * F420*
    * SKI 莽草酸
    * QRPTase_C 喹啉酸
* Oxidoreductases
    * Epimerase
    * ApbA
    * 3Beta_HSD
    * NAD_Gly3P_dh_N

* We can't interproscan all `2420143` proteins

```shell script
cd ~/data/Pseudomonas

# All proteins appeared
cat pfam_domain.tsv | cut -f 2 | sort |
    tsv-filter --not-regex 1:'^AAA' |
    tsv-filter --not-regex 1:'^GDP' |
    tsv-filter --not-regex 1:'^CBS' |
    tsv-filter --not-regex 1:'^CTP' |
    tsv-filter --not-regex 1:'^HAD' |
    tsv-filter --not-regex 1:'^F420' |
    tsv-filter --not-regex 1:'^SKI' |
    tsv-filter --not-regex 1:'^QRPTase' |
    tsv-filter --not-regex 1:'^Epimerase' |
    tsv-filter --not-regex 1:'^ApbA' |
    tsv-filter --not-regex 1:'^3Beta_HSD' |
    tsv-filter --not-regex 1:'^NAD_Gly3P' \
    > domain.lst

wc -l < domain.lst
#202

cat domain.lst |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        tsv-select -f 2 DOMAINS/{}.replace.tsv
    ' |
    sort -u \
    > DOMAINS/domains.tsv

wc -l < DOMAINS/domains.tsv
#200612

faops size PROTEINS/all.uniq.fa | wc -l
#3944568

for domain in $(cat domain.lst); do
    echo 1>&2 "==> domain [${domain}]"

    tsv-join \
        DOMAINS/domains.tsv \
        --data-fields 1 \
        -f <(
            cat DOMAINS/${domain}.replace.tsv |
                perl -nla -e 'print qq{$F[1]\tO}'
        ) \
        --key-fields 1 \
        --append-fields 2 \
        --write-all "" \
        > DOMAINS/tmp.tsv

    mv DOMAINS/tmp.tsv DOMAINS/domains.tsv
done

datamash check < DOMAINS/domains.tsv
#200612 lines, 203 fields

# Add header line
for domain in $(cat domain.lst); do
    echo "${domain}"
done |
    (echo -e "#name" && cat) |
    paste -s -d $'\t' - \
    > DOMAINS/header.tsv

cat DOMAINS/header.tsv DOMAINS/domains.tsv \
    > tmp.tsv && mv tmp.tsv DOMAINS/domains.tsv

tsv-join \
    PROTEINS/all.info.tsv \
    --data-fields 1 \
    -f DOMAINS/domains.tsv \
    --key-fields 1 \
    --append-fields 2-203 |
     keep-header -- sort -k1,1 \
    > tmp.tsv && mv tmp.tsv DOMAINS/domains.tsv

datamash check < DOMAINS/domains.tsv
#200613 lines, 206 fields

rm DOMAINS/header.tsv

```

### InterProScan

The previous `hmmsearch` steps were done to narrow down the number of target proteins so that as few
proteins as possible would be passed to the following IPS steps.

InterProScan starts slowly, so dozens of proteins from one strain are submitted at once.

```shell script
cd ~/data/Pseudomonas

mkdir -p IPS

# extract wanted sequences
cat strains.lst |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        if [ $(({#} % 10)) -eq "0" ]; then
            >&2 printf "."
        fi

        mkdir -p IPS/{}

        cat PROTEINS/all.info.tsv |
            tsv-filter --str-eq 2:{} |
            cut -f 1 |
            grep -Fx -f <(cut -f 1 DOMAINS/domains.tsv | grep "^{}") \
            > IPS/{}/wanted.lst

        faops some PROTEINS/all.replace.fa IPS/{}/wanted.lst IPS/{}/{}.fa
    '

cat strains.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        cat IPS/{}/wanted.lst | wc -l
    ' |
    tsv-summarize --quantile 1:0.25,0.5,0.75
#64      111     128

# scan proteins of each strain with InterProScan
# By default InterProScan uses 8 cpu cores
mkdir -p split
split -a 4 -l 30 -d strains.lst split/

for f in $(find split -maxdepth 1 -type f -name "[0-9]*" | sort); do
    >&2 echo "==> IPS [${f}]"
    bsub -q mpi -n 24 -J "IPS-${f}" "
        cat ${f} |
            parallel --no-run-if-empty --linebuffer -k -j 6 '
                if [ -e IPS/{}/{}.tsv ]; then
                    >&2 echo {};
                    exit;
                fi

                interproscan.sh --cpu 4 -dp -f tsv,json -i IPS/{}/{}.fa --output-file-base IPS/{}/{}
            '
        "
done

rm -fr split output.*

find IPS -type f -name "*.json" | sort |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        if [ $(({#} % 10)) -eq "0" ]; then
            >&2 printf "."
        fi
        pigz -p 3 {}
    '

# IPS family
# Some proteins belong to more than one family. Only the best one is kept here
cat strains.lst |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        if [ $(({#} % 10)) -eq "0" ]; then
            >&2 printf "."
        fi
        gzip -dcf IPS/{}/{}.json.gz |
            jq .results |
            jq -r -c '\''
                .[] |
                .xref as $name |
                .matches[] |
                .signature.entry |
                select(.type == "FAMILY") |
                [$name[0].name, .accession, .description] |
                @tsv
            '\'' |
            tsv-uniq -f 1
    ' |
    (echo -e "#name\tfamily\tdescription" && cat) \
    > IPS/family.tsv

tsv-join \
    <(cut -f 1-4 DOMAINS/domains.tsv) \
    --data-fields 1 \
    -f IPS/family.tsv \
    --key-fields 1 \
    --append-fields 2-3 \
    --write-all "" |
    tsv-join \
        --data-fields 1 \
        -f DOMAINS/domains.tsv \
        --key-fields 1 \
        --append-fields 5-206 |
     keep-header -- sort -k1,1 \
    > IPS/predicts.tsv

```

## Search for gene families with uneven members

```shell
cd ~/data/Pseudomonas

cat IPS/predicts.tsv |
    tsv-filter -H --not-empty description |
    tsv-summarize -H -g strain,family,description --count |
    keep-header -- tsv-sort -k4,4n |
    tsv-summarize -H -g family,description --values 4 |
    keep-header -- tsv-sort -k2,2 |
    tsv-filter -H --str-in-fld 3:'|' |
    tr '|' '/' |
    keep-header -- perl -nla -F"\t" -e '
        my %seen;
        my @values = split "/", $F[2];
        next if scalar @values < 200; # The family is not widely distributed
        for my $v (@values) {
            $seen{$v}++;
        }
        next if keys %seen >= 5; # Too much variation in family members
        my @valid = grep { $seen{$_} > 50 } keys %seen;
        next if scalar @valid < 2; # At least 2 categories
        print join "\t", $F[0], $F[1], join "/", sort keys %seen;
    ' |
    tsv-filter -H --istr-not-in-fld 2:"probable" |
    tsv-filter -H --istr-not-in-fld 2:"putative" |
    tsv-filter -H --istr-not-in-fld 2:"Uncharacterised" \
    > variety.tsv

cat variety.tsv |
    mlr --itsv --omd cat

```

| family    | description                                               | count_values |
|-----------|-----------------------------------------------------------|--------------|
| IPR008183 | Aldose 1-/Glucose-6-phosphate 1-epimerase                 | 1/2/3/4      |
| IPR005841 | Alpha-D-phosphohexomutase superfamily                     | 1/2/3/4      |
| IPR005888 | dTDP-glucose 4,6-dehydratase                              | 1/2/3        |
| IPR014438 | Glucan biosynthesis protein MdoG/MdoD                     | 1/2/3/4      |
| IPR025532 | Glucose-6-phosphate 1-epimerase                           | 1/2          |
| IPR005999 | Glycerol kinase                                           | 1/2          |
| IPR011837 | Glycogen debranching enzyme, GlgX type                    | 1/2/3        |
| IPR000322 | Glycoside hydrolase family 31                             | 1/2/3/4      |
| IPR000811 | Glycosyl transferase, family 35                           | 1/2/3/4      |
| IPR035461 | GmhA/DiaA                                                 | 1/2          |
| IPR026040 | Hydroxypyruvate isomerase-like                            | 1/2/3/4      |
| IPR005501 | LamB/YcsF/PxpA-like                                       | 1/2/3/4      |
| IPR037950 | Peptidoglycan deacetylase PgdA-like                       | 1/2/3/4      |
| IPR004800 | Phosphosugar isomerase, KdsD/KpsF-type                    | 1/2/3        |
| IPR023854 | Poly-beta-1,6-N-acetyl-D-glucosamine N-deacetylase PgaB   | 1/2/3        |
| IPR023853 | Poly-beta-1,6 N-acetyl-D-glucosamine synthase PgaC/IcaA   | 1/2          |
| IPR006415 | P-type ATPase, subfamily IIIB                             | 1/2/3        |
| IPR004625 | Pyridoxine kinase                                         | 1/2/3        |
| IPR002347 | Short-chain dehydrogenase/reductase SDR                   | 1/2/3        |
| IPR017583 | Tagatose/fructose phosphokinase                           | 1/2/3        |
| IPR001585 | Transaldolase/Fructose-6-phosphate aldolase               | 1/2/3/4      |
| IPR004730 | Transaldolase type 1                                      | 1/2/3/4      |
| IPR005886 | UDP-glucose 4-epimerase                                   | 1/2/3        |
| IPR005593 | Xylulose 5-phosphate/Fructose 6-phosphate phosphoketolase | 1/2/3        |

### Within species

```shell
cd ~/data/Pseudomonas

for f in $(cat variety.tsv | tsv-select -f 1 | sed '1d' | sort); do
    cat IPS/predicts.tsv |
        tsv-select -f 1-3,5,6 |
        tsv-filter -H --str-eq family:"${f}" |
        tsv-join -d 2 \
            -f strains.taxon.tsv -k 1 \
            --append-fields 4 |
        tsv-summarize -g 6,2 --count |
        keep-header -- tsv-sort -k3,3n |
        tsv-summarize -g 1 --unique-values 3 |
        tsv-filter --str-in-fld 2:'|' |
        tsv-sort -k1,1 |
        tr '|' '/' \
        > tmp.tsv

    COUNT=$(cat tmp.tsv | wc -l)
    HAS_T=$(cat tmp.tsv | grep -E "aeruginosa|fluorescens|putida|syringae" | wc -l)
    HAS_ENOUGH=$(
        cat IPS/predicts.tsv |
            tsv-select -f 1-3,5,6 |
            tsv-filter -H --str-eq family:"${f}" |
            tsv-join -d 2 \
                -f strains.taxon.tsv -k 1 \
                --append-fields 4 |
            tsv-summarize -g 6,2 --count |
            tsv-filter --ge 3:2 |
            tsv-summarize -g 1 --count |
            tsv-filter --ge 2:5 |
            wc -l
    )
    if [[ ${COUNT} -ge "1" && ${COUNT} -le "10" && ${HAS_T} -ge "1" && ${HAS_ENOUGH} -ge "1" ]]; then
        cat tmp.tsv |
            sed "1 s/^/${f}\\t/" |
            sed "2,$ s/^/\\t/"
    fi

done |
    (echo -e "#family\tspecies\tcount" && cat) |
    mlr --itsv --omd cat

```

| #family   | species                    | count |
|-----------|----------------------------|-------|
| IPR000811 | Pseudomonas balearica      | 1/2   |
|           | Pseudomonas fluorescens    | 1/2   |
|           | Pseudomonas mandelii       | 1/2   |
|           | Pseudomonas stutzeri       | 1/2   |
|           | Pseudomonas synxantha      | 1/2   |
|           | Pseudomonas veronii        | 1/2   |
| IPR004730 | Alteromonas macleodii      | 1/2   |
|           | Alteromonas mediterranea   | 1/2   |
|           | Halomonas meridiana        | 1/2   |
|           | Halomonas titanicae        | 1/2/3 |
|           | Pseudomonas azotoformans   | 1/2   |
|           | Pseudomonas fluorescens    | 1/2   |
|           | Pseudomonas synxantha      | 1/2   |
| IPR004800 | Alteromonas macleodii      | 1/2   |
|           | Halomonas piezotolerans    | 1/2   |
|           | Pseudomonas atacamensis    | 1/2   |
|           | Pseudomonas citronellolis  | 1/2   |
|           | Pseudomonas fluorescens    | 1/2   |
|           | Pseudomonas fragi          | 1/2   |
|           | Pseudomonas putida         | 1/2   |
|           | Thalassolituus oleivorans  | 1/2   |
| IPR005593 | Pseudomonas aeruginosa     | 1/2   |
|           | Pseudomonas veronii        | 1/2   |
| IPR005999 | Alteromonas australica     | 1/2   |
|           | Pseudomonas aeruginosa     | 1/2   |
| IPR035461 | Alteromonas mediterranea   | 1/2   |
|           | Pseudomonas aeruginosa     | 1/2   |
|           | Pseudomonas lalkuanensis   | 1/2   |
| IPR037950 | Pseudomonas brassicacearum | 1/3   |
|           | Pseudomonas chlororaphis   | 1/2/3 |
|           | Pseudomonas fluorescens    | 1/2/3 |

### Among species

```shell
cd ~/data/Pseudomonas

for f in $(cat variety.tsv | tsv-select -f 1 | sed '1d' | sort); do
    cat IPS/predicts.tsv |
        tsv-select -f 1-3,5,6 |
        tsv-filter -H --str-eq family:"${f}" |
        tsv-join -d 2 \
            -f strains.taxon.tsv -k 1 \
            --append-fields 4 |
        tsv-summarize -g 6,2 --count |
        keep-header -- tsv-sort -k3,3n |
        tsv-summarize -g 1 --unique-values 3 |
        tsv-filter --str-not-in-fld 2:'|' |
        tsv-sort -k1,1 |
        tr '|' '/' \
        > tmp.tsv

    IS_1=$(cat tmp.tsv | tsv-filter --eq 2:1 | wc -l)
    IS_N=$(cat tmp.tsv | tsv-filter --ne 2:1 | wc -l)
    HAS_T=$(cat tmp.tsv | tsv-filter --ne 2:1 | grep -E "aeruginosa|chlororaphis|fluorescens|protegens|putida|syringae" | wc -l)

    if [[ ${IS_1} -ge "20" && ${IS_N} -ge "1" && ${IS_N} -le "10" && ${HAS_T} -ge "1" ]]; then
        printf "%s\t%d\t%d\t%d\n" $f $IS_1 $IS_N $HAS_T
        cat tmp.tsv | tsv-filter --ne 2:1 | grep -E "aeruginosa|chlororaphis|fluorescens|protegens|putida|syringae"
    fi
done

```

### IPR005999 - Glycerol kinase

* IPR000577 - Carbohydrate kinase, FGGY
    * IPR005999 - Glycerol kinase
    * IPR006000 - Xylulokinase

```shell
cd ~/data/Pseudomonas

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR005999  |
    tsv-summarize -H -g annotation --count
#annotation      count
#glycerol kinase GlpK    1357
#glycerol kinase 7
#glycerol kinase (sn-glycerol-3-phosphate generating)    1

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR000577  |
    tsv-summarize -H -g annotation --count |
    tsv-filter -H --gt 2:5
#annotation      count
#hypothetical protein    9
#FGGY-family carbohydrate kinase 707
#glycerol kinase GlpK    839
#carbohydrate kinase     176
#xylulokinase    362
#L-fuculokinase  19
#autoinducer-2 kinase    21
#rhamnulokinase  15
#glycerol kinase 57
#pentose kinase  9
#sugar kinase    14
#ribulokinase    19
#FGGY family carbohydrate kinase 19

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR006000  |
    tsv-summarize -H -g annotation --count |
    tsv-filter -H --gt 2:5
#annotation      count
#xylulokinase    671

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR005999 |
    tsv-summarize -H -g size --count |
    keep-header -- tsv-sort -k1,1n |
    tsv-filter -H --ge count:10
#size    count
#493     10
#494     340
#495     22
#499     71
#500     32
#501     196
#502     255
#503     48
#504     38
#505     282
#506     11
#507     18

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR005999 |
    tsv-select -f 1-6 |
    tsv-join -d 2 \
        -f strains.taxon.tsv -k 1 \
        --append-fields 4 |
    tsv-filter --or --str-eq 7:"Pseudomonas aeruginosa" |
    tsv-summarize -g 7,2 --count |
    tsv-filter --gt 3:1 |
    tsv-summarize -g 1 --count
#Pseudomonas aeruginosa  195

```

All proteins have the same structure FGGY_C and FGGY_N.

```shell
cd ~/data/Pseudomonas

mkdir -p ~/data/Pseudomonas/GlpK

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR005999 |
    datamash transpose |
    perl -nl -e '
        $row = $_;
        $row =~ s/\s//g;
        length($row) > 20 and print;
    ' |
    datamash transpose \
    > GlpK/GlpK.tsv

plotr tsv GlpK/GlpK.tsv --header

cat IPS/predicts.tsv |
    tsv-filter -H --or --str-eq family:IPR005999 --str-eq family:IPR000577 --str-eq family:IPR006000 |
    datamash transpose |
    perl -nl -e '
        $row = $_;
        $row =~ s/\s//g;
        length($row) > 20 and print;
    ' |
    datamash transpose \
    > GlpK/FGGY.tsv

cat DOMAINS/domains.tsv |
    tsv-filter -H --not-empty FGGY_C --not-empty FGGY_N |
    tsv-select -H -f 1-4,FGGY_C,FGGY_N \
    > GlpK/FGGY_N_C.tsv

wc -l GlpK/*.tsv
#   5270 GlpK/FGGY_N_C.tsv
#   4329 GlpK/FGGY.tsv
#   1366 GlpK/GlpK.tsv

for f in GlpK FGGY FGGY_N_C ; do
    >&2 echo "==> ${f}"

    faops some PROTEINS/all.replace.fa <(tsv-select -f 1 GlpK/${f}.tsv) stdout \
        > GlpK/${f}.fa

    muscle -in GlpK/${f}.fa -out GlpK/${f}.aln.fa

    FastTree GlpK/${f}.aln.fa > GlpK/${f}.aln.newick

    nw_reroot GlpK/${f}.aln.newick $(nw_labels GlpK/${f}.aln.newick | grep -E "Bac_subti|Sta_aure") |
        nw_order -c n - \
        > GlpK/${f}.reoot.newick

done

```

不同基因家族之间分得很开, 应该不会相互混淆.

在 PAO1 中, 这两个基因 PA3582 (glpK) 和 PA3579 只间隔一个基因. 不能确定是 HGT.

### IPR004800 - Phosphosugar isomerase, KdsD/KpsF-type

```shell
cd ~/data/Pseudomonas

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR004800 |
    tsv-summarize -H -g annotation --count
#annotation      count
#KpsF/GutQ family sugar-phosphate isomerase      1149
#D-arabinose 5-phosphate isomerase       8
#arabinose-5-phosphate isomerase 2
#carbohydrate isomerase  1
#D-arabinose 5-phosphate isomerase GutQ  1
#D-arabinose 5-phosphate isomerase KdsD  1
#arabinose-5-phosphate isomerase KdsD    389
#putative polysialic acid capsule expression protein     1
#hypothetical protein SF2731     1

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR004800 |
    tsv-summarize -H -g size --count |
    keep-header -- tsv-sort -k1,1n |
    tsv-filter -H --ge count:10
#size    count
#310     17
#315     13
#323     40
#324     452
#325     528
#326     422
#328     18
#339     16

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR004800 |
    tsv-select -f 1-6 |
    tsv-join -d 2 \
        -f strains.taxon.tsv -k 1 \
        --append-fields 4 |
    tsv-filter --or --str-eq 7:"Pseudomonas fluorescens" --str-eq 7:"Pseudomonas putida" |
    tsv-summarize -g 7,2 --count |
    tsv-filter --gt 3:1 |
    tsv-summarize -g 1 --count
#Pseudomonas fluorescens 2
#Pseudomonas putida      24

```

All proteins have the same structure SIS.

```shell
cd ~/data/Pseudomonas

mkdir -p ~/data/Pseudomonas/KdsD

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR004800 |
    datamash transpose |
    perl -nl -e '
        $row = $_;
        $row =~ s/\s//g;
        length($row) > 20 and print;
    ' |
    datamash transpose \
    > KdsD/KdsD.tsv

plotr tsv KdsD/KdsD.tsv --header

cat DOMAINS/domains.tsv |
    tsv-filter -H --not-empty SIS |
    tsv-select -H -f 1-4,SIS \
    > KdsD/SIS.tsv

wc -l KdsD/*.tsv
#   1554 KdsD/KdsD.tsv
#   8121 KdsD/SIS.tsv

for f in KdsD ; do
    >&2 echo "==> ${f}"

    faops some PROTEINS/all.replace.fa <(tsv-select -f 1 KdsD/${f}.tsv) stdout \
        > KdsD/${f}.fa

    muscle -in KdsD/${f}.fa -out KdsD/${f}.aln.fa

    FastTree KdsD/${f}.aln.fa > KdsD/${f}.aln.newick

    nw_reroot KdsD/${f}.aln.newick $(nw_labels KdsD/${f}.aln.newick | grep -E "B_sub|St_aur") |
        nw_order -c n - \
        > KdsD/${f}.reoot.newick

done

```

### IPR005593 - Xylulose 5-phosphate/Fructose 6-phosphate phosphoketolase

```shell
cd ~/data/Pseudomonas

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR005593  |
    tsv-summarize -H -g annotation --count
#annotation      count
#D-xylulose 5-phosphate/D-fructose 6-phosphate phosphoketolase   139
#xylulose 5-phosphate 3-epimerase        82
#hypothetical protein    87
#phosphoketolase 195
#hypothetical protein PA3613     1
#phosphoketolase family protein  126

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR005593 |
    tsv-summarize -H -g size --count |
    keep-header -- tsv-sort -k1,1n |
    tsv-filter -H --ge count:10
#size    count
#788     49
#789     10
#790     13
#791     40
#801     401
#805     16
#809     54

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR005593 |
    tsv-select -f 1-6 |
    tsv-join -d 2 \
        -f strains.taxon.tsv -k 1 \
        --append-fields 4 |
    tsv-filter --or --str-eq 7:"Pseudomonas aeruginosa" |
    tsv-summarize -g 7,2 --count |
    tsv-filter --gt 3:1 |
    tsv-summarize -g 1 --count
#Pseudomonas aeruginosa  50

```

All proteins have the same structure XFP and XFP_N, some have XFP_C.

```shell
cd ~/data/Pseudomonas

mkdir -p ~/data/Pseudomonas/XFP

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR005593 |
    datamash transpose |
    perl -nl -e '
        $row = $_;
        $row =~ s/\s//g;
        length($row) > 20 and print;
    ' |
    datamash transpose \
    > XFP/XFP.tsv

plotr tsv XFP/XFP.tsv --header

cat DOMAINS/domains.tsv |
    tsv-filter -H --or --not-empty XFP --not-empty XFP_N |
    tsv-select -H -f 1-4,XFP,XFP_N \
    > XFP/XFP_N.tsv

wc -l XFP/*.tsv
#  1074 XFP/XFP_N.tsv
#   631 XFP/XFP.tsv

for f in XFP XFP_N ; do
    >&2 echo "==> ${f}"

    faops some PROTEINS/all.replace.fa <(tsv-select -f 1 XFP/${f}.tsv) stdout \
        > XFP/${f}.fa

    muscle -in XFP/${f}.fa -out XFP/${f}.aln.fa

    FastTree XFP/${f}.aln.fa > XFP/${f}.aln.newick

    nw_reroot XFP/${f}.aln.newick $(nw_labels XFP/${f}.aln.newick | grep -E "B_sub|St_aur") |
        nw_order -c n - \
        > XFP/${f}.reoot.newick

done

```

### IPR035461 - GmhA/DiaA

```shell
cd ~/data/Pseudomonas

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR035461 |
    tsv-summarize -H -g annotation --count
#annotation      count
#phosphoheptose isomerase        565
#D-sedoheptulose 7-phosphate isomerase   10
#SIS domain-containing protein   80

cat IPS/predicts.tsv |
    tsv-filter -H --str-eq family:IPR035461 |
    tsv-select -f 1-6 |
    tsv-join -d 2 \
        -f strains.taxon.tsv -k 1 \
        --append-fields 4 |
    tsv-filter --or --str-eq 7:"Pseudomonas aeruginosa" |
    tsv-summarize -g 2 --count |
    tsv-filter --gt 2:1
#Pseudom_aeru_GCF_001548335_1     2

```

## InterProScan on all proteins of typical strains

```shell
cd ~/data/Pseudomonas

faops size ASSEMBLY/Pseudom_aeru_PAO1/*_protein.faa.gz |
    wc -l
#5572

faops size ASSEMBLY/Pseudom_aeru_PAO1/*_protein.faa.gz |
    tsv-summarize --sum 2
#1858983

mkdir -p STRAINS

for S in \
    Pseudom_aeru_PAO1 \
    Pseudom_puti_KT2440_GCF_000007565_2 \
    Pseudom_chl_aureofaciens_30_84_GCF_000281915_1 \
    Pseudom_entomophi_L48_GCF_000026105_1 \
    Pseudom_fluo_SBW25_GCF_000009225_2 \
    Pseudom_prot_Pf_5_GCF_000012265_1 \
    Pseudom_sav_pv_phaseolicola_1448A_GCF_000012205_1 \
    Pseudom_stu_A1501_GCF_000013785_1 \
    Pseudom_syr_pv_syringae_B728a_GCF_000012245_1 \
    Pseudom_aeru_UCBPP_PA14_GCF_000014625_1 \
    Pseudom_aeru_PA7_GCF_000017205_1 \
    Pseudom_aeru_LESB58_GCF_000026645_1 \
    ; do
    echo ${S}
done \
    > typical.lst

for S in $(cat typical.lst); do
    mkdir -p STRAINS/${S}
    faops split-about ASSEMBLY/${S}/*_protein.faa.gz 200000 STRAINS/${S}/
done

for S in $(cat typical.lst); do
    for f in $(find STRAINS/${S}/ -maxdepth 1 -type f -name "[0-9]*.fa" | sort); do
        >&2 echo "==> ${f}"
        if [ -e ${f}.tsv ]; then
            >&2 echo ${f}
            continue
        fi

        bsub -q mpi -n 24 -J "${f}" "
            interproscan.sh --cpu 24 -dp -f tsv,json -i ${f} --output-file-base ${f}
        "
    done
done

find STRAINS -type f -name "*.json" | sort |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        if [ $(({#} % 10)) -eq "0" ]; then
            >&2 printf "."
        fi
        pigz -p 3 {}
    '

# same protein may have multiple families
for S in $(cat typical.lst); do
    for f in $(find STRAINS/${S} -maxdepth 1 -type f -name "[0-9]*.json.gz" | sort); do
        >&2 echo "==> ${f}"
        gzip -dcf ${f} |
            jq .results |
            jq -r -c '
                .[] |
                .xref as $name |
                .matches[] |
                .signature.entry |
                select(.type == "FAMILY") |
                [$name[0].name, .accession, .description] |
                @tsv
            ' |
            tsv-uniq
    done \
        > STRAINS/${S}/family.tsv
done

COUNT=
for S in $(cat typical.lst); do
    if [ ! -s STRAINS/${S}/family.tsv ]; then
        continue
    fi
    cat STRAINS/${S}/family.tsv |
        tsv-summarize -g 2,3 --count \
        > STRAINS/${S}/family-count.tsv

    COUNT=$((COUNT + 1))
done
echo $COUNT

# families in all strains
for S in $(cat typical.lst); do
    cat STRAINS/${S}/family-count.tsv
done |
    tsv-summarize -g 1,2 --count |
    tsv-filter -H --istr-not-in-fld 2:"probable" |
    tsv-filter -H --istr-not-in-fld 2:"putative" |
    tsv-filter -H --istr-not-in-fld 2:"Uncharacterised" |
    tsv-filter -H --istr-not-in-fld 2:" DUF" |
    tsv-filter --ge 3:$COUNT \
    > STRAINS/universal.tsv

# All other strains should have only 1 family member
cp STRAINS/universal.tsv STRAINS/family-1.tsv
for S in $(cat typical.lst | grep -v "_aeru_"); do
    if [ ! -s STRAINS/${S}/family-count.tsv ]; then
        continue
    fi
    cat STRAINS/${S}/family-count.tsv |
        tsv-join -k 1 -f STRAINS/family-1.tsv |
        tsv-filter --eq 3:1 \
        > STRAINS/family-tmp.tsv

    mv STRAINS/family-tmp.tsv STRAINS/family-1.tsv
done

# All P_aeru strains should have multiple family members
cp STRAINS/family-1.tsv STRAINS/family-n.tsv
for S in $(cat typical.lst | grep "_aeru_"); do
    if [ ! -s STRAINS/${S}/family-count.tsv ]; then
        continue
    fi
    cat STRAINS/${S}/family-count.tsv |
        tsv-join -k 1 -f STRAINS/family-n.tsv |
        tsv-filter --gt 3:1 \
        > STRAINS/family-tmp.tsv

    wc -l < STRAINS/family-tmp.tsv
    mv STRAINS/family-tmp.tsv STRAINS/family-n.tsv
done

wc -l STRAINS/Pseudom_aeru_PAO1/family.tsv STRAINS/universal.tsv STRAINS/family-1.tsv STRAINS/family-n.tsv
#  4084 STRAINS/Pseudom_aeru_PAO1/family.tsv
#  1567 STRAINS/universal.tsv
#   972 STRAINS/family-1.tsv
#    14 STRAINS/family-n.tsv

cat STRAINS/family-n.tsv |
    tsv-select -f 1,2 |
    (echo -e "#family\tcount" && cat) |
    mlr --itsv --omd cat

```

| #family   | count                                                         |
|-----------|---------------------------------------------------------------|
| IPR014311 | Guanine deaminase                                             |
| IPR001404 | Heat shock protein Hsp90 family                               |
| IPR005999 | Glycerol kinase                                               |
| IPR000813 | 7Fe ferredoxin                                                |
| IPR011757 | Lytic transglycosylase MltB                                   |
| IPR007416 | YggL 50S ribosome-binding protein                             |
| IPR004361 | Glyoxalase I                                                  |
| IPR024922 | Rubredoxin                                                    |
| IPR001353 | Proteasome, subunit alpha/beta                                |
| IPR002307 | Tyrosine-tRNA ligase                                          |
| IPR024088 | Tyrosine-tRNA ligase, bacterial-type                          |
| IPR037532 | Peptidoglycan D,D-transpeptidase FtsI                         |
| IPR003672 | CobN/magnesium chelatase                                      |
| IPR004685 | Branched-chain amino acid transport system II carrier protein |

### IPR007416 - YggL 50S ribosome-binding protein

* Pfam:     PF04320 - YggL_50S_bp
* PANTHER:  PTHR38778 - CYTOPLASMIC PROTEIN-RELATED (PTHR38778)

```shell
cd ~/data/Pseudomonas

cat STRAINS/Pseudom_aeru_PAO1/*.tsv |
    grep "IPR007416"

mkdir -p YggL/HMM

curl -L http://pfam.xfam.org/family/PF04320/hmm > YggL/HMM/YggL_50S_bp.hmm
curl -L www.pantherdb.org/panther/exportHmm.jsp?acc=PTHR38778 > YggL/HMM/PTHR38778.hmm

# Ribosomal protein L10 and S8
curl -L http://pfam.xfam.org/family/PF00466/hmm > YggL/HMM/Ribosomal_L10.hmm
curl -L http://pfam.xfam.org/family/PF00410/hmm > YggL/HMM/Ribosomal_S8.hmm

E_VALUE=1e-20
for domain in YggL_50S_bp PTHR38778 Ribosomal_L10 Ribosomal_S8 ; do
    >&2 echo "==> domain [${domain}]"

    if [ -e YggL/${domain}.replace.tsv ]; then
        continue;
    fi

    for ORDER in $(cat order.lst); do
        >&2 echo "==> ORDER [${ORDER}]"

        cat taxon/${ORDER} |
            parallel --no-run-if-empty --linebuffer -k -j 8 "
                gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                    hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw YggL/HMM/${domain}.hmm - |
                    grep '>>' |
                    perl -nl -e '
                        m{>>\s+(\S+)} or next;
                        \$n = \$1;
                        \$s = \$n;
                        \$s =~ s/\.\d+//;
                        printf qq{%s\t%s_%s\n}, \$n, {}, \$s;
                    '
            "
    done \
        > YggL/${domain}.replace.tsv

    >&2 echo
done

tsv-join YggL/YggL_50S_bp.replace.tsv \
    -f YggL/PTHR38778.replace.tsv \
    > YggL/YggL.replace.tsv

wc -l YggL/*.tsv
#  1559 YggL/PTHR38778.replace.tsv
#  1948 YggL/Ribosomal_L10.replace.tsv
#  1950 YggL/Ribosomal_S8.replace.tsv
#  1564 YggL/YggL_50S_bp.replace.tsv
#  1559 YggL/YggL.replace.tsv

faops some PROTEINS/all.replace.fa <(tsv-select -f 2 YggL/YggL.replace.tsv) YggL/YggL.fa

muscle -in YggL/YggL.fa -out YggL/YggL.aln.fa

FastTree YggL/YggL.aln.fa > YggL/YggL.aln.newick

nw_reroot YggL/YggL.aln.newick $(nw_labels YggL/YggL.aln.newick | grep -E "Bac_subti|Sta_aure") |
    nw_order -c n - \
    > YggL/YggL.reoot.newick

```

## Collect CDS

### `all.cds.fa`

```shell
cd ~/data/Pseudomonas

mkdir -p CDS

find ASSEMBLY -type f -name "*_cds_from_genomic.fna.gz" |
    wc -l
# 1958

# sed script converting from Contigs to Strain
for ORDER in $(cat order.lst); do
    echo 1>&2 "==> ORDER [${ORDER}]"

    for STRAIN in $(cat taxon/${ORDER}); do
        find ASSEMBLY/${STRAIN} -type f -name "*_genomic.fna.gz" |
            grep -v "_from_" |
            xargs gzip -dcf |
            grep '^>' |
            cut -d' ' -f 1 |
            sed 's/>//' |
            xargs -I{} echo -e "{}\t${STRAIN}"
    done
done \
    > CDS/contigs_to_strain.tsv

cat CDS/contigs_to_strain.tsv |
    perl -nla -e '
        print q{s/^>} . quotemeta($F[0]) . q{/>} . quotemeta($F[1]) . q{/g;};
    ' \
    > CDS/sed.script

wc -l < CDS/sed.script
# 3893

for ORDER in $(cat order.lst); do
    echo 1>&2 "==> ORDER [${ORDER}]"

    for STRAIN in $(cat taxon/${ORDER}); do
        gzip -dcf ASSEMBLY/${STRAIN}/*_cds_from_genomic.fna.gz
    done
done |
    perl -nl -e 's/^>lcl\|/>/g; print' |
    perl -nl -e 's/\s+\[.+?\]//g; print' \
    > CDS/all.cds.fa

```

### `YggL.cds.fa`

```shell
cd ~/data/Pseudomonas

for domain in YggL Ribosomal_L10 Ribosomal_S8 ; do
    cat CDS/all.cds.fa |
        grep '>' |
        grep -F -f <( cat YggL/${domain}.replace.tsv | cut -f 1 ) |
        sed 's/^>//' \
        > CDS/${domain}.lst
done

for domain in YggL Ribosomal_L10 Ribosomal_S8 ; do
    faops order CDS/all.cds.fa CDS/${domain}.lst stdout |
        sed -f CDS/sed.script \
        > CDS/${domain}.cds.fa
done

```
