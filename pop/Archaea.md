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

There are several major groups
of [Archaea](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2157), some of which are
relatively new and are not included in this study.

* Asgard group
* Candidatus Hydrothermarchaeota
* Candidatus Thermoplasmatota
* DPANN group
* Euryarchaeota
* TACK group

### List all ranks

```shell
mkdir -p ~/data/Archaea
cd ~/data/Archaea

nwr member Archaea |
    grep -v " sp." |
    tsv-summarize -H -g 3 --count |
    mlr --itsv --omd cat

```

| rank          | count |
|---------------|-------|
| superkingdom  | 1     |
| phylum        | 39    |
| no rank       | 293   |
| species       | 3320  |
| class         | 29    |
| clade         | 28    |
| order         | 49    |
| family        | 69    |
| genus         | 230   |
| strain        | 353   |
| species group | 2     |
| isolate       | 6     |

### Species with assemblies

```shell
mkdir -p ~/data/Archaea/summary
cd ~/data/Archaea/summary

# should have a valid name of genus
GENUS=$(
    nwr member Archaea -r genus |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        grep -v " sp." |
        sed '1d' |
        cut -f 1 |
        sort -n
)

for R in ${GENUS}; do
    echo "
        SELECT
            species_id,
            species,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = $R
            AND assembly_level IN ('Complete Genome', 'Chromosome')
        GROUP BY species_id
        HAVING count > 2
        " |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
    sed '1d' |
    parallel --col-sep "\t" -j 1 '
        GROUP=$(
            nwr lineage {1} |
                grep -E "\sArchaea\s" -A 1 |
                sed "1d" |
                tsv-select -f 2
        )
        printf "%s\t%s\t%s\t%s\n" {1} "$GROUP" "{2}" {3}
    ' |
    sed "s/'//g" |
    nwr append stdin -r family |
    tsv-select -f 1,2,5,3,4 |
    tsv-sort -k2,2 -k3,3 -k4,4 |
    (echo -e '#tax_id\tgroup\tfamily\tspecies\tRS' && cat) \
    > species.count.tsv

cat species.count.tsv |
    mlr --itsv --omd cat

```

| #tax_id | group         | family             | species                    | RS  |
|---------|---------------|--------------------|----------------------------|-----|
| 2841257 | Euryarchaeota | Haloarculaceae     | Halapricum desulfuricans   | 4   |
| 2242    | Euryarchaeota | Halobacteriaceae   | Halobacterium salinarum    | 3   |
| 2252    | Euryarchaeota | Haloferacaceae     | Haloferax mediterranei     | 3   |
| 2208    | Euryarchaeota | Methanosarcinaceae | Methanosarcina barkeri     | 5   |
| 2209    | Euryarchaeota | Methanosarcinaceae | Methanosarcina mazei       | 9   |
| 38027   | Euryarchaeota | Methanosarcinaceae | Methanosarcina siciliae    | 3   |
| 2261    | Euryarchaeota | Thermococcaceae    | Pyrococcus furiosus        | 3   |
| 43687   | TACK group    | Sulfolobaceae      | Metallosphaera sedula      | 6   |
| 2286    | TACK group    | Sulfolobaceae      | Saccharolobus shibatae     | 3   |
| 2287    | TACK group    | Sulfolobaceae      | Saccharolobus solfataricus | 12  |
| 2285    | TACK group    | Sulfolobaceae      | Sulfolobus acidocaldarius  | 10  |
| 43080   | TACK group    | Sulfolobaceae      | Sulfolobus islandicus      | 15  |

## Download all assemblies

```shell
cd ~/data/Archaea

echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND family IN (
            'Haloarculaceae', 'Halobacteriaceae', 'Haloferacaceae',
            'Methanosarcinaceae', 'Thermococcaceae', 'Sulfolobaceae'
        )
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND assembly_level IN ('Complete Genome', 'Chromosome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    > raw.tsv

cat raw.tsv |
    grep -v '^#' |
    tsv-uniq |
    perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "1,2,3" -s '\t' -m 3 --shortsub |
    (echo -e '#name\tftp_path\torganism\tassembly_level' && cat ) |
    perl -nl -a -F"," -e '
        BEGIN{my %seen};
        /^#/ and print and next;
        /^organism_name/i and next;
        $seen{$F[5]}++;
        $seen{$F[5]} > 1 and next;
        printf qq{%s\t%s\t%s\t%s\n}, $F[5], $F[3], $F[1], $F[4];
        ' |
    keep-header -- sort -k3,3 -k1,1 \
    > Archaea.assembly.tsv

# find potential duplicated strains or assemblies
cat Archaea.assembly.tsv |
    tsv-uniq -f 1 --repeated

# Edit .tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Archaea.assembly.tsv
# cp Archaea.assembly.tsv ~/Scripts/withncbi/pop

# Comment out unneeded strains

# Cleaning
rm raw*.*sv

```

```shell
cd ~/data/Archaea

perl ~/Scripts/withncbi/taxon/assembly_prep.pl \
    -f ~/Scripts/withncbi/pop/Archaea.assembly.tsv \
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
proxychains4 bash ASSEMBLY/Archaea.assembly.rsync.sh

bash ASSEMBLY/Archaea.assembly.collect.sh

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
