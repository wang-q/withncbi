create table strain
(
    taxonomy_id                     int,
    organism_name                   text,
    super_kingdom                   text,
    `group`                         text,
    genome_size                     double,
    gc_content                      double,
    number_of_chromosomes           int,
    number_of_plasmids              int,
    refseq_accessions               text,
    gram_stain                      text,
    shape                           text,
    arrangment                      text,
    endospores                      text,
    motility                        text,
    salinity                        text,
    oxygen_req                      text,
    habitat                         text,
    temp_range                      text,
    optimal_temp                    text,
    pathogenic_in                   text,
    disease                         text,
    released_date                   date,
    modified_date                   date,
    species                         text,
    species_id                      int,
    genus                           text,
    genus_id                        int,
    species_member                  int,
    genus_species_member            int,
    genus_strain_member             int,
    primary key (taxonomy_id)
)
ENGINE = MyISAM;

create table seq
(
    accession                       char(32),
    accession_version               text,
    genbankacc                      text,
    taxonomy_id                     int,
    organism_name                   text,
    replicon                        text,
    length                          int,
    primary key (accession),
    index(taxonomy_id)
)
ENGINE = MyISAM;

create table segment
(
    segment_id                      int             not null AUTO_INCREMENT,
    accession                       char(32),
    segment_type                    char(1),
    segment_gc_mean                 double,
    segment_gc_std                  double,
    segment_gc_cv                   double,
    segment_gc_mdcw                 double,
    primary key (segment_id),
    index(accession)
)
ENGINE = MyISAM;

create table gr
(
    taxonomy_id                     int,
    organism_name                   text,
    bioproject                      text,
    `group`                         text,
    subgroup                        text,
    genome_size                     double,
    gc_content                      double,
    chr_refseq                      text,
    plasmid_refseq                  text,
    wgs                             text,
    scaffolds                       int,
    released_date                   date,
    status                          text,
    species                         text,
    species_id                      int,
    genus                           text,
    genus_id                        int,
    species_member                  int,
    genus_species_member            int,
    genus_strain_member             int,
    primary key (taxonomy_id)
)
ENGINE = MyISAM;