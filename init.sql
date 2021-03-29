# noinspection SqlNoDataSourceInspectionForFile

create table gr
(
    taxonomy_id          int,
    organism_name        char(255),
    bioproject           char(255),
    `group`              char(255),
    subgroup             char(255),
    genome_size          double,
    gc_content           double,
    chr                  text,
    wgs                  char(255),
    scaffolds            int,
    released_date        date,
    asm_name             char(255),
    status               char(255),
    species              char(255),
    species_id           int,
    genus                char(255),
    genus_id             int,
    species_member       int,
    genus_species_member int,
    genus_strain_member  int,
    primary key (taxonomy_id)
)
ENGINE = MyISAM CHARSET=latin1;

create index gr_organism_name_K on gr
(
   organism_name
);

create table ar
(
    taxonomy_id          int,
    organism_name        char(255),
    bioproject           char(255),
    assembly_accession   char(255),
    wgs_master           char(255),
    refseq_category      char(255),
    assembly_level       char(255),
    genome_rep           char(255),
    released_date        date,
    asm_name             char(255),
    ftp_path             char(255),
    superkingdom         char(255),
    `group`              char(255),
    subgroup             char(255),
    species              char(255),
    species_id           int,
    genus                char(255),
    genus_id             int,
    species_member       int,
    genus_species_member int,
    genus_strain_member  int,
    primary key (taxonomy_id)
)
ENGINE = MyISAM CHARSET=latin1;

create index ar_organism_name_K on ar
(
   organism_name
);
