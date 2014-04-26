create table gr
(
    taxonomy_id                     int,
    organism_name                   char(255),
    bioproject                      char(255),
    `group`                         char(255),
    subgroup                        char(255),
    genome_size                     double,
    gc_content                      double,
    chr                             text,
    wgs                             char(255),
    scaffolds                       int,
    released_date                   date,
    status                          char(255),
    species                         char(255),
    species_id                      int,
    genus                           char(255),
    genus_id                        int,
    species_member                  int,
    genus_species_member            int,
    genus_strain_member             int,
    primary key (taxonomy_id)
)
ENGINE = MyISAM;

create index organism_name_K on gr
(
   organism_name
);

