
use strict;
use warnings;
use autodie;

use Bio::EnsEMBL::Utils::ConfigRegistry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

my $host = 'localhost';
my $port = 3306;
my $user = 'alignDB';
my $pass = 'alignDB';

# Arabidopsis thaliana
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Arabidopsis thaliana',
        -group   => 'core',
        -dbname  => 'arabidopsis_thaliana_core_45_98_11',
    );

    my @aliases = (
        'atha',
        'Atha',
        'Arabidopsis_thaliana',
        'atha_core_98',
        'atha_98',
        'arabidopsis_thaliana_core_45_98_11',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Arabidopsis thaliana',
        -alias   => \@aliases,
    );
}

# Aspergillus fumigatus
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Aspergillus fumigatus',
        -group   => 'core',
        -dbname  => 'aspergillus_fumigatus_core_45_98_1',
    );

    my @aliases = (
        'afum',
        'Afum',
        'Aspergillus_fumigatus',
        'afum_core_98',
        'afum_98',
        'aspergillus_fumigatus_core_45_98_1',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Aspergillus fumigatus',
        -alias   => \@aliases,
    );
}

# Caenorhabditis elegans
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Caenorhabditis elegans',
        -group   => 'core',
        -dbname  => 'caenorhabditis_elegans_core_98_269',
    );

    my @aliases = (
        'cele',
        'Cele',
        'Caenorhabditis_elegans',
        'cele_core_98',
        'cele_98',
        'caenorhabditis_elegans_core_98_269',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Caenorhabditis elegans',
        -alias   => \@aliases,
    );
}

# Dictyostelium discoideum
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Dictyostelium discoideum',
        -group   => 'core',
        -dbname  => 'dictyostelium_discoideum_core_45_98_1',
    );

    my @aliases = (
        'ddis',
        'Ddis',
        'Dictyostelium_discoideum',
        'ddis_core_98',
        'ddis_98',
        'dictyostelium_discoideum_core_45_98_1',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Dictyostelium discoideum',
        -alias   => \@aliases,
    );
}

# Drosophila melanogaster
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Drosophila melanogaster',
        -group   => 'core',
        -dbname  => 'drosophila_melanogaster_core_98_7',
    );

    my @aliases = (
        'dmel',
        'Dmel',
        'Drosophila_melanogaster',
        'dmel_core_98',
        'dmel_98',
        'drosophila_melanogaster_core_98_7',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Drosophila melanogaster',
        -alias   => \@aliases,
    );
}

# Homo sapiens
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Homo sapiens',
        -group   => 'core',
        -dbname  => 'homo_sapiens_core_98_38',
    );

    my @aliases = (
        'hsap',
        'Hsap',
        'Homo_sapiens',
        'hsap_core_98',
        'hsap_98',
        'homo_sapiens_core_98_38',
        'human',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Homo sapiens',
        -alias   => \@aliases,
    );
}

# Mus musculus
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Mus musculus',
        -group   => 'core',
        -dbname  => 'mus_musculus_core_98_38',
    );

    my @aliases = (
        'mmus',
        'Mmus',
        'Mus_musculus',
        'mmus_core_98',
        'mmus_98',
        'mus_musculus_core_98_38',
        'mouse',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Mus musculus',
        -alias   => \@aliases,
    );
}

# Oryza sativa
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Oryza sativa',
        -group   => 'core',
        -dbname  => 'oryza_sativa_core_45_98_7',
    );

    my @aliases = (
        'osat',
        'Osat',
        'Oryza_sativa',
        'osat_core_98',
        'osat_98',
        'oryza_sativa_core_45_98_7',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Oryza sativa',
        -alias   => \@aliases,
    );
}

# Plasmodium falciparum
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Plasmodium falciparum',
        -group   => 'core',
        -dbname  => 'plasmodium_falciparum_core_45_98_1',
    );

    my @aliases = (
        'pfal',
        'Pfal',
        'Plasmodium_falciparum',
        'pfal_core_98',
        'pfal_98',
        'plasmodium_falciparum_core_45_98_1',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Plasmodium falciparum',
        -alias   => \@aliases,
    );
}

# Saccharomyces cerevisiae
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Saccharomyces cerevisiae',
        -group   => 'core',
        -dbname  => 'saccharomyces_cerevisiae_core_98_4',
    );

    my @aliases = (
        'scer',
        'Scer',
        'Saccharomyces_cerevisiae',
        'scer_core_98',
        'scer_98',
        'saccharomyces_cerevisiae_core_98_4',
        'yeast',
        'S288c',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Saccharomyces cerevisiae',
        -alias   => \@aliases,
    );
}

# Schizosaccharomyces pombe
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Schizosaccharomyces pombe',
        -group   => 'core',
        -dbname  => 'schizosaccharomyces_pombe_core_45_98_2',
    );

    my @aliases = (
        'spom',
        'Spom',
        'Schizosaccharomyces_pombe',
        'spom_core_98',
        'spom_98',
        'schizosaccharomyces_pombe_core_45_98_2',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Schizosaccharomyces pombe',
        -alias   => \@aliases,
    );
}


1;
