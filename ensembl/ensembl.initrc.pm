
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
        -dbname  => 'arabidopsis_thaliana_core_29_82_10',
    );

    my @aliases = (
        'atha',
        'Atha',
        'Arabidopsis_thaliana',
        'atha_core_82',
        'atha_82',
        'arabidopsis_thaliana_core_29_82_10',
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
        -dbname  => 'aspergillus_fumigatus_core_29_82_2',
    );

    my @aliases = (
        'afum',
        'Afum',
        'Aspergillus_fumigatus',
        'afum_core_82',
        'afum_82',
        'aspergillus_fumigatus_core_29_82_2',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Aspergillus fumigatus',
        -alias   => \@aliases,
    );
}

# Brassica rapa
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Brassica rapa',
        -group   => 'core',
        -dbname  => 'brassica_rapa_core_29_82_1',
    );

    my @aliases = (
        'brap',
        'Brap',
        'Brassica_rapa',
        'brap_core_82',
        'brap_82',
        'brassica_rapa_core_29_82_1',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Brassica rapa',
        -alias   => \@aliases,
    );
}

# Caenorhabditis brenneri
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Caenorhabditis brenneri',
        -group   => 'core',
        -dbname  => 'caenorhabditis_brenneri_core_29_82_233',
    );

    my @aliases = (
        'cbre',
        'Cbre',
        'Caenorhabditis_brenneri',
        'cbre_core_82',
        'cbre_82',
        'caenorhabditis_brenneri_core_29_82_233',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Caenorhabditis brenneri',
        -alias   => \@aliases,
    );
}

# Caenorhabditis briggsae
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Caenorhabditis briggsae',
        -group   => 'core',
        -dbname  => 'caenorhabditis_briggsae_core_29_82_230',
    );

    my @aliases = (
        'cbri',
        'Cbri',
        'Caenorhabditis_briggsae',
        'cbri_core_82',
        'cbri_82',
        'caenorhabditis_briggsae_core_29_82_230',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Caenorhabditis briggsae',
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
        -dbname  => 'caenorhabditis_elegans_core_82_245',
    );

    my @aliases = (
        'cele',
        'Cele',
        'Caenorhabditis_elegans',
        'cele_core_82',
        'cele_82',
        'caenorhabditis_elegans_core_82_245',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Caenorhabditis elegans',
        -alias   => \@aliases,
    );
}

# Caenorhabditis japonica
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Caenorhabditis japonica',
        -group   => 'core',
        -dbname  => 'caenorhabditis_japonica_core_29_82_233',
    );

    my @aliases = (
        'cjap',
        'Cjap',
        'Caenorhabditis_japonica',
        'cjap_core_82',
        'cjap_82',
        'caenorhabditis_japonica_core_29_82_233',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Caenorhabditis japonica',
        -alias   => \@aliases,
    );
}

# Caenorhabditis remanei
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Caenorhabditis remanei',
        -group   => 'core',
        -dbname  => 'caenorhabditis_remanei_core_29_82_233',
    );

    my @aliases = (
        'crem',
        'Crem',
        'Caenorhabditis_remanei',
        'crem_core_82',
        'crem_82',
        'caenorhabditis_remanei_core_29_82_233',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Caenorhabditis remanei',
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
        -dbname  => 'dictyostelium_discoideum_core_29_82_1',
    );

    my @aliases = (
        'ddis',
        'Ddis',
        'Dictyostelium_discoideum',
        'ddis_core_82',
        'ddis_82',
        'dictyostelium_discoideum_core_29_82_1',
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
        -dbname  => 'drosophila_melanogaster_core_82_602',
    );

    my @aliases = (
        'dmel',
        'Dmel',
        'Drosophila_melanogaster',
        'dmel_core_82',
        'dmel_82',
        'drosophila_melanogaster_core_82_602',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Drosophila melanogaster',
        -alias   => \@aliases,
    );
}

# Drosophila sechellia
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Drosophila sechellia',
        -group   => 'core',
        -dbname  => 'drosophila_sechellia_core_29_82_1',
    );

    my @aliases = (
        'dsec',
        'Dsec',
        'Drosophila_sechellia',
        'dsec_core_82',
        'dsec_82',
        'drosophila_sechellia_core_29_82_1',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Drosophila sechellia',
        -alias   => \@aliases,
    );
}

# Drosophila simulans
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Drosophila simulans',
        -group   => 'core',
        -dbname  => 'drosophila_simulans_core_29_82_1',
    );

    my @aliases = (
        'dsim',
        'Dsim',
        'Drosophila_simulans',
        'dsim_core_82',
        'dsim_82',
        'drosophila_simulans_core_29_82_1',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Drosophila simulans',
        -alias   => \@aliases,
    );
}

# Gorilla gorilla
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Gorilla gorilla',
        -group   => 'core',
        -dbname  => 'gorilla_gorilla_core_82_31',
    );

    my @aliases = (
        'ggor',
        'Ggor',
        'Gorilla_gorilla',
        'ggor_core_82',
        'ggor_82',
        'gorilla_gorilla_core_82_31',
        'gorilla',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Gorilla gorilla',
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
        -dbname  => 'homo_sapiens_core_82_37',
    );

    my @aliases = (
        'hsap',
        'Hsap',
        'Homo_sapiens',
        'hsap_core_82',
        'hsap_82',
        'homo_sapiens_core_82_37',
        'human',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Homo sapiens',
        -alias   => \@aliases,
    );
}

# Macaca mulatta
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Macaca mulatta',
        -group   => 'core',
        -dbname  => 'macaca_mulatta_core_82_10',
    );

    my @aliases = (
        'mmul',
        'Mmul',
        'Macaca_mulatta',
        'mmul_core_82',
        'mmul_82',
        'macaca_mulatta_core_82_10',
        'rhesus',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Macaca mulatta',
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
        -dbname  => 'mus_musculus_core_82_38',
    );

    my @aliases = (
        'mmus',
        'Mmus',
        'Mus_musculus',
        'mmus_core_82',
        'mmus_82',
        'mus_musculus_core_82_38',
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
        -dbname  => 'oryza_sativa_core_29_82_7',
    );

    my @aliases = (
        'osat',
        'Osat',
        'Oryza_sativa',
        'osat_core_82',
        'osat_82',
        'oryza_sativa_core_29_82_7',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Oryza sativa',
        -alias   => \@aliases,
    );
}

# Pan troglodytes
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Pan troglodytes',
        -group   => 'core',
        -dbname  => 'pan_troglodytes_core_82_214',
    );

    my @aliases = (
        'ptro',
        'Ptro',
        'Pan_troglodytes',
        'ptro_core_82',
        'ptro_82',
        'pan_troglodytes_core_82_214',
        'chimp',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Pan troglodytes',
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
        -dbname  => 'plasmodium_falciparum_core_29_82_3',
    );

    my @aliases = (
        'pfal',
        'Pfal',
        'Plasmodium_falciparum',
        'pfal_core_82',
        'pfal_82',
        'plasmodium_falciparum_core_29_82_3',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Plasmodium falciparum',
        -alias   => \@aliases,
    );
}

# Pongo abelii
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Pongo abelii',
        -group   => 'core',
        -dbname  => 'pongo_abelii_core_82_1',
    );

    my @aliases = (
        'pabe',
        'Pabe',
        'Pongo_abelii',
        'pabe_core_82',
        'pabe_82',
        'pongo_abelii_core_82_1',
        'orangutan',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Pongo abelii',
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
        -dbname  => 'saccharomyces_cerevisiae_core_29_82_4',
    );

    my @aliases = (
        'scer',
        'Scer',
        'Saccharomyces_cerevisiae',
        'scer_core_82',
        'scer_82',
        'saccharomyces_cerevisiae_core_29_82_4',
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
        -dbname  => 'schizosaccharomyces_pombe_core_29_82_2',
    );

    my @aliases = (
        'spom',
        'Spom',
        'Schizosaccharomyces_pombe',
        'spom_core_82',
        'spom_82',
        'schizosaccharomyces_pombe_core_29_82_2',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Schizosaccharomyces pombe',
        -alias   => \@aliases,
    );
}

# Solanum lycopersicum
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Solanum lycopersicum',
        -group   => 'core',
        -dbname  => 'solanum_lycopersicum_core_29_82_250',
    );

    my @aliases = (
        'slyc',
        'Slyc',
        'Solanum_lycopersicum',
        'slyc_core_82',
        'slyc_82',
        'solanum_lycopersicum_core_29_82_250',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Solanum lycopersicum',
        -alias   => \@aliases,
    );
}

# Solanum tuberosum
{
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => 'Solanum tuberosum',
        -group   => 'core',
        -dbname  => 'solanum_tuberosum_core_29_82_4',
    );

    my @aliases = (
        'stub',
        'Stub',
        'Solanum_tuberosum',
        'stub_core_82',
        'stub_82',
        'solanum_tuberosum_core_29_82_4',
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => 'Solanum tuberosum',
        -alias   => \@aliases,
    );
}


1;
