#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use Path::Tiny;
use MCE::Flow;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/../lib";
use MyEnsembl;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::RealBin/../config.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

=head1 NAME

feature_runlists.pl - runlists of features from Ensembl

=head1 SYNOPSIS

    perl feature_runlists.pl [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port          INT     MySQL server port
        --username  -u  STR     username
        --password  -p  STR     password
        --ensembl   -e  STR     ensembl database name
        --feature   -f  STR     feature name, default is [repeat]
        --output    -o  STR     output filename
        --parallel      INT     run in parallel mode

    perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
        -e saccharomyces_cerevisiae_core_29_82_4 \
        -f repeat \
        -o repeat.yml

    perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
        -e saccharomyces_cerevisiae_core_29_82_4 \
        -f cds \
        -o cds.yml

    runlist merge repeat.yml cds.yml -o yeast.yml
    rm repeat.yml cds.yml

=cut

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'server|s=s'   => \( my $server     = $Config->{database}{server} ),
    'port=i'       => \( my $port       = $Config->{database}{port} ),
    'username|u=s' => \( my $username   = $Config->{database}{username} ),
    'password|p=s' => \( my $password   = $Config->{database}{password} ),
    'ensembl|e=s'  => \( my $ensembl_db = $Config->{database}{ensembl} ),
    'feature|f=s'  => \( my $feature    = 'repeat' ),
    'output|o=s'   => \( my $out_file ),
    'parallel=i'   => \( my $parallel   = $Config->{run}{parallel} ),
) or Getopt::Long::HelpMessage(1);

if ( !defined $out_file ) {
    $out_file = "${ensembl_db}.${feature}.yml";
}

#----------------------------------------------------------#
# Run
#----------------------------------------------------------#

# ensembl handler
my $ensembl = MyEnsembl->new(
    server => $server,
    db     => $ensembl_db,
    port   => $port,
    user   => $username,
    passwd => $password,
);

# ensembl handler
my $db_adaptor    = $ensembl->db_adaptor;
my $slice_adaptor = $db_adaptor->get_SliceAdaptor;

# get chromosome list
my @slices = @{ $slice_adaptor->fetch_all('chromosome') };
my @chrs = sort { $a->{chr_name} cmp $b->{chr_name} }
    map { { chr_name => $_->seq_region_name, chr_start => $_->start, chr_end => $_->end, } }
    @slices;

my $worker = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;
    my $chr = $chunk_ref->[0];

    my $chr_runlist = $chr->{chr_start} . "-" . $chr->{chr_end};
    printf "%s:%s\n", $chr->{chr_name}, $chr_runlist;

    eval { $ensembl->set_slice( $chr->{chr_name}, $chr->{chr_start}, $chr->{chr_end} ); };
    if ($@) {
        warn "Can't get annotation\n";
        return;
    }

    my $slice       = $ensembl->slice;
    my $ftr_chr_set = $slice->{"_$feature\_set"};

    MCE->gather( $chr->{chr_name}, $ftr_chr_set->runlist );
};

MCE::Flow::init {
    chunk_size  => 1,
    max_workers => $parallel,
};
my %feature_of = mce_flow $worker, \@chrs;
MCE::Flow::finish;

$stopwatch->block_message("Write output file [$out_file]");
YAML::Syck::DumpFile( $out_file, \%feature_of );

$stopwatch->end_message;

__END__
