#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use Path::Tiny;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::Mapper::RangeRegistry;

use AlignDB::Stopwatch;

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

chr_kary.pl - Fetch Karyotype Band from ensembl db

=head1 SYNOPSIS

    perl chr_kary.pl -e <ensembl> [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port      -P  INT     MySQL server port
        --username  -u  STR     username
        --password  -p  STR     password
        --ensembl   -e  STR     ensembl database name
        --chr           @STR    chromosome names, leave empty will write all chromosomes
        --output        STR     output file name

=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'server|s=s'   => \( my $server   = $Config->{database}{server} ),
    'port|P=i'     => \( my $port     = $Config->{database}{port} ),
    'username|u=s' => \( my $username = $Config->{database}{username} ),
    'password|p=s' => \( my $password = $Config->{database}{password} ),
    'ensembl|e=s'  => \my $ensembl_db,
    'chr=s'        => \my @chr_names,
    'output=s'     => \my $out_file,
) or HelpMessage(1);

if ( !defined $ensembl_db ) {
    die HelpMessage(1);
}

$out_file ||= "$ensembl_db.kary.tsv";

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
$stopwatch->start_message("Drawing chr bands...");

my $db_adaptor = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host   => $server,
    -dbname => $ensembl_db,
    -user   => $username,
    -pass   => $password,
) or die "Cannot connect to EnsEMBL database [$ensembl_db]\n";

my $kary_adaptor = $db_adaptor->get_KaryotypeBandAdaptor;

if (! @chr_names) {
    @chr_names = map { $_->seq_region_name } @{ $db_adaptor->get_GenomeContainer->get_karyotype };
    
    # The following codes return unsorted chromosome names
    #@chr_names = map { $_->seq_region_name } @{ $db_adaptor->get_SliceAdaptor->fetch_all('chromosome') };
    
    print "Write all available karyotypes\n";
    print Dump \@chr_names;
}

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#
my $band_fh = path($out_file)->openw;
print {$band_fh} "#chrom\tchromStart\tchromEnd\tname\tgieStain\n";
for my $chr (@chr_names) {
    print "Region [$chr]\n";
    my $band = $kary_adaptor->fetch_all_by_chr_name($chr);

    my @bands = sort { $a->start <=> $b->start } @$band;
    for (@bands) {
        printf {$band_fh} "%s\t%s\t%s\t%s\t%s\n", $chr, $_->start, $_->end, $_->name, $_->stain;
    }
}
close $band_fh;

$stopwatch->end_message;
exit;

__END__
