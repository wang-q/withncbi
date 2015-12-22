#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB::Stopwatch;
use AlignDB::Ensembl;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new();
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server     = $Config->{database}{server};
my $port       = $Config->{database}{port};
my $username   = $Config->{database}{username};
my $password   = $Config->{database}{password};
my $ensembl_db = $Config->{database}{ensembl};

my @chr_names;
my $out_file;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    's|server=s'   => \$server,
    'P|port=i'     => \$port,
    'u|username=s' => \$username,
    'p|password=s' => \$password,
    'e|ensembl=s'  => \$ensembl_db,
    'chr=s'        => \@chr_names,
    'output=s'     => \$out_file,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

$out_file ||= "chromosome.band.$ensembl_db.txt";
@chr_names = qw{ 1 } unless @chr_names;

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Drawing chr bands...");

# ensembl object
my $ensembl = AlignDB::Ensembl->new(
    server => $server,
    db     => $ensembl_db,
    user   => $username,
    passwd => $password,
);

my $db_adaptor   = $ensembl->db_adaptor;
my $kary_adaptor = $db_adaptor->get_KaryotypeBandAdaptor;

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#
open my $band_fh, '>', $out_file;
print {$band_fh} "#chrom\tchromStart\tchromEnd\tname\tgieStain\n";
for my $chr (@chr_names) {
    print "chr$chr\n";
    my $band = $kary_adaptor->fetch_all_by_chr_name($chr);

    my @bands = sort { $a->start <=> $b->start } @$band;
    for (@bands) {
        print {$band_fh} "chr$chr\t@{[$_->start]}\t@{[$_->end]}\t";
        print {$band_fh} "@{[$_->name]}\t@{[$_->stain]}\n";
    }
}
close $band_fh;

$stopwatch->end_message;
exit;

__END__

=head1 NAME

    chr_kary.pl - Fetch Karyotype Band from ensembl db

=head1 SYNOPSIS

    chr_kary.pl [options]
        Options:
            --help              brief help message
            --man               full documentation
            --server            MySQL server IP/Domain name
            --username          username
            --password          password
            --ensembl           ensembl db name
            --chr_name          chromosome name
            --output            output file name

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do someting
useful with the contents thereof.

=cut
