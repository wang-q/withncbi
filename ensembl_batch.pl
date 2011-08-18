#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Find::Rule;
use File::Spec;
use Path::Class;

use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};

# write_axt parameter
my $yml_file = "$FindBin::Bin/aspergillus.yml";

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    'server=s'     => \$server,
    'port=i'       => \$port,
    'username=s'   => \$username,
    'password=s'   => \$password,
    'y|yml_file=s' => \$yml_file,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Ensembl batching for $yml_file");

my $dispatch = LoadFile($yml_file);

my $mysql_dir = dir( $dispatch->{dir}{mysql} );
my $fasta_dir = dir( $dispatch->{dir}{fasta} );

my $species_ref = $dispatch->{species};
my @species     = sort keys %{$species_ref};

my $version = $dispatch->{meta}{version};

#----------------------------------------------------------#
# Write .axt files from alignDB
#----------------------------------------------------------#
for my $sp (@species) {
    print "$sp\n";
    if ( $species_ref->{$sp}{core} ) {
        my ( $ensembl_dir, $ensembl_db ) = match_ensembl( $sp, 'core' );
        if ( $species_ref->{$sp}{db} ) {
            $ensembl_db = $species_ref->{$sp}{db};
        }

        my $cmd
            = "perl $FindBin::Bin/build_ensembl.pl"
            . " -s=$server"
            . " --port=$port"
            . " -u=$username"
            . " --password=$password"
            . " --initdb"
            . " --db=$ensembl_db"
            . " --ensembl=$ensembl_dir";

        $stopwatch->block_message("Build ensembl for $sp");
        $stopwatch->block_message($cmd);
        system $cmd;
        $stopwatch->block_message("Finish build");
    }
}

$stopwatch->end_message;

sub match_ensembl {
    my $sp   = shift;
    my $type = shift;

    my ( $genus, $species ) = split " ", $sp;

    my $str = join "_", ( $genus, $species, $type );
    $str = lc $str;

    my $ensembl_dir;
    while ( my $child = $mysql_dir->next ) {
        next unless -d $child;
        my $dir = $child->stringify;
        my $flag = index $dir, $str;
        next if $flag == -1;
        $ensembl_dir = $dir;
    }

    my $ensembl_db
        = lc( substr( $genus, 0, 1 ) )
        . substr( $species, 0, 3 )
        . "_$type"
        . "_$version";

    return ( $ensembl_dir, $ensembl_db );
}

__END__

=head1 NAME

    write_axt_slice.pl - extract alignment slices from alignDB

=head1 SYNOPSIS

    write_axt_slice.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --db                database name
        --username          username
        --password          password
        -y, --yaml_dir      dir of yaml

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

format 1: chr1.yml
  --- 1-25744,815056-817137
  
  output axt: chr1/chr1.axt
format 2: bin1.yml
  ---
  chr1: '1-25744,815056-817137'

  output axt: bin1/chr1.axt


=cut

