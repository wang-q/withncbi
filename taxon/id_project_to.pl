#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use Scalar::Util::Numeric qw(isint);
use Bio::DB::Taxonomy;

use FindBin;
use lib "$FindBin::Bin/../lib";
use MyUtil qw(replace_home);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::Bin/../config.ini");

my $column    = 1;
my $seperator = '\s+';

my $rank;
my $rankid;

# running options
my $td_dir = replace_home( $Config->{path}{td} );    # taxdmp

my $man  = 0;
my $help = 0;

GetOptions(
    'help'          => \$help,
    'man'           => \$man,
    'c|column=i'    => \$column,
    's|seperator=s' => \$seperator,
    'rank=s'        => \$rank,
    'rankid'        => \$rankid,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$|++;

# check $rank
my @ranks = qw{species genus family order class phylum superkingdom};
my %valid = map { $_ => 1 } @ranks;

if ( $rank and !$valid{$rank} ) {
    die "--rank [$rank} is invalid, it should be one of [@ranks]\n";
}

my $taxon_db;
if ( -e "$td_dir/nodes.dmp" ) {
    $taxon_db = Bio::DB::Taxonomy->new(
        -source    => 'flatfile',
        -directory => $td_dir,
        -nodesfile => "$td_dir/nodes.dmp",
        -namesfile => "$td_dir/names.dmp",
    );
}
else {
    $taxon_db = Bio::DB::Taxonomy->new( -source => 'entrez', );
}

#----------------------------------------------------------#
# start
#----------------------------------------------------------#
while ( my $line = <> ) {
    chomp $line;
    my @row = split /$seperator/, $line;

    my $taxon_id = $row[ $column - 1 ];

    if ( isint($taxon_id) ) {
        my $strain = $taxon_db->get_taxon( -taxonid => $taxon_id );
        if ( !$strain ) {
            warn "Can't find taxon for [$taxon_id]\n";
            push @row, 'NA';
            push @row, '0' if $rankid;
        }
        else {
            if ($rank) {
                my $taxon_obj = find_ancestor( $strain, $rank );
                if ( !$taxon_obj ) {
                    warn "Can't find [$rank] for [$taxon_id]\n";
                    push @row, 'NA';
                    push @row, '0' if $rankid;
                }
                else {
                    push @row, $taxon_obj->scientific_name;
                    push @row, $taxon_obj->id if $rankid;
                }
            }
            else {
                push @row, $strain->scientific_name;
                push @row, $strain->id if $rankid;
            }
        }
    }
    else {
        warn "[$taxon_id] doesn't look like a taxonmy id.\n";
        push @row, 'NA';
        push @row, '0' if $rankid;
    }

    print join( ",", @row ), "\n";
}

exit;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#
sub find_ancestor {
    my $taxon = shift;
    my $rank = shift || 'species';

    return $taxon if $taxon->rank eq $rank;

RANK: while (1) {
        $taxon = $taxon->ancestor;
        last RANK unless defined $taxon;
        return $taxon if $taxon->rank eq $rank;
    }

    return;
}

__END__

=head1 NAME

id_project_to.pl - Project taxonomy ids to names, sepcies, genus or higher ranks.

=head1 SYNOPSIS

    cat <file> | perl id_project_to.pl [options]
      Options:
        --help                  brief help message
        --man                   full documentation
        -c, --column INT        column order where ids are, start from 1
        -s, --seperator STR     seperator of the line, default is "\s+"
        --rank STR              Project to which rank, default is scientific name.
        --rankid                Also append rank id

=head1 EXAMPLE

    $ echo 9606 | perl taxon/id_project_to.pl 
    9606,Homo sapiens

    $ echo 9606 | perl taxon/id_project_to.pl --rank class 
    9606,Mammalia

    $ echo 9606 Human | perl taxon/id_project_to.pl --rank class --rankid
    9606,Human,Mammalia,40674

    $ echo Human,9606  | perl taxon/id_project_to.pl -c 2 -s ","
    Human,9606,Homo sapiens

=cut
