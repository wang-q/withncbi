#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use DBI;
use Text::CSV_XS;
use DateTime::Format::Natural;
use List::MoreUtils qw(any all uniq);

use Bio::DB::Taxonomy;

use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use MyUtil qw(replace_home);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../config.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

my @ids;
my %name_of;

# for arbitrary ids
my %species_of;

# for unrecorded strains, give them arbitrary ids
my $arbitrary = 10_000_000;

# running options
my $td_dir = replace_home( $Config->{path}{td} );    # taxdmp

my $filename = "strains_taxon_info.csv";

# simplify strain name
my $simple;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'    => \$help,
    'man'       => \$man,
    'id=i'      => \@ids,
    'name=s'    => \%name_of,
    'species=s' => \%species_of,
    'o|file=s'    => \$filename,
    'simple'    => \$simple,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Writing strains summary...");

my $taxon_db;
if ( -e "$td_dir/nodes.dmp" ) {
    $stopwatch->block_message("Load ncbi taxdmp.");
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

my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n" } );
open my $csv_fh, ">", $filename;

# headers
$csv->print( $csv_fh,
    [ map { ( $_, $_ . "_id" ) } qw{strain species genus family order} ] );

ID: for my $taxon_id (@ids) {
    my @row;

    my $strain;
    if ( exists $species_of{$taxon_id} ) {

        $strain = $taxon_db->get_taxon( -taxonid => $species_of{$taxon_id} );
        if ( !$strain ) {
            warn
                "Can't find taxon for $species_of{$taxon_id}. Fake id is $taxon_id\n";
            next;
        }

        if ( !exists $name_of{$taxon_id} ) {
            warn "Fake id should has its own name.\n";
            next;
        }
        push @row, ( "SHOULD BE REPLACED", $taxon_id, );
    }
    else {
        $strain = $taxon_db->get_taxon( -taxonid => $taxon_id );
        if ( !$strain ) {
            warn "Can't find taxon for $taxon_id\n";
            next;
        }
        push @row, ( $strain->scientific_name, $strain->id, );
    }

    for my $level (qw{species}) {
        my $taxon_obj = find_ancestor( $strain, $level );
        if ( !$taxon_obj ) {
            warn "Can't find $level for $taxon_id\n";
            next ID;
        }
        push @row, ( $taxon_obj->scientific_name, $taxon_obj->id, );
    }

    if ( exists $name_of{$taxon_id} ) {
        $row[0] = $name_of{$taxon_id};
    }
    elsif ($simple) {
        my $sub_name = $row[0];
        $sub_name =~ s/^$row[2]\s*//;
        $sub_name =~ s/\W/_/g;
        $sub_name =~ s/_+/_/g;
        $row[0] = $sub_name;
    }

    for my $level (qw{genus family order}) {
        my $taxon_obj = find_ancestor( $strain, $level );
        if ( !$taxon_obj ) {
            warn "Can't find $level for $taxon_id\n";
            push @row, ( undef, undef, );
        }
        else {
            push @row, ( $taxon_obj->scientific_name, $taxon_obj->id, );
        }
    }

    $csv->print( $csv_fh, \@row );
}
close $csv_fh;

#----------------------------#
# Finish
#----------------------------#
$stopwatch->end_message;
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

strain_info.pl - generate a csv file for taxonomy info

=head1 SYNOPSIS

    perl strain_info.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        -o, --file          output filename
        --id                @ids
        --name              %name_of
        --species           %species_of, fake id need this
        --simple            means use subspecies strain name as name
        -a, --aria2         url file is for aria2
        --fix               sometimes WGS records miss assigning strain id

    perl strain_info.pl \
        --file   yeast_ncbi.csv \
        --simple \
        --id     559292         \
        --id     285006         \
        --id     307796         \
        --id     226125         \
        --name   226125=Spar

=cut
