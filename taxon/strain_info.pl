#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use Path::Tiny;
use Text::CSV_XS;
use Bio::DB::Taxonomy;

use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/../lib";
use MyUtil qw(find_ancestor);

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

strain_info.pl - generate a csv file for taxonomy info

=head1 SYNOPSIS

    perl strain_info.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        -o, --file STR      output filename
        --id @INT           ids
        --stdin             Read ids from stdin, --id is still working
        --withname          stdin is id,name, 9606,Human
        --name INT=STR      Assign name to id
        --species INT=STR   fake id need this
        --simple            means use subspecies strain name as name
        --entrez            don't use local taxdmp

=head1 EXAMPLE

    perl strain_info.pl \
        --file   yeast_ncbi.csv \
        --simple \
        --id     559292         \
        --id     285006         \
        --id     307796         \
        --id     226125         \
        --name   226125=Spar

=cut

# running options
my $td_dir = path( $Config->{path}{td} )->stringify;    # taxdmp

GetOptions(
    'help|?'    => sub { HelpMessage(0) },
    'id=i'      => \my @ids,
    'stdin'     => \my $stdin,
    'withname'  => \my $stdin_with_name,
    'name=s'    => \my %name_of,
    'species=s' => \my %species_of,
    'file|o=s' => \( my $filename = "strains_taxon_info.csv" ),
    'simple'   => \my $simple,
    'entrez'   => \my $entrez,
) or HelpMessage(1);

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Writing strains summary...");

my $taxon_db;
if ( !$entrez and -e "$td_dir/nodes.dmp" ) {
    $stopwatch->block_message("Load ncbi taxdmp.");
    $taxon_db = Bio::DB::Taxonomy->new(
        -source    => 'flatfile',
        -directory => $td_dir,
        -nodesfile => "$td_dir/nodes.dmp",
        -namesfile => "$td_dir/names.dmp",
    );
}
else {
    $stopwatch->block_message("Use online ncbi entrez.");
    $taxon_db = Bio::DB::Taxonomy->new( -source => 'entrez', );
}

if ($stdin) {
    while (<>) {
        chomp;
        if ($stdin_with_name) {
            my ( $in_id, $in_name ) = split /,/, $_;
            push @ids, $in_id;
            $name_of{$in_id} = $in_name;
        }
        else {
            push @ids, $_;
        }
    }
}

$stopwatch->block_message("Processing...");

my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n" } );
open my $csv_fh, ">", $filename;

# headers
$csv->print( $csv_fh, [ map { ( $_, $_ . "_id" ) } qw{strain species genus family order} ] );

ID: for my $taxon_id (@ids) {
    my @row;

    my $strain;
    if ( exists $species_of{$taxon_id} ) {

        $strain = $taxon_db->get_taxon( -taxonid => $species_of{$taxon_id} );
        if ( !$strain ) {
            warn "Can't find taxon for $species_of{$taxon_id}. Fake id is $taxon_id\n";
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

__END__
