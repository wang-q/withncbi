#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use List::MoreUtils::PP;

use lib "$FindBin::RealBin/../lib";
use MyUtil;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

abbr_name.pl - Abbreviate strain scientific names.

=head1 SYNOPSIS

    cat <file> | perl abbr_name.pl [options]
      Options:
        --help              brief help message
        --column    -c  STR Columns of strain, species, genus, default is 1,2,3.
                            If there's no strain, use 1,1,2.
                            Don't need the strain part, use 2,2,3
                            When there's only strain, use 1,1,1
        --separator -s  STR separator of the line, default is "\s+"
        --min INT           mininal length for abbreviation of species
        --tight             no underscore between Genus and species
        --shortsub          clean subspecies parts

=head1 EXAMPLE

    $ echo -e 'Homo sapiens,Homo\nHomo erectus,Homo\n' \
        | perl abbr_name.pl -s ',' -c "1,1,2"
    H_sap
    H_ere
    
    $ echo -e 'Homo sapiens,Homo\nHomo erectus,Homo\n' \
        | perl abbr_name.pl -s ',' -c "1,1,2" --tight
    Hsap
    Here
    
    $ echo -e 'Homo sapiens sapiens,Homo sapiens,Homo\nHomo erectus,Homo erectus,Homo\n' \
        | perl abbr_name.pl -s ',' -c "1,2,3" --tight
    Hsap_sapiens
    Here
    
=cut

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'column|c=s'    => \( my $column      = '1,2,3' ),
    'separator|s=s' => \( my $separator   = '\s+' ),
    'min|m=i'       => \( my $min_species = 3 ),
    'tight'         => \my $tight,
    'shortsub'      => \my $shortsub,
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$|++;

#----------------------------------------------------------#
# start
#----------------------------------------------------------#
my @columns = map { $_ - 1 } split( /,/, $column );
my @fields;
my @rows;
while ( my $line = <> ) {
    chomp $line;
    next unless $line;
    my @row = split /$separator/, $line;
    s/"|'//g for @row;

    my ( $strain, $species, $genus ) = @row[@columns];
    if ( $genus eq $species ) {
        my @O = split /\s+/, $strain;    # Organism
        if ( scalar @O >= 2 ) {
            $genus  = shift @O;
            $strain = shift @O;
            $strain = join "_", @O;
        }
        else {
            warn "Parsing strain [$strain] error.\n";
        }
    }
    else {
        $strain =~ s/^$species ?//;
        $species =~ s/^$genus //;
    }

    # Remove `Candidatus`
    $genus =~ s/\bCandidatus \b/C/g;

    # Clean long subspecies names
    if ($shortsub) {
        $strain =~ s/\bsubsp\b//g;
        $strain =~ s/\bserovar\b//g;
        $strain =~ s/\bstr\b//g;
        $strain =~ s/\bstrain\b//g;
        $strain =~ s/\bserotype\b//g;
        $strain =~ s/\bbiovar\b//g;
    }

    s/\W+/_/g for ( $strain, $species, $genus );
    s/_+/_/g  for ( $strain, $species, $genus );
    s/_$//    for ( $strain, $species, $genus );
    s/^_//    for ( $strain, $species, $genus );

    push @fields, [ $strain, $species, $genus ];

    push @rows, \@row;
}

my $count = scalar @fields;

my @ge = map { $_->[2] } @fields;
my @sp = map { $_->[1] } @fields;
my @st = map { $_->[0] } @fields;

my $ge_of = MyUtil::abbr_most( [ List::MoreUtils::PP::uniq(@ge) ], 1,            "Yes" );
my $sp_of = MyUtil::abbr_most( [ List::MoreUtils::PP::uniq(@sp) ], $min_species, "Yes" );

for my $i ( 0 .. $count - 1 ) {
    my $spacer = $tight ? '' : '_';
    my $ge_sp = join $spacer, grep { defined and length } $ge_of->{ $ge[$i] }, $sp_of->{ $sp[$i] };
    my $organism = join "_", grep { defined and length } $ge_sp, $st[$i];

    print join( ",", @{ $rows[$i] }, $organism ), "\n";
}

exit;

__END__
