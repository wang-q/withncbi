#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use List::MoreUtils qw(uniq);

use FindBin;
use lib "$FindBin::Bin/../lib";
use MyUtil qw(abbr_most);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $column    = '1,2,3';
my $seperator = '\s+';

my $min_species = 3;

my $tight;

my $man  = 0;
my $help = 0;

GetOptions(
    'help'          => \$help,
    'man'           => \$man,
    'c|column=s'    => \$column,
    's|seperator=s' => \$seperator,
    'm|min=i'       => \$min_species,
    'tight'         => \$tight,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$|++;

#----------------------------------------------------------#
# start
#----------------------------------------------------------#
my @columns = map { $_ - 1 } split( /,/, $column );
my @rows;
while ( my $line = <> ) {
    chomp $line;
    next unless $line;
    my @row = split /$seperator/, $line;
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

    s/\W+/_/g for ( $strain, $species, $genus );
    push @rows, [ $strain, $species, $genus ];
}

my $count = scalar @rows;

my @ge = map { $_->[2] } @rows;
my @sp = map { $_->[1] } @rows;
my @st = map { $_->[0] } @rows;

my $ge_of = abbr_most( [ uniq(@ge) ], 1,            "Yes" );
my $sp_of = abbr_most( [ uniq(@sp) ], $min_species, "Yes" );

for my $i ( 0 .. $count - 1 ) {
    my $spacer = $tight ? '' : '_';
    my $ge_sp = join $spacer, grep { defined and length } $ge_of->{ $ge[$i] },
        $sp_of->{ $sp[$i] };
    my $organism = join "_", grep { defined and length } $ge_sp, $st[$i];

    print $organism, "\n";
}

# TODO: append abbr to existing columns, just like id_project_to.pl

exit;

__END__

=head1 NAME

abbr_name.pl - Abbreviate strain scientific names.

=head1 SYNOPSIS

    cat <file> | perl id_project_to.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        -c, --column STR    Columns of strain, species, genus, default is 1,2,3.
                            If there's no strain, use 1,1,2.
                            Don't need the strain part, use 2,2,3
                            When there's only strain, use 1,1,1
                            
        -s, --seperator STR seperator of the line, default is "\s+"
        -m, --min INT       mininal length for abbreviation of species
        --tight             No underscore between Genus and species

=head1 EXAMPLE

    $ echo -e 'Homo sapiens,Homo\nHomo erectus,Homo\n' | perl abbr_name.pl -s ',' -c "1,1,2"
    H_sap
    H_ere
    
    $ echo -e 'Homo sapiens,Homo\nHomo erectus,Homo\n' | perl abbr_name.pl -s ',' -c "1,1,2" --tight
    Hsap
    Here
    
    $ echo -e 'Homo sapiens sapiens,Homo sapiens,Homo\nHomo erectus,Homo erectus,Homo\n' \
        | perl abbr_name.pl -s ',' -c "1,2,3" --tight
    Hsap_sapiens
    Here
    
=cut
