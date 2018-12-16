#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Path::Tiny qw();
use Mojo::DOM;
use Tie::IxHash;
use YAML::Syck qw();

# http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=2759&opt=plastid
# http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=33090&opt=organelle

my $file = shift;
die "Provide a valid html file!\n" unless $file;

my $html = Path::Tiny::path($file)->slurp;

my $dom = Mojo::DOM->new($html);

my $ref = $dom->find('#Tbl1 tr[bgcolor=#F0F0F0] > td > a')->map(
    sub {
        my $el = shift;
        if ( $el->text =~ /NC_\d+/ ) {
            return $el->text;
        }
        elsif ( $el->{href} =~ /id\=(\d+)/ ) {
            if ( $el->text !~ /not exist/ ) {
                return $1;
            }
        }
    }
);

#348535
#NC_014345
#484906
#NC_011395
#505693
#NC_014340
#400756
#NC_014287
#480043
#NC_004823
#160619
#NC_014267
#137071
#NC_017928

my @els = grep {$_} @{$ref};
tie my %acc_of, "Tie::IxHash";
{
    my $cur_id;
    while (@els) {
        my $el = shift @els;
        last unless $el;
        if ( $el =~ /^\d+$/ ) {
            $cur_id = $el;
        }
        elsif ( $el =~ /NC_\d+/ ) {
            if ( !exists $acc_of{$cur_id} ) {
                $acc_of{$cur_id} = [];
            }
            push @{ $acc_of{$cur_id} }, $el;
        }
    }
}

print "#strain_taxon_id,accession\n";
for my $key ( keys %acc_of ) {
    print join ",", ( $key, @{ $acc_of{$key} } );
    print "\n";
}
