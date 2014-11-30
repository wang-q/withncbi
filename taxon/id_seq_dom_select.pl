#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use File::Slurp;
use Mojo::DOM;
use YAML qw(Dump Load DumpFile LoadFile);

# 'http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=2759&opt=plastid';
# 'http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=33090&opt=organelle';

my $file = shift;
die "Provide a valid html file!\n" unless $file;

my $html = read_file($file);

my $dom = Mojo::DOM->new($html);

my $string = $dom->find('#Tbl1 tr[bgcolor=#F0F0F0] > td > a')->map(
    sub {
        my $el = shift;
        if ( $el->text =~ /NC_\d+/ ) {
            return $el->text;
        }
        elsif ( $el->{href} =~ /id\=(\d+)/ ) {
            return $1;
        }
        else {
            return $el;
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

my $acc_of = {};
{
    my @els = split /\n/, $string;
    my $cur_id;
    while (@els) {
        my $el = shift @els;
        last unless $el;
        if ( $el =~ /^\d+$/ ) {
            $cur_id = $el;
        }
        elsif ( $el =~ /NC_\d+/ ) {
            if ( !exists $acc_of->{$cur_id} ) {
                $acc_of->{$cur_id} = [];
            }
            push @{ $acc_of->{$cur_id} }, $el;
        }
    }
}

print "strain_taxon_id,accession\n";
for my $key ( sort keys %{$acc_of} ) {
    print join ",", ( $key, @{ $acc_of->{$key} } );
    print "\n";
}

