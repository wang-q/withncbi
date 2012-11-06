#!/usr/bin/perl
use strict;
use warnings;

use Template;
use File::Basename;
use File::Find::Rule;
use File::Remove qw(remove);
use File::Spec;
use String::Compare;
use YAML qw(Dump Load DumpFile LoadFile);

my $store_dir = shift
    || File::Spec->catdir( $ENV{HOME}, "data/alignment/trichoderma" );
my $parallel = 12;

{    # on linux
    my $data_dir
        = File::Spec->catdir( $ENV{HOME}, "data/alignment/trichoderma" );
    my $pl_dir      = File::Spec->catdir( $ENV{HOME}, "Scripts" );
    my $kentbin_dir = File::Spec->catdir( $ENV{HOME}, "bin/x86_64" );

    # ensembl genomes 65
    my $fasta_dir
        = File::Spec->catdir( $ENV{HOME}, "data/alignment/trichoderma/WGS" );

    my $tt = Template->new;

    my @data = (
        {   taxon    => 452589,
            name     => "Tatr",
            sciname  => "Trichoderma atroviride",
            prefix   => "ABDG02",
            coverage => "8.26x Sanger",
        },
        {   taxon    => 431241,
            name     => "Tree",
            sciname  => "Trichoderma reesei",
            prefix   => "AAIL02",
            coverage => "9x Sanger",
        },
        {   taxon    => 413071,
            name     => "Tvir",
            sciname  => "Trichoderma virens",
            prefix   => "ABDF02",
            coverage => "8.05x Sanger",
        },
    );

    my @files_fasta = File::Find::Rule->file->name('*.fasta.gz')->in($fasta_dir);

    for my $item (@data) {

        # match the most similar name
        my ($fasta) = map { $_->[0] }
            sort { $b->[1] <=> $a->[1] }
            map { [ $_, compare( basename($_), $item->{prefix} ) ] }
            @files_fasta;
        $item->{fasta} = $fasta;

        # prepare working dir
        my $dir = File::Spec->catdir( $data_dir, $item->{name} );
        mkdir $dir if !-e $dir;
        $item->{dir} = $dir;
    }
    
    print Dump \@data;

    # taxon.csv
    my $text = <<'EOF';
[% FOREACH item IN data -%]
[% item.taxon %],[% item.sciname FILTER replace(' ', ',') %],[% item.name %],,
[% END -%]
EOF
    $tt->process(
        \$text,
        { data => \@data, },
        File::Spec->catfile( $store_dir, "taxon.csv" )
    ) or die Template->error;

    # chr_length.csv
    $text = <<'EOF';
[% FOREACH item IN data -%]
[% item.taxon %],chrUn,999999999,[% item.name %]/WGS
[% END -%]
EOF
    $tt->process(
        \$text,
        { data => \@data, },
        File::Spec->catfile( $store_dir, "chr_length.csv" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# unzip, filter and split
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %] 
echo [% item.name %]

cd [% item.dir %]
gzip -d -c [% item.fasta %] > toplevel.fa
perl -p -i -e '/>/ and s/\>gi\|(\d+).*/\>gi_$1/' toplevel.fa
[% kentbin_dir %]/faCount toplevel.fa | perl -aln -e 'next if $F[0] eq 'total'; print $F[0] if $F[1] > 100000; print $F[0] if $F[1] > 5000 and $F[6]/$F[1] < 0.05' | uniq > listFile
[% kentbin_dir %]/faSomeRecords toplevel.fa listFile toplevel.filtered.fa
[% kentbin_dir %]/faSplit byname toplevel.filtered.fa .
rm toplevel.fa toplevel.filtered.fa listFile

rename 's/fa$/fasta/' *.fa;
    
[% END -%]

EOF

    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "file.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# RepeatMasker
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
echo [% item.name %]

cd [% item.dir %]
RepeatMasker [% item.dir %]/*.fasta -species Fungi -xsmall --parallel [% parallel %]

[% END -%]

EOF

    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir,
            parallel    => $parallel,
        },
        File::Spec->catfile( $store_dir, "rm.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# RepeatMasker
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
echo [% item.name %]

for i in [% item.dir %]/*.fasta;
do
    if [ -f $i.masked ];
    then
        rename 's/fasta.masked$/fa/' $i.masked;
        find [% item.dir %] -type f -name "`basename $i`*" | xargs rm 
    fi;
done;

[% END -%]

EOF

    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir,
            parallel    => $parallel,
        },
        File::Spec->catfile( $store_dir, "clean-rm.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# blastz 
#----------------------------#
[% FOREACH i IN [ 0 .. data.max ] -%]
[% FOREACH j IN [ i .. data.max ] -%]
[% NEXT IF i == j -%]
# [% data.$i.name %] versus [% data.$j.name %]
perl [% pl_dir %]/blastz/bz.pl \
    -dt [% data_dir %]/[% data.$i.name %] \
    -dq [% data_dir %]/[% data.$j.name %] \
    -dl [% data_dir %]/[% data.$i.name %]vs[% data.$j.name %] \
    -s set01 -p [% parallel %] --noaxt -pb lastz --lastz

[% END -%]
[% END -%]

EOF
    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir,
            parallel    => $parallel,
        },
        File::Spec->catfile( $store_dir, "bz.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# lpcna 
#----------------------------#
[% FOREACH i IN [ 0 .. data.max ] -%]
[% FOREACH j IN [ i .. data.max ] -%]
[% NEXT IF i == j -%]
# [% data.$i.name %] versus [% data.$j.name %]
perl [% pl_dir %]/blastz/lpcna.pl \
    -dt [% data_dir %]/[% data.$i.name %] \
    -dq [% data_dir %]/[% data.$j.name %] \
    -dl [% data_dir %]/[% data.$i.name %]vs[% data.$j.name %] \
    -p [% parallel %]

[% END -%]
[% END -%]

EOF
    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir,
            parallel    => $parallel,
        },
        File::Spec->catfile( $store_dir, "lpcna.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# amp 
#----------------------------#
[% FOREACH i IN [ 0 .. data.max ] -%]
[% FOREACH j IN [ i .. data.max ] -%]
[% NEXT IF i == j -%]
# [% data.$i.name %] versus [% data.$j.name %]
perl [% pl_dir %]/blastz/amp.pl -syn \
    -dt [% data_dir %]/[% data.$i.name %] \
    -dq [% data_dir %]/[% data.$j.name %] \
    -dl [% data_dir %]/[% data.$i.name %]vs[% data.$j.name %] \
    -p [% parallel %]

[% END -%]
[% END -%]

EOF
    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir,
            parallel    => $parallel,
        },
        File::Spec->catfile( $store_dir, "amp.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# stat
#----------------------------#
[% FOREACH i IN [ 0 .. data.max ] -%]
[% FOREACH j IN [ i .. data.max ] -%]
[% NEXT IF i == j -%]
# [% data.$i.name %] versus [% data.$j.name %]
perl [% pl_dir %]/alignDB/extra/two_way_batch.pl \
    -d [% data.$i.name %]vs[% data.$j.name %] \
    -t "[% data.$i.taxon %],[% data.$i.name %]" \
    -q "[% data.$j.taxon %],[% data.$j.name %]" \
    -a [% data_dir %]/[% data.$i.name %]vs[% data.$j.name %] \
    -at 5000 -st 0 -ct 0 --parallel [% parallel %] --run 1-3,21,40

[% END -%]
[% END -%]


EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
            parallel => $parallel,
        },
        File::Spec->catfile( $store_dir, "pair_stat.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash

#----------------------------#
# only keeps chr.2bit files
#----------------------------#
# find [% data_dir %] -name "*.fa" | xargs rm

#----------------------------#
# clean pairwise maf
#----------------------------#
# find [% data_dir %] -name "mafSynNet" | xargs rm -fr
# find [% data_dir %] -name "mafNet" | xargs rm -fr

#----------------------------#
# gzip maf, fas
#----------------------------#
find [% data_dir %] -name "*.maf" | parallel gzip
find [% data_dir %] -name "*.maf.fas" | parallel gzip

#----------------------------#
# clean maf-fasta
#----------------------------#
# rm -fr [% data_dir %]/*_fasta

EOF
    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir,
        },
        File::Spec->catfile( $store_dir, "clean.sh" )
    ) or die Template->error;
}

#{    # multiz
#    my $data_dir
#        = File::Spec->catdir( $ENV{HOME}, "data/alignment/aspergillus" );
#    my $pl_dir = File::Spec->catdir( $ENV{HOME}, "Scripts" );
#
#    my $tt = Template->new;
#    my $strains_of
#        = { AfumvsVII => [qw{ Acla Afla Anid Anig Aory Ater Nfis }], };
#
#    my @data;
#    for my $key ( sort keys %{$strains_of} ) {
#        my @strains = @{ $strains_of->{$key} };
#        push @data,
#            {
#            out_dir => $key,
#            strains => \@strains,
#            };
#    }
#
#    my $text = <<'EOF';
##!/bin/bash
#
##----------------------------#
## mz
##----------------------------#
## find . -name "*MT.synNet*" | xargs rm
#
#[% FOREACH item IN data -%]
## [% item.out_dir %]
#perl [% pl_dir %]/blastz/mz.pl \
#    [% FOREACH st IN item.strains -%]
#    -d [% data_dir %]/Afumvs[% st %] \
#    [% END -%]
#    --tree [% data_dir %]/8way.nwk \
#    --out [% data_dir %]/[% item.out_dir %] \
#    -syn -p [% parallel %]
#
#[% END -%]
#
#EOF
#    $tt->process(
#        \$text,
#        {   data     => \@data,
#            data_dir => $data_dir,
#            pl_dir   => $pl_dir,
#            parallel => $parallel,
#        },
#        File::Spec->catfile( $store_dir, "mz.sh" )
#    ) or die Template->error;
#
#    $text = <<'EOF';
##----------------------------#
## maf2fasta
##----------------------------#
#[% FOREACH item IN data -%]
## [% item.out_dir %]
#perl [% pl_dir %]/alignDB/util/maf2fasta.pl \
#    --has_outgroup --id 330879 -p [% parallel %] --block \
#    -i [% data_dir %]/[% item.out_dir %] \
#    -o [% data_dir %]/[% item.out_dir %]_fasta
#
#[% END -%]
#
##----------------------------#
## mafft
##----------------------------#
#[% FOREACH item IN data -%]
## [% item.out_dir %]
#perl [% pl_dir %]/alignDB/util/refine_fasta.pl \
#    --msa mafft --block -p [% parallel %] \
#    -i [% data_dir %]/[% item.out_dir %]_fasta \
#    -o [% data_dir %]/[% item.out_dir %]_mft
#
#[% END -%]
#
##----------------------------#
## muscle-quick
##----------------------------#
##[% FOREACH item IN data -%]
### [% item.out_dir %]
##perl [% pl_dir %]/alignDB/util/refine_fasta.pl \
##    --msa muscle --quick --block -p [% parallel %] \
##    -i [% data_dir %]/[% item.out_dir %]_fasta \
##    -o [% data_dir %]/[% item.out_dir %]_mslq
##
##[% END -%]
#
#EOF
#
#    $tt->process(
#        \$text,
#        {   data     => \@data,
#            data_dir => $data_dir,
#            pl_dir   => $pl_dir,
#            parallel => $parallel,
#        },
#        File::Spec->catfile( $store_dir, "maf_fasta.sh" )
#    ) or die Template->error;
#
#    $text = <<'EOF';
##!/bin/bash
#
##----------------------------#
## multi_way_batch
##----------------------------#
#[% FOREACH item IN data -%]
## [% item.out_dir %]
## mafft
#perl [% pl_dir %]/alignDB/extra/multi_way_batch.pl \
#    -d [% item.out_dir %] -e Afum_65 \
#    --block --id 330879 \
#    -f [% data_dir %]/[% item.out_dir %]_mft  \
#    -lt 5000 -st 0 -ct 0 --parallel [% parallel %] --run 1-3,21,40
#
#[% END -%]
#
#EOF
#    $tt->process(
#        \$text,
#        {   data     => \@data,
#            data_dir => $data_dir,
#            pl_dir   => $pl_dir,
#            parallel => $parallel,
#        },
#        File::Spec->catfile( $store_dir, "multi.sh" )
#    ) or die Template->error;
#}
