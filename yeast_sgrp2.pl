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

#----------------------------------------------------------#
# RepeatMasker has been done
#----------------------------------------------------------#
my $store_dir = shift
    || File::Spec->catdir( $ENV{HOME}, "data/alignment/yeast_sgrp2" );

{    # on linux
    my $data_dir
        = File::Spec->catdir( $ENV{HOME}, "data/alignment/yeast_sgrp2" );
    my $pl_dir      = File::Spec->catdir( $ENV{HOME}, "Scripts" );
    my $kentbin_dir = File::Spec->catdir( $ENV{HOME}, "bin/x86_64" );

    # sgrp2 SGA assembly
    my $seq_dir = File::Spec->catdir( $ENV{HOME}, "data/SGRP2/assembly/" );

    my $tt   = Template->new;
    my @data = (

        #{ taxon => 226125, name => 'Spar',   coverage => '7x', },
        #{ taxon => 285006, name => 'RM11',   coverage => '10x', },
        #{ taxon => 307796, name => 'YJM789', coverage => '10x', },

        # sgrp2
        { taxon => 901501, name => "YPS128",    coverage => "56x illumina", },
        { taxon => 901502, name => "DBVPG1106", coverage => "34x illumina", },
        { taxon => 901503, name => "SK1",       coverage => "40x illumina", },
        { taxon => 901504, name => "L1528",     coverage => "34x illumina", },
        { taxon => 901505, name => "Y12",       coverage => "29x illumina", },
        { taxon => 901506, name => "UWOPS83",   coverage => "638x illumina", },
        { taxon => 901507, name => "UWOPS87",   coverage => "834x illumina", },
        { taxon => 901508, name => "UWOPS03",   coverage => "24x illumina", },
        { taxon => 901509, name => "W303",      coverage => "33x illumina", },
        { taxon => 901510, name => "YJM975",    coverage => "22x illumina", },
        { taxon => 901511, name => "DBVPG6765", coverage => "37x illumina", },
        { taxon => 901512, name => "DBVPG6044", coverage => "32x illumina", },
        { taxon => 901513, name => "DBVPG1373", coverage => "26x illumina", },
        { taxon => 901514, name => "Y55",       coverage => "28x illumina", },
    );

    my @files = File::Find::Rule->file->name('*.fa')->in($seq_dir);

    for my $item ( sort @data ) {
        my $name = $item->{name};

        # match the most similar name
        my ($file) = map { $_->[0] }
            sort { $b->[1] <=> $a->[1] }
            map {
            [   $_,
                compare(
                    lc basename( $_, ".scaffolds.fa" ),
                    lc( "sc_" . $item->{name} )
                )
            ]
            }
            grep {/$item->{name}/i} @files;
        $item->{seq} = $file;

        # prepare working dir
        my $dir = File::Spec->catdir( $data_dir, $name );
        mkdir $dir if !-e $dir;
        $item->{dir} = $dir;
    }

    push @data,
        {
        taxon    => 226125,
        name     => 'Spar',
        coverage => '7x',
        dir      => File::Spec->catdir( $data_dir, 'Spar' ),
        seq      => File::Spec->catfile( $data_dir, 'Spar', 'Spar.fa' ),
        };

    print Dump \@data;

    # taxon.csv
    my $text = <<'EOF';
[% FOREACH item IN data -%]
[% item.taxon %],Saccharomyces,cerevisiae,[% item.name %],,
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
[% item.taxon %],chrUn,999999999,[% item.name %]/sgrp2
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
find [% item.dir %] -name "*.fa" -o -name "*.fasta" | xargs rm
[% kentbin_dir %]/faCount [% item.seq %] | perl -aln -e 'next if $F[0] eq 'total'; print $F[0] if $F[1] > 10000; print $F[0] if $F[1] > 1000  and $F[6]/$F[1] < 0.05' | uniq > listFile
[% kentbin_dir %]/faSomeRecords [% item.seq %] listFile [% item.name %].fasta
rm listFile

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
# repeatmasker on all fasta
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
RepeatMasker [% item.dir %]/*.fasta -species "yeast" -xsmall --parallel 8

[% END -%]

# find [% data_dir %] -name "*.fasta.masked" | sed "s/\.fasta\.masked$//" | xargs -i echo mv {}.fasta.masked {}.fa | sh

EOF

    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "rm.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# blastz
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
perl [% pl_dir %]/blastz/bz.pl -dt [% data_dir %]/S288C -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/S288Cvs[% item.name %] -s set01 -p 8 --noaxt -pb lastz --lastz

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "bz.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# lpcna
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
perl [% pl_dir %]/blastz/lpcna.pl -dt [% data_dir %]/S288C -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/S288Cvs[% item.name %] -p 8

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "lpcna.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
    
#----------------------------#
# amp
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
perl [% pl_dir %]/blastz/amp.pl -syn -dt [% data_dir %]/S288C -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/S288Cvs[% item.name %] -p 8

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data        => [@data],
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "amp.sh" )
    ) or die Template->error;

    $text = <<'EOF';
cd [% data_dir %]

#----------------------------#
# stat
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
find [% data_dir %]/S288Cvs[% item.name %]/axtNet -name "*.axt.gz" | xargs gzip -d
perl [% pl_dir %]/alignDB/extra/two_way_batch.pl -d S288Cvs[% item.name %] -t="4932,S288C" -q "[% item.taxon %],[% item.name %]" -a [% data_dir %]/S288Cvs[% item.name %] -lt 10000 -st 0 --parallel 8 --run 1-3,21,40
gzip [% data_dir %]/S288Cvs[% item.name %]/axtNet/*.axt

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => [@data],
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
        },
        File::Spec->catfile( $store_dir, "pair_stat.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
    
#----------------------------#
# tar-gzip
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
cd [% data_dir %]/S288Cvs[% item.name %]/

tar -czvf lav.tar.gz   [*.lav   --remove-files
tar -czvf psl.tar.gz   [*.psl   --remove-files
tar -czvf chain.tar.gz [*.chain --remove-files
gzip *.chain
gzip net/*
gzip axtNet/*.axt

[% END -%]

#----------------------------#
# clean RepeatMasker outputs
#----------------------------#
# find [% data_dir %] -name "*.fasta*" | xargs rm

#----------------------------#
# only keeps chr.2bit files
#----------------------------#
# find [% data_dir %] -name "*.fa" | xargs rm

#----------------------------#
# clean pairwise maf
#----------------------------#
find [% data_dir %] -name "mafSynNet" | xargs rm -fr
find [% data_dir %] -name "mafNet" | xargs rm -fr

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
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "clean.sh" )
    ) or die Template->error;
}

{    # multiz
    my $data_dir
        = File::Spec->catdir( $ENV{HOME}, "data/alignment/yeast_sgrp2" );
    my $pl_dir = File::Spec->catdir( $ENV{HOME}, "Scripts" );

    my $tt         = Template->new;
    my $strains_of = {

        # sum-align > 10M
        S288CvsXIS2 => [
            qw{ Spar DBVPG1373 DBVPG6044 DBVPG6765 L1528 SK1 UWOPS83 UWOPS87
                W303 Y55 YPS128 }
        ],

        ## sum-align > 9M
        #S288CvsXIIIS2 => [
        #    qw{ Spar DBVPG1373 DBVPG6044 DBVPG6765 L1528 SK1 UWOPS83 UWOPS87
        #    W303 Y55 YPS128 DBVPG1106 Y12 }
        #],

        # sum-align > 3M
        S288CvsXVS2 => [
            qw{ Spar DBVPG1373 DBVPG6044 DBVPG6765 L1528 SK1 UWOPS83 UWOPS87
                W303 Y55 YPS128 DBVPG1106 Y12 UWOPS03 YJM975}
        ],
    };

    my @data;
    for my $key ( sort keys %{$strains_of} ) {
        my @strains = @{ $strains_of->{$key} };
        push @data,
            {
            out_dir => $key,
            strains => \@strains,
            };
    }

    my $text = <<'EOF';
#!/bin/bash

#----------------------------#
# mz
#----------------------------#
# find [% data_dir %] -name "chrMito*maf" | xargs rm

[% FOREACH item IN data -%]
# [% item.out_dir %]
perl [% pl_dir %]/blastz/mz.pl \
    [% FOREACH st IN item.strains -%]
    -d [% data_dir %]/S288Cvs[% st %] \
    [% END -%]
    --tree [% data_dir %]/15way.nwk \
    --out [% data_dir %]/[% item.out_dir %] \
    -syn -p 8

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
        },
        File::Spec->catfile( $store_dir, "mz.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#----------------------------#
# maf2fasta
#----------------------------#
[% FOREACH item IN data -%]
# [% item.out_dir %]
perl [% pl_dir %]/blastz/maf2fasta.pl \
    --has_outgroup --id 4932 -p 8 --block \
    -i [% data_dir %]/[% item.out_dir %] \
    -o [% data_dir %]/[% item.out_dir %]_fasta

[% END -%]

#----------------------------#
# mafft
#----------------------------#
[% FOREACH item IN data -%]
# [% item.out_dir %]
perl [% pl_dir %]/blastz/refine_fasta.pl \
    --msa mafft --block -p 8 \
    -i [% data_dir %]/[% item.out_dir %]_fasta \
    -o [% data_dir %]/[% item.out_dir %]_mafft

[% END -%]

#----------------------------#
# muscle-quick
#----------------------------#
[% FOREACH item IN data -%]
# [% item.out_dir %]
perl [% pl_dir %]/blastz/refine_fasta.pl \
    --msa muscle --quick --block -p 8 \
    -i [% data_dir %]/[% item.out_dir %]_fasta \
    -o [% data_dir %]/[% item.out_dir %]_muscle

[% END -%]

#----------------------------#
# clean
#----------------------------#
[% FOREACH item IN data -%]
# [% item.out_dir %]
cd [% data_dir %]
rm -fr [% item.out_dir %]_fasta

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
        },
        File::Spec->catfile( $store_dir, "maf_fasta.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash

#----------------------------#
# multi_way_batch
#----------------------------#
[% FOREACH item IN data -%]
# [% item.out_dir %]
# mafft
perl [% pl_dir %]/alignDB/extra/multi_way_batch.pl \
    -d [% item.out_dir %] -e yeast_65 \
    --block --id 4932 \
    -f [% data_dir %]/[% item.out_dir %]_mafft  \
    -lt 5000 -st 1000000 --parallel 8 --run 1-3,21,40

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
        },
        File::Spec->catfile( $store_dir, "multi.sh" )
    ) or die Template->error;

}
