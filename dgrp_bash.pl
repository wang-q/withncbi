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
    || File::Spec->catdir( $ENV{HOME}, "data/alignment/dgrp" );

{    # on linux
    my $data_dir    = File::Spec->catdir( $ENV{HOME}, "data/alignment/dgrp" );
    my $pl_dir      = File::Spec->catdir( $ENV{HOME}, "Scripts" );
    my $kentbin_dir = File::Spec->catdir( $ENV{HOME}, "bin/x86_64" );

    # nature 2012
    my $seq_dir = File::Spec->catdir( $ENV{HOME},
        "data/DGRP/projects/dgrp/freeze1_July_2010/sequences" );

    my $tt = Template->new;

    my @data = (
        { taxon => 900501, name => "Line138", coverage => 34, },
        { taxon => 900502, name => "Line176", coverage => 38.48, },
        { taxon => 900503, name => "Line181", coverage => 30.66, },
        { taxon => 900504, name => "Line208", coverage => 34.51, },
        { taxon => 900505, name => "Line321", coverage => 41.75, },
        { taxon => 900506, name => "Line332", coverage => 31.53, },
        { taxon => 900507, name => "Line375", coverage => 41.91, },
        { taxon => 900508, name => "Line38",  coverage => 34.1, },
        { taxon => 900509, name => "Line380", coverage => 36.73, },
        { taxon => 900510, name => "Line391", coverage => 47.62, },
        { taxon => 900511, name => "Line40",  coverage => 41.27, },
        { taxon => 900512, name => "Line406", coverage => 30.3, },
        { taxon => 900513, name => "Line443", coverage => 33.35, },
        { taxon => 900514, name => "Line517", coverage => 45.97, },
        { taxon => 900515, name => "Line57",  coverage => 38.28, },
        { taxon => 900516, name => "Line727", coverage => 33.64, },
        { taxon => 900517, name => "Line738", coverage => 31.74, },
        { taxon => 900518, name => "Line757", coverage => 32.87, },
        { taxon => 900519, name => "Line852", coverage => 40.42, },
        { taxon => 900520, name => "Line897", coverage => 32.81, },
    );

    my @files = File::Find::Rule->file->name('*.gz')->in($seq_dir);

    for my $item ( sort @data ) {
        my $name = $item->{name};

        # match the most similar name
        my $regex = $name . "_Chr";
        my @line_files = grep {/$regex/} @files;
        $item->{seqs} = \@line_files;

        # prepare working dir
        my $dir = File::Spec->catdir( $data_dir, $name );
        mkdir $dir if !-e $dir;
        $item->{dir} = $dir;
    }

    # taxon.csv
    my $text = <<'EOF';
[% FOREACH item IN data -%]
[% item.taxon %],Drosophila,melanogaster,[% item.name %],,
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
[% item.taxon %],chrUn,999999999,[% item.name %]/DGRP
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
if [ -f toplevel.fa ]; then rm toplevel.fa; fi;
[% FOREACH seq IN item.seqs -%]
gzip -d -c [% seq %] >> toplevel.fa
echo >> toplevel.fa
[% END -%]
perl -p -i -e '/>/ and s/\>[% item.name %]_/\>/' toplevel.fa
[% kentbin_dir %]/faCount toplevel.fa | perl -aln -e 'next if $F[0] eq 'total'; print $F[0] if $F[1] > 100000; print $F[0] if $F[1] > 10000  and $F[6]/$F[1] < 0.05' | uniq > listFile
[% kentbin_dir %]/faSomeRecords toplevel.fa listFile toplevel.filtered.fa
[% kentbin_dir %]/faSplit byname toplevel.filtered.fa .
rm toplevel.fa toplevel.filtered.fa listFile

~/perl5/bin/rename 's/fa$/fasta/' *.fa

[% END -%]

EOF

    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "auto_dgrp_file.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# repeatmasker on all fasta
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
bsub -q mpi_2 -n 8 -J [% item.name %]-rm RepeatMasker [% item.dir %]/*.fasta -species Flies -xsmall --parallel 8

[% END -%]

#----------------------------#
# find failed rm jobs
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
find ~/data/alignment/dgrp/[% item.name %] -name "*fasta" \
    | perl -e \
    'while(<>) {chomp; s/\.fasta$//; next if -e qq{$_.fasta.masked}; next if -e qq{$_.fa}; print qq{ bsub -n 8 -J [% item.name %]_ RepeatMasker $_.fasta -species Flies -xsmall --parallel 8 \n};}' >> catchup.txt

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
        File::Spec->catfile( $store_dir, "auto_dgrp_rm.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# blastz
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
bsub -q mpi_2 -n 8 -J [% item.name %]-bz perl [% pl_dir %]/blastz/bz.pl -dt [% data_dir %]/Dmel_65 -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/Dmelvs[% item.name %] -s set01 -p 8 --noaxt -pb lastz --lastz

[% END -%]

#----------------------------#
# lpcna
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
perl [% pl_dir %]/blastz/lpcna.pl -dt [% data_dir %]/Dmel_65 -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/Dmelvs[% item.name %] -p 8

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "auto_dgrp_bz.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
    
#----------------------------#
# amp
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
perl [% pl_dir %]/blastz/amp.pl -syn -dt [% data_dir %]/Dmel_65 -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/Dmelvs[% item.name %] -p 8

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data        => [ { name => "Dsim_65" }, @data ],
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "auto_dgrp_amp.sh" )
    ) or die Template->error;

    $text = <<'EOF';
cd [% data_dir %]

#----------------------------#
# stat
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
find [% data_dir %]/Dmelvs[% item.name %]/axtNet -name "*.axt.gz" | xargs gzip -d
perl [% pl_dir %]/alignDB/extra/two_way_batch.pl -d Dmelvs[% item.name %] -t="7227,Dmel" -q "[% item.taxon %],[% item.name %]" -a [% data_dir %]/Dmelvs[% item.name %] -at 10000 -st 10000000 --parallel 6 --run 1-3,21,40
gzip [% data_dir %]/Dmelvs[% item.name %]/axtNet/*.axt

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => [ { name => "Dsim_65", taxon => 7240 }, @data ],
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
        },
        File::Spec->catfile( $store_dir, "auto_dgrp_stat.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
    
#----------------------------#
# tar-gzip
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
cd [% data_dir %]/Dmelvs[% item.name %]/

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
find [% data_dir %] -name "*.maf" | xargs gzip
find [% data_dir %] -name "*.maf.fas" | xargs gzip

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
        File::Spec->catfile( $store_dir, "auto_dgrp_clean.sh" )
    ) or die Template->error;
}

{    # on linux
    my $data_dir = File::Spec->catdir( $ENV{HOME}, "data/alignment/dgrp" );
    my $pl_dir   = File::Spec->catdir( $ENV{HOME}, "Scripts" );

    my $tt = Template->new;

    my $strains_of = {
        Dmelvs21 => [
            qw{ Dsim_65 Line138 Line176 Line181 Line208 Line321 Line332 Line375
                Line38 Line380 Line391 Line40 Line406 Line443 Line517 Line57
                Line727 Line738 Line757 Line852 Line897 }
        ],
    };

    my @group;
    for my $dbname ( sort keys %{$strains_of} ) {
        my @strains = @{ $strains_of->{$dbname} };
        my $dbs     = join ',', map { "Dmelvs" . $_ } @strains;
        my $queries = join ',',
            map { $_ . "query" } ( 1 .. scalar @strains - 1 );
        push @group,
            {
            goal_db  => $dbname,
            outgroup => '0query',
            target   => '0target',
            dbs      => $dbs,
            queries  => $queries,
            };
    }

    my $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

[% FOREACH item IN data -%]
# [% item.goal_db %]
perl [% pl_dir %]/alignDB/extra/join_dbs.pl --crude_only \
    --dbs [% item.dbs %] \
    --goal_db [% item.goal_db %] --outgroup [% item.outgroup %] \
    --target [% item.target %] \
    --queries [% item.queries %] \
    --no_insert=1 --trimmed_fasta=1 --length 10000

perl [% pl_dir %]/alignDB/util/refine_fasta.pl \
    --msa mafft -p 4 \
    -i [% data_dir %]/[% item.goal_db %].crude \
    -o [% data_dir %]/[% item.goal_db %]_mft

perl [% pl_dir %]/tool/catfasta2phyml.pl -f [% data_dir %]/[% item.goal_db %]_mft/*.fas > [% data_dir %]/all.fasta

perl [% pl_dir %]/alignDB/extra/multi_way_batch.pl -d [% item.goal_db %] -e fly_65 -f [% data_dir %]/[% item.goal_db %]_mft  -lt 1000 -st 10000000 --parallel 8 --run all

[% END -%]
EOF

    $tt->process(
        \$text,
        { data => \@group, data_dir => $data_dir, pl_dir => $pl_dir, },
        File::Spec->catfile( $store_dir, "auto_dgrp_joins.sh" )
    ) or die Template->error;
}

{    # multiz
    my $data_dir = File::Spec->catdir( $ENV{HOME}, "data/alignment/dgrp" );
    my $pl_dir   = File::Spec->catdir( $ENV{HOME}, "Scripts" );

    my $tt         = Template->new;
    my $strains_of = {
        DmelvsVIIGE40x =>
            [ qw{ Dsim_65 Line321 Line375 Line391 Line40 Line517 Line852 } ],
        DmelvsXXI => [
            qw{ Dsim_65 Line138 Line176 Line181 Line208 Line321 Line332 Line375
                Line38 Line380 Line391 Line40 Line406 Line443 Line517 Line57
                Line727 Line738 Line757 Line852 Line897 }
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
[% FOREACH item IN data -%]
# [% item.out_dir %]
bsub -q mpi_2 -n 8 -J [% item.out_dir %]-mz perl [% pl_dir %]/blastz/mz.pl \
    [% FOREACH st IN item.strains -%]
    -d [% data_dir %]/Dmelvs[% st FILTER ucfirst %] \
    [% END -%]
    --tree [% data_dir %]/22way.nwk \
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
        File::Spec->catfile( $store_dir, "auto_dgrp_mz.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#----------------------------#
# maf2fasta
#----------------------------#
[% FOREACH item IN data -%]
# [% item.out_dir %]
perl [% pl_dir %]/alignDB/util/maf2fasta.pl \
    --has_outgroup --id 7227 -p 8 --block \
    -i [% data_dir %]/[% item.out_dir %] \
    -o [% data_dir %]/[% item.out_dir %]_fasta

[% END -%]

#----------------------------#
# mafft
#----------------------------#
[% FOREACH item IN data -%]
# [% item.out_dir %]
bsub -q mpi_2 -n 8 -J [% item.out_dir %]-mft perl [% pl_dir %]/alignDB/util/refine_fasta.pl \
    --msa mafft --block -p 8 \
    -i [% data_dir %]/[% item.out_dir %]_fasta \
    -o [% data_dir %]/[% item.out_dir %]_mft

[% END -%]

#----------------------------#
# muscle
#----------------------------#
#[% FOREACH item IN data -%]
## [% item.out_dir %]
#bsub -q mpi_2 -n 8 -J [% item.out_dir %]-msl perl [% pl_dir %]/alignDB/util/refine_fasta.pl \
#    --msa muscle --block -p 8 \
#    -i [% data_dir %]/[% item.out_dir %]_fasta \
#    -o [% data_dir %]/[% item.out_dir %]_msl
#
#[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
        },
        File::Spec->catfile( $store_dir, "auto_dgrp_maf_fasta.sh" )
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
    -d [% item.out_dir %] -e fly_65 \
    --block --id 7227 \
    -f [% data_dir %]/[% item.out_dir %]_mft  \
    -lt 5000 -st 1000000 --parallel 4 --run 1-3,21,40

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
        },
        File::Spec->catfile( $store_dir, "auto_dgrp_multi.sh" )
    ) or die Template->error;
}
