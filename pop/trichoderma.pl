#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Template;
use File::Basename;
use File::Find::Rule;
use File::Remove qw(remove);
use File::Spec;
use String::Compare;
use YAML qw(Dump Load DumpFile LoadFile);

my $parallel = 8;

my $group_name = 'trichoderma';
my $base_dir   = File::Spec->catdir( $ENV{HOME}, "data/alignment" );
my $data_dir   = File::Spec->catdir( $base_dir, $group_name );

#my $phylo_tree = File::Spec->catfile( $data_dir, "primates_13way.nwk" );
my $pl_dir = File::Spec->catdir( $ENV{HOME}, "Scripts" );

# NCBI WGS
my $fasta_dir
    = File::Spec->catdir( $ENV{HOME}, "data/alignment/trichoderma/WGS" );

my $tt = Template->new;

my @data = (
    {   taxon    => 452589,
        name     => "Tatr",
        sciname  => "Trichoderma atroviride IMI 206040",
        prefix   => "ABDG02",
        coverage => "8.26x Sanger",
    },
    {   taxon    => 5544,
        name     => "Thar",
        sciname  => "Trichoderma harzianum",
        prefix   => "JNNP01",
        coverage => "20.0x Illumina HiSeq",
    },
    {   taxon    => 1234776,
        name     => "Tpse",
        sciname  => "Trichoderma longibrachiatum SMF2",
        prefix   => "ANBJ01",
        coverage => "69x 454; Illumina HiSeq",
    },
    {   taxon    => 431241,
        name     => "Tree_QM6a",
        sciname  => "Trichoderma reesei QM6a",
        prefix   => "AAIL02",
        coverage => "9x Sanger",
    },
    {   taxon    => 1344414,
        name     => "Tree_RUT_C_30",
        sciname  => "Trichoderma reesei RUT C-30",
        prefix   => "JABP01",
        coverage => "47.6x Illumina",
    },
    {   taxon    => 1331945,
        name     => "Tvir_FT_333",
        sciname  => "Trichoderma virens FT-333",
        prefix   => "JTGJ01",
        coverage => "51.0x SOLiD",
    },
    {   taxon    => 413071,
        name     => "Tvir_Gv29_8",
        sciname  => "Trichoderma virens Gv29-8",
        prefix   => "ABDF02",
        coverage => "8.05x Sanger",
    },

    # contigs are too short
    #{   taxon    => 1247866,
    #    name     => "Tham",
    #    sciname  => "Trichoderma hamatum GD12",
    #    prefix   => "ANCB01",
    #    coverage => "40.0x Illumina HiSeq",
    #},
);

my @subdirs_fasta = File::Find::Rule->file->name('*.fsa_nt.gz')->in($fasta_dir);

for my $item (@data) {

    # match the most similar name
    my ($fasta) = map { $_->[0] }
        sort { $b->[1] <=> $a->[1] }
        map { [ $_, compare( basename($_), $item->{prefix} ) ] } @subdirs_fasta;
    $item->{fasta} = $fasta;

    # prepare working dir
    my $dir = File::Spec->catdir( $data_dir, $item->{name} );
    mkdir $dir if !-e $dir;
    $item->{dir} = $dir;
}

my $text;

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
faops count toplevel.fa | perl -aln -e 'next if $F[0] eq 'total'; print $F[0] if $F[1] > 100000; print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05' | uniq > listFile
faops some toplevel.fa listFile toplevel.filtered.fa
[% IF item.name == 'Tatr' or item.name == 'Tvir_Gv29_8' or item.name == 'Tree_QM6a' -%]
faops split-name toplevel.filtered.fa .
[% ELSE -%]
faops split-about toplevel.filtered.fa 10000000 .
[% END -%]
rm toplevel.fa toplevel.filtered.fa listFile

rename 's/fa$/fasta/' *.fa;
    
[% END -%]

EOF

$tt->process(
    \$text,
    {   data     => \@data,
        data_dir => $data_dir,
        pl_dir   => $pl_dir,
    },
    File::Spec->catfile( $data_dir, "01_file.sh" )
) or die Template->error;

$text = <<'EOF';
#!/bin/bash

#----------------------------------------------------------#
# RepeatMasker
#----------------------------------------------------------#
cd [% data_dir %]
echo Doing RepeatMasker

[% FOREACH item IN data -%]
#----------------------------#
# [% item.name %] [% item.coverage %]
#----------------------------#
echo [% item.name %]

cd [% item.dir %]
RepeatMasker [% item.dir %]/*.fasta -species Fungi -xsmall --parallel [% parallel %]

[% END -%]

#----------------------------------------------------------#
# Clean RepeatMasker
#----------------------------------------------------------#
cd [% data_dir %]
echo Cleaning RepeatMasker

[% FOREACH item IN data -%]
#----------------------------#
# [% item.name %] [% item.coverage %]
#----------------------------#
echo [% item.name %]

cd [% item.dir %]
for i in *.fasta;
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
    {   data     => \@data,
        data_dir => $data_dir,
        pl_dir   => $pl_dir,
        parallel => $parallel,
    },
    File::Spec->catfile( $data_dir, "02_rm.sh" )
) or die Template->error;

$text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# generate taxon file
#----------------------------#
perl [% pl_dir %]/withncbi/taxon/strain_info.pl \
[% FOREACH item IN data -%]
    --id   [% item.taxon %] \
    --name [% item.taxon %]=[% item.name %] \
[% END -%]
    --file [% data_dir %]/[% group_name %].csv

#----------------------------#
# multi genome alignment plan
#----------------------------#
# don't copy sequences (RepeatMasker done)
# Execute the following lines by copy & paste.

# Trichoderma_7way
cd [% data_dir %]
perl [% pl_dir %]/withncbi/taxon/strain_bz.pl \
    --file [% data_dir %]/[% group_name %].csv \
    -w     [% base_dir %] \
    --name [% group_name %] \
    --multi_name Trichoderma_7way \
    --use_name \
    --parallel [% parallel %]\
    --norm \
    -t Tatr \
    -q Thar \
    -q Tpse \
    -q Tree_QM6a \
    -q Tree_RUT_C_30 \
    -q Tvir_FT_333 \
    -q Tvir_Gv29_8

EOF

$tt->process(
    \$text,
    {   data       => \@data,
        group_name => $group_name,
        data_dir   => $data_dir,
        base_dir   => $base_dir,
        pl_dir     => $pl_dir,
        parallel   => $parallel,
    },
    File::Spec->catfile( $data_dir, "03_prepare.sh" )
) or die Template->error;

__END__

# create withncbi/doc/trichoderma.tsv manually

mkdir -p ~/data/alignment/trichoderma
cd ~/data/alignment/trichoderma

perl ~/Scripts/withncbi/util/wgs_prep.pl \
    -f ~/Scripts/withncbi/doc/trichoderma.tsv \
    -o WGS \
    -a 

aria2c -x 6 -s 3 -c -i WGS/trichoderma.url.txt

find WGS -name "*.gz" | xargs gzip -t 

# edit ~/Scripts/withncbi/pop/trichoderma.pl, add contents from trichoderma.data.txt

perl ~/Scripts/withncbi/pop/trichoderma.pl
sh 01_file.sh
sh 02_rm.sh

# execute 03_prepare.sh by copy & paste  

# for each multi_name, execute the following bash file
sh 1_real_chr.sh
sh 3_pair_cmd.sh
sh 4_rawphylo.sh
sh 5_multi_cmd.sh
sh 6_multi_db_only.sh
