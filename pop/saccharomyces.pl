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

my $group_name = 'saccharomyces';
my $base_dir   = File::Spec->catdir( $ENV{HOME}, "data/alignment" );
my $data_dir   = File::Spec->catdir( $base_dir, $group_name );

#my $phylo_tree = File::Spec->catfile( $data_dir, "primates_13way.nwk" );
my $pl_dir = File::Spec->catdir( $ENV{HOME}, "Scripts" );

# NCBI WGS
my $fasta_dir
    = File::Spec->catdir( $ENV{HOME}, "data/alignment/saccharomyces/WGS" );

my $tt   = Template->new;
my @data = (
    {   taxon    => 1160507,
        name     => "Sarb",
        sciname  => "Saccharomyces arboricola H-6",
        prefix   => "ALIE01",
        coverage => "50.0x 454; SOLiD",
    },
    {   taxon    => 226231,
        name     => "Sbay",
        sciname  => "Saccharomyces bayanus 623-6C",
        prefix   => "AACG02",
        coverage => " ",
    },
    {   taxon    => 252598,
        name     => "Sbou",
        sciname  => "Saccharomyces sp. 'boulardii'",
        prefix   => "JPJH02",
        coverage => "308x SOLiD; Illumina MiSeq",
    },
    {   taxon    => 1073566,
        name     => "Scar",
        sciname  => "Saccharomyces carlsbergensis CBS 1513",
        prefix   => "AZCJ01",
        coverage => "18.0x 454",
    },
    {   taxon    => 1163645,
        name     => "Skud",
        sciname  => "Saccharomyces kudriavzevii FM1066",
        prefix   => "AJHS01",
        coverage => "8x Illumina GA",
    },
    {   taxon    => 226126,
        name     => "Smik",
        sciname  => "Saccharomyces mikatae IFO 1815",
        prefix   => "AABZ01",
        coverage => " ",
    },
    {   taxon    => 226125,
        name     => "Spar",
        sciname  => "Saccharomyces paradoxus NRRL Y-17217",
        prefix   => "AABY01",
        coverage => " ",
    },
    {   taxon    => 520522,
        name     => "Spas",
        sciname  => "Saccharomyces pastorianus Weihenstephan 34/70",
        prefix   => "AZAA01",
        coverage => "18x Illumina MiSeq",
    },
    {   taxon    => 226127,
        name     => "Suva",
        sciname  => "Saccharomyces uvarum MCYC 623",
        prefix   => "AACA01",
        coverage => " ",
    },
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

# Add pre-masked S288c here
unshift @data,
    {
    taxon   => 559292,
    name    => "Scer",
    sciname => "Saccharomyces cerevisiae S288c",
    };

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

# Saccharomyces_10way
cd [% data_dir %]
perl [% pl_dir %]/withncbi/taxon/strain_bz.pl \
    --file [% data_dir %]/[% group_name %].csv \
    -w     [% base_dir %] \
    --name [% group_name %] \
    --multi_name Saccharomyces_10way \
    --use_name \
    --parallel [% parallel %] \
    --norm \
    -t Scer \
    -q Sarb \
    -q Sbay \
    -q Sbou \
    -q Scar \
    -q Skud \
    -q Smik \
    -q Spar \
    -q Spas \
    -q Suva

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

# create withncbi/doc/saccharomyces.tsv manually

mkdir -p ~/data/alignment/saccharomyces
cd ~/data/alignment/saccharomyces

perl ~/Scripts/withncbi/util/wgs_prep.pl \
    -f ~/Scripts/withncbi/doc/saccharomyces.tsv \
    -o WGS \
    -a 

aria2c -x 6 -s 3 -c -i WGS/saccharomyces.url.txt

find WGS -name "*.gz" | xargs gzip -t 

# edit ~/Scripts/withncbi/pop/saccharomyces.pl, add contents from saccharomyces.data.txt

perl ~/Scripts/withncbi/pop/saccharomyces.pl
sh 01_file.sh
sh 02_rm.sh

# copy S288c sequences
mkdir Scer
cp ~/data/alignment/yeast_genome/S288c/* Scer

# execute 03_prepare.sh by copy & paste  

# for each multi_name, execute the following bash file
sh 1_real_chr.sh
sh 3_pair_cmd.sh
sh 4_rawphylo.sh
sh 5_multi_cmd.sh
sh 6_multi_db_only.sh
