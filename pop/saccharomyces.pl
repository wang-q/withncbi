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
my $pl_dir     = File::Spec->catdir( $ENV{HOME}, "Scripts" );

# NCBI WGS
my $fasta_dir
    = File::Spec->catdir( $ENV{HOME}, "data/alignment/saccharomyces/WGS" );

my @data = (
    {   taxon    => 1160507,
        name     => "Sarb_H_6",
        sciname  => "Saccharomyces arboricola H-6",
        prefix   => "ALIE01",
        coverage => "50.0x 454; SOLiD",
    },
    {   taxon    => 226231,
        name     => "Sbay_623_6C",
        sciname  => "Saccharomyces bayanus 623-6C",
        prefix   => "AACG02",
        coverage => " ",
    },
    {   taxon       => 100000001,
        name        => "Sbou_17",
        sciname     => "Saccharomyces sp. 'boulardii' 17",
        prefix      => "JPJH02",
        coverage    => "308x SOLiD; Illumina MiSeq",
        original_id => 252598,
    },
    {   taxon       => 100000002,
        name        => "Sbou_ATCC_MYA_796",
        sciname     => "Saccharomyces sp. 'boulardii' ATCC MYA-796",
        prefix      => "JRHY01",
        coverage    => "403.0x Illumina HiSeq",
        original_id => 252598,
    },
    {   taxon    => 1343043,
        name     => "Sbou_EDRL",
        sciname  => "Saccharomyces sp. EDRL",
        prefix   => "ATCS01",
        coverage => "50.0x 454",
    },
    {   taxon    => 1073566,
        name     => "Scar_CBS_1513",
        sciname  => "Saccharomyces carlsbergensis CBS 1513",
        prefix   => "AZCJ01",
        coverage => "18.0x 454",
    },
    {   taxon   => 1095631,
        name    => "ScerSkud_VIN7",
        sciname => "Saccharomyces cerevisiae x Saccharomyces kudriavzevii VIN7",
        prefix  => "AGVY01",
        coverage => "20x 454 GS FLX titanium",
    },
    {   taxon    => 1163633,
        name     => "Skud_FM1056",
        sciname  => "Saccharomyces kudriavzevii FM1056",
        prefix   => "AJHP01",
        coverage => "9x Illumina GA",
    },
    {   taxon    => 226230,
        name     => "Skud_IFO_1802",
        sciname  => "Saccharomyces kudriavzevii IFO 1802",
        prefix   => "AACI03",
        coverage => "3.4x Sanger",
    },
    {   taxon    => 1165456,
        name     => "Skud_ZP591",
        sciname  => "Saccharomyces kudriavzevii ZP591",
        prefix   => "AJIH01",
        coverage => "36x Illumina GA",
    },
    {   taxon    => 226126,
        name     => "Smik_IFO_1815_1",
        sciname  => "Saccharomyces mikatae IFO 1815",
        prefix   => "AABZ01",
        coverage => " ",
    },

    # same strain and same taxon_id, keep one based on filtered sequence length
    #{   taxon    => 226126,
    #    name     => "Smik_IFO_1815_2",
    #    sciname  => "Saccharomyces mikatae IFO 1815",
    #    prefix   => "AACH01",
    #    coverage => " ",
    #},
    {   taxon    => 226125,
        name     => "Spar_NRRL_Y_17217",
        sciname  => "Saccharomyces paradoxus NRRL Y-17217",
        prefix   => "AABY01",
        coverage => " ",
    },
    {   taxon       => 100000003,
        name        => "Spas_CBS_1483",
        sciname     => "Saccharomyces pastorianus CBS 1483",
        prefix      => "JTFI01",
        coverage    => "52.2x Illumina HiSeq",
        original_id => 27292,
    },
    {   taxon    => 1214527,
        name     => "Spas_CCY48_91",
        sciname  => "Saccharomyces pastorianus CCY48 - 91",
        prefix   => "ALJS01",
        coverage => "30.0x 454 FLX",
    },

    # same strain and same taxon_id, keep one based on filtered sequence length
    #{   taxon    => 520522,
    #    name     => "Spas_Weihenstephan_34_70_1",
    #    sciname  => "Saccharomyces pastorianus Weihenstephan 34/70",
    #    prefix   => "ABPO01",
    #    coverage => " ",
    #},
    {   taxon    => 520522,
        name     => "Spas_Weihenstephan_34_70_2",
        sciname  => "Saccharomyces pastorianus Weihenstephan 34/70",
        prefix   => "AZAA01",
        coverage => "18x Illumina MiSeq",
    },
    {   taxon       => 100000004,
        name        => "Sunv_A9",
        sciname     => "Saccharomyces uvarum A9",
        prefix      => "JNVO01",
        coverage    => "30.0x IonTorrent",
        original_id => 230603,
    },
    {   taxon    => 226127,
        name     => "Suva_MCYC_623",
        sciname  => "Saccharomyces uvarum MCYC 623",
        prefix   => "AACA01",
        coverage => " ",
    },

    # short contigs
    #{   taxon    => 1163634,
    #    name     => "Skud_FM1057",
    #    sciname  => "Saccharomyces kudriavzevii FM1057",
    #    prefix   => "AJHQ01",
    #    coverage => "9x Illumina GA",
    #},
    #{   taxon    => 1163641,
    #    name     => "Skud_FM1062",
    #    sciname  => "Saccharomyces kudriavzevii FM1062",
    #    prefix   => "AJHR01",
    #    coverage => "8x Illumina GA",
    #},
    #{   taxon    => 1163645,
    #    name     => "Skud_FM1066",
    #    sciname  => "Saccharomyces kudriavzevii FM1066",
    #    prefix   => "AJHS01",
    #    coverage => "8x Illumina GA",
    #},
    #{   taxon    => 1163632,
    #    name     => "Skud_FM1069",
    #    sciname  => "Saccharomyces kudriavzevii FM1069",
    #    prefix   => "AJHT01",
    #    coverage => "20.2897236279x Illumina GA",
    #},
    #{   taxon    => 1165457,
    #    name     => "Skud_FM1071",
    #    sciname  => "Saccharomyces kudriavzevii FM1071",
    #    prefix   => "AJHU01",
    #    coverage => "13x Illumina GA",
    #},
    #{   taxon    => 1163638,
    #    name     => "Skud_FM1072",
    #    sciname  => "Saccharomyces kudriavzevii FM1072",
    #    prefix   => "AJHV01",
    #    coverage => "9x Illumina GA",
    #},
    #{   taxon    => 1163639,
    #    name     => "Skud_FM1073",
    #    sciname  => "Saccharomyces kudriavzevii FM1073",
    #    prefix   => "AJHW01",
    #    coverage => "11x Illumina GA",
    #},
    #{   taxon    => 1163640,
    #    name     => "Skud_FM1074",
    #    sciname  => "Saccharomyces kudriavzevii FM1074",
    #    prefix   => "AJHX01",
    #    coverage => "16x Illumina GA",
    #},
    #{   taxon    => 1163642,
    #    name     => "Skud_FM1075",
    #    sciname  => "Saccharomyces kudriavzevii FM1075",
    #    prefix   => "AJHY01",
    #    coverage => "9x Illumina GA",
    #},
    #{   taxon    => 1163643,
    #    name     => "Skud_FM1076",
    #    sciname  => "Saccharomyces kudriavzevii FM1076",
    #    prefix   => "AJHZ01",
    #    coverage => "8x Illumina GA",
    #},
    #{   taxon    => 1163644,
    #    name     => "Skud_FM1077",
    #    sciname  => "Saccharomyces kudriavzevii FM1077",
    #    prefix   => "AJIA01",
    #    coverage => "13x Illumina GA",
    #},
    #{   taxon    => 1163636,
    #    name     => "Skud_FM1078",
    #    sciname  => "Saccharomyces kudriavzevii FM1078",
    #    prefix   => "AJIB01",
    #    coverage => "9x Illumina GA",
    #},
    #{   taxon    => 1163637,
    #    name     => "Skud_FM1079",
    #    sciname  => "Saccharomyces kudriavzevii FM1079",
    #    prefix   => "AJIC01",
    #    coverage => "15x Illumina GA",
    #},
    #{   taxon    => 1165458,
    #    name     => "Skud_FM1094",
    #    sciname  => "Saccharomyces kudriavzevii FM1094",
    #    prefix   => "AJID01",
    #    coverage => "8x Illumina GA",
    #},
    #{   taxon    => 1163646,
    #    name     => "Skud_IFO10990",
    #    sciname  => "Saccharomyces kudriavzevii IFO10990",
    #    prefix   => "AJIE01",
    #    coverage => "9x Illumina GA",
    #},
    #{   taxon    => 1163647,
    #    name     => "Skud_IFO10991",
    #    sciname  => "Saccharomyces kudriavzevii IFO10991",
    #    prefix   => "AJIF01",
    #    coverage => "9x Illumina GA",
    #},
    #{   taxon    => 1163648,
    #    name     => "Skud_IFO1803",
    #    sciname  => "Saccharomyces kudriavzevii IFO1803",
    #    prefix   => "AJIG01",
    #    coverage => "24x Illumina GA",
    #},
);

my @subdirs_fasta = File::Find::Rule->file->name('*.fsa_nt.gz')->in($fasta_dir);

for my $item (@data) {

    # match the most similar name
    my ($fasta) = map { $_->[0] }
        sort { $b->[1] <=> $a->[1] }
        map { [ $_, compare( basename($_), $item->{prefix} ) ] }
        grep {/$item->{prefix}/} @subdirs_fasta;
    $item->{fasta} = $fasta;

    if ( index( $item->{fasta}, $item->{name} ) == -1 ) {
        printf "[%s] with [%s] matches to [%s]\n", $item->{name},
            $item->{prefix}, $item->{fasta};
        die "Match errors. Please check.\n";
    }

    # prepare working dir
    my $dir = File::Spec->catdir( $data_dir, $item->{name} );
    mkdir $dir if !-e $dir;
    $item->{dir} = $dir;
}

my $tt = Template->new;
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
faops split-about toplevel.filtered.fa 10000000 .
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
    name    => "Scer_S288c",
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
[% IF item.original_id -%]
    --species [% item.taxon %]=[% item.original_id %] \
[% END -%]
[% END -%]
    --file [% data_dir %]/[% group_name %].csv

#----------------------------#
# multi genome alignment plan
#----------------------------#
# don't copy sequences (RepeatMasker done)
# Execute the following lines by copy & paste.

# Saccharomyces_18way
cd [% data_dir %]
perl [% pl_dir %]/withncbi/taxon/strain_bz.pl \
    --file [% data_dir %]/[% group_name %].csv \
    -w     [% base_dir %] \
    --name [% group_name %] \
    --multi_name Saccharomyces_18way \
    --use_name \
    --parallel [% parallel %] \
    --norm \
    -t Scer_S288c \
    -q Sarb_H_6 \
    -q Sbay_623_6C \
    -q Sbou_17 \
    -q Sbou_ATCC_MYA_796 \
    -q Sbou_EDRL \
    -q Scar_CBS_1513 \
    -q ScerSkud_VIN7 \
    -q Skud_FM1056 \
    -q Skud_IFO_1802 \
    -q Skud_ZP591 \
    -q Smik_IFO_1815_1 \
    -q Spar_NRRL_Y_17217 \
    -q Spas_CBS_1483 \
    -q Spas_CCY48_91 \
    -q Spas_Weihenstephan_34_70_2 \
    -q Sunv_A9 \
    -q Suva_MCYC_623

# exclude Skud_FM1056
# Saccharomyces_17way
cd [% data_dir %]
perl [% pl_dir %]/withncbi/taxon/strain_bz.pl \
    --file [% data_dir %]/[% group_name %].csv \
    -w     [% base_dir %] \
    --name [% group_name %] \
    --multi_name Saccharomyces_17way \
    --use_name \
    --parallel [% parallel %] \
    --norm \
    --phylo_tree [% data_dir %]/phylo/Saccharomyces_18way.nwk \
    -t Scer_S288c \
    -q Sarb_H_6 \
    -q Sbay_623_6C \
    -q Sbou_17 \
    -q Sbou_ATCC_MYA_796 \
    -q Sbou_EDRL \
    -q Scar_CBS_1513 \
    -q ScerSkud_VIN7 \
    -q Skud_IFO_1802 \
    -q Skud_ZP591 \
    -q Smik_IFO_1815_1 \
    -q Spar_NRRL_Y_17217 \
    -q Spas_CBS_1483 \
    -q Spas_CCY48_91 \
    -q Spas_Weihenstephan_34_70_2 \
    -q Sunv_A9 \
    -q Suva_MCYC_623

# Use for outgroup
# Saccharomyces_5way_outgroup
cd [% data_dir %]
perl [% pl_dir %]/withncbi/taxon/strain_bz.pl \
    --file [% data_dir %]/[% group_name %].csv \
    -w     [% base_dir %] \
    --name [% group_name %] \
    --multi_name Saccharomyces_5way_outgroup \
    --use_name \
    --parallel [% parallel %] \
    --norm \
    --phylo_tree [% data_dir %]/phylo/Saccharomyces_18way.nwk \
    -t Scer_S288c \
    -q Sbou_ATCC_MYA_796 \
    -q Scar_CBS_1513 \
    -q Spar_NRRL_Y_17217 \
    -q Spas_CBS_1483

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

$text = <<'EOF';
#!/bin/bash

cd [% data_dir %]

# Copy & paste what you want to execute.

#----------------------------#
# clean RepeatMasker outputs
#----------------------------#
# find [% data_dir %] -type f -name "*.fasta*" | xargs rm

#----------------------------#
# only keeps chr.2bit files
#----------------------------#
# find [% data_dir %] -type f -name "*.fa" | xargs rm

#----------------------------#
# clean pairwise maf
#----------------------------#
# find [% data_dir %] -type d -name "mafSynNet" | xargs rm -fr
# find [% data_dir %] -type d -name "mafNet" | xargs rm -fr

#----------------------------#
# clean maf-fasta
#----------------------------#
# rm -fr [% data_dir %]/*_fasta

#----------------------------#
# Restore to the beginning
#----------------------------#
# find . -maxdepth 1 -type d -not -path "*WGS" | xargs rm -fr
# rm *.xlsx *.csv *.sh *.bat *.nwk

EOF
$tt->process(
    \$text,
    {   data     => \@data,
        data_dir => $data_dir,
        pl_dir   => $pl_dir,
    },
    File::Spec->catfile( $data_dir, "XXX_clean.sh" )
) or die Template->error;

__END__

# create withncbi/doc/saccharomyces.tsv manually

mkdir -p ~/data/alignment/saccharomyces
cd ~/data/alignment/saccharomyces

perl ~/Scripts/withncbi/util/wgs_prep.pl \
    -f ~/Scripts/withncbi/doc/saccharomyces.tsv \
    --fix \
    -o WGS \
    -a 

aria2c -x 6 -s 3 -c -i WGS/saccharomyces.url.txt

find WGS -name "*.gz" | xargs gzip -t 

# edit ~/Scripts/withncbi/pop/saccharomyces.pl, add contents from saccharomyces.data.txt

perl ~/Scripts/withncbi/pop/saccharomyces.pl
sh 01_file.sh
sh 02_rm.sh

# copy S288c sequences
mkdir Scer_S288c
cp ~/data/alignment/yeast_genome/S288c/* Scer_S288c
mv Scer_S288c/chrMito.fa Scer_S288c/chrMito.fa.skip

### execute 03_prepare.sh by copy & paste
# perl /home/wangq/Scripts/withncbi/taxon/strain_info.pl \
    ...

# for Saccharomyces_18way, execute the following bash file
# perl /home/wangq/Scripts/withncbi/taxon/strain_bz.pl \
    ...

sh 1_real_chr.sh
sh 3_pair_cmd.sh
sh 4_rawphylo.sh
sh 5_multi_cmd.sh
sh 6_multi_db_only.sh

# for other multi names, execute the following bash file
# perl /home/wangq/Scripts/withncbi/taxon/strain_bz.pl \
    ...

sh 5_multi_cmd.sh
sh 6_multi_db_only.sh
