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

my $group_name = 'scer_new';
my $base_dir   = File::Spec->catdir( $ENV{HOME}, "data/alignment" );
my $data_dir   = File::Spec->catdir( $base_dir, $group_name );
my $pl_dir     = File::Spec->catdir( $ENV{HOME}, "Scripts" );

# NCBI WGS
my $fasta_dir = File::Spec->catdir( $ENV{HOME}, "data/alignment/scer_new/WGS" );

my @data = (
    {   taxon       => 100000001,
        name        => "10560_6B",
        sciname     => "Saccharomyces cerevisiae 10560-6B",
        prefix      => "JRIQ01",
        coverage    => "191.0x Illumina HiSeq",
        original_id => 4932,
    },
    {   taxon       => 100000002,
        name        => "9464",
        sciname     => "Saccharomyces cerevisiae 9464",
        prefix      => "JSAC01",
        coverage    => "117.0x PacBio",
        original_id => 4932,
    },
    {   taxon       => 545124,
        name        => "AWRI1631",
        sciname     => "Saccharomyces cerevisiae AWRI1631",
        prefix      => "ABSV01",
        coverage    => " ",
    },
    {   taxon       => 764097,
        name        => "AWRI796",
        sciname     => "Saccharomyces cerevisiae AWRI796",
        prefix      => "ADVS01",
        coverage    => "20x 454",
    },
    {   taxon       => 100000003,
        name        => "BC187",
        sciname     => "Saccharomyces cerevisiae BC187",
        prefix      => "JRII01",
        coverage    => "177.0x Illumina HiSeq",
        original_id => 4932,
    },
    {   taxon       => 1247190,
        name        => "BY4741",
        sciname     => "Saccharomyces cerevisiae BY4741",
        prefix      => "JRIS01",
        coverage    => "209.0x Illumina HiSeq",
    },
    {   taxon       => 100000004,
        name        => "BY4742",
        sciname     => "Saccharomyces cerevisiae BY4742",
        prefix      => "JRIR01",
        coverage    => "103.0x Illumina HiSeq",
        original_id => 4932,
    },
    {   taxon       => 929587,
        name        => "CBS_7960",
        sciname     => "Saccharomyces cerevisiae CBS 7960",
        prefix      => "AEWL01",
        coverage    => "17.0X 454; ABI 3730",
    },
    {   taxon       => 889517,
        name        => "CEN_PK113_7D",
        sciname     => "Saccharomyces cerevisiae CEN.PK113-7D",
        prefix      => "AEHG01",
        coverage    => "18x 454; llumina",
    },
    {   taxon       => 100000005,
        name        => "CEN_PK2_1Ca",
        sciname     => "Saccharomyces cerevisiae CEN.PK2-1Ca",
        prefix      => "JRIV01",
        coverage    => "89.0x Illumina HiSeq",
        original_id => 4932,
    },
    {   taxon       => 464025,
        name        => "CLIB215",
        sciname     => "Saccharomyces cerevisiae CLIB215",
        prefix      => "AEWP01",
        coverage    => "16.9X 454; ABI 3730",
    },
    {   taxon       => 100000006,
        name        => "D273_10B",
        sciname     => "Saccharomyces cerevisiae D273-10B",
        prefix      => "JRIY01",
        coverage    => "112.0x Illumina HiSeq",
        original_id => 4932,
    },
    {   taxon       => 100000007,
        name        => "DBVPG6044",
        sciname     => "Saccharomyces cerevisiae DBVPG6044",
        prefix      => "JRIG01",
        coverage    => "176.0x Illumina HiSeq",
        original_id => 4932,
    },
    {   taxon       => 100000008,
        name        => "EBY_VW4000",
        sciname     => "Saccharomyces cerevisiae EBY.VW4000",
        prefix      => "JSFO01",
        coverage    => "90.0x Illumina HiSeq",
        original_id => 4932,
    },
    {   taxon       => 1095001,
        name        => "EC9_8",
        sciname     => "Saccharomyces cerevisiae EC9-8",
        prefix      => "AGSJ01",
        coverage    => "30x 454 Roche GS Junior",
    },
    {   taxon       => 947036,
        name        => "FL100",
        sciname     => "Saccharomyces cerevisiae FL100",
        prefix      => "JRIT01",
        coverage    => "184.0x Illumina HiSeq",
    },
    {   taxon       => 100000009,
        name        => "FY1679",
        sciname     => "Saccharomyces cerevisiae FY1679",
        prefix      => "JRIN01",
        coverage    => "329.0x Illumina HiSeq",
        original_id => 4932,
    },
    {   taxon       => 764102,
        name        => "FostersB",
        sciname     => "Saccharomyces cerevisiae FostersB",
        prefix      => "AEHH01",
        coverage    => "20x 454",
    },
    {   taxon       => 764101,
        name        => "FostersO",
        sciname     => "Saccharomyces cerevisiae FostersO",
        prefix      => "AEEZ01",
        coverage    => "20x 454",
    },
    {   taxon       => 1352824,
        name        => "IR_2",
        sciname     => "Saccharomyces cerevisiae IR-2",
        prefix      => "BAUI01",
        coverage    => "26.1x 454 GS FLX; SOLiD 3; ABI 3730xl",
    },
    {   taxon       => 574961,
        name        => "JAY291",
        sciname     => "Saccharomyces cerevisiae JAY291",
        prefix      => "ACFL01",
        coverage    => "12x 454; 58x Solexa single-end reads; 95x Solexa paired-end 454; Solexa",
    },
    {   taxon       => 100000010,
        name        => "JK9_3d",
        sciname     => "Saccharomyces cerevisiae JK9-3d",
        prefix      => "JRIZ01",
        coverage    => "154.0x Illumina HiSeq",
        original_id => 4932,
    },
    {   taxon       => 100000011,
        name        => "K11",
        sciname     => "Saccharomyces cerevisiae K11",
        prefix      => "JRIJ01",
        coverage    => "189.0x Illumina HiSeq",
        original_id => 4932,
    },
    {   taxon       => 721032,
        name        => "Kyokai_no_7",
        sciname     => "Saccharomyces cerevisiae Kyokai no. 7",
        prefix      => "BABQ01",
        coverage    => "9.1x ABI 3730xl",
    },
    {   taxon       => 100000012,
        name        => "L1528",
        sciname     => "Saccharomyces cerevisiae L1528",
        prefix      => "JRIK01",
        coverage    => "186.0x Illumina HiSeq",
        original_id => 4932,
    },
    {   taxon       => 764098,
        name        => "Lalvin_QA23",
        sciname     => "Saccharomyces cerevisiae Lalvin QA23",
        prefix      => "ADVV01",
        coverage    => "20x 454",
    },
    {   taxon       => 1149757,
        name        => "M3707",
        sciname     => "Saccharomyces cerevisiae M3707",
        prefix      => "AMQB01",
        coverage    => "125.6x Illumina",
    },
    {   taxon       => 1162671,
        name        => "M3836",
        sciname     => "Saccharomyces cerevisiae M3836",
        prefix      => "AMQC01",
        coverage    => "172.6x Illumina",
    },
    {   taxon       => 1162672,
        name        => "M3837",
        sciname     => "Saccharomyces cerevisiae M3837",
        prefix      => "AMQD01",
        coverage    => "185.6x Illumina",
    },
    {   taxon       => 1162673,
        name        => "M3838",
        sciname     => "Saccharomyces cerevisiae M3838",
        prefix      => "AMQE01",
        coverage    => "95.6x Illumina",
    },
    {   taxon       => 1162674,
        name        => "M3839",
        sciname     => "Saccharomyces cerevisiae M3839",
        prefix      => "AMQF01",
        coverage    => "165.5x Illumina",
    },
    {   taxon       => 1201112,
        name        => "N85",
        sciname     => "Saccharomyces cerevisiae N85",
        prefix      => "CBYJ01",
        coverage    => " ",
    },
    {   taxon       => 1352823,
        name        => "NAM34_4C",
        sciname     => "Saccharomyces cerevisiae NAM34-4C",
        prefix      => "BAUH01",
        coverage    => "31.2x 454 GS FLX; SOLiD 3; ABI 3730xl",
    },
    {   taxon       => 1331972,
        name        => "NY1308",
        sciname     => "Saccharomyces cerevisiae NY1308",
        prefix      => "ASJZ01",
        coverage    => "70.0x 454; Illumina",
    },
    {   taxon       => 1177187,
        name        => "P283",
        sciname     => "Saccharomyces cerevisiae P283",
        prefix      => "AOJC01",
        coverage    => "14.0x 454",
    },
    {   taxon       => 1182968,
        name        => "P301",
        sciname     => "Saccharomyces cerevisiae P301",
        prefix      => "APIQ01",
        coverage    => "14.0x 454",
    },
    {   taxon       => 947039,
        name        => "PW5",
        sciname     => "Saccharomyces cerevisiae PW5",
        prefix      => "AFDC01",
        coverage    => "16.10x 454",
    },
    {   taxon       => 1182966,
        name        => "R008",
        sciname     => "Saccharomyces cerevisiae R008",
        prefix      => "APIP01",
        coverage    => "14.0x 454",
    },
    {   taxon       => 1182967,
        name        => "R103",
        sciname     => "Saccharomyces cerevisiae R103",
        prefix      => "APIR01",
        coverage    => "14.0x 454",
    },
    {   taxon       => 285006,
        name        => "RM11_1a",
        sciname     => "Saccharomyces cerevisiae RM11-1a",
        prefix      => "AAEG01",
        coverage    => " ",
    },
    {   taxon       => 100000013,
        name        => "RedStar",
        sciname     => "Saccharomyces cerevisiae RedStar",
        prefix      => "JRIL01",
        coverage    => "180.0x Illumina HiSeq",
        original_id => 4932,
    },
    {   taxon       => 100000014,
        name        => "SEY6210",
        sciname     => "Saccharomyces cerevisiae SEY6210",
        prefix      => "JRIW01",
        coverage    => "106.0x Illumina HiSeq",
        original_id => 4932,
    },
    {   taxon       => 580239,
        name        => "SK1",
        sciname     => "Saccharomyces cerevisiae SK1",
        prefix      => "JRIH01",
        coverage    => "261.0x Illumina HiSeq",
    },
    {   taxon       => 100000015,
        name        => "Sbou",
        sciname     => "Saccharomyces sp. 'boulardii' ATCC MYA-796",
        prefix      => "JRHY01",
        coverage    => "403.0x Illumina HiSeq",
        original_id => 252598,
    },
    {   taxon       => 1073566,
        name        => "Scar",
        sciname     => "Saccharomyces carlsbergensis CBS 1513",
        prefix      => "AZCJ01",
        coverage    => "18.0x 454",
    },
    {   taxon       => 658763,
        name        => "Sigma1278b",
        sciname     => "Saccharomyces cerevisiae Sigma1278b",
        prefix      => "ACVY01",
        coverage    => " ",
    },
    {   taxon       => 226125,
        name        => "Spar",
        sciname     => "Saccharomyces paradoxus NRRL Y-17217",
        prefix      => "AABY01",
        coverage    => " ",
    },
    {   taxon       => 100000016,
        name        => "Spas",
        sciname     => "Saccharomyces pastorianus CBS 1483",
        prefix      => "JTFI01",
        coverage    => "52.2x Illumina HiSeq",
        original_id => 27292,
    },
    {   taxon       => 929585,
        name        => "T7",
        sciname     => "Saccharomyces cerevisiae T7",
        prefix      => "AFDE01",
        coverage    => "25.4x 454; ABI 3730",
    },
    {   taxon       => 947040,
        name        => "UC5",
        sciname     => "Saccharomyces cerevisiae UC5",
        prefix      => "AFDD01",
        coverage    => "15.7x 454",
    },
    {   taxon       => 1434269,
        name        => "UFMG_A_905",
        sciname     => "Saccharomyces cerevisiae UFMG A-905",
        prefix      => "JACO02",
        coverage    => "202x Illumina MiSeq",
    },
    {   taxon       => 100000018,
        name        => "UWOPS05_217_3",
        sciname     => "Saccharomyces cerevisiae UWOPS05_217_3",
        prefix      => "JRIM01",
        coverage    => "57.0x Illumina HiSeq",
        original_id => 4932,
    },
    {   taxon       => 764100,
        name        => "VL3",
        sciname     => "Saccharomyces cerevisiae VL3",
        prefix      => "AEJS01",
        coverage    => "20x 454",
    },
    {   taxon       => 764099,
        name        => "Vin13",
        sciname     => "Saccharomyces cerevisiae Vin13",
        prefix      => "ADXC01",
        coverage    => "20x 454",
    },
    {   taxon       => 580240,
        name        => "W303",
        sciname     => "Saccharomyces cerevisiae W303",
        prefix      => "JRIU01",
        coverage    => "301.0x Illumina HiSeq",
    },
    {   taxon       => 100000019,
        name        => "X2180_1A",
        sciname     => "Saccharomyces cerevisiae X2180-1A",
        prefix      => "JRIX01",
        coverage    => "112.0x Illumina HiSeq",
        original_id => 4932,
    },
    {   taxon       => 100000020,
        name        => "Y55",
        sciname     => "Saccharomyces cerevisiae Y55",
        prefix      => "JRIF01",
        coverage    => "112.0x Illumina HiSeq",
        original_id => 4932,
    },
    {   taxon       => 929586,
        name        => "YJM269",
        sciname     => "Saccharomyces cerevisiae YJM269",
        prefix      => "AEWN01",
        coverage    => "16.7X 454; ABI 3730",
    },
    {   taxon       => 1337529,
        name        => "YJM339",
        sciname     => "Saccharomyces cerevisiae YJM339",
        prefix      => "JRIE01",
        coverage    => "102.0x Illumina HiSeq",
    },
    {   taxon       => 307796,
        name        => "YJM789",
        sciname     => "Saccharomyces cerevisiae YJM789",
        prefix      => "AAFW02",
        coverage    => " ",
    },
    {   taxon       => 1087981,
        name        => "YJSH1",
        sciname     => "Saccharomyces cerevisiae YJSH1",
        prefix      => "AGAW01",
        coverage    => "24x 454 GS FLX Titanium; Sanger",
    },
    {   taxon       => 100000021,
        name        => "YPH499",
        sciname     => "Saccharomyces cerevisiae YPH499",
        prefix      => "JRIO01",
        coverage    => "69.0x Illumina HiSeq",
        original_id => 4932,
    },
    {   taxon       => 100000022,
        name        => "YPS128",
        sciname     => "Saccharomyces cerevisiae YPS128",
        prefix      => "JRID01",
        coverage    => "95.0x Illumina HiSeq",
        original_id => 4932,
    },
    {   taxon       => 538976,
        name        => "YPS163",
        sciname     => "Saccharomyces cerevisiae YPS163",
        prefix      => "JRIC01",
        coverage    => "96.0x Illumina HiSeq",
    },
    {   taxon       => 100000023,
        name        => "YS9",
        sciname     => "Saccharomyces cerevisiae YS9",
        prefix      => "JRIB01",
        coverage    => "100.0x Illumina HiSeq",
        original_id => 4932,
    },
    {   taxon       => 1227742,
        name        => "ZTW1",
        sciname     => "Saccharomyces cerevisiae ZTW1",
        prefix      => "AMDD01",
        coverage    => "20.0x 454; Sanger",
    },
    
    # short contigs
    #{   taxon       => 929629,
    #    name        => "CLIB324",
    #    sciname     => "Saccharomyces cerevisiae CLIB324",
    #    prefix      => "AEWM01",
    #    coverage    => "7.14X 454; ABI 3730",
    #},
    #{   taxon       => 947035,
    #    name        => "CLIB382",
    #    sciname     => "Saccharomyces cerevisiae CLIB382",
    #    prefix      => "AFDG01",
    #    coverage    => "5.96x 454",
    #},
    #{   taxon       => 538975,
    #    name        => "M22",
    #    sciname     => "Saccharomyces cerevisiae M22",
    #    prefix      => "ABPC01",
    #    coverage    => " ",
    #},
    #{   taxon       => 1337645,
    #    name        => "M5",
    #    sciname     => "Saccharomyces cerevisiae M5",
    #    prefix      => "JPXA01",
    #    coverage    => "31.0x Illumina HiSeq",
    #},
    #{   taxon       => 471859,
    #    name        => "T73",
    #    sciname     => "Saccharomyces cerevisiae T73",
    #    prefix      => "AFDF01",
    #    coverage    => "13.9x 454",
    #},
    #{   taxon       => 100000017,
    #    name        => "UCD51",
    #    sciname     => "Saccharomyces cerevisiae UCD51",
    #    prefix      => "JPXB01",
    #    coverage    => "34.0x Illumina HiSeq",
    #    original_id => 4932,
    #},
    #{   taxon       => 462210,
    #    name        => "Y10",
    #    sciname     => "Saccharomyces cerevisiae Y10",
    #    prefix      => "AEWK01",
    #    coverage    => "6.6X 454; ABI 3730",
    #},
);

# Add fake WGS records here
# http://www.pnas.org/content/suppl/2009/09/09/0904673106.DCSupplemental/ST1_PDF.pdf
push @data,
    {
    taxon    => 643680,
    name     => "EC1118",
    sciname  => "Saccharomyces cerevisiae EC1118",
    prefix   => "EC1118",
    coverage => "17.6x 454; 6x Sanger",
    };

# http://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=189879
push @data,
    {
    taxon    => 1294331,
    name     => "YJM993",
    sciname  => "Saccharomyces cerevisiae YJM993",
    prefix   => "YJM993",
    coverage => "~175x Illumina HiSeq",
    };

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
    name    => "S288c",
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

# Use this as guild tree
# Scer_69way
cd [% data_dir %]
perl [% pl_dir %]/withncbi/taxon/strain_bz.pl \
    --file [% data_dir %]/[% group_name %].csv \
    -w     [% base_dir %] \
    --name [% group_name %] \
    --multi_name Scer_69way \
    --use_name \
    --parallel [% parallel %] \
    --norm \
    -t S288c \
    -q=10560_6B    -q=9464          -q=AWRI1631  -q=AWRI796      -q=BC187       \
    -q=BY4741      -q=BY4742        -q=CBS_7960  -q=CEN_PK113_7D -q=CEN_PK2_1Ca \
    -q=CLIB215     -q=D273_10B      -q=DBVPG6044 -q=EBY_VW4000   -q=EC9_8       \
    -q=FL100       -q=FY1679        -q=FostersB  -q=FostersO     -q=IR_2        \
    -q=JAY291      -q=JK9_3d        -q=K11       -q=Kyokai_no_7  -q=L1528       \
    -q=Lalvin_QA23 -q=M3707         -q=M3836     -q=M3837        -q=M3838       \
    -q=M3839       -q=N85           -q=NAM34_4C  -q=NY1308       -q=P283        \
    -q=P301        -q=PW5           -q=R008      -q=R103         -q=RM11_1a     \
    -q=RedStar     -q=SEY6210       -q=SK1       -q=Sbou         -q=Scar        \
    -q=Sigma1278b  -q=Spar          -q=Spas      -q=T7           -q=UC5         \
    -q=UFMG_A_905  -q=UWOPS05_217_3 -q=VL3       -q=Vin13        -q=W303        \
    -q=X2180_1A    -q=Y55           -q=YJM269    -q=YJM339       -q=YJM789      \
    -q=YJSH1       -q=YPH499        -q=YPS128    -q=YPS163       -q=YS9         \
    -q=ZTW1        -q=EC1118        -q=YJM993

# Use this to find the bottleneck
# Get rid of bottlenecks
# Scer_test_WoO
cd [% data_dir %]
perl [% pl_dir %]/withncbi/taxon/strain_bz.pl \
    --file [% data_dir %]/[% group_name %].csv \
    -w     [% base_dir %] \
    --name [% group_name %] \
    --multi_name Scer_test_WoO \
    --use_name \
    --parallel [% parallel %] \
    --norm \
    --phylo_tree [% data_dir %]/phylo/Scer_69way.nwk \
    -t S288c \
    -q=10560_6B    -q=9464                       -q=AWRI796      -q=BC187       \
    -q=BY4741      -q=BY4742                     -q=CEN_PK113_7D -q=CEN_PK2_1Ca \
                   -q=D273_10B      -q=DBVPG6044 -q=EBY_VW4000   -q=EC9_8       \
    -q=FL100       -q=FY1679        -q=FostersB  -q=FostersO     -q=IR_2        \
    -q=JAY291      -q=JK9_3d        -q=K11       -q=Kyokai_no_7  -q=L1528       \
    -q=Lalvin_QA23 -q=M3707         -q=M3836     -q=M3837        -q=M3838       \
    -q=M3839       -q=N85           -q=NAM34_4C  -q=NY1308       -q=P283        \
    -q=P301                         -q=R008      -q=R103         -q=RM11_1a     \
    -q=RedStar     -q=SEY6210       -q=SK1                                      \
    -q=Sigma1278b  -q=Spar                       -q=T7           -q=UC5         \
    -q=UFMG_A_905                   -q=VL3       -q=Vin13        -q=W303        \
    -q=X2180_1A    -q=Y55           -q=YJM269    -q=YJM339       -q=YJM789      \
    -q=YJSH1       -q=YPH499        -q=YPS128    -q=YPS163       -q=YS9         \
    -q=ZTW1        -q=EC1118        -q=YJM993                                   \
    -o Spar

# Scer_61way_Spar
cd [% data_dir %]
perl [% pl_dir %]/withncbi/taxon/strain_bz.pl \
    --file [% data_dir %]/[% group_name %].csv \
    -w     [% base_dir %] \
    --name [% group_name %] \
    --multi_name Scer_61way_Spar \
    --use_name \
    --parallel [% parallel %] \
    --norm \
    --phylo_tree [% data_dir %]/phylo/Scer_69way.nwk \
    -t S288c \
    -q=10560_6B    -q=9464                       -q=AWRI796      -q=BC187       \
    -q=BY4741      -q=BY4742                     -q=CEN_PK113_7D -q=CEN_PK2_1Ca \
                   -q=D273_10B      -q=DBVPG6044 -q=EBY_VW4000   -q=EC9_8       \
    -q=FL100       -q=FY1679        -q=FostersB  -q=FostersO     -q=IR_2        \
    -q=JAY291      -q=JK9_3d        -q=K11       -q=Kyokai_no_7  -q=L1528       \
    -q=Lalvin_QA23 -q=M3707         -q=M3836     -q=M3837        -q=M3838       \
    -q=M3839       -q=N85           -q=NAM34_4C  -q=NY1308       -q=P283        \
    -q=P301                         -q=R008      -q=R103         -q=RM11_1a     \
    -q=RedStar     -q=SEY6210       -q=SK1                                      \
    -q=Sigma1278b                                -q=T7           -q=UC5         \
    -q=UFMG_A_905                   -q=VL3       -q=Vin13        -q=W303        \
    -q=X2180_1A    -q=Y55           -q=YJM269    -q=YJM339       -q=YJM789      \
    -q=YJSH1       -q=YPH499        -q=YPS128    -q=YPS163       -q=YS9         \
    -q=ZTW1        -q=EC1118        -q=YJM993

# Scer_66way_refSbou
cd [% data_dir %]
perl [% pl_dir %]/withncbi/taxon/strain_bz.pl \
    --file [% data_dir %]/[% group_name %].csv \
    -w     [% base_dir %] \
    --name [% group_name %] \
    --multi_name Scer_70way_GT10M \
    --use_name \
    --parallel [% parallel %] \
    --norm \
    --phylo_tree [% data_dir %]/phylo/Scer_69way.nwk \
    -t S288c \
    -q=10560_6B    -q=9464          -q=AWRI1631  -q=AWRI796      -q=BC187       \
    -q=BY4741      -q=BY4742        -q=CBS_7960  -q=CEN_PK113_7D -q=CEN_PK2_1Ca \
    -q=CLIB215     -q=D273_10B      -q=DBVPG6044 -q=EBY_VW4000   -q=EC9_8       \
    -q=FL100       -q=FY1679        -q=FostersB  -q=FostersO     -q=IR_2        \
    -q=JAY291      -q=JK9_3d        -q=K11       -q=Kyokai_no_7  -q=L1528       \
    -q=Lalvin_QA23 -q=M3707         -q=M3836     -q=M3837        -q=M3838       \
    -q=M3839       -q=N85           -q=NAM34_4C  -q=NY1308       -q=P283        \
    -q=P301        -q=PW5           -q=R008      -q=R103         -q=RM11_1a     \
    -q=RedStar     -q=SEY6210       -q=SK1       -q=Sbou                        \
    -q=Sigma1278b                                -q=T7           -q=UC5         \
    -q=UFMG_A_905  -q=UWOPS05_217_3 -q=VL3       -q=Vin13        -q=W303        \
    -q=X2180_1A    -q=Y55           -q=YJM269    -q=YJM339       -q=YJM789      \
    -q=YJSH1       -q=YPH499        -q=YPS128    -q=YPS163       -q=YS9         \
    -q=ZTW1        -q=EC1118        -q=YJM993                                   \
    -o Sbou

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

# create withncbi/doc/scer_new.tsv manually

mkdir -p ~/data/alignment/scer_new
cd ~/data/alignment/scer_new

perl ~/Scripts/withncbi/util/wgs_prep.pl \
    -f ~/Scripts/withncbi/doc/scer_new.tsv \
    --fix \
    -o WGS \
    -a 

aria2c -x 6 -s 3 -c -i WGS/scer_new.url.txt

find WGS -name "*.gz" | xargs gzip -t 

# edit ~/Scripts/withncbi/pop/scer_new.pl, add contents to @data from saccharomyces.data.txt

# EC1118 and YJM993 are not in WGS
# EC1118 http://www.ncbi.nlm.nih.gov/assembly/GCA_000218975.1/
# YJM993 http://www.ncbi.nlm.nih.gov/assembly/GCA_000662435.1/
# Then create "yeast_name_seq.csv" manually and put it in "~/data/alignment/yeast_genome"
cd ~/data/alignment/yeast_genome
perl ~/Scripts/withncbi/util/batch_get_seq.pl yeast_name_seq.csv  2>&1 | tee yeast_name_seq.log

# pretend to be WGS
mkdir ~/data/alignment/scer_new/WGS/EC1118
cat ~/data/alignment/yeast_genome/EC1118/*.fasta > ~/data/alignment/scer_new/WGS/EC1118/EC1118.fsa_nt
gzip ~/data/alignment/scer_new/WGS/EC1118/EC1118.fsa_nt

mkdir ~/data/alignment/scer_new/WGS/YJM993
cat ~/data/alignment/yeast_genome/YJM993/*.fasta > ~/data/alignment/scer_new/WGS/YJM993/YJM993.fsa_nt
gzip ~/data/alignment/scer_new/WGS/YJM993/YJM993.fsa_nt

# add EC1118 and YJM993 to @data

# Run the perl script again if you changed it.
perl ~/Scripts/withncbi/pop/scer_new.pl
sh 01_file.sh
sh 02_rm.sh

# copy S288c sequences
mkdir S288c
cp ~/data/alignment/yeast_genome/S288c/* S288c
mv S288c/chrMito.fa S288c/chrMito.fa.skip

### execute 03_prepare.sh by copy & paste
# perl /home/wangq/Scripts/withncbi/taxon/strain_info.pl \
    ...

# for Scer_69way, execute the following bash file
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
