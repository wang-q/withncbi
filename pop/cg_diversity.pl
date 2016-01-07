#!/usr/bin/perl
use strict;
use warnings;

use Template;
use File::Basename;
use File::Find::Rule;
use File::Spec;
use String::Compare;
use YAML qw(Dump Load DumpFile LoadFile);

my $store_dir = shift
    || File::Spec->catdir( $ENV{HOME}, "data/alignment/human_cg" );
my $parallel = 12;

{    # on linux
    my $data_dir = File::Spec->catdir( $ENV{HOME}, "data/alignment/human_cg" );
    my $pl_dir   = File::Spec->catdir( $ENV{HOME}, "Scripts" );
    my $kentbin_dir = File::Spec->catdir( $ENV{HOME}, "bin/x86_64" );
    my $vcfbin_dir  = File::Spec->catdir( $ENV{HOME}, "share/vcftools" );

    # CG 69 genomes
    my $all_vcf = File::Spec->catdir( $ENV{HOME},
        "data/CG/Public_Genome_Summary_Analysis/Complete_Public_Genomes_69genomes_all_VCF.txt"
    );

    my @data = (
        {   taxon  => 902001,
            name   => "NA19700",
            pop    => "ASW",
            column => "NA19700-200-37-ASM",
        },

 #{taxon =>902002,name => "NA19701",pop =>"ASW",column =>"NA19701-200-37-ASM",},
 #{taxon =>902003,name => "NA19703",pop =>"ASW",column =>"NA19703-200-37-ASM",},
 #{taxon =>902004,name => "NA19704",pop =>"ASW",column =>"NA19704-200-37-ASM",},
 #{taxon =>902005,name => "NA19834",pop =>"ASW",column =>"NA19834-200-37-ASM",},
        {   taxon  => 902006,
            name   => "NA06985",
            pop    => "CEU",
            column => "NA06985-200-37-ASM",
        },

 #{taxon =>902007,name => "NA06994",pop =>"CEU",column =>"NA06994-200-37-ASM",},
 #{taxon =>902008,name => "NA07357",pop =>"CEU",column =>"NA07357-200-37-ASM",},
 #{taxon =>902009,name => "NA10851",pop =>"CEU",column =>"NA10851-200-37-ASM",},
 #{taxon =>902010,name => "NA12004",pop =>"CEU",column =>"NA12004-200-37-ASM",},
        {   taxon  => 902011,
            name   => "NA18526",
            pop    => "CHB",
            column => "NA18526-200-37-ASM",
        },

 #{taxon =>902012,name => "NA18537",pop =>"CHB",column =>"NA18537-200-37-ASM",},
 #{taxon =>902013,name => "NA18555",pop =>"CHB",column =>"NA18555-200-37-ASM",},
 #{taxon =>902014,name => "NA18558",pop =>"CHB",column =>"NA18558-200-37-ASM",},
        {   taxon  => 902015,
            name   => "NA20845",
            pop    => "GIH",
            column => "NA20845-200-37-ASM",
        },

 #{taxon =>902016,name => "NA20846",pop =>"GIH",column =>"NA20846-200-37-ASM",},
 #{taxon =>902017,name => "NA20847",pop =>"GIH",column =>"NA20847-200-37-ASM",},
 #{taxon =>902018,name => "NA20850",pop =>"GIH",column =>"NA20850-200-37-ASM",},
        {   taxon  => 902019,
            name   => "NA18940",
            pop    => "JPT",
            column => "NA18940-200-37-ASM",
        },

 #{taxon =>902020,name => "NA18942",pop =>"JPT",column =>"NA18942-200-37-ASM",},
 #{taxon =>902021,name => "NA18947",pop =>"JPT",column =>"NA18947-200-37-ASM",},
 #{taxon =>902022,name => "NA18956",pop =>"JPT",column =>"NA18956-200-37-ASM",},
        {   taxon  => 902023,
            name   => "NA19017",
            pop    => "LWK",
            column => "NA19017-200-37-ASM",
        },

 #{taxon =>902024,name => "NA19020",pop =>"LWK",column =>"NA19020-200-37-ASM",},
 #{taxon =>902025,name => "NA19025",pop =>"LWK",column =>"NA19025-200-37-ASM",},
 #{taxon =>902026,name => "NA19026",pop =>"LWK",column =>"NA19026-200-37-ASM",},
        {   taxon  => 902027,
            name   => "NA21732",
            pop    => "MKK",
            column => "NA21732-200-37-ASM",
        },

 #{taxon =>902028,name => "NA21733",pop =>"MKK",column =>"NA21733-200-37-ASM",},
 #{taxon =>902029,name => "NA21737",pop =>"MKK",column =>"NA21737-200-37-ASM",},
 #{taxon =>902030,name => "NA21767",pop =>"MKK",column =>"NA21767-200-37-ASM",},
        {   taxon  => 902031,
            name   => "NA19735",
            pop    => "MXL",
            column => "NA19735-200-37-ASM",
        },

 #{taxon =>902032,name => "NA19648",pop =>"MXL",column =>"NA19648-200-37-ASM",},
 #{taxon =>902033,name => "NA19649",pop =>"MXL",column =>"NA19649-200-37-ASM",},
 #{taxon =>902034,name => "NA19669",pop =>"MXL",column =>"NA19669-200-37-ASM",},
 #{taxon =>902035,name => "NA19670",pop =>"MXL",column =>"NA19670-200-37-ASM",},
        {   taxon  => 902036,
            name   => "NA20502",
            pop    => "TSI",
            column => "NA20502-200-37-ASM",
        },

 #{taxon =>902037,name => "NA20509",pop =>"TSI",column =>"NA20509-200-37-ASM",},
 #{taxon =>902038,name => "NA20510",pop =>"TSI",column =>"NA20510-200-37-ASM",},
 #{taxon =>902039,name => "NA20511",pop =>"TSI",column =>"NA20511-200-37-ASM",},
        {   taxon  => 902040,
            name   => "NA18501",
            pop    => "YRI",
            column => "NA18501-200-37-ASM",
        },

 #{taxon =>902041,name => "NA18502",pop =>"YRI",column =>"NA18502-200-37-ASM",},
 #{taxon =>902042,name => "NA18504",pop =>"YRI",column =>"NA18504-200-37-ASM",},
 #{taxon =>902043,name => "NA18505",pop =>"YRI",column =>"NA18505-200-37-ASM",},
 #{taxon =>902044,name => "NA18508",pop =>"YRI",column =>"NA18508-200-37-ASM",},
 #{taxon =>902045,name => "NA18517",pop =>"YRI",column =>"NA18517-200-37-ASM",},
 #{taxon =>902046,name => "NA19129",pop =>"YRI",column =>"NA19129-200-37-ASM",},

    );

    for my $item ( sort @data ) {
        my $name = $item->{name};

        ## match the most similar name
       #my ($file) = map { $_->[0] }
       #    sort { $b->[1] <=> $a->[1] }
       #    map { [ $_, compare( lc basename($_), lc $item->{name} ) ] } @files;
       #$item->{seq} = $file;

        # prepare working dir
        my $dir = File::Spec->catdir( $data_dir, $name );
        mkdir $dir if !-e $dir;
        $item->{dir} = $dir;
    }

    #for my $d ( qw{ maf } ) {
    #    my $dir = File::Spec->catdir( $data_dir, $d );
    #    mkdir $dir if !-e $dir;
    #}

    print Dump \@data;

    my $tt = Template->new;

    # taxon.csv
    my $text = <<'EOF';
[% FOREACH item IN data -%]
[% item.taxon %],Homo,sapiens,[% item.name %],,
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
[% item.taxon %],chrUn,999999999,[% item.name %]/cg69
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

perl [% pl_dir %]/alignDB/util/merge_csv.pl -t [% pl_dir %]/alignDB/init/taxon.csv -m [% data_dir %]/taxon.csv

perl [% pl_dir %]/alignDB/util/merge_csv.pl -t [% pl_dir %]/alignDB/init/chr_length.csv -m [% data_dir %]/chr_length.csv

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
        },
        File::Spec->catfile( $store_dir, "merge_chr.sh" )
    ) or die Template->error;

    # id2name.csv
    $text = <<'EOF';
[% FOREACH item IN data -%]
[% item.taxon %],[% item.name %]
[% END -%]
EOF
    $tt->process(
        \$text,
        {   data => [
                @data,
                {   name  => "human",
                    taxon => 9606,
                },
                {   name  => "chimp",
                    taxon => 9598,
                },
            ],
        },
        File::Spec->catfile( $store_dir, "id2name.csv" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# vcf
#----------------------------#
[% FOREACH item IN data -%]
echo [% item.name %] [% item.pop %] 
cd [% item.dir %]
[% vcfbin_dir %]/vcftools --vcf [% all_vcf %] \
    --out [% item.name %] \
    --recode \
    --indv [% item.column %]

perl -nl -e "/#/ and print; next if /\.\/\./ or /\.\|\./; next if /0\/0/ or /0\|0/; print;" \
    [% item.name %].recode.vcf > [% item.name %].vcf

[% END -%]

EOF

    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir,
            vcfbin_dir  => $vcfbin_dir,
            all_vcf     => $all_vcf,
        },
        File::Spec->catfile( $store_dir, "vcf.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# axt to maf
#----------------------------#
if [ ! -d [% data_dir %]/maf ]
then
    mkdir [% data_dir %]/maf
fi

echo uncompress axt
find [% data_dir %]/HumanvsChimp/axtNet -name "*.axt.gz" | parallel gzip -d

echo axt to maf
for f in [% data_dir %]/HumanvsChimp/axtNet/*.axt ; do
    [% kentbin_dir %]/axtToMaf -tPrefix=human. -qPrefix=chimp. \
    "$f" \
    [% data_dir %]/human/chr.sizes \
    [% data_dir %]/chimp/chr.sizes \
    [% data_dir %]/maf/`basename $f .axt`.maf
done

echo compress axt
find [% data_dir %]/HumanvsChimp/axtNet -name "*.axt" | parallel gzip

if [ -d [% data_dir %]/fas ]
then
    rm -fr [% data_dir %]/fas
fi

#----------------------------#
# maf to fas
#----------------------------#
echo maf to fas
perl [% pl_dir %]/blastz/maf2fasta.pl \
    --parallel [% parallel %] --block --length 5000 \
    -i [% data_dir %]/maf \
    -o [% data_dir %]/fas

EOF

    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir,
            parallel    => $parallel,
        },
        File::Spec->catfile( $store_dir, "axt_maf_fas.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

if [ -d [% data_dir %]/fas_1 ]
then
    rm -fr [% data_dir %]/fas_1
fi

if [ -d [% data_dir %]/fas_2 ]
then
    rm -fr [% data_dir %]/fas_2
fi

cp -R [% data_dir %]/fas [% data_dir %]/fas_1

#----------------------------#
# add vcf to fas
#----------------------------#
[% FOREACH item IN data -%]
echo [% item.name %] [% item.pop %]
cd [% item.dir %]
perl [% pl_dir %]/alignDB/util/add_vcf2fa.pl --parallel [% parallel %] \
    -f [% item.dir %]/[% item.name %].vcf \
    -i [% data_dir %]/fas_1 \
    -o [% data_dir %]/fas_2

rm -fr [% data_dir %]/fas_1
mv [% data_dir %]/fas_2 [% data_dir %]/fas_1

[% END -%]

cd [% data_dir %]
mv [% data_dir %]/fas_1 [% data_dir %]/fas_final

EOF

    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir,
            parallel    => $parallel,
        },
        File::Spec->catfile( $store_dir, "add_vcf.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
    
#----------------------------#
# tar-gzip
#----------------------------#

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

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# mafft
#----------------------------#
# [% item.out_dir %]
perl [% pl_dir %]/blastz/refine_fasta.pl \
    --msa mafft --parallel [% parallel %] \
    --block --outgroup \
    -i [% data_dir %]/fas_final \
    -o [% data_dir %]/fas_mft

EOF
    $tt->process(
        \$text,
        {    #data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
            parallel => $parallel,
        },
        File::Spec->catfile( $store_dir, "refine.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]
    
#----------------------------#
# multi_way_batch
#----------------------------#
# mafft
perl [% pl_dir %]/alignDB/extra/multi_way_batch.pl \
    -d HumanvsXI -e human_65 \
    --block --outgroup \
    --id [% data_dir %]/id2name.csv \
    -da [% data_dir %]/fas_mft  \
    -lt 5000 --parallel [% parallel %] --run common

EOF
    $tt->process(
        \$text,
        {    #data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
            parallel => $parallel,
        },
        File::Spec->catfile( $store_dir, "multi.sh" )
    ) or die Template->error;
}
