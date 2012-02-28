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
my $base_dir = shift
    || File::Spec->catdir( $ENV{HOME}, "data/DGRP" );

{    # sra to fa
    my $bin_dir = {
        stk  => File::Spec->catdir( $ENV{HOME}, "share/sratoolkit" ),
        gatk => File::Spec->catdir( $ENV{HOME}, "share/GenomeAnalysisTK" ),
        pcd  => File::Spec->catdir( $ENV{HOME}, "share/picard" ),
    };
    my $data_dir = {
        sra  => File::Spec->catdir( $ENV{HOME}, "data/DGRP/SRP000694" ),
        proc => File::Spec->catdir( $ENV{HOME}, "data/DGRP/process" ),
        bash => File::Spec->catdir( $ENV{HOME}, "data/DGRP/bash" ),
    };
    my $ref_file = {
        seq => File::Spec->catfile( $ENV{HOME}, "data/DGRP/ref", "Dmel_65.fa" ),
        vcf =>
            File::Spec->catfile( $ENV{HOME}, "data/DGRP/ref", "Dmel_65.vcf" ),
    };

    my $parallel = 8;
    my $memory   = 4;

    my @ids = qw{ 101 177 138 176 181 208 321 332 375 38 380 391 40 406 443 517
        57 727 738 757 852 897};
    my @data = map { { name => "DGRP-$_" } } @ids;

    my $yml = LoadFile("DGRP.yml");

ITEM: for my $item (@data) {
        my $name = $item->{name};
        print "$name\n";
        my $dir = File::Spec->catdir( $data_dir->{proc}, $item->{name} );
        $item->{dir}   = $dir;
        $item->{lanes} = [];

        my @srxs = sort keys %{ $yml->{$name} };
        if ( !scalar @srxs ) {
            print "There are no srx for $name\n";
            $item = undef;
            next ITEM;
        }
        for my $srx (@srxs) {
            my $info = $yml->{$name}{$srx};
            print " " x 4, "$srx\n";

            my $platform = $info->{platform};
            my $layout   = $info->{layout};
            next unless $platform =~ /illumina|solexa/i;
            next unless $layout   =~ /pair/i;

            for my $i ( 0 .. scalar @{ $info->{srr} } - 1 ) {
                my $srr = $info->{srr}[$i];
                print " " x 8, "$srr\n";
                my $file = File::Spec->catfile( $data_dir->{sra}, "$srr.sra" );
                if ( !-e $file ) {
                    print "Can't find $srr.sra for $name\n";
                    $item = undef;
                    next ITEM;
                }
                my $rg_str
                    = '@RG'
                    . "\\tID:$srr"
                    . "\\tLB:$srx"
                    . "\\tPL:$platform"
                    . "\\tSM:$name";
                my $lane = {
                    srr    => $srr,
                    file   => $file,
                    rg_str => $rg_str,
                };

                push @{ $item->{lanes} }, $lane;
            }
        }
    }
    @data = grep { defined $_ } @data;

    #print Dump \@data;

    my $tt = Template->new;

    my $text = <<'EOF';
#!/bin/bash
start_time=`date +%s`

cd [% base_dir %]

# prepare reference
# this should be done at the first time
# ~/bin/x86_64/twoBitToFa ~/data/alignment/dgrp/Dmel_65/chr.2bit Dmel_65.fa.masked
# perl -p -e '/>/ and next; $_ = uc' Dmel_65.fa.masked > ref/Dmel_65.fa
# rm Dmel_65.fa.masked

# index reference genome
# bwa index -a bwtsw ref/Dmel_65.fa
# samtools faidx ref/Dmel_65.fa

# prepare vcf
# perl ~/Scripts/alignDB/util/gvf2vcf.pl ~/data/ensembl65/variation/gvf/drosophila_melanogaster/Drosophila_melanogaster.gvf > Dmel_65.unsort.vcf
# ~/share/vcftools/vcf-sort Dmel_65.unsort.vcf > ref/Dmel_65.vcf
# rm Dmel_65.unsort.vcf

if [ -d [% item.dir %] ];
then
    rm -fr [% item.dir %] ;
fi;
mkdir [% item.dir %]


[% FOREACH lane IN item.lanes -%]
#
# [% lane.srr %]
#
mkdir [% item.dir %]/[% lane.srr %]

# sra to fastq (pair end)
[% bin_dir.stk %]/fastq-dump [% lane.file %] \
    --split-files --gzip -O [% item.dir %]/[% lane.srr %]
[ $? -ne 0 ] && echo `date` [% item.name %] [% lane.srr %] [fastq dump] failed >> [% base_dir %]/fail.log && exit 255

# align pair reads to reference genome
bwa aln -q 15 -t [% parallel %] [% ref_file.seq %] [% item.dir %]/[% lane.srr %]/[% lane.srr %]_1.fastq.gz \
    > [% item.dir %]/[% lane.srr %]/[% lane.srr %]_1.sai
[ $? -ne 0 ] && echo `date` [% item.name %] [% lane.srr %] [bwa aln] failed >> [% base_dir %]/fail.log && exit 255

# align pair reads to reference genome
bwa aln -q 15 -t [% parallel %] [% ref_file.seq %] [% item.dir %]/[% lane.srr %]/[% lane.srr %]_2.fastq.gz \
    > [% item.dir %]/[% lane.srr %]/[% lane.srr %]_2.sai
[ $? -ne 0 ] && echo `date` [% item.name %] [% lane.srr %] [bwa aln] failed >> [% base_dir %]/fail.log && exit 255

# convert sai to sam
# add read groups info
bwa sampe -r "[% lane.rg_str %]" \
    [% ref_file.seq %] \
    [% item.dir %]/[% lane.srr %]/[% lane.srr %]*.sai \
    [% item.dir %]/[% lane.srr %]/[% lane.srr %]*.fastq.gz \
    | gzip > [% item.dir %]/[% lane.srr %]/[% lane.srr %].sam.gz
[ $? -ne 0 ] && echo `date` [% item.name %] [bwa sampe] failed >> [% base_dir %]/fail.log && exit 255

# convert sam to bam
samtools view -uS [% item.dir %]/[% lane.srr %]/[% lane.srr %].sam.gz \
    | samtools sort - [% item.dir %]/[% lane.srr %]/[% lane.srr %].tmp1
samtools fixmate [% item.dir %]/[% lane.srr %]/[% lane.srr %].tmp1.bam - \
    | samtools sort - [% item.dir %]/[% lane.srr %]/[% lane.srr %]

# clean
mv [% item.dir %]/[% lane.srr %]/[% lane.srr %].sam.gz [% item.dir %]/[% lane.srr %].sam.gz
mv [% item.dir %]/[% lane.srr %]/[% lane.srr %].bam    [% item.dir %]/[% lane.srr %].bam
rm -fr [% item.dir %]/[% lane.srr %]/

[% END -%]

[% IF item.lanes.size > 1 -%]
# merge with samtools
perl -e 'print "[% FOREACH lane IN item.lanes %]\[% lane.rg_str %]\n[% END %]"' > [% item.dir %]/rg.txt
samtools merge -rh [% item.dir %]/rg.txt [% item.dir %]/[% item.name %].bam [% FOREACH lane IN item.lanes %] [% item.dir %]/[% lane.srr %].bam [% END %]
rm [% FOREACH lane IN item.lanes %] [% item.dir %]/[% lane.srr %].bam [% END %]

# sort bam
samtools sort [% item.dir %]/[% item.name %].bam [% item.dir %]/[% item.name %].sort

[% ELSE -%]
# rename bam
mv [% item.dir %]/[% lane.srr %].bam [% item.dir %]/[% item.name %].sort.bam

[% END -%]
# index bam
samtools index [% item.dir %]/[% item.name %].sort.bam

# index regions for realignment
java -Xmx[% memory %]g -jar [% bin_dir.gatk %]/GenomeAnalysisTK.jar -nt [% parallel %] \
    -T RealignerTargetCreator \
    -R [% ref_file.seq %] \
    -I [% item.dir %]/[% item.name %].sort.bam  \
    --out [% item.dir %]/[% item.name %].intervals
[ $? -ne 0 ] && echo `date` [% item.name %] [gatk target] failed >> [% base_dir %]/fail.log && exit 255

# realign bam to get better Indel calling
java -Xmx[% memory %]g -jar [% bin_dir.gatk %]/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R [% ref_file.seq %] \
    -I [% item.dir %]/[% item.name %].sort.bam \
    -targetIntervals [% item.dir %]/[% item.name %].intervals \
    --out [% item.dir %]/[% item.name %].realign.bam
[ $? -ne 0 ] && echo `date` [% item.name %] [gatk realign] failed >> [% base_dir %]/fail.log && exit 255

# dup marking
java -Xmx[% memory %]g -jar [% bin_dir.pcd %]/MarkDuplicates.jar \
    INPUT=[% item.dir %]/[% item.name %].realign.bam \
    OUTPUT=[% item.dir %]/[% item.name %].dedup.bam \
    METRICS_FILE=[% item.dir %]/output.metrics \
    ASSUME_SORTED=true \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=LENIENT
[ $? -ne 0 ] && echo `date` [% item.name %] [picard dedup] failed >> [% base_dir %]/fail.log && exit 255

# reindex the realigned dedup BAM
samtools index [% item.dir %]/[% item.name %].dedup.bam

# recalibration - Count covariates
java -Xmx[% memory %]g -jar [% bin_dir.gatk %]/GenomeAnalysisTK.jar -nt [% parallel %] \
    -T CountCovariates  \
    -R [% ref_file.seq %] \
    -I [% item.dir %]/[% item.name %].dedup.bam \
    -knownSites [% ref_file.vcf %] \
    -recalFile [% item.dir %]/recal_data.csv \
    -cov ReadGroupCovariate \
    -cov QualityScoreCovariate \
    -cov CycleCovariate \
    -cov DinucCovariate
[ $? -ne 0 ] && echo `date` [% item.name %] [gatk covariates] failed >> [% base_dir %]/fail.log && exit 255

# recalibration - Tabulate recalibration
java -Xmx[% memory %]g -jar [% bin_dir.gatk %]/GenomeAnalysisTK.jar \
    -T TableRecalibration  \
    -R [% ref_file.seq %] \
    -I [% item.dir %]/[% item.name %].dedup.bam \
    -o [% item.dir %]/[% item.name %].recal.bam \
    -recalFile [% item.dir %]/recal_data.csv
[ $? -ne 0 ] && echo `date` [% item.name %] [gatk recal] failed >> [% base_dir %]/fail.log && exit 255

# generate fastq from bam
samtools mpileup -uf [% ref_file.seq %] [% item.dir %]/[% item.name %].recal.bam \
    | bcftools view -cg - \
    | vcfutils.pl vcf2fq > [% item.dir %]/[% item.name %].fq

# convert fastq to fasta 
# mask bases with quality lower than 20 to lowercases
seqtk fq2fa [% item.dir %]/[% item.name %].fq 20 > [% item.dir %]/[% item.name %].fa

# let's clean up
find [% item.dir %] -type f \
    -name "*.sai"    -o -name "*.fastq.gz" \
    -o -name "*.bam" -o -name "*.bai" \
    -o -name "*.csv" -o -name "*.intervals" \
    | grep -v "recal" | xargs rm

echo run time is $(expr `date +%s` - $start_time) s

EOF
    for my $item (@data) {
        $tt->process(
            \$text,
            {   base_dir => $base_dir,
                item     => $item,
                bin_dir  => $bin_dir,
                data_dir => $data_dir,
                ref_file => $ref_file,
                parallel => $parallel,
                memory   => $memory,
            },
            File::Spec->catfile( $data_dir->{bash}, $item->{name} . "_sra.sh" )
        ) or die Template->error;
    }

    $text = <<'EOF';
#!/bin/bash
cd [% base_dir %]
rm [% base_dir %]/fail.log

[% FOREACH item IN data -%]
# [% item.name %]
bsub -q mpi_2 -n [% parallel %] -J [% item.name %] "sh [% data_dir.bash %]/[% item.name %]_sra.sh"
# sh [% data_dir.bash %]/[% item.name %]_sra.sh

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            base_dir => $base_dir,
            data_dir => $data_dir,
            parallel => $parallel,
        },
        File::Spec->catfile( $base_dir, "master.sh" )
    ) or die Template->error;
}
