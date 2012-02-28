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
    || File::Spec->catdir( $ENV{HOME}, "data/yeast10000" );

{    # sra to fa
    my $stkbin_dir = File::Spec->catdir( $ENV{HOME}, "share/sratoolkit" );
    my $gatkbin_dir
        = File::Spec->catdir( $ENV{HOME}, "share/GenomeAnalysisTK" );
    my $pcdbin_dir
        = File::Spec->catdir( $ENV{HOME}, "share/picard" );

    my $sra_dir = File::Spec->catdir( $ENV{HOME}, "data/yeast10000/ERP000547" );
    my $ref_seq
        = File::Spec->catfile( $ENV{HOME}, "data/yeast10000/ref", "S288C.fa" );
    my $ref_vcf
        = File::Spec->catfile( $ENV{HOME}, "data/yeast10000/ref", "S288C.vcf" );
    my $proc_dir = File::Spec->catdir( $ENV{HOME}, "data/yeast10000/process" );
    my $bash_dir = File::Spec->catdir( $ENV{HOME}, "data/yeast10000/bash" );

    my $parallel = 8;
    my $memory = 4;

    my @files = File::Find::Rule->file->name('*.sra')->in($sra_dir);

    my @data;
    for my $file ( sort @files ) {

        my $item = {
            name => basename( $file, ".sra" ),
            file => $file,
        };

        # prepare working dir
        my $dir = File::Spec->catdir( $proc_dir, $item->{name} );
        $item->{dir} = $dir;

        push @data, $item;
    }

    my $tt = Template->new;

    my $text = <<'EOF';
#!/bin/bash
start_time=`date +%s`

cd [% base_dir %]

# prepare reference
# this should be done at the first time
# cat ~/data/alignment/yeast65/S288C/*.fa > ref/S288C.fa.masked
# perl -p -e '/>/ and next; $_ = uc' ref/S288C.fa.masked > ref/S288C.fa
# rm ref/S288C.fa.masked

# index reference genome
# bwa index -a bwtsw ref/S288C.fa
# samtools faidx ref/S288C.fa

# prepare vcf
# perl ~/Scripts/alignDB/gvf2vcf.pl ~/data/ensembl65/variation/gvf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.gvf > S288C.unsort.vcf
# ~/share/vcftools/vcf-sort S288C.unsort.vcf > ref/S288C.vcf
# rm S288C.unsort.vcf

if [ -d [% item.dir %] ];
then
    rm -fr [% item.dir %] ;
fi;
mkdir [% item.dir %]

# sra to fastq (pair end)
[% stkbin_dir %]/fastq-dump [% item.file %] \
    --split-files --gzip -O [% item.dir %]
[ $? -ne 0 ] && echo `date` [% item.name %] [fastq dump] failed >> [% base_dir %]/fail.log && exit 255

# align pair reads to reference genome
bwa aln -q 15 -t [% parallel %] [% ref_seq %] [% item.dir %]/[% item.name %]_1.fastq.gz \
    > [% item.dir %]/[% item.name %]_1.sai
[ $? -ne 0 ] && echo `date` [% item.name %] [bwa aln] failed >> [% base_dir %]/fail.log && exit 255

bwa aln -q 15 -t [% parallel %] [% ref_seq %] [% item.dir %]/[% item.name %]_2.fastq.gz \
    > [% item.dir %]/[% item.name %]_2.sai
[ $? -ne 0 ] && echo `date` [% item.name %] [bwa aln] failed >> [% base_dir %]/fail.log && exit 255

# convert sai to sam
# add read groups info
bwa sampe -r "@RG\tID:[% item.name %]\tLB:[% item.name %]\tPL:ILLUMINA\tSM:[% item.name %]" \
    [% ref_seq %] \
    [% item.dir %]/[% item.name %]*.sai [% item.dir %]/[% item.name %]*.fastq.gz \
    | gzip > [% item.dir %]/[% item.name %].sam.gz
[ $? -ne 0 ] && echo `date` [% item.name %] [bwa sampe] failed >> [% base_dir %]/fail.log && exit 255

# convert sam to bam and sort
samtools view -uS [% item.dir %]/[% item.name %].sam.gz \
    | samtools sort -n - [% item.dir %]/[% item.name %].tmp1
samtools fixmate [% item.dir %]/[% item.name %].tmp1.bam - \
    | samtools sort - [% item.dir %]/[% item.name %].sort
rm [% item.dir %]/[% item.name %].tmp1.bam

# index bam
samtools index [% item.dir %]/[% item.name %].sort.bam

# index regions for realignment
java -Xmx[% memory %]g -jar [% gatkbin_dir %]/GenomeAnalysisTK.jar -nt [% parallel %] \
    -T RealignerTargetCreator \
    -R [% ref_seq %] \
    -I [% item.dir %]/[% item.name %].sort.bam  \
    --out [% item.dir %]/[% item.name %].intervals
[ $? -ne 0 ] && echo `date` [% item.name %] [gatk target] failed >> [% base_dir %]/fail.log && exit 255

# realign bam to get better Indel calling
java -Xmx[% memory %]g -jar [% gatkbin_dir %]/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R [% ref_seq %] \
    -I [% item.dir %]/[% item.name %].sort.bam \
    -targetIntervals [% item.dir %]/[% item.name %].intervals \
    --out [% item.dir %]/[% item.name %].realign.bam
[ $? -ne 0 ] && echo `date` [% item.name %] [gatk realign] failed >> [% base_dir %]/fail.log && exit 255

# dup marking
java -Xmx[% memory %]g -jar [% pcdbin_dir %]/MarkDuplicates.jar \
    INPUT=[% item.dir %]/[% item.name %].realign.bam \
    OUTPUT=[% item.dir %]/[% item.name %].dedup.bam \
    METRICS_FILE=[% item.dir %]/output.metrics \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=LENIENT
[ $? -ne 0 ] && echo `date` [% item.name %] [picard dedup] failed >> [% base_dir %]/fail.log && exit 255

# reindex the realigned dedup BAM
samtools index [% item.dir %]/[% item.name %].dedup.bam

# recalibration - Count covariates
java -Xmx[% memory %]g -jar [% gatkbin_dir %]/GenomeAnalysisTK.jar -nt [% parallel %] \
    -T CountCovariates  \
    -R [% ref_seq %] \
    -I [% item.dir %]/[% item.name %].dedup.bam \
    -knownSites [% ref_vcf %] \
    -recalFile [% item.dir %]/recal_data.csv \
    -cov ReadGroupCovariate \
    -cov QualityScoreCovariate \
    -cov CycleCovariate \
    -cov DinucCovariate
[ $? -ne 0 ] && echo `date` [% item.name %] [gatk covariates] failed >> [% base_dir %]/fail.log && exit 255

# recalibration - Tabulate recalibration
java -Xmx[% memory %]g -jar [% gatkbin_dir %]/GenomeAnalysisTK.jar \
    -T TableRecalibration  \
    -R [% ref_seq %] \
    -I [% item.dir %]/[% item.name %].dedup.bam \
    -o [% item.dir %]/[% item.name %].recal.bam \
    -recalFile [% item.dir %]/recal_data.csv
[ $? -ne 0 ] && echo `date` [% item.name %] [gatk recal] failed >> [% base_dir %]/fail.log && exit 255

# generate fastq from bam
samtools mpileup -uf [% ref_seq %] [% item.dir %]/[% item.name %].recal.bam \
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
            {   base_dir    => $base_dir,
                item        => $item,
                stkbin_dir  => $stkbin_dir,
                gatkbin_dir => $gatkbin_dir,
                pcdbin_dir => $pcdbin_dir,
                sra_dir     => $sra_dir,
                ref_seq     => $ref_seq,
                ref_vcf     => $ref_vcf,
                proc_dir    => $proc_dir,
                parallel    => $parallel,
                memory    => $memory,
            },
            File::Spec->catfile( $bash_dir, $item->{name} . "_sra.sh" )
        ) or die Template->error;
    }

    my @jobs;
    while ( scalar @data ) {
        my @batching = splice @data, 0, 5;
        push @jobs, [@batching];
    }

    $text = <<'EOF';
#!/bin/bash
cd [% base_dir %]
rm [% base_dir %]/fail.log

[% FOREACH job IN jobs -%]
# [% job.0.name %]
bsub -q mpi_2 -n [% parallel %] -J [% job.0.name %] "[% FOREACH item IN job -%] sh [% bash_dir %]/[% item.name %]_sra.sh && [% END -%] sleep 1"
# sh [% bash_dir %]/[% item.name %]_sra.sh

[% END -%]

EOF
    $tt->process(
        \$text,
        {   jobs     => \@jobs,
            base_dir => $base_dir,
            bash_dir => $bash_dir,
            parallel => $parallel,
        },
        File::Spec->catfile( $base_dir, "master.sh" )
    ) or die Template->error;
}
