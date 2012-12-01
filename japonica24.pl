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
    || File::Spec->catdir( $ENV{HOME}, "data/alignment/rice" );
my $parallel = 12;
{    # on linux
    my $data_dir    = File::Spec->catdir( $ENV{HOME}, "data/alignment/rice" );
    my $pl_dir      = File::Spec->catdir( $ENV{HOME}, "Scripts" );
    my $kentbin_dir = File::Spec->catdir( $ENV{HOME}, "bin/x86_64" );

    # rice50
    my $seq_dir = File::Spec->catdir( $ENV{HOME}, "data/japonica24/process/" );

    my @data = (

        #TEJ
        { taxon => 900901, name => "IRGC1107",  coverage => 11.33, },
        { taxon => 900902, name => "IRGC2540",  coverage => 17.49, },
        { taxon => 900903, name => "IRGC27630", coverage => 9.6, },
        { taxon => 900904, name => "IRGC32399", coverage => 13.21, },
        { taxon => 900905, name => "IRGC418",   coverage => 13.27, },
        { taxon => 900906, name => "IRGC55471", coverage => 15.99, },
        { taxon => 900907, name => "IRGC8191",  coverage => 13.75, },
        { taxon => 900908, name => "IRGC38698", coverage => 16.07, },

        #TRJ
        { taxon => 900909, name => "IRGC11010", coverage => 11.51, },
        { taxon => 900910, name => "IRGC17757", coverage => 10.74, },
        { taxon => 900911, name => "IRGC328",   coverage => 12.5, },
        { taxon => 900912, name => "IRGC43325", coverage => 11.79, },
        { taxon => 900913, name => "IRGC43675", coverage => 10.97, },
        { taxon => 900914, name => "IRGC50448", coverage => 15.35, },
        { taxon => 900915, name => "IRGC66756", coverage => 11.83, },
        { taxon => 900916, name => "IRGC8244",  coverage => 11.27, },
        { taxon => 900917, name => "IRGC26872", coverage => 12.78, },

        #ARO
        { taxon => 900918, name => "IRGC12793", coverage => 14.64, },
        { taxon => 900919, name => "IRGC38994", coverage => 13.6, },
        { taxon => 900920, name => "IRGC9060",  coverage => 12.13, },
        { taxon => 900921, name => "IRGC9062",  coverage => 13.11, },
        { taxon => 900922, name => "RA4952",    coverage => 13.98, },
        { taxon => 900923, name => "IRGC31856", coverage => 12.45, },
    );

    my @files = File::Find::Rule->file->name('*.vcf.fasta')->in($seq_dir);

    for my $item ( sort @data ) {
        my $name = $item->{name};

        # match the most similar name
        my ($file) = map { $_->[0] }
            sort { $b->[1] <=> $a->[1] }
            map { [ $_, compare( lc basename($_), lc $item->{name} ) ] } @files;
        $item->{seq} = $file;

        # prepare working dir
        my $dir = File::Spec->catdir( $data_dir, $name );
        mkdir $dir if !-e $dir;
        $item->{dir} = $dir;
    }

    print Dump \@data;

    my @data_with_exists = (
        {   name  => "nip_65",
            taxon => 39947,
            dir   => File::Spec->catdir( $data_dir, "nip_65" )
        },
        {   name  => "9311_65",
            taxon => 39946,
            dir   => File::Spec->catdir( $data_dir, "9311_65" )
        },
        @data
    );

    my $tt = Template->new;

    # taxon.csv
    my $text = <<'EOF';
[% FOREACH item IN data -%]
[% item.taxon %],Oryza,sativa,[% item.name %],[% item.name %],
[% END -%]
EOF
    $tt->process(
        \$text,
        { data => \@data, },
        File::Spec->catfile( $store_dir, "taxon.csv" )
    ) or die Template->error;

    # chr_length_chrUn.csv
    $text = <<'EOF';
[% FOREACH item IN data -%]
[% item.taxon %],chrUn,999999999,[% item.name %]/rice50
[% END -%]
EOF
    $tt->process(
        \$text,
        { data => \@data_with_exists, },
        File::Spec->catfile( $store_dir, "chr_length_chrUn.csv" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

if [ -f real_chr.csv ]; then
    rm real_chr.csv;
fi;

[% FOREACH item IN data -%]
perl -aln -F"\t" -e 'print qq{[% item.taxon %],$F[0],$F[1],[% item.name %]/rice50}' [% item.dir %]/chr.sizes >> real_chr.csv
[% END -%]

cat chr_length_chrUn.csv real_chr.csv > chr_length.csv
rm real_chr.csv

echo Run the following cmds to merge csv files
echo
echo perl [% pl_dir %]/alignDB/util/merge_csv.pl -t [% pl_dir %]/alignDB/init/taxon.csv -m [% data_dir %]/taxon.csv
echo
echo perl [% pl_dir %]/alignDB/util/merge_csv.pl -t [% pl_dir %]/alignDB/init/chr_length.csv -m [% data_dir %]/chr_length.csv
echo

EOF
    $tt->process(
        \$text,
        {   data     => \@data_with_exists,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
        },
        File::Spec->catfile( $store_dir, "real_chr.sh" )
    ) or die Template->error;

    # id2name.csv
    $text = <<'EOF';
[% FOREACH item IN data -%]
[% item.taxon %],[% item.name %]
[% END -%]
EOF
    $tt->process(
        \$text,
        { data => \@data_with_exists, },
        File::Spec->catfile( $store_dir, "id2name.csv" )
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
[% kentbin_dir %]/faSplit byname [% item.seq %] .
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
bsub -q mpi_2 -n 8 -J [% item.name %]-rm RepeatMasker [% item.dir %]/*.fasta -species rice -xsmall --parallel 8

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
# find failed rm jobs
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
find [% item.dir %] -name "*fasta" \
    | perl -e \
    'while(<>) {chomp; s/\.fasta$//; next if -e qq{$_.fasta.masked}; next if -e qq{$_.fa}; print qq{ bsub -n 8 -J [% item.name %]_ RepeatMasker $_.fasta -species rice -xsmall --parallel 8 \n};}' >> catchup.txt

[% END -%]

EOF

    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "rm_failed.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# blastz
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
bsub -q mpi_2 -n 8 -J [% item.name %]-bz perl [% pl_dir %]/blastz/bz.pl -dt [% data_dir %]/nip_65 -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/Nipvs[% item.name %] -s set01 -p 8 --noaxt -pb lastz --lastz

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
perl [% pl_dir %]/blastz/lpcna.pl -dt [% data_dir %]/nip_65 -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/Nipvs[% item.name %] -p 8

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
perl [% pl_dir %]/blastz/amp.pl -syn -dt [% data_dir %]/nip_65 -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/Nipvs[% item.name %] -p 8

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data => [ { name => "9311_65", taxon => 39946 }, @data ],
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
find [% data_dir %]/Nipvs[% item.name %]/axtNet -name "*.axt.gz" | xargs gzip -d
perl [% pl_dir %]/alignDB/extra/two_way_batch.pl -d Nipvs[% item.name %] -t="39947,Nip" -q "[% item.taxon %],[% item.name %]" -a [% data_dir %]/Nipvs[% item.name %] -at 10000 -st 0 --parallel 8 --run 1-3,21,40
gzip [% data_dir %]/Nipvs[% item.name %]/axtNet/*.axt

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data => [ { name => "9311_65", taxon => 39946 }, @data ],
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
cd [% data_dir %]/Nipvs[% item.name %]/

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
    my $data_dir = File::Spec->catdir( $ENV{HOME}, "data/alignment/rice" );
    my $pl_dir   = File::Spec->catdir( $ENV{HOME}, "Scripts" );

    my $tt         = Template->new;
    my $strains_of = {
        NipvsVII => [
            qw{ 9311_65 IRGC2540 IRGC38698 IRGC50448 IRGC26872 IRGC12793 RA4952 }
        ],
        NipvsXXIV => [
            qw{ 9311_65 IRGC1107 IRGC2540 IRGC27630 IRGC32399 IRGC418 IRGC55471
                IRGC8191 IRGC38698 IRGC11010 IRGC17757 IRGC328 IRGC43325
                IRGC43675 IRGC50448 IRGC66756 IRGC8244 IRGC26872 IRGC12793
                IRGC38994 IRGC9060 IRGC9062 RA4952 IRGC31856 }
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
    -d [% data_dir %]/Nipvs[% st FILTER ucfirst %] \
    [% END -%]
    --tree [% data_dir %]/25way.nwk \
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
perl [% pl_dir %]/alignDB/util/maf2fasta.pl \
    --has_outgroup --id 39947 -p 8 --block \
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
    -d [% item.out_dir %] -e nip_65 \
    --block --id [% data_dir %]/id2name.csv \
    -f [% data_dir %]/[% item.out_dir %]_mft  \
    -lt 5000 -st 0 -ct 0 --parallel [% parallel %] --run common

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
            parallel => $parallel,
        },
        File::Spec->catfile( $store_dir, "multi.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

if [ ! -d [% data_dir %]/phylo ]
then
    mkdir [% data_dir %]/phylo
fi

#----------------------------#
# concat
#----------------------------#
[% FOREACH item IN data -%]
[% IF item.strains.size > 3 -%]
# [% item.out_dir %]
# concat mafft fas to relaxed phylip
if [ ! -f [% data_dir %]/phylo/[% item.out_dir %].phy ]
then
    perl [% pl_dir %]/alignDB/util/concat_fasta.pl \
        -i [% data_dir %]/[% item.out_dir %]_mft  \
        -o [% data_dir %]/phylo/[% item.out_dir %].phy \
        -p
fi

[% END -%]
[% END -%]

cd [% data_dir %]/phylo

#----------------------------#
# phylo with raxml (ML + rapid bootstrap)
#----------------------------#
[% FOREACH item IN data -%]
[% IF item.strains.size > 3 -%]
# [% item.out_dir %]
if [ -f [% data_dir %]/phylo/[% item.out_dir %].phy.reduced ]
then
    raxml -T 6 -f a -m GTRGAMMA -p $RANDOM -N 100 -x $RANDOM -O \
        -o 9311_65 -n [% item.out_dir %] \
        -s [% data_dir %]/phylo/[% item.out_dir %].phy.reduced
elif [ -f [% data_dir %]/phylo/[% item.out_dir %].phy ]
then
    raxml -T 6 -f a -m GTRGAMMA -p $RANDOM -N 100 -x $RANDOM \
        -o 9311_65 -n [% item.out_dir %] \
        -s [% data_dir %]/phylo/[% item.out_dir %].phy
fi

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
        File::Spec->catfile( $store_dir, "phylo.sh" )
    ) or die Template->error;

}
