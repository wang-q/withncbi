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
    || File::Spec->catdir( $ENV{HOME}, "data/alignment/arabidopsis19" );
my $parallel = 12;
{    # on linux
    my $data_dir
        = File::Spec->catdir( $ENV{HOME}, "data/alignment/arabidopsis19" );
    my $pl_dir      = File::Spec->catdir( $ENV{HOME}, "Scripts" );
    my $kentbin_dir = File::Spec->catdir( $ENV{HOME}, "bin/x86_64" );

    # nature 2011
    my $seq_dir
        = File::Spec->catdir( $ENV{HOME}, "data/1001/19genomes/fasta/MASKED" );

    my @data = (
        {   taxon    => 900201,
            name     => "Bur_0",
            coverage => 25,
            origin   => "Ireland"
        },
        {   taxon    => 900202,
            name     => "Can_0",
            coverage => 47,
            origin   => "Canary Isles"
        },
        { taxon => 900203, name => "Ct_1", coverage => 50, origin => "Italy" },
        {   taxon    => 900204,
            name     => "Edi_0",
            coverage => 52,
            origin   => "Scotland"
        },
        {   taxon    => 900205,
            name     => "Hi_0",
            coverage => 33,
            origin   => "Netherlands"
        },
        {   taxon    => 900206,
            name     => "Kn_0",
            coverage => 28,
            origin   => "Lithuania"
        },
        {   taxon    => 900207,
            name     => "Ler_0",
            coverage => 27,
            origin   => "Poland"
        },
        { taxon => 900208, name => "Mt_0", coverage => 30, origin => "Libya" },
        {   taxon    => 900209,
            name     => "No_0",
            coverage => 38,
            origin   => "Germany"
        },
        { taxon => 900210, name => "Oy_0", coverage => 54, origin => "Norway" },
        {   taxon    => 900211,
            name     => "Po_0",
            coverage => 41,
            origin   => "Germany"
        },
        {   taxon    => 900212,
            name     => "Rsch_4",
            coverage => 38,
            origin   => "Russia"
        },
        { taxon => 900213, name => "Sf_2",  coverage => 40, origin => "Spain" },
        { taxon => 900214, name => "Tsu_0", coverage => 48, origin => "Japan" },
        {   taxon    => 900215,
            name     => "Wil_2",
            coverage => 40,
            origin   => "Russia"
        },
        { taxon => 900216, name => "Ws_0", coverage => 33, origin => "Russia" },
        {   taxon    => 900217,
            name     => "Wu_0",
            coverage => 26,
            origin   => "Germany"
        },
        {   taxon    => 900218,
            name     => "Zu_0",
            coverage => 31,
            origin   => "Germany"
        },
    );

    my @files = File::Find::Rule->file->name('*.fas')->in($seq_dir);

    for my $item ( sort @data ) {

        # match the most similar name
        my ($file) = map { $_->[0] }
            sort { $b->[1] <=> $a->[1] }
            map { [ $_, compare( lc basename($_), lc $item->{name} ) ] } @files;
        $item->{file} = $file;

        # prepare working dir
        my $dir = File::Spec->catdir( $data_dir, $item->{name} );
        mkdir $dir if !-e $dir;
        $item->{dir} = $dir;
    }

    my $basecount = File::Spec->catfile( $data_dir, "basecount.txt" );
    remove( \1, $basecount ) if -e $basecount;

    print Dump \@data;

    my @data_with_exists = (
        {   name  => "ath_65",
            taxon => 3702,
            dir   => File::Spec->catdir( $data_dir, "ath_65" )
        },
        {   name  => "lyrata_65",
            taxon => 59689,
            dir   => File::Spec->catdir( $data_dir, "lyrata_65" )
        },
        @data
    );

    my $tt = Template->new;

    # taxon.csv
    my $text = <<'EOF';
[% FOREACH item IN data -%]
[% item.taxon %],Arabidopsis,thaliana,[% item.name %],[% item.name %],
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
[% item.taxon %],chrUn,999999999,[% item.name %]/arabidopsis19/nature2011
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
perl -aln -F"\t" -e 'print qq{[% item.taxon %],$F[0],$F[1],[% item.name %]/arabidopsis19/nature2011}' [% item.dir %]/chr.sizes >> real_chr.csv
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

    #
    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# basecount and split
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.origin %]
echo [% item.name %]
[% kentbin_dir %]/faSplit byname [% item.file %] [% item.dir %]/

# uncovered regions are masked by lowercase
perl -p -i -e '/>/ and next; s/[a-z]/n/g' [% item.dir %]/*.fa

echo [% item.name %] >> [% data_dir %]/basecount.txt
[% kentbin_dir %]/faCount [% item.dir %]/*.fa >> [% data_dir %]/basecount.txt
echo >> [% data_dir %]/basecount.txt

rename 's/fa$/fasta/' [% item.dir %]/*.fa

[% END -%]

EOF

    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir,
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
RepeatMasker [% item.dir %]/*.fasta -species arabidopsis -xsmall --parallel [% parallel %]

[% END -%]

EOF

    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir,
            parallel    => $parallel,
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
    'while(<>) {chomp; s/\.fasta$//; next if -e qq{$_.fasta.masked}; next if -e qq{$_.fa}; print qq{ RepeatMasker $_.fasta -species arabidopsis -xsmall --parallel [% parallel %] \n};}' >> catchup.txt

[% END -%]

EOF

    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir,
            parallel    => $parallel,
        },
        File::Spec->catfile( $store_dir, "rm_failed.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# RepeatMasker
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
echo [% item.name %]

for i in [% item.dir %]/*.fasta;
do
    if [ -f $i.masked ];
    then
        rename 's/fasta.masked$/fa/' $i.masked;
        find [% item.dir %] -type f -name "`basename $i`*" | xargs rm 
    fi;
done;

[% END -%]

echo Please check the following files
find [% data_dir %] -name "*.fasta"

EOF

    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir,
            parallel    => $parallel,
        },
        File::Spec->catfile( $store_dir, "clean-rm.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# blastz
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
perl [% pl_dir %]/blastz/bz.pl \
    -dt [% data_dir %]/ath_65 \
    -dq [% data_dir %]/[% item.name %] \
    -dl [% data_dir %]/Athvs[% item.name FILTER ucfirst %] \
    -s set01 -p [% parallel %] --noaxt -pb lastz --lastz

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data        => [ { name => "lyrata_65" }, @data ],
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir,
            parallel    => $parallel,
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
perl [% pl_dir %]/blastz/lpcna.pl \
    -dt [% data_dir %]/ath_65 \
    -dq [% data_dir %]/[% item.name %] \
    -dl [% data_dir %]/Athvs[% item.name FILTER ucfirst %] \
    -p [% parallel %]

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data        => [ { name => "lyrata_65" }, @data ],
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir,
            parallel    => $parallel,
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
perl [% pl_dir %]/blastz/amp.pl -syn \
    -dt [% data_dir %]/ath_65 \
    -dq [% data_dir %]/[% item.name %] \
    -dl [% data_dir %]/Athvs[% item.name FILTER ucfirst %] \
    -p [% parallel %]

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data        => [ { name => "lyrata_65" }, @data ],
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir,
            parallel    => $parallel,
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
perl [% pl_dir %]/alignDB/extra/two_way_batch.pl \
    -d Athvs[% item.name FILTER ucfirst %] \
    -t="3702,Ath" -q "[% item.taxon %],[% item.name %]" \
    -da [% data_dir %]/Athvs[% item.name FILTER ucfirst %] \
    -lt 10000 -st 0 -ct 0 --parallel [% parallel %] \
    --run 1-5,21,40

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data => [ { name => "lyrata_65", taxon => 59689 }, @data ],
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
            parallel => $parallel,
        },
        File::Spec->catfile( $store_dir, "pair_stat.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash

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
# find [% data_dir %] -name "mafSynNet" | xargs rm -fr
# find [% data_dir %] -name "mafNet" | xargs rm -fr

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

{    # on windows
    my $data_dir = "d:/data/alignment/arabidopsis19";
    my $pl_dir   = "d:/wq/Scripts";

    my $tt = Template->new;

    my @data = (
        { taxon => 900201, name => "Bur_0",  coverage => 25, },
        { taxon => 900202, name => "Can_0",  coverage => 47, },
        { taxon => 900203, name => "Ct_1",   coverage => 50, },
        { taxon => 900204, name => "Edi_0",  coverage => 52, },
        { taxon => 900205, name => "Hi_0",   coverage => 33, },
        { taxon => 900206, name => "Kn_0",   coverage => 28, },
        { taxon => 900207, name => "Ler_0",  coverage => 27, },
        { taxon => 900208, name => "Mt_0",   coverage => 30, },
        { taxon => 900209, name => "No_0",   coverage => 38, },
        { taxon => 900210, name => "Oy_0",   coverage => 54, },
        { taxon => 900211, name => "Po_0",   coverage => 41, },
        { taxon => 900212, name => "Rsch_4", coverage => 38, },
        { taxon => 900213, name => "Sf_2",   coverage => 40, },
        { taxon => 900214, name => "Tsu_0",  coverage => 48, },
        { taxon => 900215, name => "Wil_2",  coverage => 40, },
        { taxon => 900216, name => "Ws_0",   coverage => 33, },
        { taxon => 900217, name => "Wu_0",   coverage => 26, },
        { taxon => 900218, name => "Zu_0",   coverage => 31, },
    );

    my $text = <<'EOF';
cd /d [% data_dir %]

REM #----------------------------#
REM # multi
REM #----------------------------#
perl [% pl_dir %]/alignDB/extra/join_dbs.pl --dbs [% dbs %] --goal_db [% goal_db %] --outgroup [% outgroup %] --target [% target %] --queries [% queries %] --no_insert=1 --trimmed_fasta=1 --length 1000

perl [% pl_dir %]/alignDB/extra/multi_way_batch.pl -d [% goal_db %] -e ath_65 -da [% data_dir %]/[% goal_db %]  -lt 1000 -st 100000 --parallel 4 --run all

EOF

    my @names = ( "Lyrata", map { $_->{name} } @data );
    my $dbs = join ',', map { "Athvs" . $_ } @names;
    my $queries = join ',', map { $_ . "query" } ( 1 .. scalar @names - 1 );
    $tt->process(
        \$text,
        {   goal_db  => "Athvs19",
            outgroup => '0query',
            target   => '0target',
            dbs      => $dbs,
            queries  => $queries,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
        },
        File::Spec->catfile( $store_dir, "joins.bat" )
    ) or die Template->error;
}

{    # multiz
    my $data_dir
        = File::Spec->catdir( $ENV{HOME}, "data/alignment/arabidopsis19" );
    my $pl_dir = File::Spec->catdir( $ENV{HOME}, "Scripts" );

    my $tt         = Template->new;
    my $strains_of = {
        AthvsV   => [qw{ lyrata_65 Bur_0 Zu_0 No_0 Ler_0  }],
        AthvsXIX => [
            qw{ lyrata_65 Bur_0 Can_0 Ct_1 Edi_0 Hi_0 Kn_0 Ler_0 Mt_0 No_0 Oy_0 Po_0
                Rsch_4 Sf_2 Tsu_0 Wil_2 Ws_0 Wu_0 Zu_0 }
        ],
        AthvsXVIII_pop => [
            qw{ Bur_0 Can_0 Ct_1 Edi_0 Hi_0 Kn_0 Ler_0 Mt_0 No_0 Oy_0 Po_0
                Rsch_4 Sf_2 Tsu_0 Wil_2 Ws_0 Wu_0 Zu_0 }
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
perl [% pl_dir %]/blastz/mz.pl \
[% FOREACH st IN item.strains -%]
    -d [% data_dir %]/Athvs[% st FILTER ucfirst %] \
[% END -%]
    --tree [% data_dir %]/20way.nwk \
    --out [% data_dir %]/[% item.out_dir %] \
    -syn -p [% parallel %]

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
            parallel => $parallel,
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
    -p [% parallel %] --block \
    -i [% data_dir %]/[% item.out_dir %] \
    -o [% data_dir %]/[% item.out_dir %]_fasta

[% END -%]

#----------------------------#
# mafft
#----------------------------#
[% FOREACH item IN data -%]
# [% item.out_dir %]
perl [% pl_dir %]/blastz/refine_fasta.pl \
    --msa mafft --block -p [% parallel %] \
    --block [% UNLESS item.out_dir.match('_pop$') %]--outgroup[% END %] \
    -i [% data_dir %]/[% item.out_dir %]_fasta \
    -o [% data_dir %]/[% item.out_dir %]_mft

[% END -%]

#----------------------------#
# muscle
#----------------------------#
#[% FOREACH item IN data -%]
## [% item.out_dir %]
#perl [% pl_dir %]/blastz/refine_fasta.pl \
#    --msa muscle --block -p [% parallel %] \
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
            parallel => $parallel,
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
    -d [% item.out_dir %] -e ath_65 \
    --block [% UNLESS item.out_dir.match('_pop$') %]--outgroup[% END %] \
    --id [% data_dir %]/id2name.csv \
    -da [% data_dir %]/[% item.out_dir %]_mft  \
    -lt 5000 --parallel [% parallel %] --run common

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
    perl [% pl_dir %]/blastz/concat_fasta.pl \
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
        -o lyrata_65 -n [% item.out_dir %] \
        -s [% data_dir %]/phylo/[% item.out_dir %].phy.reduced
elif [ -f [% data_dir %]/phylo/[% item.out_dir %].phy ]
then
    raxml -T 6 -f a -m GTRGAMMA -p $RANDOM -N 100 -x $RANDOM \
        -o lyrata_65 -n [% item.out_dir %] \
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
