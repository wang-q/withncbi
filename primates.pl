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
    || File::Spec->catdir( $ENV{HOME}, "data/alignment/primates" );
my $parallel = 12;
{    # on linux
    my $data_dir    = File::Spec->catdir( $ENV{HOME}, "data/alignment/primates" );
    my $pl_dir      = File::Spec->catdir( $ENV{HOME}, "Scripts" );
    my $kentbin_dir = File::Spec->catdir( $ENV{HOME}, "bin/x86_64" );

    # ensembl 65
    my $fasta_dir = File::Spec->catdir( $ENV{HOME}, "data/ensembl65" );
    my $mysql_dir = File::Spec->catdir( $ENV{HOME}, "data/ensembl65" );

    my $tt = Template->new;

    my @data = (
        {   taxon    => 9606,
            name     => "human",
            sciname  => "homo_sapiens",
            coverage => "best",
        },
        {   taxon    => 9598,
            name     => "chimp",
            sciname  => "pan_troglodytes",
            coverage => "6x sanger",
        },
        {   taxon    => 9595,
            name     => "gorilla",
            sciname  => "gorilla_gorilla",
            coverage => "2.1x sanger; 35x solexa",
        },
        {   taxon    => 9601,
            name     => "orangutan",
            sciname  => "pongo_abelii",
            coverage => "6x sanger",
        },
        {   taxon    => 61853,
            name     => "gibbon",
            sciname  => "nomascus_leucogenys",
            coverage => "5.6x sanger",
        },
        {   taxon    => 9544,
            name     => "rhesus",
            sciname  => "macaca_mulatta",
            coverage => "6.1x sanger",
        },
        {   taxon    => 9483,
            name     => "marmoset",
            sciname  => "callithrix_jacchus",
            coverage => "6x sanger",
        },
        {   taxon    => 9478,
            name     => "tarsier",
            sciname  => "tarsius_syrichta",
            coverage => "1.82x sanger",
        },
        {   taxon    => 30608,
            name     => "lemur",
            sciname  => "microcebus_murinus",
            coverage => "1.93x sanger",
        },
        {   taxon    => 30611,
            name     => "bushbaby",
            sciname  => "otolemur_garnettii",
            coverage => "2x sanger",
        },
    );

    my @subdirs_fasta = File::Find::Rule->directory->in($fasta_dir);
    my @subdirs_mysql = File::Find::Rule->directory->in($mysql_dir);

    for my $item (@data) {

        # match the most similar name
        my ($fasta) = map { $_->[0] }
            sort { $b->[1] <=> $a->[1] }
            map {
            [ $_, compare( lc basename($_), $item->{sciname} . "/dna" ) ]
            } @subdirs_fasta;
        $item->{fasta} = $fasta;

        my ($mysql) = map { $_->[0] }
            sort { $b->[1] <=> $a->[1] }
            map {
            [ $_, compare( lc basename($_), $item->{sciname} . "_core" ) ]
            } @subdirs_mysql;
        $item->{mysql} = $mysql;

        $item->{db} = $item->{name} . "_65";

        # prepare working dir
        my $dir = File::Spec->catdir( $data_dir, $item->{name} );
        mkdir $dir if !-e $dir;
        $item->{dir} = $dir;
    }

    #my $basecount = File::Spec->catfile( $data_dir, "basecount.txt" );
    #remove( \1, $basecount ) if -e $basecount;

    # taxon.csv
    my $text = <<'EOF';
[% FOREACH item IN data -%]
[% item.taxon %],[% item.sciname FILTER ucfirst FILTER replace('_', ',') %],[% item.name %],[% item.name %],
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
[% item.taxon %],chrUn,999999999,[% item.name %]/ensembl65
[% END -%]
EOF
    $tt->process(
        \$text,
        { data => \@data, },
        File::Spec->catfile( $store_dir, "chr_length_chrUn.csv" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

if [ -f real_chr.csv ]; then
    rm real_chr.csv;
fi;

[% FOREACH item IN data -%]
perl -aln -F"\t" -e 'print qq{[% item.taxon %],$F[0],$F[1],[% item.name %]/ensembl65}' [% item.dir %]/chr.sizes >> real_chr.csv
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
        {   data     => \@data,
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
        { data => \@data, },
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
find [% item.fasta %] -name "*dna.toplevel*" | xargs gzip -d -c > toplevel.fa
[% kentbin_dir %]/faCount toplevel.fa | perl -aln -e 'next if $F[0] eq 'total'; print $F[0] if $F[1] > 100000; print $F[0] if $F[1] > 10000  and $F[6]/$F[1] < 0.05' | uniq > listFile
[% kentbin_dir %]/faSomeRecords toplevel.fa listFile toplevel.filtered.fa
[% kentbin_dir %]/faSplit byname toplevel.filtered.fa .
rm toplevel.fa toplevel.filtered.fa listFile

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
# Ensembl annotation or RepeatMasker
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
echo [% item.name %]

cd [% item.dir %]

if [ ! -f [% item.db %]_repeat.yml ]; then perl [% pl_dir %]/alignDB/util/build_ensembl.pl --initdb --db [% item.db %] --ensembl [% item.mysql %];  fi;
if [ ! -f [% item.db %]_repeat.yml ]; then perl [% pl_dir %]/alignDB/util/write_masked_chr.pl -e [% item.db %]; fi;
perl [% pl_dir %]/alignDB/util/write_masked_chr.pl -y [% item.db %]_repeat.yml --dir [% item.dir %]

find . -name "*fa" | xargs rm
~/perl5/bin/rename 's/\.masked//' *.fa.masked
~/perl5/bin/rename 's/^/chr/' *.fa

if [ -f chrUn.fasta ];
then
    [% kentbin_dir %]/faSplit about [% item.dir %]/chrUn.fasta 100000000 [% item.dir %]/;
    rm [% item.dir %]/chrUn.fasta;    
    ~/perl5/bin/rename 's/fa$/fasta/' [0-9][0-9].fa;
fi;

bsub -n 8 -J [% item.name %]-RM RepeatMasker [% item.dir %]/*.fasta -species Primates -xsmall --parallel 8;
if [ -f *.fasta.masked ];
then
    ~/perl5/bin/rename 's/fasta.masked$/fa/' *.fasta.masked;
fi;
find [% item.dir %] -type f -name "*fasta*" | xargs rm 

[% END -%]

EOF

    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "ensemblrm.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# repeatmasker on all fasta
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
# for i in [% item.dir %]/*.fasta; do bsub -n 8 -J [% item.name %]_`basename $i .fasta` RepeatMasker $i -species Primates -xsmall --parallel 8; done;
bsub -n 8 -J [% item.name %]-rm RepeatMasker [% item.dir %]/*.fasta -species Primates -xsmall --parallel 8

[% END -%]

#----------------------------#
# find failed rm jobs
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
perl -e 'for $i (0..30) { $i = sprintf qq{%02d}, $i; $str = qq{[% item.dir %]/$i}; next if ! -e qq{$str.fasta}; next if -e qq{$str.fasta.masked}; next if -e qq{$str.fa}; print qq{ bsub -n 8 -J [% item.name %]_$i RepeatMasker $str.fasta -species Primates -xsmall --parallel 8 \n};}' >> catchup.txt

[% END -%]

# find [% data_dir %] -name "*.fasta.masked" | sed "s/\.fasta\.masked$//" | xargs -i echo mv {}.fasta.masked {}.fa | sh
# find [% data_dir %] -type f | grep -v fa$ | xargs rm -fr

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
# blastz
#----------------------------#
[% FOREACH item IN data -%]
[% IF item.name != 'human' -%]
# [% item.name %] [% item.coverage %]
# To avoid memory overflow, set parallel to 6
bsub -n 8 -J [% item.name %]-bz perl [% pl_dir %]/blastz/bz.pl -dt [% data_dir %]/human -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/Humanvs[% item.name FILTER ucfirst %] -s set01 -p 6 --noaxt -pb lastz --lastz

[% END -%]
[% END -%]

#----------------------------#
# lpcna
#----------------------------#
[% FOREACH item IN data -%]
[% IF item.name != 'human' -%]
# [% item.name %] [% item.coverage %]
bsub -n 8 -J [% item.name %]-axt perl [% pl_dir %]/blastz/lpcna.pl -dt [% data_dir %]/human -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/Humanvs[% item.name FILTER ucfirst %] -p 8

[% END -%]
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

#----------------------------#
# clean pairwise maf
#----------------------------#
find [% data_dir %] -name "mafSynNet" | xargs rm -fr
find [% data_dir %] -name "mafNet" | xargs rm -fr
    
#----------------------------#
# only keeps chr.2bit files
#----------------------------#
# find [% data_dir %] -name "*.fa" | xargs rm
# find [% data_dir %] -name "*.fasta*" | xargs rm

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

    $text = <<'EOF';
#!/bin/bash
    
#----------------------------#
# amp
#----------------------------#
[% FOREACH item IN data -%]
[% IF item.name != 'human' -%]
# [% item.name %] [% item.coverage %]
perl [% pl_dir %]/blastz/amp.pl -syn -dt [% data_dir %]/human -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/Humanvs[% item.name FILTER ucfirst %] -p 8

[% END -%]
[% END -%]

EOF
    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "amp.sh" )
    ) or die Template->error;
}

{    # on windows
    my $data_dir = File::Spec->catdir("d:/data/alignment/primates");
    my $pl_dir   = File::Spec->catdir("d:/wq/Scripts");

    my $tt = Template->new;

    my @data = (
        { taxon => 9606,  name => "human", },
        { taxon => 9598,  name => "chimp", },
        { taxon => 9595,  name => "gorilla", },
        { taxon => 9601,  name => "orangutan", },
        { taxon => 61853, name => "gibbon", },
        { taxon => 9544,  name => "rhesus", },
        { taxon => 9483,  name => "marmoset", },
        { taxon => 9478,  name => "tarsier", },
        { taxon => 30608, name => "lemur", },
        { taxon => 30611, name => "bushbaby", },
    );

    my $text = <<'EOF';
cd /d [% data_dir %]

#----------------------------#
# stat
#----------------------------#
[% FOREACH item IN data -%]
[% IF item.name != 'human' -%]
# [% item.name %]
perl [% pl_dir %]/alignDB/extra/two_way_batch.pl -d Humanvs[% item.name FILTER ucfirst %] -t="9606,human" -q "[% item.taxon %],[% item.name %]" -a [% data_dir %]/Humanvs[% item.name FILTER ucfirst %] -at 10000 -st 10000000 --parallel 4 --run 1-3,21,40

[% END -%]
[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
        },
        File::Spec->catfile( $store_dir, "stat.bat" )
    ) or die Template->error;
}

{    # multiz
    my $data_dir = File::Spec->catdir( $ENV{HOME}, "data/alignment/primates" );
    my $pl_dir   = File::Spec->catdir( $ENV{HOME}, "Scripts" );

    my $tt         = Template->new;
    my $strains_of = {
        HumanvsCR   => [qw{ chimp rhesus }],
        HumanvsCGOR => [qw{ chimp gorilla orangutan rhesus }],
        HumanvsVI   => [qw{ chimp gorilla orangutan gibbon rhesus marmoset }],
        HumanvsVIII => [
            qw{ chimp gorilla orangutan gibbon rhesus marmoset lemur bushbaby}],
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
    -d [% data_dir %]/Humanvs[% st FILTER ucfirst %] \
    [% END -%]
    --tree [% data_dir %]/primates_11way.nwk \
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
    --has_outgroup --id 9606 -p 8 --block \
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
# muscle-quick
#----------------------------#
#[% FOREACH item IN data -%]
## [% item.out_dir %]
#bsub -q mpi_2 -n 8 -J [% item.out_dir %]-msl perl [% pl_dir %]/alignDB/util/refine_fasta.pl \
#    --msa muscle --quick --block -p 8 \
#    -i [% data_dir %]/[% item.out_dir %]_fasta \
#    -o [% data_dir %]/[% item.out_dir %]_mslq
#
#[% END -%]

#----------------------------#
# clean
#----------------------------#
[% FOREACH item IN data -%]
# [% item.out_dir %]
cd [% data_dir %]
rm -fr [% item.out_dir %]_fasta

[% END -%]

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
    -d [% item.out_dir %] -e human_65 \
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
}
