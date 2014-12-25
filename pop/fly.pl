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

my $group_name = 'fly65';
my $base_dir   = File::Spec->catdir( $ENV{HOME}, "data/alignment" );
my $data_dir   = File::Spec->catdir( $base_dir, $group_name );
my $phylo_tree = File::Spec->catfile( $data_dir, "fly_7way.nwk" );
my $pl_dir     = File::Spec->catdir( $ENV{HOME}, "Scripts" );

# ensembl genomes 12
my $fasta_dir = File::Spec->catdir( $ENV{HOME},
    "data/ensemblgenomes12_65/metazoa/fasta" );
my $mysql_dir = File::Spec->catdir( $ENV{HOME},
    "data/ensemblgenomes12_65/metazoa/mysql" );

my $tt = Template->new;

my @data = (
    {   taxon    => 7227,
        name     => "Dmel",
        sciname  => "Drosophila_melanogaster",
        coverage => "",
    },
    {   taxon    => 7240,
        name     => "Dsim",
        sciname  => "Drosophila_simulans",
        coverage => "",
    },
    {   taxon    => 7238,
        name     => "Dsech",
        sciname  => "Drosophila_sechellia",
        coverage => "",
    },
    {   taxon    => 7245,
        name     => "Dyak",
        sciname  => "Drosophila_yakuba",
        coverage => "",
    },
    {   taxon    => 7220,
        name     => "Dere",
        sciname  => "Drosophila_erecta",
        coverage => "",
    },
    {   taxon    => 7237,
        name     => "Dpse",
        sciname  => "Drosophila_pseudoobscura",
        coverage => "",
    },
    {   taxon    => 7234,
        name     => "Dper",
        sciname  => "Drosophila_persimilis",
        coverage => "",
    },
);

my @subdirs_fasta = File::Find::Rule->directory->in($fasta_dir);
my @subdirs_mysql = File::Find::Rule->directory->in($mysql_dir);

for my $item (@data) {

    # match the most similar name
    my ($fasta) = map { $_->[0] }
        sort { $b->[1] <=> $a->[1] }
        map { [ $_, compare( lc basename($_), $item->{sciname} . "/dna" ) ] }
        @subdirs_fasta;
    $item->{fasta} = $fasta;

    my ($mysql) = map { $_->[0] }
        sort { $b->[1] <=> $a->[1] }
        map { [ $_, compare( lc basename($_), $item->{sciname} . "_core" ) ] }
        @subdirs_mysql;
    $item->{mysql} = $mysql;

    $item->{db} = lc( $item->{name} ) . "_65";

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
find [% item.fasta %] -name "*dna.toplevel*" | xargs gzip -d -c > toplevel.fa
faops count toplevel.fa | perl -aln -e 'next if $F[0] eq 'total'; print $F[0] if $F[1] > 100000; print $F[0] if $F[1] > 10000  and $F[6]/$F[1] < 0.05' | uniq > listFile
faops some toplevel.fa listFile toplevel.filtered.fa
faops split-name toplevel.filtered.fa .
rm toplevel.fa toplevel.filtered.fa listFile

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
cd [% data_dir %]

#----------------------------------------------------------#
# Ensembl annotation or RepeatMasker
#----------------------------------------------------------#

[% FOREACH item IN data -%]
#----------------------------#
# [% item.name %] [% item.coverage %]
#----------------------------#
echo [% item.name %]

cd [% item.dir %]

if [ ! -f [% item.db %]_repeat.yml ]; then
    # check ensembl database existence
    result=$(mysql -s -N -ualignDB -palignDB -e "SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA WHERE SCHEMA_NAME='[% item.db %]'");
    if [ -z "$result" ]; then
        echo "[% item.db %] does not exists";
        perl [% pl_dir %]/alignDB/util/build_ensembl.pl --initdb --db [% item.db %] --ensembl [% item.mysql %];
    fi
fi;
if [ ! -f [% item.db %]_repeat.yml ]; then
    perl [% pl_dir %]/alignDB/util/write_masked_chr.pl -e [% item.db %] --parallel [% parallel %];
fi;
perl [% pl_dir %]/alignDB/util/write_masked_chr.pl -y [% item.db %]_repeat.yml --dir [% item.dir %] --parallel [% parallel %]

find . -name "*fa" | xargs rm
rename 's/\.masked//' *.fa.masked
rename 's/^/chr/' *.fa

if [ -f chrUn.fasta ]; then
    faops split-about [% item.dir %]/chrUn.fasta 100000000 [% item.dir %]/;
    rm [% item.dir %]/chrUn.fasta;    
    rename 's/fa$/fasta/' [0-9]*.fa;
fi;

if [ -f *.fasta ]; then
    RepeatMasker [% item.dir %]/*.fasta -species Drosophila -xsmall --parallel [% parallel %];
fi;
if [ -f *.fasta.masked ]; then
    rename 's/fasta.masked$/fa/' *.fasta.masked;
fi;
find [% item.dir %] -type f -name "*fasta*" | xargs rm

[% END -%]

EOF

$tt->process(
    \$text,
    {   data     => \@data,
        data_dir => $data_dir,
        pl_dir   => $pl_dir,
        parallel => $parallel,
    },
    File::Spec->catfile( $data_dir, "02_ensemblrm.sh" )
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

# Dmel_8way
cd [% data_dir %]
perl [% pl_dir %]/withncbi/taxon/strain_bz.pl \
    --file [% data_dir %]/[% group_name %].csv \
    -w     [% base_dir %] \
    --name [% group_name %] \
    --multi_name Dmel_8way \
    --use_name \
    --phylo_tree [% phylo_tree %] \
    --parallel [% parallel %]\
    -t Dmel \
    -q Dsim \
    -q Dsech \
    -q Dyak \
    -q Dere \
    -q Dpse \
    -q Dper

# Dmel_3way
cd [% data_dir %]
perl [% pl_dir %]/withncbi/taxon/strain_bz.pl \
    --file [% data_dir %]/[% group_name %].csv \
    -w     [% base_dir %] \
    --name [% group_name %] \
    --multi_name HumanvsCGOR \
    --use_name \
    --phylo_tree [% phylo_tree %] \
    --parallel [% parallel %]\
    -t Human \
    -q Chimp \
    -q Gorilla \
    -q Orangutan \
    -q Rhesus \
    -o Rhesus
    
EOF

$tt->process(
    \$text,
    {   data       => \@data,
        group_name => $group_name,
        data_dir   => $data_dir,
        base_dir   => $base_dir,
        pl_dir     => $pl_dir,
        phylo_tree => $phylo_tree,
        parallel   => $parallel,
    },
    File::Spec->catfile( $data_dir, "03_prepare.sh" )
) or die Template->error;

__END__

mkdir -p cd ~/data/alignment/fly65
cd ~/data/alignment/fly65

# download from http://hgdownload.soe.ucsc.edu/goldenPath/dm6/multiz27way/
~/share/phast/bin/tree_doctor dm6.27way.scientificNames.nh --newick \
    --rename "Drosophila_melanogaster -> Dmel ; Drosophila_simulans -> Dsim ; Drosophila_sechellia -> Dsech" \
    > temp1.nwk

~/share/phast/bin/tree_doctor temp1.nwk --newick \
    --rename "Drosophila_yakuba -> Dyak ; Drosophila_erecta -> Dere" \
    > temp2.nwk

~/share/phast/bin/tree_doctor temp2.nwk --newick \
    --rename "Drosophila_pseudoobscura_pseudoobscura -> Dpse ; Drosophila_persimilis -> Dper" \
    > temp3.nwk

~/share/phast/bin/tree_doctor temp3.nwk --newick \
    --prune-all-but Dmel,Dsim,Dsech,Dyak,Dere,Dpse,Dper \
    > fly_7way.nwk

rm temp[0-9].nwk

perl ~/Scripts/withncbi/pop/fly.pl 
sh 01_file.sh
sh 02_ensemblrm.sh

# execute 03_prepare.sh by copy & paste  

