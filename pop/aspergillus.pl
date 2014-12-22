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

my $group_name = 'aspergillus65';
my $base_dir   = File::Spec->catdir( $ENV{HOME}, "data/alignment" );
my $data_dir   = File::Spec->catdir( $base_dir, $group_name );
my $phylo_tree = File::Spec->catfile( $data_dir, "aspergillus_8way.nwk" );
my $pl_dir     = File::Spec->catdir( $ENV{HOME}, "Scripts" );

# ensembl genomes 12
my $fasta_dir
    = File::Spec->catdir( $ENV{HOME}, "data/ensemblgenomes12_65/fungi/fasta" );
my $mysql_dir
    = File::Spec->catdir( $ENV{HOME}, "data/ensemblgenomes12_65/fungi/mysql" );

my $tt = Template->new;

my @data = (
    {   taxon     => 330879,
        name      => "Afum",
        sciname   => "Aspergillus_fumigatus",
        othername => "Aspergillus fumigatus Af293",
        coverage  => "10.5x sanger",
    },
    {   taxon     => 5062,
        name      => "Aory",
        sciname   => "Aspergillus_oryzae",
        othername => "Eurotium nidulans",
        coverage  => "9x sanger",
    },
    {   taxon    => 5057,
        name     => "Acla",
        sciname  => "Aspergillus_clavatus",
        coverage => "11.4x sanger",
    },
    {   taxon    => 5059,
        name     => "Afla",
        sciname  => "Aspergillus_flavus",
        coverage => "5x sanger",
    },
    {   taxon     => 162425,
        name      => "Anid",
        sciname   => "Aspergillus_nidulans",
        othername => "Emericella nidulans",
        coverage  => "13x sanger",
    },
    {   taxon    => 5061,
        name     => "Anig",
        sciname  => "Aspergillus_niger",
        coverage => "7.5x sanger",
    },
    {   taxon    => 33178,
        name     => "Ater",
        sciname  => "Aspergillus_terreus",
        coverage => "11.05x sanger",
    },
    {   taxon     => 36630,
        name      => "Nfis",
        sciname   => "Neosartorya_fischeri",
        othername => "Aspergillus fischeri",
        coverage  => "11.0x sanger",
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

[% IF item.name != 'Afum' and item.name != 'Aory' -%]
rename 's/fa$/fasta/' *.fa
[% END -%]

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

#----------------------------#
# Ensembl annotation or RepeatMasker
#----------------------------#
[% FOREACH item IN data -%]
[% IF item.name == 'Afum' or item.name == 'Aory' -%]
# [% item.name %] [% item.coverage %]
echo [% item.name %]

cd [% item.dir %]

if [ ! -f [% item.db %]_repeat.yml ]; then perl [% pl_dir %]/alignDB/util/build_ensembl.pl --initdb --db [% item.db %] --ensembl [% item.mysql %];  fi;
if [ ! -f [% item.db %]_repeat.yml ]; then perl [% pl_dir %]/alignDB/util/write_masked_chr.pl -e [% item.db %]; fi;
perl [% pl_dir %]/alignDB/util/write_masked_chr.pl -y [% item.db %]_repeat.yml --dir [% item.dir %]

find . -name "*fa" | xargs rm
rename 's/\.masked//' *.fa.masked
rename 's/^/chr/' *.fa

if [ -f chrUn.fasta ];
then
    [% kentbin_dir %]/faSplit about [% item.dir %]/chrUn.fasta 100000000 [% item.dir %]/;
    rm [% item.dir %]/chrUn.fasta;    
    rename 's/fa$/fasta/' [0-9][0-9].fa;
fi;

RepeatMasker [% item.dir %]/*.fasta -species Fungi -xsmall --parallel [% parallel %]
if [ -f *.fasta.masked ];
then
    rename 's/fasta.masked$/fa/' *.fasta.masked;
fi;
find [% item.dir %] -type f -name "*fasta*" | xargs rm 

[% ELSE %]
# [% item.name %] [% item.coverage %]
echo [% item.name %]

cd [% item.dir %]
RepeatMasker [% item.dir %]/*.fasta -species Fungi -xsmall --parallel [% parallel %]
rename 's/fasta.masked$/fa/' *.fasta.masked
find [% item.dir %]  -type f -name "*fasta*" | xargs rm

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
    File::Spec->catfile( $data_dir, "02_ensemblrm.sh" )
) or die Template->error;

$text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# blastz Afum
#----------------------------#
[% FOREACH item IN data -%]
[% IF item.name != 'Afum' -%]
# [% item.name %] [% item.coverage %]
perl [% pl_dir %]/blastz/bz.pl \
    -dt [% data_dir %]/Afum -dq [% data_dir %]/[% item.name %] \
    -dl [% data_dir %]/Afumvs[% item.name %] \
    -s set01 -p [% parallel %] --noaxt -pb lastz --lastz

[% END -%]
[% END -%]

#----------------------------#
# blastz Aory
#----------------------------#
[% FOREACH item IN data -%]
[% IF item.name != 'Aory' -%]
# [% item.name %] [% item.coverage %]
perl [% pl_dir %]/blastz/bz.pl \
    -dt [% data_dir %]/Aory -dq [% data_dir %]/[% item.name %] \
    -dl [% data_dir %]/Aoryvs[% item.name %] \
    -s set01 -p [% parallel %] --noaxt -pb lastz --lastz

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
    File::Spec->catfile( $data_dir, "bz.sh" )
) or die Template->error;

{    # multiz
    my $data_dir
        = File::Spec->catdir( $ENV{HOME}, "data/alignment/aspergillus" );
    my $pl_dir = File::Spec->catdir( $ENV{HOME}, "Scripts" );

    my $tt = Template->new;
    my $strains_of
        = { AfumvsVII => [qw{ Acla Afla Anid Anig Aory Ater Nfis }], };

}
