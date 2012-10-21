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
    || File::Spec->catdir( $ENV{HOME}, "data/alignment/aspergillus" );

{    # on linux
    my $data_dir
        = File::Spec->catdir( $ENV{HOME}, "data/alignment/aspergillus" );
    my $pl_dir      = File::Spec->catdir( $ENV{HOME}, "Scripts" );
    my $kentbin_dir = File::Spec->catdir( $ENV{HOME}, "bin/x86_64" );

    # ensembl genomes 65
    my $fasta_dir = File::Spec->catdir( $ENV{HOME},
        "data/ensemblgenomes12_65/fungi/fasta" );

    my $tt = Template->new;

    my @data = (
        {   taxon    => 5057,
            name     => "Acla",
            sciname  => "Aspergillus clavatus",
            coverage => "11.4x sanger",
        },
        {   taxon    => 5059,
            name     => "Afla",
            sciname  => "Aspergillus flavus",
            coverage => "5x sanger",
        },
        {   taxon     => 330879,
            name      => "Afum",
            sciname   => "Aspergillus fumigatus",
            othername => "Aspergillus fumigatus Af293",
            coverage  => "10.5x sanger",
        },
        {   taxon     => 162425,
            name      => "Anid",
            sciname   => "Aspergillus nidulans",
            othername => "Emericella nidulans",
            coverage  => "13x sanger",
        },
        {   taxon    => 5061,
            name     => "Anig",
            sciname  => "Aspergillus niger",
            coverage => "7.5x sanger",
        },
        {   taxon     => 5062,
            name      => "Aory",
            sciname   => "Aspergillus oryzae",
            othername => "Eurotium nidulans",
            coverage  => "9x sanger",
        },
        {   taxon    => 33178,
            name     => "Ater",
            sciname  => "Aspergillus terreus",
            coverage => "11.05x sanger",
        },
        {   taxon     => 36630,
            name      => "Nfis",
            sciname   => "Neosartorya fischeri",
            othername => "Aspergillus fischeri",
            coverage  => "11.0x sanger",
        },
    );

    my @subdirs_fasta = File::Find::Rule->directory->in($fasta_dir);

    for my $item (@data) {
        my $folder = $item->{sciname};
        $folder =~ s/ /_/g;
        $folder = lc $folder;

        #$folder .= "/dna";

        # match the most similar name
        my ($fasta) = map { $_->[0] }
            sort { $b->[1] <=> $a->[1] }
            map { [ $_, compare( lc basename($_), $folder ) ] } @subdirs_fasta;
        $item->{fasta} = $fasta;

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
[% item.taxon %],[% item.sciname FILTER replace(' ', ',') %],[% item.name %],,
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
[% item.taxon %],chrUn,999999999,[% item.name %]/ensemblgenomes65
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

~/perl5/bin/rename 's/fa$/fasta/' *.fa

[% IF item.name = 'Aory' -%]
~/perl5/bin/rename 's/^/chr/' *.fasta
perl -p -i -e '/>/ and s/\>/\>chr/' *.fasta
[% IF item.name != 'Aory' -%]

[% END -%]

EOF

    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "auto_aspergillus_file.sh" )
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
bsub -n 8 -J [% item.name %]-RM RepeatMasker [% item.dir %]/*.fasta -species Fungi -xsmall --parallel 8;

[% END -%]

# find [% data_dir %] -name "*.fasta.masked" | sed "s/\.fasta\.masked$//" | xargs -i echo mv {}.fasta.masked {}.fa | sh
# find [% data_dir %] -type f -name "*fasta*" | xargs rm 

EOF

    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "auto_aspergillus_rm.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# blastz
#----------------------------#
[% FOREACH item IN data -%]
[% IF item.name != 'Aory' -%]
# [% item.name %] [% item.coverage %]
bsub -n 8 -J [% item.name %]-bz perl [% pl_dir %]/blastz/bz.pl -dt [% data_dir %]/Aory -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/Aoryvs[% item.name %] -s set01 -p 8 --noaxt -pb lastz --lastz

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
        File::Spec->catfile( $store_dir, "auto_aspergillus_bz_aory.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# lpcna
#----------------------------#
[% FOREACH item IN data -%]
[% IF item.name != 'Aory' -%]
# [% item.name %] [% item.coverage %]
perl [% pl_dir %]/blastz/lpcna.pl -dt [% data_dir %]/Aory -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/Aoryvs[% item.name %] -p 8

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
        File::Spec->catfile( $store_dir, "auto_aspergillus_lpcna_aory.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
    
#----------------------------#
# amp
#----------------------------#
[% FOREACH item IN data -%]
[% IF item.name != 'Aory' -%]
# [% item.name %] [% item.coverage %]
perl [% pl_dir %]/blastz/amp.pl -syn -dt [% data_dir %]/Aory -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/Aoryvs[% item.name %] -p 8

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
        File::Spec->catfile( $store_dir, "auto_aspergillus_amp_aory.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# stat
#----------------------------#
[% FOREACH item IN data -%]
[% IF item.name != 'Aory' -%]
# [% item.name %]
find [% data_dir %]/Aoryvs[% item.name %]/axtNet -name "*.axt.gz" | xargs gzip -d
perl [% pl_dir %]/alignDB/extra/two_way_batch.pl -d Aoryvs[% item.name %] -t="5062,Aory" -q "[% item.taxon %],[% item.name %]" -a [% data_dir %]/Aoryvs[% item.name %] -at 10000 -st 1000000 --parallel 8 --run 1-3,21,40
gzip [% data_dir %]/Aoryvs[% item.name %]/axtNet/*.axt

[% END -%]
[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
        },
        File::Spec->catfile( $store_dir, "auto_aspergillus_stat_aory.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
    
#----------------------------#
# 2bit to fa
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
cd [% data_dir %]/[% item.name %]

[% kentbin_dir %]/twoBitToFa [% data_dir %]/[% item.name %]/chr.2bit [% data_dir %]/[% item.name %]/[% item.name %].fa
[% kentbin_dir %]/faSplit byname [% data_dir %]/[% item.name %]/[% item.name %].fa .
rm [% data_dir %]/[% item.name %]/[% item.name %].fa

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
# find [% data_dir %] -name "mafSynNet" | xargs rm -fr
# find [% data_dir %] -name "mafNet" | xargs rm -fr

#----------------------------#
# gzip maf, fas
#----------------------------#
# find [% data_dir %] -name "*.maf" | parallel gzip
# find [% data_dir %] -name "*.maf.fas" | parallel gzip

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
        File::Spec->catfile( $store_dir, "auto_aspergillus_clean_aory.sh" )
    ) or die Template->error;
}

{    # multiz
    my $data_dir
        = File::Spec->catdir( $ENV{HOME}, "data/alignment/aspergillus" );
    my $pl_dir = File::Spec->catdir( $ENV{HOME}, "Scripts" );

    my $tt = Template->new;
    my $strains_of
        = { AoryvsVII => [qw{ Acla Afla Afum Anid Anig Aory Ater Nfis }], };

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
find . -name "*MT.synNet*" | xargs rm
find . -name "47Q*.synNet*" | xargs rm

[% FOREACH item IN data -%]
# [% item.out_dir %]
bsub -q mpi_2 -n 8 -J [% item.out_dir %]-mz perl [% pl_dir %]/blastz/mz.pl \
    [% FOREACH st IN item.strains -%]
    -d [% data_dir %]/Aoryvs[% st %] \
    [% END -%]
    --tree [% data_dir %]/8way.nwk \
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
        File::Spec->catfile( $store_dir, "auto_aspergillus_mz_aory.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#----------------------------#
# maf2fasta
#----------------------------#
[% FOREACH item IN data -%]
# [% item.out_dir %]
perl [% pl_dir %]/alignDB/util/maf2fasta.pl \
    --has_outgroup --id 5062 -p 8 --block \
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
        File::Spec->catfile( $store_dir, "auto_aspergillus_maf_fasta_aory.sh" )
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
    -d [% item.out_dir %] -e Aory_65 \
    --block --id 5062 \
    -f [% data_dir %]/[% item.out_dir %]_mft  \
    -lt 5000 -st 1000000 --parallel 8 --run 1-3,21,40

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
        },
        File::Spec->catfile( $store_dir, "auto_aspergillus_multi_aory.sh" )
    ) or die Template->error;
}
