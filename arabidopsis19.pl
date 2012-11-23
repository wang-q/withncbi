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
        { taxon => 900201, name => "Bur_0",  origin => "Ireland" },
        { taxon => 900202, name => "Can_0",  origin => "Canary Isles" },
        { taxon => 900203, name => "Ct_1",   origin => "Italy" },
        { taxon => 900204, name => "Edi_0",  origin => "Scotland" },
        { taxon => 900205, name => "Hi_0",   origin => "Netherlands" },
        { taxon => 900206, name => "Kn_0",   origin => "Lithuania" },
        { taxon => 900207, name => "Ler_0",  origin => "Poland" },
        { taxon => 900208, name => "Mt_0",   origin => "Libya" },
        { taxon => 900209, name => "No_0",   origin => "Germany" },
        { taxon => 900210, name => "Oy_0",   origin => "Norway" },
        { taxon => 900211, name => "Po_0",   origin => "Germany" },
        { taxon => 900212, name => "Rsch_4", origin => "Russia" },
        { taxon => 900213, name => "Sf_2",   origin => "Spain" },
        { taxon => 900214, name => "Tsu_0",  origin => "Japan" },
        { taxon => 900215, name => "Wil_2",  origin => "Russia" },
        { taxon => 900216, name => "Ws_0",   origin => "Russia" },
        { taxon => 900217, name => "Wu_0",   origin => "Germany" },
        { taxon => 900218, name => "Zu_0",   origin => "Germany" },
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
    my $tt = Template->new;

    # taxon.csv
    my $text = <<'EOF';
[% FOREACH item IN data -%]
[% item.taxon %],Arabidopsis,thaliana,[% item.name %],,
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
[% item.taxon %],chrUn,999999999,[% item.name %]/arabidopsis19/1001
[% END -%]
EOF
    $tt->process(
        \$text,
        { data => \@data, },
        File::Spec->catfile( $store_dir, "chr_length.csv" )
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

~/perl5/bin/rename 's/fa$/fasta/' [% item.dir %]/*.fa
# find [% item.dir %] -name "*.fa" | sed "s/\.fa$//" | xargs -i echo mv {}.fa {}.fasta | sh

[% END -%]

EOF

    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "auto_ath19_file.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# repeatmasker on all fasta
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
bsub -q mpi_2 -n 8 -J [% item.name %]-rm RepeatMasker [% item.dir %]/*.fasta -species arabidopsis -xsmall --parallel 8

[% END -%]

#----------------------------#
# find failed rm jobs
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
find [% data_dir %]/[% item.name %] -name "*fasta" \
    | perl -e \
    'while(<>) {chomp; s/\.fasta$//; next if -e qq{$_.fasta.masked}; next if -e qq{$_.fa}; print qq{ bsub -n 8 -J [% item.name %]_ RepeatMasker $_.fasta -species arabidopsis -xsmall --parallel 8 \n};}' >> catchup.txt

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
        File::Spec->catfile( $store_dir, "auto_ath19_rm.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# blastz
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.origin %]
bsub -q mpi_2 -n 8 -J [% item.name %]-bz perl [% pl_dir %]/blastz/bz.pl \
    -dt [% data_dir %]/ath_65 \
    -dq [% data_dir %]/[% item.name %] \
    -dl [% data_dir %]/Athvs[% item.name %] \
    -s set01 -p 8 --noaxt -pb lastz --lastz 

[% END -%]

#----------------------------#
# lpcna
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.origin %]
perl [% pl_dir %]/blastz/lpcna.pl \
    -dt [% data_dir %]/ath_65 \
    -dq [% data_dir %]/[% item.name %] \
    -dl [% data_dir %]/Athvs[% item.name %] \
    -p 8

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "auto_ath19_bz.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
    
#----------------------------#
# amp
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
perl [% pl_dir %]/blastz/amp.pl -syn -dt [% data_dir %]/ath_65 -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/Athvs[% item.name FILTER ucfirst %] -p 8

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data        => [ { name => "lyrata_65" }, @data ],
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "auto_ath19_amp.sh" )
    ) or die Template->error;

    $text = <<'EOF';
cd [% data_dir %]

#----------------------------#
# stat
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
perl [% pl_dir %]/alignDB/extra/two_way_batch.pl -d Athvs[% item.name %] -t="3702,Ath" -q "[% item.taxon %],[% item.name %]" -a [% data_dir %]/Athvs[% item.name %] -at 10000 -st 1000000 --parallel 4 --run 1-3,21,40

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data => [ { name => "lyrata_65", taxon => 59689 }, @data ],
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
        },
        File::Spec->catfile( $store_dir, "auto_ath19_stat.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
    
#----------------------------#
# tar-gzip
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
cd [% data_dir %]/Athvs[% item.name %]/

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
find [% data_dir %] -name "*.maf" | xargs gzip
find [% data_dir %] -name "*.maf.fas" | xargs gzip

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
        File::Spec->catfile( $store_dir, "auto_ath19_clean.sh" )
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

perl [% pl_dir %]/alignDB/extra/multi_way_batch.pl -d [% goal_db %] -e ath_65 -f [% data_dir %]/[% goal_db %]  -lt 1000 -st 100000 --parallel 4 --run all

EOF

    my @names = ( "Lyrata", map { $_->{name} } @data );
    my $dbs = join ',', map { "Athvs" . $_ } @names;
    my $queries = join ',', map { $_ . "query" } ( 1 .. scalar @names - 1 );
    $tt->process(
        \$text,
        {   goal_db  => "AthvsNineteen",
            outgroup => '0query',
            target   => '0target',
            dbs      => $dbs,
            queries  => $queries,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
        },
        File::Spec->catfile( $store_dir, "auto_ath19_joins.bat" )
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
    -d [% data_dir %]/Athvs[% st FILTER ucfirst %] \
    [% END -%]
    --tree [% data_dir %]/20way.nwk \
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
        File::Spec->catfile( $store_dir, "auto_ath19_mz.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#----------------------------#
# maf2fasta
#----------------------------#
[% FOREACH item IN data -%]
# [% item.out_dir %]
perl [% pl_dir %]/alignDB/util/maf2fasta.pl \
    --has_outgroup --id 3702 -p 8 --block \
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
        File::Spec->catfile( $store_dir, "auto_ath19_maf_fasta.sh" )
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
    --block --id 3702 \
    -f [% data_dir %]/[% item.out_dir %]_mft  \
    -lt 5000 -st 1000000 --parallel 4 --run 1-3,21,40

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
        },
        File::Spec->catfile( $store_dir, "auto_ath19_multi.sh" )
    ) or die Template->error;
}
