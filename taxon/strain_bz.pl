#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use DBI;
use Text::CSV_XS;
use DateTime::Format::Natural;
use List::MoreUtils qw(any all uniq);
use Template;

use IPC::Cmd qw(can_run);
use Path::Tiny;
use File::Find::Rule;

use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::RealBin/../config.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

=head1 NAME

strain_bz.pl - Full procedure for multiple genome alignments.

=head1 SYNOPSIS

    perl strain_bz.pl [options]
      Options:
        --help          -?          brief help message
        --working_dir   -w  STR     Default is [.]
        --seq_dir       -s  STR     Will do prep_fa() from this dir or use seqs store in $working_dir
        --keep              @STR    don't touch anything inside these fasta files
        --target        -t  STR
        --queries       -q  @STR
        --outgroup      -o  STR
        --csv_taxon     -c  STR     All taxons in this project (may also contain unused taxons)
        --name_str      -n  STR     Default is []
        --phylo_tree    -p  STR     Predefined phylogenetic tree
        --multi_name    -m  STR     Naming multiply alignment, the default value is $name_str
                                    This option is for more than one align combination.
        --msa               STR     Aligning program for refine. Default is [mafft]
        --norm                      RepeatMasker has been done.
        --nostat                    Don't do stat stuffs
        --norawphylo                Skip rawphylo
        --parallel          INT     number of child processes

=cut

my $aligndb   = path( $Config->{run}{aligndb} )->stringify;
my $egaz      = path( $Config->{run}{egaz} )->stringify;
my $fig_table = path( $Config->{run}{fig_table} )->stringify;

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'working_dir|w=s' => \( my $working_dir = "." ),
    'seq_dir|s=s'     => \my $seq_dir,
    'keep=s'          => \my @keep,
    'target|t=s'      => \my $target,
    'queries|q=s'     => \my @queries,
    'outgroup|o|r=s'  => \my $outgroup,
    'csv_taxon|c=s'   => \( my $csv_taxon   = "strains_taxon_info.csv" ),
    'name_str|n=s'    => \( my $name_str    = "working" ),
    'phylo_tree|p=s'  => \my $phylo_tree,
    'multi_name|m=s'  => \my $multi_name,
    'msa=s'           => \( my $msa         = 'mafft' ),
    'norm'            => \my $norm,
    'nostat'          => \my $nostat,
    'norawphylo'      => \my $norawphylo,
    'parallel=i'      => \( my $parallel    = $Config->{run}{parallel} ),
) or Getopt::Long::HelpMessage(1);

if ( defined $phylo_tree ) {
    if ( !-e $phylo_tree ) {
        warn "$phylo_tree does not exists. Unset it.\n";
        undef $phylo_tree;
    }
}

if ( !defined $multi_name ) {
    $multi_name = $name_str;
}

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Writing strains summary...");

# prepare working dir
{
    print "Working on $name_str\n";
    $working_dir = path( $working_dir, $name_str )->absolute;
    $working_dir->mkpath;
    $working_dir = $working_dir->stringify;
    print " " x 4, "Working dir is $working_dir\n";

    path( $working_dir, 'Genomes' )->mkpath;
    path( $working_dir, 'Pairwise' )->mkpath;
}

# move $outgroup to last
if ($outgroup) {
    my ($exist) = grep { $_ eq $outgroup } @queries;
    if ( !defined $exist ) {
        die "outgroup does not exist!\n";
    }

    @queries = grep { $_ ne $outgroup } @queries;
    push @queries, $outgroup;
}

# build basic information
my @data;
{
    my %taxon_of;
    if ($csv_taxon) {
        for my $line ( path($csv_taxon)->lines ) {
            my @fields = split /,/, $line;
            if ( $#fields >= 2 ) {
                $taxon_of{ $fields[0] } = $fields[1];
            }
        }
    }
    @data = map {
        {   name  => $_,
            taxon => exists $taxon_of{$_} ? $taxon_of{$_} : 0,
            dir   => path( $working_dir, 'Genomes', $_ )->stringify,
        }
    } ( $target, @queries );
}

# if seqs is not in working dir, copy them from seq_dir
if ($seq_dir) {
    print "Get seqs from [$seq_dir]\n";

    for my $id ( $target, @queries ) {
        my ($keep) = grep { $_ eq $id } @keep;
        print " " x 4 . "Copy seq of [$id]\n";
        if ( defined $keep ) {
            print " " x 8 . "Don't change fasta header for [$id]\n";
        }

        my $original_dir = path( $seq_dir, $id )->stringify;
        my $cur_dir = path( $working_dir, 'Genomes', $id );
        $cur_dir->mkpath;
        $cur_dir = $cur_dir->stringify;

        my @fa_files
            = File::Find::Rule->file->name( '*.fna', '*.fa', '*.fas',
            '*.fasta' )->in($original_dir);

        printf " " x 8 . "Total %d fasta file(s)\n", scalar @fa_files;

        for my $fa_file (@fa_files) {
            my $basename;
            if ( defined $keep ) {
                $basename = prep_fa( $fa_file, $cur_dir, 1 );
            }
            else {
                $basename = prep_fa( $fa_file, $cur_dir );
            }

            my $gff_file = path( $original_dir, "$basename.gff" );
            if ( $gff_file->is_file ) {
                $gff_file->copy($cur_dir);
            }
            my $rm_gff_file = path( $original_dir, "$basename.rm.gff" );
            if ( $rm_gff_file->is_file ) {
                $rm_gff_file->copy($cur_dir);
            }
        }
    }
}

# If there's no phylo tree, generate a fake one.
if ( !defined $phylo_tree ) {
    print "Create fake_tree.nwk\n";
    my $fh = path( $working_dir, "fake_tree.nwk" )->openw;
    print {$fh} "(" x scalar(@queries) . "$target";
    for my $id (@queries) {
        print {$fh} ",$id)";
    }
    print {$fh} ";\n";
    close $fh;
}

{
    my $tt = Template->new;
    my $text;
    my $sh_name;

    #----------------------------#
    # all *.sh files
    #----------------------------#

    # real_chr.sh
    $sh_name = "1_real_chr.sh";
    print "Create $sh_name\n";
    $text = <<'EOF';
#!/bin/bash
cd [% working_dir %]

sleep 1;

cat << DELIMITER > chrUn.csv
common_name,taxon_id,chr,length,assembly
[% FOREACH item IN data -%]
[% item.name %],[% item.taxon %],chrUn,999999999,
[% END -%]
DELIMITER

if [ -f real_chr.csv ]; then
    rm real_chr.csv;
fi;

[% FOREACH item IN data -%]
# [% item.name %]
faops size [% item.dir %]/*.fa > [% item.dir %]/chr.sizes;
perl -aln -F"\t" -e 'print qq{[% item.name %],[% item.taxon %],$F[0],$F[1],}' [% item.dir %]/chr.sizes >> real_chr.csv;

[% END -%]

cat chrUn.csv real_chr.csv > chr_length.csv
echo "chr_length.csv generated."

rm chrUn.csv
rm real_chr.csv

EOF
    $tt->process(
        \$text,
        {   data        => \@data,
            working_dir => $working_dir,
        },
        path( $working_dir, $sh_name )->stringify
    ) or die Template->error;

    # file-rm.sh
    if ( !$norm ) {
        $sh_name = "2_file_rm.sh";
        print "Create $sh_name\n";
        $text = <<'EOF';
#!/bin/bash
cd [% working_dir %]

sleep 1;

#----------------------------#
# repeatmasker on all fasta
#----------------------------#
[% FOREACH item IN data -%]

for f in `find [% item.dir%] -name "*.fa"` ; do
    rename 's/fa$/fasta/' $f ;
done

for f in `find [% item.dir%] -name "*.fasta"` ; do
    RepeatMasker $f -xsmall --parallel [% parallel %] ;
done

for f in `find [% item.dir%] -name "*.fasta.out"` ; do
    rmOutToGFF3.pl $f > `dirname $f`/`basename $f .fasta.out`.rm.gff;
done

for f in `find [% item.dir%] -name "*.fasta"` ; do
    if [ -f $f.masked ];
    then
        rename 's/fasta.masked$/fa/' $f.masked;
        find [% item.dir%] -type f -name "`basename $f`*" | parallel --no-run-if-empty rm;
    else
        rename 's/fasta$/fa/' $f;
        echo `date` "RepeatMasker on $f failed.\n" >> RepeatMasker.log
        find [% item.dir%] -type f -name "`basename $f`*" | parallel --no-run-if-empty rm;
    fi;
done;

[% END %]

EOF

        $tt->process(
            \$text,
            {   data        => \@data,
                parallel    => $parallel,
                working_dir => $working_dir,
                target      => $target,
                queries     => \@queries,
            },
            path( $working_dir, $sh_name )->stringify
        ) or die Template->error;
    }

    # pair_cmd.sh
    $sh_name = "3_pair_cmd.sh";
    print "Create $sh_name\n";
    $text = <<'EOF';
#!/bin/bash
# strain_bz.pl
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]

sleep 1;

#----------------------------#
# z_batch
#----------------------------#
[% FOREACH q IN queries -%]
perl [% egaz %]/z_batch.pl \
    -dt [% working_dir %]/Genomes/[% target %] \
    -dq [% working_dir %]/Genomes/[% q %] \
    -dw [% working_dir %]/Pairwise \
    -r 2-4 \
    --clean \
    --parallel [% parallel %]

[% END -%]

EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            egaz        => $egaz,
            nostat      => $nostat,
            target      => $target,
            outgroup    => $outgroup,
            queries     => \@queries,
        },
        path( $working_dir, $sh_name )->stringify
    ) or die Template->error;

    # rawphylo.sh
    if ( !$norawphylo and !defined $phylo_tree ) {
        $sh_name = "4_rawphylo.sh";
        print "Create $sh_name\n";
        $text = <<'EOF';
#!/bin/bash
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]

sleep 1;

#----------------------------#
# Clean previous directories
#----------------------------#
if [ -d [% working_dir %]/[% multi_name %]_raw ]; then
    rm -fr [% working_dir %]/[% multi_name %]_raw;
fi;
mkdir -p [% working_dir %]/[% multi_name %]_raw;

#----------------------------#
# maf2fas
#----------------------------#
echo "==> Convert maf to fas"

[% FOREACH q IN queries -%]
mkdir -p [% working_dir %]/[% multi_name %]_raw/[% target %]vs[% q %]
find [% working_dir %]/Pairwise/[% target %]vs[% q %] -name "*.maf" -or -name "*.maf.gz" \
    | parallel --no-run-if-empty -j 1 \
        fasops maf2fas {} -o [% working_dir %]/[% multi_name %]_raw/[% target %]vs[% q %]/{/}.fas
sleep 1;
fasops covers \
    [% working_dir %]/[% multi_name %]_raw/[% target %]vs[% q %]/*.fas \
    -n [% target %] -l 1000 -t 10 \
    -o [% working_dir %]/[% multi_name %]_raw/[% target %]vs[% q %].yml

[% END -%]

#----------------------------#
# intersect
#----------------------------#
echo "==> Intersect"

runlist compare --op intersect \
[% FOREACH q IN queries -%]
    [% working_dir %]/[% multi_name %]_raw/[% target %]vs[% q %].yml \
[% END -%]
    -o [% working_dir %]/[% multi_name %]_raw/intersect.raw.yml

runlist span --op excise -n 1000 \
    [% working_dir %]/[% multi_name %]_raw/intersect.raw.yml \
    -o [% working_dir %]/[% multi_name %]_raw/intersect.filter.yml

#----------------------------#
# slicing
#----------------------------#
echo "==> Slicing with intersect"

[% FOREACH q IN queries -%]
if [ -e [% working_dir %]/[% multi_name %]_raw/[% target %]vs[% q %].slice.fas ];
then
    rm [% working_dir %]/[% multi_name %]_raw/[% target %]vs[% q %].slice.fas
fi
find [% working_dir %]/[% multi_name %]_raw/[% target %]vs[% q %]/ -name "*.fas" -or -name "*.fas.gz" \
    | sort \
    | parallel --no-run-if-empty --keep-order -j 1 " \
        fasops slice {} \
            [% working_dir %]/[% multi_name %]_raw/intersect.filter.yml \
            -n [% target %] -l 1000 -o stdout \
            >> [% working_dir %]/[% multi_name %]_raw/[% target %]vs[% q %].slice.fas
    "

[% END -%]

#----------------------------#
# join
#----------------------------#
echo "==> Join intersects"

fasops join \
[% FOREACH q IN queries -%]
    [% working_dir %]/[% multi_name %]_raw/[% target %]vs[% q %].slice.fas \
[% END -%]
    -n [% target %] \
    -o [% working_dir %]/[% multi_name %]_raw/join.raw.fas

echo [% target %] > [% working_dir %]/[% multi_name %]_raw/names.list
[% FOREACH q IN queries -%]
echo [% q %] >> [% working_dir %]/[% multi_name %]_raw/names.list
[% END -%]

fasops subset \
    [% working_dir %]/[% multi_name %]_raw/join.raw.fas \
    [% working_dir %]/[% multi_name %]_raw/names.list \
    --required \
    -o [% working_dir %]/[% multi_name %]_raw/join.filter.fas

fasops refine \
    --msa mafft --parallel [% parallel %] \
    [% working_dir %]/[% multi_name %]_raw/join.filter.fas \
    -o [% working_dir %]/[% multi_name %]_raw/join.refine.fas

#----------------------------#
# RAxML
#----------------------------#
# raw phylo guiding tree
if [ -d [% working_dir %]/[% multi_name %]_rawphylo ]; then
    rm -fr [% working_dir %]/[% multi_name %]_rawphylo;
fi;
mkdir -p [% working_dir %]/[% multi_name %]_rawphylo;

cd [% working_dir %]/[% multi_name %]_rawphylo

[% IF queries.size > 2 -%]
perl [% egaz%]/concat_fasta.pl \
    -i [% working_dir %]/[% multi_name %]_raw/join.refine.fas \
    -o [% working_dir %]/[% multi_name %]_rawphylo/[% multi_name %].phy \
    --sampling --total 10_000_000 --relaxed

[% IF avx -%]
raxmlHPC-PTHREADS-AVX -T [% IF parallel > 8 %] 8 [% ELSIF parallel > 3 %] [% parallel - 1 %] [% ELSE %] 2 [% END %] \
    -f a -m GTRGAMMA -p $(openssl rand 3 | od -DAn) -N 100 -x $(openssl rand 3 | od -DAn) \
[% IF outgroup -%]
    -o [% outgroup %] \
[% END -%]
    --no-bfgs -n [% multi_name %] -s [% working_dir %]/[% multi_name %]_rawphylo/[% multi_name %].phy
[% ELSE -%]
raxmlHPC-PTHREADS -T [% IF parallel > 8 %] 8 [% ELSIF parallel > 3 %] [% parallel - 1 %] [% ELSE %] 2 [% END %] \
    -f a -m GTRGAMMA -p $(openssl rand 3 | od -DAn) -N 100 -x $(openssl rand 3 | od -DAn) \
[% IF outgroup -%]
    -o [% outgroup %] \
[% END -%]
    --no-bfgs -n [% multi_name %] -s [% working_dir %]/[% multi_name %]_rawphylo/[% multi_name %].phy
[% END -%]

cp [% working_dir %]/[% multi_name %]_rawphylo/RAxML_best* [% working_dir %]/[% multi_name %]_rawphylo/[% multi_name %].nwk

[% ELSIF queries.size == 2 -%]
echo "(([% target %],[% queries.0 %]),[% queries.1 %]);" > [% working_dir %]/[% multi_name %]_rawphylo/[% multi_name %].nwk

[% ELSE -%]

echo "([% target %],[% queries.0 %]);" > [% working_dir %]/[% multi_name %]_rawphylo/[% multi_name %].nwk

[% END -%]

EOF
        $tt->process(
            \$text,
            {   stopwatch   => $stopwatch,
                parallel    => $parallel,
                working_dir => $working_dir,
                aligndb     => $aligndb,
                egaz        => $egaz,
                target      => $target,
                outgroup    => $outgroup,
                queries     => \@queries,
                multi_name  => $multi_name,
                nostat      => $nostat,
                avx         => can_run('raxmlHPC-PTHREADS-AVX'),
            },
            path( $working_dir, $sh_name )->stringify
        ) or die Template->error;
    }

    # multi_cmd.sh
    $sh_name = "5_multi_cmd.sh";
    print "Create $sh_name\n";
    $text = <<'EOF';
#!/bin/bash
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]

sleep 1;

#----------------------------#
# Clean previous directories
#----------------------------#
if [ -d [% working_dir %]/[% multi_name %]_mz ]; then
    rm -fr [% working_dir %]/[% multi_name %]_mz;
fi;
mkdir -p [% working_dir %]/[% multi_name %]_mz;

if [ -d [% working_dir %]/[% multi_name %]_fasta ]; then
    rm -fr [% working_dir %]/[% multi_name %]_fasta;
fi;
mkdir -p [% working_dir %]/[% multi_name %]_fasta;

if [ -d [% working_dir %]/[% multi_name %]_refined ]; then
    rm -fr [% working_dir %]/[% multi_name %]_refined;
fi;
mkdir -p [% working_dir %]/[% multi_name %]_refined;

if [ -d [% working_dir %]/[% multi_name %]_phylo ]; then
    rm -fr [% working_dir %]/[% multi_name %]_phylo;
fi;
mkdir -p [% working_dir %]/[% multi_name %]_phylo;

#----------------------------#
# mz
#----------------------------#
[% IF phylo_tree -%]
perl [% egaz %]/mz.pl \
    [% FOREACH id IN queries -%]
    -d [% working_dir %]/Pairwise/[% target %]vs[% id %] \
    [% END -%]
    --tree [% phylo_tree %] \
    --out [% working_dir %]/[% multi_name %]_mz \
    -p [% parallel %]
[% ELSE %]
if [ -f [% working_dir %]/[% multi_name %]_rawphylo/[% multi_name %].nwk ]
then
    perl [% egaz %]/mz.pl \
        [% FOREACH id IN queries -%]
        -d [% working_dir %]/Pairwise/[% target %]vs[% id %] \
        [% END -%]
        --tree [% working_dir %]/[% multi_name %]_rawphylo/[% multi_name %].nwk \
        --out [% working_dir %]/[% multi_name %]_mz \
        -p [% parallel %]
else
    perl [% egaz %]/mz.pl \
        [% FOREACH id IN queries -%]
        -d [% working_dir %]/Pairwise/[% target %]vs[% id %] \
        [% END -%]
        --tree [% working_dir %]/fake_tree.nwk \
        --out [% working_dir %]/[% multi_name %]_mz \
        -p [% parallel %]
fi
[% END -%]

find [% working_dir %]/[% multi_name %]_mz -type f -name "*.maf" | parallel --no-run-if-empty -j [% parallel %] gzip

#----------------------------#
# maf2fas
#----------------------------#
echo "Convert maf to fas"
find [% working_dir %]/[% multi_name %]_mz -name "*.maf" -or -name "*.maf.gz" \
    | parallel --no-run-if-empty -j [% parallel %] \
        fasops maf2fas {} -o [% working_dir %]/[% multi_name %]_fasta/{/}.fas

#----------------------------#
# refine fasta
#----------------------------#
echo "Refine fasta"
find [% working_dir %]/[% multi_name %]_fasta -name "*.fas" -or -name "*.fas.gz" \
    | parallel --no-run-if-empty -j [% parallel %] \
        fasops refine {} \
        --msa [% msa %] --parallel [% parallel %] \
        --quick --expand 100 --join 100 \
[% IF outgroup -%]
        --outgroup \
[% END -%]
        -o [% working_dir %]/[% multi_name %]_refined/{/}

find [% working_dir %]/[% multi_name %]_refined -type f -name "*.fas" | parallel -j [% parallel %] gzip

#----------------------------#
# RAxML
#----------------------------#
[% IF queries.size > 2 -%]
cd [% working_dir %]/[% multi_name %]_phylo

perl [% egaz %]/concat_fasta.pl \
    -i [% working_dir %]/[% multi_name %]_refined \
    -o [% working_dir %]/[% multi_name %]_phylo/[% multi_name %].phy \
    --sampling --total 10_000_000 --relaxed

find [% working_dir %]/[% multi_name %]_phylo -type f -name "RAxML*" | parallel --no-run-if-empty rm

[% IF avx -%]
raxmlHPC-PTHREADS-AVX -T [% IF parallel > 8 %] 8 [% ELSIF parallel > 3 %] [% parallel - 1 %] [% ELSE %] 2 [% END %] \
    -f a -m GTRGAMMA -p $(openssl rand 3 | od -DAn) -N 100 -x $(openssl rand 3 | od -DAn) \
[% IF outgroup -%]
    -o [% outgroup %] \
[% END -%]
    --no-bfgs -n [% multi_name %] -s [% working_dir %]/[% multi_name %]_phylo/[% multi_name %].phy
[% ELSE -%]
raxmlHPC-PTHREADS -T [% IF parallel > 8 %] 8 [% ELSIF parallel > 3 %] [% parallel - 1 %] [% ELSE %] 2 [% END %] \
    -f a -m GTRGAMMA -p $(openssl rand 3 | od -DAn) -N 100 -x $(openssl rand 3 | od -DAn) \
[% IF outgroup -%]
    -o [% outgroup %] \
[% END -%]
    --no-bfgs -n [% multi_name %] -s [% working_dir %]/[% multi_name %]_phylo/[% multi_name %].phy
[% END -%]

cp [% working_dir %]/[% multi_name %]_phylo/RAxML_bipartitions.* [% working_dir %]/[% multi_name %]_phylo/[% multi_name %].nwk

nw_display -s -b 'visibility:hidden' [% working_dir %]/[% multi_name %]_phylo/[% multi_name %].nwk > [% working_dir %]/[% multi_name %]_phylo/[% multi_name %].svg

[% END -%]

EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            egaz        => $egaz,
            target      => $target,
            outgroup    => $outgroup,
            queries     => \@queries,
            phylo_tree  => $phylo_tree,
            multi_name  => $multi_name,
            msa         => $msa,
            avx         => can_run('raxmlHPC-PTHREADS-AVX'),
        },
        path( $working_dir, $sh_name )->stringify
    ) or die Template->error;

    # var_list.sh
    $sh_name = "6_var_list.sh";
    print "Create $sh_name\n";
    $text = <<'EOF';
#!/bin/bash
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]

sleep 1;

if [ -d [% working_dir %]/[% multi_name %]_vcf ]; then
    rm -fr [% working_dir %]/[% multi_name %]_vcf;
fi;

mkdir -p [% working_dir %]/[% multi_name %]_vcf

#----------------------------#
# var_list
#----------------------------#
find [% working_dir %]/[% multi_name %]_refined -type f -name "*.fas" -or -name "*.fas.gz" \
    | parallel --no-run-if-empty basename {} \
    | parallel --no-run-if-empty -j [% parallel %] \
        perl [% egaz %]/fas2vcf.pl \
            -s [% working_dir %]/Genomes/[% target %]/chr.sizes \
            -i [% working_dir %]/[% multi_name %]_refined/{} \
            -o [% working_dir %]/[% multi_name %]_vcf/{}.vcf

EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            egaz        => $egaz,
            target      => $target,
            multi_name  => $multi_name,
        },
        path( $working_dir, $sh_name )->stringify
    ) or die Template->error;

    if ( !$nostat ) {
        $sh_name = "7_multi_db_only.sh";
        print "Create $sh_name\n";
        $text = <<'EOF';
#!/bin/bash
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]

sleep 1;

#----------------------------#
# Create anno.yml
#----------------------------#
perl [% aligndb %]/util/gff2anno.pl \
    --type CDS --remove \
    [% working_dir %]/Genomes/[% target %]/*.gff \
    > [% working_dir %]/cds.yml

perl [% aligndb %]/util/gff2anno.pl \
    --remove \
    [% working_dir %]/Genomes/[% target %]/*.rm.gff \
    > [% working_dir %]/repeat.yml

runlist merge [% working_dir %]/repeat.yml [% working_dir %]/cds.yml -o [% working_dir %]/anno.yml
rm [% working_dir %]/repeat.yml [% working_dir %]/cds.yml

#----------------------------#
# multi_way_batch
#----------------------------#
perl [% aligndb %]/util/multi_way_batch.pl \
    -d [% multi_name %] \
    -da [% working_dir %]/[% multi_name %]_refined \
    -a [% working_dir %]/anno.yml \
[% IF outgroup -%]
    --outgroup \
[% END -%]
    -chr [% working_dir %]/chr_length.csv \
    -lt 1000 --parallel [% parallel %] --batch 5 \
    --run 1,2,5,10,21,30-32,40-42,44

EOF
        $tt->process(
            \$text,
            {   stopwatch   => $stopwatch,
                parallel    => $parallel,
                working_dir => $working_dir,
                aligndb     => $aligndb,
                target      => $target,
                outgroup    => $outgroup,
                queries     => \@queries,
                multi_name  => $multi_name,
            },
            path( $working_dir, $sh_name )->stringify
        ) or die Template->error;

        # chart.sh
        print "Create chart.sh\n";
        $text = <<'EOF';
# basicstat
perl [% fig_table %]/collect_common_basic.pl -d .

EOF
        $tt->process(
            \$text,
            {   stopwatch   => $stopwatch,
                working_dir => $working_dir,
                fig_table   => $fig_table,
            },
            path( $working_dir, "chart.sh" )->stringify
        ) or die Template->error;
    }

    # pack_it_up.sh
    $sh_name = "9_pack_it_up.sh";
    print "Create $sh_name\n";
    $text = <<'EOF';
#!/bin/bash
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]

sleep 1;

find . -type f \
    | grep -v -E "\.(sh|2bit|bat)$" \
    | grep -v -E "(chr_length|id2name)\.csv$" \
    | grep -v -F "fake_tree.nwk" \
    > file_list.txt

tar -czvf [% multi_name %].tar.gz -T file_list.txt

EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            multi_name  => $multi_name,
        },
        path( $working_dir, $sh_name )->stringify
    ) or die Template->error;

    # message
    $stopwatch->block_message("Execute *.sh files in order.");
}

#----------------------------#
# Finish
#----------------------------#
$stopwatch->end_message;
exit;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#
sub prep_fa {
    my $infile = shift;
    my $dir    = shift;
    my $keep   = shift;

    my $basename = path($infile)->basename( '.fna', '.fa', '.fas', '.fasta' );
    my $in_fh    = path($infile)->openr;
    my $out_fh   = path( $dir, "$basename.fa" )->openw;
    while (<$in_fh>) {
        if ($keep) {
            print {$out_fh} $_;
        }
        else {
            if (/>/) {
                print {$out_fh} ">$basename\n";
            }
            else {
                print {$out_fh} $_;
            }
        }
    }
    close $out_fh;
    close $in_fh;

    return $basename;
}

__END__
