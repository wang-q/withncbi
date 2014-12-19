#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use DBI;
use Text::CSV_XS;
use DateTime::Format::Natural;
use List::MoreUtils qw(any all uniq);
use Template;

use File::Copy::Recursive qw(fcopy);
use File::Spec;
use File::Find::Rule;
use File::Basename;
use Cwd qw(realpath);

use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use MyUtil qw(replace_home);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../config.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

my $working_dir = ".";
my $seq_dir;    #  will prep_fa from this dir ~/Scripts/alignDB/taxon
                #  or use seqs store in $working_dir
my @keep;       # don't touch anything inside fasta files

my $target_id;
my @query_ids;
my $outgroup_id;

my $name_str = "working";

# All taxons in this project (may also contain unused taxons)
my $taxon_file = "strains_taxon_info.csv";

# predefined phylogenetic tree
my $phylo_tree;

# Naming multiply alignment, the default value is $name_str
# This option is for more than one align combination.
my $multi_name;

my $clustalw;

my $aligndb  = replace_home( $Config->{run}{aligndb} );     # alignDB path
my $egaz     = replace_home( $Config->{run}{egaz} );        # egaz path
my $kent_bin = replace_home( $Config->{run}{kent_bin} );    # exes of Jim Kent
my $bat_dir  = $Config->{run}{bat};                         # Windows scripts

# Use name instead of taxon_id as identifier. These names should only contain
# alphanumeric value and match with sequence directory names.
# For strains not recorded in NCBI taxonomy, you should assign them fake ids.
# If this option set to be true, all $target_id, @query_ids is actually names.
my $use_name;

# Use Ensembl database as annotation source
my $ensembl;

# run in parallel mode
my $parallel = $Config->{run}{parallel};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'          => \$help,
    'man'             => \$man,
    'file=s'          => \$taxon_file,
    'w|working_dir=s' => \$working_dir,
    's|seq_dir=s'     => \$seq_dir,
    'keep=s'          => \@keep,
    't|target_id=s'   => \$target_id,
    'q|query_ids=s'   => \@query_ids,
    'o|r|outgroup=s'  => \$outgroup_id,
    'n|name_str=s'    => \$name_str,
    'un|use_name'     => \$use_name,
    'p|phylo_tree=s'  => \$phylo_tree,
    'm|multi_name=s'  => \$multi_name,
    'clustalw'        => \$clustalw,
    'e|ensembl=s'     => \$ensembl,
    'parallel=i'      => \$parallel,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

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

die "$taxon_file doesn't exist\n" unless -e $taxon_file;

# prepare working dir
{
    print "Working on $name_str\n";
    $working_dir = File::Spec->catdir( $working_dir, $name_str );
    $working_dir = File::Spec->rel2abs($working_dir);
    mkdir $working_dir unless -d $working_dir;
    $working_dir = realpath($working_dir);
    print " " x 4, "Working dir is $working_dir\n";
}

# move $outgroup_id to last
if ($outgroup_id) {
    my ($exist) = grep { $_ eq $outgroup_id } @query_ids;
    if ( !defined $exist ) {
        die "outgroup does not exist!\n";
    }

    @query_ids = grep { $_ ne $outgroup_id } @query_ids;
    push @query_ids, $outgroup_id;
}

my @data;
{
    if ( !$use_name ) {
        @data
            = map { taxon_info( $_, $working_dir ) } ( $target_id, @query_ids );
    }
    else {
        @data
            = map { taxon_info_name( $_, $working_dir ) }
            ( $target_id, @query_ids );
    }
}

# if seqs is not in working dir, copy them from seq_dir
my @target_seqs;    # for gff and rm-gff files
if ($seq_dir) {
    print "Get seqs from [$seq_dir]\n";

    for my $id ( $target_id, @query_ids ) {
        my ($keep) = grep { $_ eq $id } @keep;
        print " " x 4 . "Copy seq of [$id]\n";
        if ( defined $keep ) {
            print " " x 8 . "Don't change fasta header for [$id]\n";
        }

        my $original_dir = File::Spec->catdir( $seq_dir,     $id );
        my $cur_dir      = File::Spec->catdir( $working_dir, $id );
        mkdir $cur_dir unless -d $cur_dir;

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
            if ( $id eq $target_id ) {
                push @target_seqs, $basename;
            }
            if ( !$ensembl ) {
                my $gff_file
                    = File::Spec->catdir( $original_dir, "$basename.gff" );
                if ( -e $gff_file ) {
                    fcopy( $gff_file, $cur_dir );
                }

                my $rm_gff_file
                    = File::Spec->catdir( $original_dir, "$basename.rm.gff" );
                if ( -e $rm_gff_file ) {
                    fcopy( $rm_gff_file, $cur_dir );
                }
            }
        }
    }
}

{
    # write seq_pair.csv and left seq_pair_batch.pl to handle other things
    print "Create seq_pair.csv\n";
    my $seq_pair_file = File::Spec->catfile( $working_dir, "seq_pair.csv" );
    open my $fh, '>', $seq_pair_file;
    for my $query_id (@query_ids) {
        print {$fh} File::Spec->catdir( $working_dir, $target_id ), ",",
            File::Spec->catdir( $working_dir, $query_id ), "\n";
    }
    close $fh;
}

{
    print "Create id2name.csv\n";
    my $id2name_file = File::Spec->catfile( $working_dir, "id2name.csv" );
    open my $fh, '>', $id2name_file;
    if ( !$use_name ) {

        # if not set --use_name, use id as strain name
        for my $id ( $target_id, @query_ids ) {
            print {$fh} "$id,$id\n";
        }
    }
    else {
        for my $name ( $target_id, @query_ids ) {
            my ($id)
                = map { $_->{taxon} } grep { $_->{name} eq $name } @data;
            print {$fh} "$id,$name\n";
        }
    }
    close $fh;
}

# If there's no phylo tree, generate a fake one.
if ( !defined $phylo_tree ) {
    print "Create fake_tree.nwk\n";
    my $fake_tree_file = File::Spec->catfile( $working_dir, "fake_tree.nwk" );
    open my $fh, '>', $fake_tree_file;
    print {$fh} "(" x scalar(@query_ids) . "$target_id";
    for my $id (@query_ids) {
        print {$fh} ",$id)";
    }
    print {$fh} ";\n";
    close $fh;
}

{
    my $tt = Template->new;
    my $text;
    my $sh_name;

    # taxon.csv
    print "Create taxon.csv\n";
    $text = <<'EOF';
taxon_id,genus,species,sub_species,common_name
[% FOREACH item IN data -%]
[% item.taxon %],[% item.genus %],[% item.species %],[% item.subname %],[% item.name %]
[% END -%]
EOF
    $tt->process(
        \$text,
        { data => \@data, },
        File::Spec->catfile( $working_dir, "taxon.csv" )
    ) or die Template->error;

    #----------------------------#
    # all *.sh files
    #----------------------------#

    # real_chr.sh
    $sh_name = "1_real_chr.sh";
    print "Create $sh_name\n";
    $text = <<'EOF';
#!/bin/bash
cd [% working_dir %]

cat << DELIMITER > chrUn.csv
taxon_id,chr,length,name,assembly
[% FOREACH item IN data -%]
[% item.taxon %],chrUn,999999999,[% item.name %]
[% END -%]
DELIMITER

if [ -f real_chr.csv ]; then
    rm real_chr.csv;
fi;

[% FOREACH item IN data -%]
if [ ! -f [% item.dir%]/chr.sizes ]; then 
    [% kent_bin %]/faSize -detailed [% item.dir%]/*.fa > [% item.dir%]/chr.sizes;
fi;
perl -aln -F"\t" -e 'print qq{[% item.taxon %],$F[0],$F[1],[% item.name %]}' [% item.dir %]/chr.sizes >> real_chr.csv;
[% END -%]

cat chrUn.csv real_chr.csv > chr_length.csv
echo "chr_length.csv generated."

rm chrUn.csv
rm real_chr.csv

echo '# If you want, run the following cmds to merge csv files'
echo
echo perl [% aligndb %]/util/merge_csv.pl -t [% aligndb %]/data/taxon.csv -m [% working_dir %]/taxon.csv -f 0 -f 1
echo
echo perl [% aligndb %]/util/merge_csv.pl -t [% aligndb %]/data/chr_length.csv -m [% working_dir %]/chr_length.csv -f 0 -f 1
echo

EOF
    $tt->process(
        \$text,
        {   data        => \@data,
            working_dir => $working_dir,
            aligndb     => $aligndb,
            kent_bin    => $kent_bin,
        },
        File::Spec->catfile( $working_dir, $sh_name )
    ) or die Template->error;

    # file-rm.sh
    $sh_name = "2_file_rm.sh";
    print "Create $sh_name\n";
    $text = <<'EOF';
#!/bin/bash

cd [% working_dir %]

#----------------------------#
# repeatmasker on all fasta
#----------------------------#
for f in `find . -name "*.fa"` ; do
    rename 's/fa$/fasta/' $f ;
done

for f in `find . -name "*.fasta"` ; do
    RepeatMasker $f -xsmall --parallel [% parallel %] ;
done

for f in `find . -name "*.fasta.out"` ; do
    rmOutToGFF3.pl $f > `dirname $f`/`basename $f .fasta.out`.rm.gff;
done

for f in `find . -name "*.fasta"` ; do
    if [ -f $f.masked ];
    then
        rename 's/fasta.masked$/fa/' $f.masked;
        find . -type f -name "`basename $f`*" | xargs rm;
    fi;
done;

EOF

    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            target_id   => $target_id,
            query_ids   => \@query_ids,
        },
        File::Spec->catfile( $working_dir, $sh_name )
    ) or die Template->error;

    # pair_cmd.sh
    $sh_name = "3_pair_cmd.sh";
    print "Create $sh_name\n";
    $text = <<'EOF';
#!/bin/bash
# strain_bz.pl
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]

#----------------------------#
# seq_pair
#----------------------------#
perl [% aligndb %]/extra/seq_pair_batch.pl \
    --bin [% kent_bin %] \
[% IF use_name -%]
    --id [% working_dir %]/id2name.csv \
[% ELSE -%]
    --dir_as_taxon \
[% END -%]
    --parallel [% parallel %] \
    -f [% working_dir %]/seq_pair.csv \
    -taxon [% working_dir %]/taxon.csv \
    -chr [% working_dir %]/chr_length.csv \
    -lt 1000 -r 100-102

perl [% aligndb %]/extra/seq_pair_batch.pl \
    --bin [% kent_bin %] \
[% IF use_name -%]
    --id [% working_dir %]/id2name.csv \
[% ELSE -%]
    --dir_as_taxon \
[% END -%]
    --parallel [% parallel %] \
    -f [% working_dir %]/seq_pair.csv \
    -taxon [% working_dir %]/taxon.csv \
    -chr [% working_dir %]/chr_length.csv \
    -lt 1000 -r 1,2,5,21,40

EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            aligndb     => $aligndb,
            kent_bin    => $kent_bin,
            use_name    => $use_name,
            target_id   => $target_id,
            outgroup_id => $outgroup_id,
            query_ids   => \@query_ids,
        },
        File::Spec->catfile( $working_dir, $sh_name )
    ) or die Template->error;

    # rawphylo.sh
    if ( !defined $phylo_tree ) {
        $sh_name = "4_rawphylo.sh";
        print "Create $sh_name\n";
        $text = <<'EOF';
#!/bin/bash
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]

# join_dbs.pl
perl [% aligndb %]/extra/join_dbs.pl \
    --no_insert --block --trimmed_fasta --length 1000 \
    --goal_db [% multi_name %]_raw --target 0target \
[% IF outgroup_id -%]
    --outgroup [% query_ids.max %]query \
    --queries [% maxq = query_ids.max - 1 %][% FOREACH i IN [ 0 .. maxq ] %][% i %]query,[% END %] \
[% ELSE -%]
    --queries [% FOREACH i IN [ 0 .. query_ids.max ] %][% i %]query,[% END %] \
[% END -%]
    --dbs [% FOREACH id IN query_ids %][% target_id %]vs[% id %],[% END %]

#----------------------------#
# RAxML
#----------------------------#
# raw phylo guiding tree
if [ ! -d [% working_dir %]/rawphylo ]
then
    mkdir [% working_dir %]/rawphylo
fi

cd [% working_dir %]/rawphylo

rm [% working_dir %]/rawphylo/RAxML*

[% IF query_ids.size > 2 -%]
perl [% egaz%]/concat_fasta.pl \
    -i [% working_dir %]/[% multi_name %]_raw \
    -o [% working_dir %]/rawphylo/[% multi_name %].phy \
    -p

raxml -T 5 -f a -m GTRGAMMA -p $RANDOM -N 100 -x $RANDOM \
[% IF outgroup_id -%]
    -o [% outgroup_id %] \
[% END -%]
    -n [% multi_name %] -s [% working_dir %]/rawphylo/[% multi_name %].phy

cp [% working_dir %]/rawphylo/RAxML_best* [% working_dir %]/rawphylo/[% multi_name %].nwk

[% ELSIF query_ids.size == 2 -%]
echo "(([% target_id %],[% query_ids.0 %]),[% query_ids.1 %]);" > [% working_dir %]/rawphylo/[% multi_name %].nwk

[% ELSE -%]

echo "([% target_id %],[% query_ids.0 %]);" > [% working_dir %]/rawphylo/[% multi_name %].nwk

[% END -%]

EOF
        $tt->process(
            \$text,
            {   stopwatch   => $stopwatch,
                parallel    => $parallel,
                working_dir => $working_dir,
                aligndb     => $aligndb,
                egaz        => $egaz,
                target_id   => $target_id,
                outgroup_id => $outgroup_id,
                query_ids   => \@query_ids,
                multi_name  => $multi_name,
            },
            File::Spec->catfile( $working_dir, $sh_name )
        ) or die Template->error;
    }

    # multi_cmd.sh
    $sh_name = "5_multi_cmd.sh";
    print "Create $sh_name\n";
    $text = <<'EOF';
#!/bin/bash
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]
mkdir [% working_dir %]/[% multi_name %]

if [ -d [% working_dir %]/[% multi_name %]_fasta ]
then
    rm -fr [% working_dir %]/[% multi_name %]_fasta
fi

if [ -d [% working_dir %]/[% multi_name %]_mft ]
then
    rm -fr [% working_dir %]/[% multi_name %]_mft
fi

if [ -d [% working_dir %]/[% multi_name %]_clw ]
then
    rm -fr [% working_dir %]/[% multi_name %]_clw
fi

if [ -d [% working_dir %]/phylo ]
then
    rm -fr [% working_dir %]/phylo
fi

mkdir [% working_dir %]/phylo

#----------------------------#
# mz
#----------------------------#
[% IF phylo_tree -%]
perl [% egaz %]/mz.pl \
    [% FOREACH id IN query_ids -%]
    -d [% working_dir %]/[% target_id %]vs[% id %] \
    [% END -%]
    -bin [% kent_bin %] \
    --tree [% phylo_tree %] \
    --out [% working_dir %]/[% multi_name %] \
    -syn -p [% parallel %]
[% ELSE %]
if [ -f [% working_dir %]/rawphylo/[% multi_name %].nwk ]
then
    perl [% egaz %]/mz.pl \
        [% FOREACH id IN query_ids -%]
        -d [% working_dir %]/[% target_id %]vs[% id %] \
        [% END -%]
        -bin [% kent_bin %] \
        --tree [% working_dir %]/rawphylo/[% multi_name %].nwk \
        --out [% working_dir %]/[% multi_name %] \
        -syn -p [% parallel %]
else
    perl [% egaz %]/mz.pl \
        [% FOREACH id IN query_ids -%]
        -d [% working_dir %]/[% target_id %]vs[% id %] \
        [% END -%]
        -bin [% kent_bin %] \
        --tree [% working_dir %]/fake_tree.nwk \
        --out [% working_dir %]/[% multi_name %] \
        -syn -p [% parallel %]
fi
[% END -%]

find [% working_dir %]/[% multi_name %] -type f -name "*.maf" | parallel -j [% parallel %] gzip

#----------------------------#
# maf2fasta
#----------------------------#
perl [% egaz %]/maf2fasta.pl \
    -p [% parallel %] --block \
    -i [% working_dir %]/[% multi_name %] \
    -o [% working_dir %]/[% multi_name %]_fasta

#----------------------------#
# mafft
#----------------------------#
perl [% egaz %]/refine_fasta.pl \
    --msa mafft --block -p [% parallel %] \
[% IF outgroup_id -%]
    --outgroup \
[% END -%]
    -i [% working_dir %]/[% multi_name %]_fasta \
    -o [% working_dir %]/[% multi_name %]_mft

find [% working_dir %]/[% multi_name %]_mft -type f -name "*.fas" | parallel -j [% parallel %] gzip

[% IF clustalw -%]
#----------------------------#
# clustalw
#----------------------------#
perl [% egaz %]/refine_fasta.pl \
    --msa clustalw --block -p [% parallel %] \
[% IF outgroup_id -%]
    --outgroup \
[% END -%]
    -i [% working_dir %]/[% multi_name %]_fasta \
    -o [% working_dir %]/[% multi_name %]_clw

find [% working_dir %]/[% multi_name %]_clw -type f -name "*.fas" | parallel -j [% parallel %] gzip
[% END -%]

#----------------------------#
# RAxML
#----------------------------#
[% IF query_ids.size > 2 -%]
cd [% working_dir %]/phylo

perl [% egaz %]/concat_fasta.pl \
[% IF clustalw -%]
    -i [% working_dir %]/[% multi_name %]_clw \
[% ELSE -%]
    -i [% working_dir %]/[% multi_name %]_mft  \
[% END -%] 
    -o [% working_dir %]/phylo/[% multi_name %].phy \
    -p

rm [% working_dir %]/phylo/RAxML*

raxml -T 5 -f a -m GTRGAMMA -p $RANDOM -N 100 -x $RANDOM \
[% IF outgroup_id -%]
    -o [% outgroup_id %] \
[% END -%]
    -n [% multi_name %] -s [% working_dir %]/phylo/[% multi_name %].phy

cp [% working_dir %]/phylo/RAxML_best* [% working_dir %]/phylo/[% multi_name %].nwk
[% END -%]

EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            aligndb     => $aligndb,
            egaz        => $egaz,
            kent_bin    => $kent_bin,
            target_id   => $target_id,
            outgroup_id => $outgroup_id,
            query_ids   => \@query_ids,
            target_seqs => \@target_seqs,
            phylo_tree  => $phylo_tree,
            multi_name  => $multi_name,
            clustalw    => $clustalw,
        },
        File::Spec->catfile( $working_dir, $sh_name )
    ) or die Template->error;

    # multi_db_only.sh
    $sh_name = "6_multi_db_only.sh";
    print "Create $sh_name\n";
    $text = <<'EOF';
#!/bin/bash
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]

#----------------------------#
# multi_way_batch
#----------------------------#
perl [% aligndb %]/extra/multi_way_batch.pl \
    -d [% multi_name %] \
[% IF clustalw -%]
    -da [% working_dir %]/[% multi_name %]_clw \
[% ELSE -%]
    -da [% working_dir %]/[% multi_name %]_mft \
[% END -%]
    --gff_file [% FOREACH seq IN target_seqs %][% working_dir %]/[% target_id %]/[% seq %].gff,[% END %] \
    --rm_gff_file [% FOREACH seq IN target_seqs %][% working_dir %]/[% target_id %]/[% seq %].rm.gff,[% END %] \
    --block \
    --id [% working_dir %]/id2name.csv \
[% IF outgroup_id -%]
    --outgroup \
[% END -%]
    -taxon [% working_dir %]/taxon.csv \
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
            target_id   => $target_id,
            outgroup_id => $outgroup_id,
            query_ids   => \@query_ids,
            target_seqs => \@target_seqs,
            multi_name  => $multi_name,
            clustalw    => $clustalw,
        },
        File::Spec->catfile( $working_dir, $sh_name )
    ) or die Template->error;

    # chart.bat
    print "Create chart.bat\n";
    $text = <<'EOF';
REM strain_bz.pl
REM perl [% stopwatch.cmd_line %]

REM basicstat
perl [% bat_dir %]/fig_table/collect_common_basic.pl -d .

REM common chart
if exist [% multi_name %].common.xlsx perl [% bat_dir %]/alignDB/stat/common_chart_factory.pl -i [% multi_name %].common.xlsx

REM multi chart
if exist [% multi_name %].multi.xlsx  perl [% bat_dir %]/alignDB/stat/multi_chart_factory.pl -i [% multi_name %].multi.xlsx

REM gc chart
if exist [% multi_name %].gc.xlsx     perl [% bat_dir %]/alignDB/stat/gc_chart_factory.pl --add_trend 1 -i [% multi_name %].gc.xlsx

EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            bat_dir     => $bat_dir,
            multi_name  => $multi_name,
        },
        File::Spec->catfile( $working_dir, "chart.bat" )
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
sub taxon_info {
    my $taxon_id = shift;
    my $dir      = shift;

    my $dbh = DBI->connect("DBI:CSV:");

    $dbh->{csv_tables}->{t0} = {
        eol       => "\n",
        sep_char  => ",",
        file      => $taxon_file,
        col_names => [
            map { ( $_, $_ . "_id" ) } qw{strain species genus family order}
        ],
    };

    my $query
        = qq{ SELECT strain_id, strain, genus, species FROM t0 WHERE strain_id = ? };
    my $sth = $dbh->prepare($query);
    $sth->execute($taxon_id);
    my ( $taxonomy_id, $organism_name, $genus, $species )
        = $sth->fetchrow_array;
    $species =~ s/^$genus\s+//;
    my $sub_name = $organism_name;
    $sub_name =~ s/^$genus\s+//;
    $sub_name =~ s/^$species\s*//;
    $organism_name =~ s/\W/_/g;
    $organism_name =~ s/_+/_/g;

    return {
        taxon   => $taxonomy_id,
        name    => $organism_name,
        genus   => $genus,
        species => $species,
        subname => $sub_name,
        dir     => File::Spec->catdir( $working_dir, $taxon_id ),
    };
}

sub taxon_info_name {
    my $name = shift;
    my $dir  = shift;

    my $dbh = DBI->connect("DBI:CSV:");

    $dbh->{csv_tables}->{t0} = {
        eol       => "\n",
        sep_char  => ",",
        file      => $taxon_file,
        col_names => [
            map { ( $_, $_ . "_id" ) } qw{strain species genus family order}
        ],
    };

    my $query
        = qq{ SELECT strain_id, strain, genus, species FROM t0 WHERE strain = ? };
    my $sth = $dbh->prepare($query);
    $sth->execute($name);
    my ( $taxonomy_id, $organism_name, $genus, $species )
        = $sth->fetchrow_array;
    $species =~ s/^$genus\s+//;
    my $sub_name = $organism_name;
    $sub_name =~ s/^$genus[\s_]+//;
    $sub_name =~ s/^$species[\s_]*//;

    return {
        taxon   => $taxonomy_id,
        name    => $name,
        genus   => $genus,
        species => $species,
        subname => $sub_name,
        dir     => File::Spec->catdir( $working_dir, $name ),
    };
}

sub prep_fa {
    my $infile = shift;
    my $dir    = shift;
    my $keep   = shift;

    my $basename = basename( $infile, '.fna', '.fa', '.fas', '.fasta' );

    my $outfile = File::Spec->catfile( $dir, "$basename.fa" );
    open my $in_fh,  '<', $infile;
    open my $out_fh, '>', $outfile;
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
perl strain_bz.pl 
