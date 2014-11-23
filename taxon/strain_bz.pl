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

use File::Copy::Recursive qw(fcopy);
use File::Spec;
use File::Find::Rule;
use File::Basename;
use Cwd qw(realpath);

use Template;

use FindBin;

use AlignDB::Stopwatch;

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

my $target_id;
my $outgroup_id;
my @query_ids;

my $clustalw;

my $name_str = "working";

my $filename = "strains_taxon_info.csv";

my $aligndb = $Config->{run}{aligndb};    # current alignDB path
my $bat_dir = $Config->{run}{bat};        # Windows alignDB path

# run in parallel mode
my $parallel = $Config->{run}{parallel};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'          => \$help,
    'man'             => \$man,
    'file=s'          => \$filename,
    'w|working_dir=s' => \$working_dir,
    's|seq_dir=s'     => \$seq_dir,
    'b|bat_dir=s'     => \$bat_dir,
    't|target_id=i'   => \$target_id,
    'o|r|outgroup=i'  => \$outgroup_id,
    'q|query_ids=i'   => \@query_ids,
    'n|name_str=s'    => \$name_str,
    'clustalw'        => \$clustalw,
    'parallel=i'      => \$parallel,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Writing strains summary...");

die "$filename doesn't exist\n" unless -e $filename;

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

# if seqs is not in working dir, copy them from seq_dir
my @target_accs;
if ($seq_dir) {
    print "Get seqs from [$seq_dir]\n";

    for my $id ( $target_id, @query_ids ) {
        print " " x 4, "Copy seq of $id\n";

        my $original_dir = File::Spec->catdir( $seq_dir,     $id );
        my $cur_dir      = File::Spec->catdir( $working_dir, $id );
        mkdir $cur_dir unless -d $cur_dir;

        my @fa_files
            = File::Find::Rule->file->name( '*.fna', '*.fa', '*.fas',
            '*.fasta' )->in($original_dir);

        for my $fa_file (@fa_files) {
            my $basename = prep_fa( $fa_file, $cur_dir );
            if ( $id eq $target_id ) {
                push @target_accs, $basename;
            }
            my $gff_file = File::Spec->catdir( $original_dir, "$basename.gff" );
            if ( -e $gff_file ) {
                fcopy( $gff_file, $cur_dir );
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
    for my $id ( $target_id, @query_ids ) {
        print {$fh} "$id,$id\n";
    }
    close $fh;
}

{
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
    my @data
        = map { taxon_info( $_, $working_dir ) } ( $target_id, @query_ids );

    # taxon.csv
    print "Create taxon.csv\n";
    $text = <<'EOF';
[% FOREACH item IN data -%]
[% item.taxon %],[% item.genus %],[% item.species %],[% item.subname %],[% item.name %],
[% END -%]
EOF
    $tt->process(
        \$text,
        { data => \@data, },
        File::Spec->catfile( $working_dir, "taxon.csv" )
    ) or die Template->error;

    # chr_length.csv
    print "Create chr_length.csv\n";
    $text = <<'EOF';
[% FOREACH item IN data -%]
[% item.taxon %],chrUn,999999999,[% item.name %]
[% END -%]
EOF
    $tt->process(
        \$text,
        { data => \@data, },
        File::Spec->catfile( $working_dir, "chr_length_chrUn.csv" )
    ) or die Template->error;

    # file-rm.sh
    print "Create file-rm.sh\n";
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
        File::Spec->catfile( $working_dir, "file-rm.sh" )
    ) or die Template->error;

    # real_chr.sh
    print "Create real_chr.sh\n";
    $text = <<'EOF';
#!/bin/bash
cd [% working_dir %]

if [ -f real_chr.csv ]; then
    rm real_chr.csv;
fi;

[% FOREACH item IN data -%]
faSize -detailed [% item.dir%]/*.fa > [% item.dir%]/chr.sizes
perl -aln -F"\t" -e 'print qq{[% item.taxon %],$F[0],$F[1],[% item.name %]}' [% item.dir %]/chr.sizes >> real_chr.csv
[% END -%]

cat chr_length_chrUn.csv real_chr.csv > chr_length.csv
rm real_chr.csv

echo '# Run the following cmds to merge csv files'
echo
echo perl [% aligndb %]/util/merge_csv.pl -t [% aligndb %]/init/taxon.csv -m [% working_dir %]/taxon.csv -f 0 -f 1
echo
echo perl [% aligndb %]/util/merge_csv.pl -t [% aligndb %]/init/chr_length.csv -m [% working_dir %]/chr_length.csv -f 0 -f 1
echo

EOF
    $tt->process(
        \$text,
        {   data        => \@data,
            working_dir => $working_dir,
            aligndb     => $aligndb,
        },
        File::Spec->catfile( $working_dir, "real_chr.sh" )
    ) or die Template->error;

    # pair_cmd.sh
    print "Create pair_cmd.sh\n";
    $text = <<'EOF';
#!/bin/bash
# strain_bz.pl
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]

#----------------------------#
# seq_pair
#----------------------------#
perl [% aligndb %]/extra/seq_pair_batch.pl \
    --dir_as_taxon \
    --parallel [% parallel %] \
    -f [% working_dir %]/seq_pair.csv \
    -lt 1000 -r 100-102

perl [% aligndb %]/extra/seq_pair_batch.pl \
    --dir_as_taxon \
    --parallel [% parallel %] \
    -f [% working_dir %]/seq_pair.csv \
    -lt 1000 -r 1,2,5,21,40

EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            aligndb     => $aligndb,
            name_str    => $name_str,
            target_id   => $target_id,
            outgroup_id => $outgroup_id,
            query_ids   => \@query_ids,
        },
        File::Spec->catfile( $working_dir, "pair_cmd.sh" )
    ) or die Template->error;

    # rawphylo.sh
    print "Create rawphylo.sh\n";
    $text = <<'EOF';
#!/bin/bash
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]

# join_dbs.pl
perl [% aligndb %]/extra/join_dbs.pl \
    --no_insert --block --trimmed_fasta --length 1000 \
    --goal_db [% name_str %]_raw --target 0target \
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
perl [% aligndb %]/../blastz/concat_fasta.pl \
    -i [% working_dir %]/[% name_str %]_raw \
    -o [% working_dir %]/rawphylo/[% name_str %].phy \
    -p

raxml -T 5 -f a -m GTRGAMMA -p $RANDOM -N 100 -x $RANDOM \
[% IF outgroup_id -%]
    -o [% outgroup_id %] \
[% END -%]
    -n [% name_str %] -s [% working_dir %]/rawphylo/[% name_str %].phy

cp [% working_dir %]/rawphylo/RAxML_best* [% working_dir %]/rawphylo/[% name_str %].nwk

[% ELSIF query_ids.size == 2 -%]
echo "(([% target_id %],[% query_ids.0 %]),[% query_ids.1 %]);" > [% working_dir %]/rawphylo/[% name_str %].nwk

[% ELSE -%]

echo "([% target_id %],[% query_ids.0 %]);" > [% working_dir %]/rawphylo/[% name_str %].nwk

[% END -%]

EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            aligndb     => $aligndb,
            name_str    => $name_str,
            target_id   => $target_id,
            outgroup_id => $outgroup_id,
            query_ids   => \@query_ids,
        },
        File::Spec->catfile( $working_dir, "rawphylo.sh" )
    ) or die Template->error;

    # cmd.bat
    print "Create cmd.bat\n";
    $text = <<'EOF';
REM strain_bz.pl
REM perl [% stopwatch.cmd_line %]

REM basicstat
perl [% bat_dir %]/fig/collect_common_basic.pl -d .

REM common chart
if exist [% name_str %].common.xlsx perl [% bat_dir %]/stat/common_chart_factory.pl -i [% name_str %].common.xlsx

REM multi chart
if exist [% name_str %].multi.xlsx  perl [% bat_dir %]/stat/multi_chart_factory.pl -i [% name_str %].multi.xlsx

REM gc chart
if exist [% name_str %].gc.xlsx     perl [% bat_dir %]/stat/gc_chart_factory.pl --add_trend 1 -i [% name_str %].gc.xlsx

EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            bat_dir     => $bat_dir,
            name_str    => $name_str,
        },
        File::Spec->catfile( $working_dir, "chart.bat" )
    ) or die Template->error;

    # multi_cmd.sh
    print "Create multi_cmd.sh\n";
    $text = <<'EOF';
#!/bin/bash
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]
mkdir [% working_dir %]/[% name_str %]

if [ -d [% working_dir %]/[% name_str %]_fasta ]
then
    rm -fr [% working_dir %]/[% name_str %]_fasta
fi

if [ -d [% working_dir %]/[% name_str %]_mft ]
then
    rm -fr [% working_dir %]/[% name_str %]_mft
fi

if [ -d [% working_dir %]/[% name_str %]_clw ]
then
    rm -fr [% working_dir %]/[% name_str %]_clw
fi

if [ -d [% working_dir %]/phylo ]
then
    rm -fr [% working_dir %]/phylo
fi

mkdir [% working_dir %]/phylo

#----------------------------#
# mz
#----------------------------#
if [ -f [% working_dir %]/rawphylo/[% name_str %].nwk ]
then
    perl [% aligndb %]/../blastz/mz.pl \
        [% FOREACH id IN query_ids -%]
        -d [% working_dir %]/[% target_id %]vs[% id %] \
        [% END -%]
        --tree [% working_dir %]/rawphylo/[% name_str %].nwk \
        --out [% working_dir %]/[% name_str %] \
        -syn -p [% parallel %]
else
    perl [% aligndb %]/../blastz/mz.pl \
        [% FOREACH id IN query_ids -%]
        -d [% working_dir %]/[% target_id %]vs[% id %] \
        [% END -%]
        --tree [% working_dir %]/fake_tree.nwk \
        --out [% working_dir %]/[% name_str %] \
        -syn -p [% parallel %]
fi

find [% working_dir %]/[% name_str %] -type f -name "*.maf" | parallel -j [% parallel %] gzip

#----------------------------#
# maf2fasta
#----------------------------#
perl [% aligndb %]/../blastz/maf2fasta.pl \
    -p [% parallel %] --block \
    -i [% working_dir %]/[% name_str %] \
    -o [% working_dir %]/[% name_str %]_fasta

#----------------------------#
# mafft
#----------------------------#
perl [% aligndb %]/../blastz/refine_fasta.pl \
    --msa mafft --block -p [% parallel %] \
[% IF outgroup_id -%]
    --outgroup \
[% END -%]
    -i [% working_dir %]/[% name_str %]_fasta \
    -o [% working_dir %]/[% name_str %]_mft

find [% working_dir %]/[% name_str %]_mft -type f -name "*.fas" | parallel -j [% parallel %] gzip

[% IF clustalw -%]
#----------------------------#
# clustalw
#----------------------------#
perl [% aligndb %]/../blastz/refine_fasta.pl \
    --msa clustalw --block -p [% parallel %] \
[% IF outgroup_id -%]
    --outgroup \
[% END -%]
    -i [% working_dir %]/[% name_str %]_fasta \
    -o [% working_dir %]/[% name_str %]_clw

find [% working_dir %]/[% name_str %]_clw -type f -name "*.fas" | parallel -j [% parallel %] gzip
[% END -%]

#----------------------------#
# multi_way_batch
#----------------------------#
perl [% aligndb %]/extra/multi_way_batch.pl \
    -d [% name_str %] \
[% IF clustalw -%]
    -da [% working_dir %]/[% name_str %]_clw \
[% ELSE -%]
    -da [% working_dir %]/[% name_str %]_mft \
[% END -%]
    --gff_file [% FOREACH acc IN target_accs %][% working_dir %]/[% target_id %]/[% acc %].gff,[% END %] \
    --rm_gff_file [% FOREACH acc IN target_accs %][% working_dir %]/[% target_id %]/[% acc %].rm.gff,[% END %] \
    --block --id [% working_dir %]/id2name.csv \
[% IF outgroup_id -%]
    --outgroup \
[% END -%]
    -lt 1000 --parallel [% parallel %] --batch 5 \
    --run 1,2,5,10,21,30-32,40-42,44

#----------------------------#
# RAxML
#----------------------------#
[% IF query_ids.size > 2 -%]
cd [% working_dir %]/phylo

perl [% aligndb %]/../blastz/concat_fasta.pl \
    -i [% working_dir %]/[% name_str %]_mft  \
    -o [% working_dir %]/phylo/[% name_str %].phy \
    -p

rm [% working_dir %]/phylo/RAxML*

raxml -T 5 -f a -m GTRGAMMA -p $RANDOM -N 100 -x $RANDOM \
[% IF outgroup_id -%]
    -o [% outgroup_id %] \
[% END -%]
    -n [% name_str %] -s [% working_dir %]/phylo/[% name_str %].phy

cp [% working_dir %]/phylo/RAxML_best* [% working_dir %]/phylo/[% name_str %].nwk
[% END -%]

EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            aligndb     => $aligndb,
            name_str    => $name_str,
            target_id   => $target_id,
            outgroup_id => $outgroup_id,
            query_ids   => \@query_ids,
            target_accs => \@target_accs,
            clustalw    => $clustalw,
        },
        File::Spec->catfile( $working_dir, "multi_cmd.sh" )
    ) or die Template->error;

    # multi_db_only.sh
    print "Create multi_db_only.sh\n";
    $text = <<'EOF';
#!/bin/bash
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]

#----------------------------#
# multi_way_batch
#----------------------------#
perl [% aligndb %]/extra/multi_way_batch.pl \
    -d [% name_str %] \
[% IF clustalw -%]
    -da [% working_dir %]/[% name_str %]_clw \
[% ELSE -%]
    -da [% working_dir %]/[% name_str %]_mft \
[% END -%]
    --gff_file [% FOREACH acc IN target_accs %][% working_dir %]/[% target_id %]/[% acc %].gff,[% END %] \
    --rm_gff_file [% FOREACH acc IN target_accs %][% working_dir %]/[% target_id %]/[% acc %].rm.gff,[% END %] \
    --block --id [% working_dir %]/id2name.csv \
[% IF outgroup_id -%]
    --outgroup \
[% END -%]
    -lt 1000 --parallel [% parallel %] --batch 5 \
    --run 1,2,5,10,21,30-32,40-42,44

    
EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            aligndb     => $aligndb,
            name_str    => $name_str,
            target_id   => $target_id,
            outgroup_id => $outgroup_id,
            query_ids   => \@query_ids,
            target_accs => \@target_accs,
            clustalw    => $clustalw,
        },
        File::Spec->catfile( $working_dir, "multi_db_only.sh" )
    ) or die Template->error;
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
        file      => $filename,
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

sub prep_fa {
    my $infile = shift;
    my $dir    = shift;

    my $basename = basename( $infile, '.fna', '.fa', '.fas', '.fasta' );

    my $outfile = File::Spec->catfile( $dir, "$basename.fa" );
    open my $in_fh,  '<', $infile;
    open my $out_fh, '>', $outfile;
    while (<$in_fh>) {
        if (/>/) {
            print {$out_fh} ">$basename\n";
        }
        else {
            print {$out_fh} $_;
        }
    }
    close $out_fh;
    close $in_fh;

    return $basename;
}

__END__

perl bac_strains.pl 
