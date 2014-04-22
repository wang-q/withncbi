#!/usr/bin/perl
use strict;
use warnings;

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
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

my $working_dir = ".";
my $seq_dir;    #  will prep_fa from this dir ~/Scripts/alignDB/taxon
                #  or use seqs store in $working_dir

my $bat_dir = "d:/wq/Scripts/alignDB";    # Windows alignDB path

my $target_id;
my @query_ids;

my $clustalw;

my $length = 1000;

my $name_str = "working";

my $filename = "strains_taxon_info.csv";

# run in parallel mode
my $parallel = $Config->{generate}{parallel};

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
    'q|query_ids=i'   => \@query_ids,
    'n|name_str=s'    => \$name_str,
    'clustalw'        => \$clustalw,
    'length=i'        => \$length,
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

# if seqs is not in working dir, copy them from seq_dir
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
            my $gff_file = File::Spec->catdir( $original_dir, "$basename.gff" );
            if ( -e $gff_file ) {
                fcopy( $gff_file, $cur_dir );
            }
        }
    }
}

my $acc_of = {};
{
    for my $id ( $target_id, @query_ids ) {
        my $cur_dir = File::Spec->catdir( $working_dir, $id );
        my @fa_files
            = File::Find::Rule->file->name( '*.fna', '*.fa', '*.fas',
            '*.fasta' )->in($cur_dir);

        $acc_of->{$id} = [];
        for my $fa_file (@fa_files) {
            my $basename
                = basename( $fa_file, '.fna', '.fa', '.fas', '.fasta' );
            push @{ $acc_of->{$id} }, $basename;
        }
    }
}

{    # use id as species name
    print "Create id2name.csv\n";
    my $id2name_file = File::Spec->catfile( $working_dir, "id2name.csv" );
    open my $fh, '>', $id2name_file;
    for my $id ( $target_id, @query_ids ) {
        print {$fh} "$id,$id\n";
    }
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
echo perl [% findbin %]/../util/merge_csv.pl -t [% findbin %]/../init/taxon.csv -m [% working_dir %]/taxon.csv -f 0 -f 1
echo
echo perl [% findbin %]/../util/merge_csv.pl -t [% findbin %]/../init/chr_length.csv -m [% working_dir %]/chr_length.csv -f 0 -f 1
echo

EOF
    $tt->process(
        \$text,
        {   data        => \@data,
            working_dir => $working_dir,
            findbin     => $FindBin::Bin,
        },
        File::Spec->catfile( $working_dir, "real_chr.sh" )
    ) or die Template->error;

    # self_cmd.sh
    print "Create self_cmd.sh\n";
    $text = <<'EOF';
#!/bin/bash
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]

[% FOREACH item IN data -%]
#----------------------------------------------------------#
# [% item.taxon %] [% item.name %]
#----------------------------------------------------------#
if [ -d [% working_dir %]/[% item.taxon %]vsselfalign ]
then
    rm -fr [% working_dir %]/[% item.taxon %]vsselfalign
fi

#----------------------------#
# self bz
#----------------------------#
perl [% findbin %]/../../blastz/bz.pl \
    --is_self \
    -s set01 -C 0 --noaxt -pb lastz --lastz \
    -dt [% working_dir %]/[% item.taxon %] \
    -dq [% working_dir %]/[% item.taxon %] \
    -dl [% working_dir %]/[% item.taxon %]vsselfalign \
    --parallel [% parallel %]

#----------------------------#
# lpcna
#----------------------------#
perl [% findbin %]/../../blastz/lpcna.pl \
    -dt [% working_dir %]/[% item.taxon %] \
    -dq [% working_dir %]/[% item.taxon %] \
    -dl [% working_dir %]/[% item.taxon %]vsselfalign \
    --parallel [% parallel %]

[% END -%]

EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            findbin     => $FindBin::Bin,
            name_str    => $name_str,
            target_id   => $target_id,
            query_ids   => \@query_ids,
            data        => \@data,
        },
        File::Spec->catfile( $working_dir, "self_cmd.sh" )
    ) or die Template->error;

    # proc_cmd.sh
    print "Create proc_cmd.sh\n";
    $text = <<'EOF';
#!/bin/bash
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]

[% FOREACH item IN data -%]
#----------------------------------------------------------#
# [% item.taxon %] [% item.name %]
#----------------------------------------------------------#
if [ -d [% working_dir %]/[% item.taxon %]_proc ]
then
    find [% working_dir %]/[% item.taxon %]_proc -type f -not -name "circos.conf" | xargs rm
else
    mkdir [% working_dir %]/[% item.taxon %]_proc
fi

if [ ! -d [% working_dir %]/[% item.taxon %]_result ]
then
    mkdir [% working_dir %]/[% item.taxon %]_result
fi

cd [% working_dir %]/[% item.taxon %]_proc

#----------------------------#
# quick and dirty coverage
#----------------------------#
perl [% findbin %]/../slice/gather_info_axt.pl -l [% length %] -d [% working_dir %]/[% item.taxon %]vsselfalign --nomatch
perl [% findbin %]/../slice/compare_runlist.pl -op intersect -f1 [% item.taxon %]vsselfalign.runlist.0.yml -f2 [% item.taxon %]vsselfalign.runlist.1.yml -o [% item.taxon %]vsselfalign.runlist.i.yml
perl [% findbin %]/../slice/compare_runlist.pl -op xor -f1 [% item.taxon %]vsselfalign.runlist.0.yml -f2 [% item.taxon %]vsselfalign.runlist.1.yml -o [% item.taxon %]vsselfalign.runlist.x.yml
perl [% findbin %]/../slice/compare_runlist.pl -op union -f1 [% item.taxon %]vsselfalign.runlist.0.yml -f2 [% item.taxon %]vsselfalign.runlist.1.yml -o [% item.taxon %]vsselfalign.runlist.u.yml

for op in 0 1 i x u
do
    perl [% findbin %]/../slice/stat_runlist.pl --size [% working_dir %]/[% item.taxon %]/chr.sizes -f [% item.taxon %]vsselfalign.runlist.$op.yml;
done

echo "group,name,length,size,coverage" > [% working_dir %]/[% item.taxon %]_result/[% item.taxon %].[% length %].csv
for op in 0 1 i x u
do
    OP=$op perl -nl -e '/^name/ and next; print qq{$ENV{OP},$_};' [% item.taxon %]vsselfalign.runlist.$op.yml.csv;
done >> [% working_dir %]/[% item.taxon %]_result/[% item.taxon %].[% length %].csv

for op in 0 1 i x u
do
    rm [% item.taxon %]vsselfalign.runlist.$op.yml;
    rm [% item.taxon %]vsselfalign.runlist.$op.yml.csv;
done

#----------------------------#
# genome blast
#----------------------------#
find [% working_dir %]/[% item.taxon %] -type f -name "*.fa" \
    | sort | xargs cat \
    | perl -nl -e '/^>/ or uc; print' \
    > [% item.taxon %].genome.fasta

echo "* build genome blast db [% item.taxon %].genome.fasta"
formatdb -p F -o T -i [% item.taxon %].genome.fasta

perl [% findbin %]/../slice/gather_seq_axt.pl -l [% length %] -d [% working_dir %]/[% item.taxon %]vsselfalign

echo "* blast [% item.taxon %]vsselfalign.axt.fasta"
blastall -p blastn -F "m D" -m 0 -b 10 -v 10 -e 1e-3 -a 8 -i [% item.taxon %]vsselfalign.axt.fasta -d [% item.taxon %].genome.fasta -o [% item.taxon %]vsselfalign.axt.blast

#----------------------------#
# paralog sequences
#----------------------------#
# Omit genome locations in .axt
# There are errors, especially for queries
perl [% findbin %]/../slice/blastn_genome_location.pl -f [% item.taxon %]vsselfalign.axt.blast -m 0 -i 90 -c 0.95

#----------------------------#
# paralog blast
#----------------------------#
echo "* build paralog blast db [% item.taxon %]vsselfalign.gl.fasta"
formatdb -p F -o T -i [% item.taxon %]vsselfalign.gl.fasta

echo "* blast [% item.taxon %]vsselfalign.gl.fasta"
blastall -p blastn -F "m D" -m 0 -b 10 -v 10 -e 1e-3 -a 8 -i [% item.taxon %]vsselfalign.gl.fasta -d [% item.taxon %]vsselfalign.gl.fasta -o [% item.taxon %]vsselfalign.gl.blast

#----------------------------#
# merge
#----------------------------#
perl [% findbin %]/../slice/blastn_paralog.pl -f [% item.taxon %]vsselfalign.gl.blast -m 0 -i 90 -c 0.9

perl [% findbin %]/../slice/merge_node.pl    -v -f [% item.taxon %]vsselfalign.blast.tsv -o [% item.taxon %]vsselfalign.merge.yml -c 0.9
perl [% findbin %]/../slice/paralog_graph.pl -v -f [% item.taxon %]vsselfalign.blast.tsv -m [% item.taxon %]vsselfalign.merge.yml --nonself -o [% item.taxon %]vsselfalign.merge.graph.yml
perl [% findbin %]/../slice/cc.pl               -f [% item.taxon %]vsselfalign.merge.graph.yml
perl [% findbin %]/../slice/proc_cc_chop.pl     -f [% item.taxon %]vsselfalign.cc.yml --size [% working_dir %]/[% item.taxon %]/chr.sizes
perl [% findbin %]/../slice/proc_cc_stat.pl     -f [% item.taxon %]vsselfalign.cc.yml --size [% working_dir %]/[% item.taxon %]/chr.sizes

#----------------------------#
# result
#----------------------------#
cp [% item.taxon %]vsselfalign.cc.yml [% working_dir %]/[% item.taxon %]_result
mv [% item.taxon %]vsselfalign.cc.csv [% working_dir %]/[% item.taxon %]_result

#----------------------------#
# clean
#----------------------------#
find [% working_dir %]/[% item.taxon %]_proc -type f -name "*genome.fasta*" | xargs rm
find [% working_dir %]/[% item.taxon %]_proc -type f -name "*gl.fasta*" | xargs rm
find [% working_dir %]/[% item.taxon %]_proc -type f -name "*.blast" | xargs rm

[% END -%]
EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            findbin     => $FindBin::Bin,
            name_str    => $name_str,
            target_id   => $target_id,
            query_ids   => \@query_ids,
            data        => \@data,
            length      => $length,
        },
        File::Spec->catfile( $working_dir, "proc_cmd.sh" )
    ) or die Template->error;

    # circos.conf
    $text = <<'EOF';
<image>
dir*   = [% working_dir %]/[% taxon_id %]_result
file*  = [% taxon_id %].circos.png
background*     = white

# radius of inscribed circle in image
radius         = 1500p
background     = white

# by default angle=0 is at 3 o'clock position
angle_offset   = -90

24bit             = yes
auto_alpha_colors = yes
auto_alpha_steps  = 5
</image>

karyotype = karyotype.[% taxon_id %].txt
chromosomes_units = 1000

chromosomes_display_default = yes

<links>

<link>
file          = [% taxon_id %]vsselfalign.cc.link4.txt
radius        = 0.88r
bezier_radius = 0.2r
color         = purple
thickness     = 2
ribbon        = yes
stroke_color  = purple
stroke_thickness = 2
</link>

<link>
file          = [% taxon_id %]vsselfalign.cc.link3.txt
radius        = 0.88r
bezier_radius = 0.1r
color         = dgreen
thickness     = 3
ribbon        = yes
stroke_color  = dgreen
stroke_thickness = 2
</link>

<link>
file          = [% taxon_id %]vsselfalign.cc.link2.txt
radius        = 0.88r
bezier_radius = 0r
color         = dorange
thickness     = 3
ribbon        = yes
stroke_color  = dorange
stroke_thickness = 2
</link>

</links>

<highlights>

<highlight>
file = highlight.features.[% taxon_id %].txt
r0 = 0.95r
r1 = 0.98r
</highlight>

<highlight>
file = highlight.repeats.[% taxon_id %].txt
r0 = 0.93r
r1 = 0.98r
</highlight>

<highlight>
file = [% taxon_id %]vsselfalign.cc.linkN.txt
r0 = 0.89r
r1 = 0.92r
stroke_thickness = 2
stroke_color = grey
</highlight>

</highlights>

<ideogram>

<spacing>
    default = 0u
    break   = 0u
</spacing>

# thickness (px) of chromosome ideogram
thickness        = 20p
stroke_thickness = 2p

# ideogram border color
stroke_color     = dgrey
fill             = yes

# the default chromosome color is set here and any value
# defined in the karyotype file overrides it
fill_color       = black

# fractional radius position of chromosome ideogram within image
radius         = 0.85r
show_label     = no
label_font     = condensedbold
label_radius   = dims(ideogram,radius) + 0.05r
label_size     = 36

label_parallel   = yes

show_bands            = yes
fill_bands            = yes
band_stroke_thickness = 2
band_stroke_color     = white
band_transparency     = 1

</ideogram>

show_ticks       = yes
show_tick_labels = yes

show_grid        = no
grid_start       = dims(ideogram,radius_inner)-0.5r
grid_end         = dims(ideogram,radius_inner)

<ticks>
    skip_first_label           = yes
    skip_last_label            = no
    radius                     = dims(ideogram,radius_outer)
    tick_separation            = 2p
    min_label_distance_to_edge = 0p
    label_separation           = 5p
    label_offset               = 2p
    label_size                 = 8p
    multiplier                 = 0.001
    color                      = black

<tick>
    spacing        = 10u
    size           = 8p
    thickness      = 2p
    color          = black
    show_label     = no
    grid           = yes
    grid_color     = grey
    grid_thickness = 1p
</tick>

<tick>
    spacing        = 100u
    size           = 8p
    thickness      = 2p
    color          = black
    show_label     = yes
    suffix = " kb"
    label_size     = 36p
    label_offset   = 5p
    format         = %s
    grid           = yes
    grid_color     = dgrey
    grid_thickness = 1p
</tick>

</ticks>

<<include etc/colors_fonts_patterns.conf>>

<<include etc/housekeeping.conf>>

EOF
    for my $taxon_id ( $target_id, @query_ids ) {
        print "    Create circos.conf for $taxon_id\n";
        $tt->process(
            \$text,
            {   stopwatch   => $stopwatch,
                parallel    => $parallel,
                working_dir => $working_dir,
                findbin     => $FindBin::Bin,
                name_str    => $name_str,
                taxon_id    => $taxon_id,
            },
            File::Spec->catfile(
                $working_dir, "${taxon_id}_proc", "circos.conf"
            )
        ) or die Template->error;
    }

    # circos_cmd.sh
    print "Create circos_cmd.sh\n";
    $text = <<'EOF';
#!/bin/bash
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]

[% FOREACH item IN data -%]
#----------------------------------------------------------#
# [% item.taxon %] [% item.name %]
#----------------------------------------------------------#

#----------------------------#
# karyotype
#----------------------------#
cd [% working_dir %]/[% item.taxon %]_proc

# generate karyotype file
perl -anl -e '$i++; print qq{chr - $F[0] $F[0] 0 $F[1] chr$i}' [% working_dir %]/[% item.taxon %]/chr.sizes > karyotype.[% item.taxon %].txt

#----------------------------#
# gff to highlight
#----------------------------#
# coding and other features
perl -anl -e '/^#/ and next; $F[0] =~ s/\.\d+//; $color = q{}; $F[2] eq q{CDS} and $color = q{chr9}; $F[2] eq q{ncRNA} and $color = q{dark2-8-qual-1}; $F[2] eq q{rRNA} and $color = q{dark2-8-qual-2}; $F[2] eq q{tRNA} and $color = q{dark2-8-qual-3}; $F[2] eq q{tmRNA} and $color = q{dark2-8-qual-4}; $color and ($F[4] - $F[3] > 49) and print qq{$F[0] $F[3] $F[4] fill_color=$color};' [% working_dir %]/[% item.taxon %]/*.gff > highlight.features.[% item.taxon %].txt

# repeats
perl -anl -e '/^#/ and next; $F[0] =~ s/\.\d+//; $color = q{}; $F[2] eq q{region} and $F[8] =~ /mobile_element|Transposon/i and $color = q{chr15}; $F[2] =~ /repeat/ and $F[8] !~ /RNA/ and $color = q{chr15}; $color and ($F[4] - $F[3] > 49) and print qq{$F[0] $F[3] $F[4] fill_color=$color};' [% working_dir %]/[% item.taxon %]/*.gff > highlight.repeats.[% item.taxon %].txt

#----------------------------#
# run circos
#----------------------------#
perl ~/share/circos/bin/circos -noparanoid -conf circos.conf

[% END -%]

EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            findbin     => $FindBin::Bin,
            name_str    => $name_str,
            target_id   => $target_id,
            query_ids   => \@query_ids,
            data        => \@data,
        },
        File::Spec->catfile( $working_dir, "circos_cmd.sh" )
    ) or die Template->error;

    # feature_cmd.sh
    print "Create feature_cmd.sh\n";
    $text = <<'EOF';
#!/bin/bash
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]

[% FOREACH item IN data -%]
#----------------------------------------------------------#
# [% item.taxon %] [% item.name %]
#----------------------------------------------------------#

#----------------------------#
# gff to feature
#----------------------------#
cd [% working_dir %]/[% item.taxon %]_proc

# coding
perl -anl -e '/^#/ and next; $F[0] =~ s/\.\d+//; $F[2] eq q{CDS} and print qq{$F[0]\t$F[3]\t$F[4]};' [% working_dir %]/[% item.taxon %]/*.gff > feature.coding.[% item.taxon %].bed

# repeats
perl -anl -e '/^#/ and next; $F[0] =~ s/\.\d+//; $F[2] eq q{region} and $F[8] =~ /mobile_element|Transposon/i and print qq{$F[0]\t$F[3]\t$F[4]};' [% working_dir %]/[% item.taxon %]/*.gff > feature.repeats.[% item.taxon %].bed
perl -anl -e '/^#/ and next; $F[0] =~ s/\.\d+//; $F[2] =~ /repeat/ and $F[8] !~ /RNA/ and print qq{$F[0]\t$F[3]\t$F[4]};' [% working_dir %]/[% item.taxon %]/*.gff >> feature.repeats.[% item.taxon %].bed

# others
perl -anl -e '/^#/ and next; $F[0] =~ s/\.\d+//; $F[2] eq q{ncRNA} and print qq{$F[0]\t$F[3]\t$F[4]};' [% working_dir %]/[% item.taxon %]/*.gff > feature.ncRNA.[% item.taxon %].bed
perl -anl -e '/^#/ and next; $F[0] =~ s/\.\d+//; $F[2] eq q{rRNA} and print qq{$F[0]\t$F[3]\t$F[4]};' [% working_dir %]/[% item.taxon %]/*.gff > feature.rRNA.[% item.taxon %].bed
perl -anl -e '/^#/ and next; $F[0] =~ s/\.\d+//; $F[2] eq q{tRNA} and print qq{$F[0]\t$F[3]\t$F[4]};' [% working_dir %]/[% item.taxon %]/*.gff > feature.tRNA.[% item.taxon %].bed

#----------------------------#
# merge bed and stat
#----------------------------#
for ftr in coding repeats ncRNA rRNA tRNA
do
    if [ -s feature.$ftr.[% item.taxon %].bed ]
    then
        # there are some data in .bed file
        perl [% findbin %]/../ofg/bed_op.pl --op merge_to_runlist --file feature.$ftr.[% item.taxon %].bed --name feature.$ftr.[% item.taxon %].yml;
    else
        # .bed file is empty
        # create empty runlists from chr.sizes
        perl -ane'BEGIN { print qq{---\n} }; print qq{$F[0]: "-"\n}; END {print qq{\n}};' [% working_dir %]/[% item.taxon %]/chr.sizes > feature.$ftr.[% item.taxon %].yml;
    fi;
    perl [% findbin %]/../slice/stat_runlist.pl --size [% working_dir %]/[% item.taxon %]/chr.sizes -f feature.$ftr.[% item.taxon %].yml;
done

echo "feature,name,length,size,coverage" > [% working_dir %]/[% item.taxon %]_result/[% item.taxon %].feature.csv
for ftr in coding repeats ncRNA rRNA tRNA
do
    FTR=$ftr perl -nl -e '/^name/ and next; print qq{$ENV{FTR},$_};' feature.$ftr.[% item.taxon %].yml.csv;
done >> [% working_dir %]/[% item.taxon %]_result/[% item.taxon %].feature.csv

for ftr in coding repeats ncRNA rRNA tRNA
do
    perl [% findbin %]/../slice/compare_runlist.pl -op intersect --mk -f1 [% item.taxon %]vsselfalign.cc.runlist.yml -f2 feature.$ftr.[% item.taxon %].yml -o [% item.taxon %]vsselfalign.cc.runlist.$ftr.yml
done

for ftr in coding repeats ncRNA rRNA tRNA
do
    perl [% findbin %]/../slice/stat_runlist.pl --mk --size [% working_dir %]/[% item.taxon %]/chr.sizes -f [% item.taxon %]vsselfalign.cc.runlist.$ftr.yml;
done

echo "feature,copy,name,length,size,coverage" > [% working_dir %]/[% item.taxon %]_result/[% item.taxon %]vsselfalign.feature.copies.csv
for ftr in coding repeats ncRNA rRNA tRNA
do
    FTR=$ftr perl -nl -e '/^key/ and next; /\,all\,/ or next; print qq{$ENV{FTR},$_};' [% item.taxon %]vsselfalign.cc.runlist.$ftr.yml.csv;
done >> [% working_dir %]/[% item.taxon %]_result/[% item.taxon %]vsselfalign.feature.copies.csv

for ftr in coding repeats ncRNA rRNA tRNA
do
    rm feature.$ftr.[% item.taxon %].bed;
    rm feature.$ftr.[% item.taxon %].yml;
    rm feature.$ftr.[% item.taxon %].yml.csv;
    rm [% item.taxon %]vsselfalign.cc.runlist.$ftr.yml;
    rm [% item.taxon %]vsselfalign.cc.runlist.$ftr.yml.csv;
done

[% END -%]

EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            findbin     => $FindBin::Bin,
            name_str    => $name_str,
            target_id   => $target_id,
            query_ids   => \@query_ids,
            data        => \@data,
        },
        File::Spec->catfile( $working_dir, "feature_cmd.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]

[% FOREACH item IN data -%]
#----------------------------------------------------------#
# [% item.taxon %] [% item.name %]
#----------------------------------------------------------#
# gen_alignDB
perl [% findbin %]/../extra/two_way_batch.pl \
    -d [% item.taxon %]vs[% item.taxon %] \
    -t [% item.taxon %],[% item.taxon %] \
    -q [% item.taxon %],[% item.taxon %] \
    -da [% working_dir %]/[% item.taxon %]vsselfalign \
    --gff_files [% FOREACH acc IN acc_of.${item.taxon} %][% working_dir %]/[% item.taxon %]/[% acc %].gff,[% END %] \
    --rm_gff_files [% FOREACH acc IN acc_of.${item.taxon} %][% working_dir %]/[% item.taxon %]/[% acc %].rm.gff,[% END %] \
    -lt 1000 --parallel [% parallel %] --run 1-5,21,40

[% END -%]

#----------------------------------------------------------#
# [% name_str %]
#----------------------------------------------------------#
# init db
perl [% findbin %]/../extra/two_way_batch.pl \
    -d [% name_str %]_paralog \
    -r 1

#----------------------------#
# gen_alignDB.pl
#----------------------------#

[% FOREACH item IN data -%]
# [% item.taxon %]
# gen_alignDB to existing database
perl [% findbin %]/../extra/two_way_batch.pl \
    -d [% name_str %]_paralog \
    -t [% item.taxon %],[% item.name %] \
    -q [% item.taxon %],[% item.name %] \
    -da [% working_dir %]/[% item.taxon %]vsselfalign \
    -lt 1000 --parallel [% parallel %] --run 2

[% END -%]

#----------------------------#
# rest steps
#----------------------------#
perl [% findbin %]/../extra/two_way_batch.pl \
    -d [% name_str %]_paralog \
[% FOREACH item IN data -%]
    --gff_files [% FOREACH acc IN acc_of.${item.taxon} %][% working_dir %]/[% item.taxon %]/[% acc %].gff,[% END %] \
    --rm_gff_files [% FOREACH acc IN acc_of.${item.taxon} %][% working_dir %]/[% item.taxon %]/[% acc %].rm.gff,[% END %] \
[% END -%]
    -lt 1000 --parallel [% parallel %] --batch 5 \
    --run 5,10,21,30-32,40-42,44

EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            findbin     => $FindBin::Bin,
            name_str    => $name_str,
            data        => \@data,
            acc_of      => $acc_of,
        },
        File::Spec->catfile( $working_dir, "pair_stat.sh" )
    ) or die Template->error;

    # cmd.bat
    print "Create cmd.bat\n";
    $text = <<'EOF';
REM strain_bz_self.pl
REM perl [% stopwatch.cmd_line %]

REM basicstat
perl [% bat_dir %]/fig/collect_common_basic.pl -d .

REM common chart
if exist [% name_str %]_paralog.common.xlsx perl [% bat_dir %]/stat/common_chart_factory.pl -i [% name_str %]_paralog.common.xlsx

REM multi chart
if exist [% name_str %]_paralog.multi.xlsx  perl [% bat_dir %]/stat/multi_chart_factory.pl -i [% name_str %]_paralog.multi.xlsx

REM gc chart
if exist [% name_str %]_paralog.gc.xlsx     perl [% bat_dir %]/stat/gc_chart_factory.pl --add_trend 1 -i [% name_str %]_paralog.gc.xlsx

EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            findbin     => $FindBin::Bin,
            bat_dir     => $bat_dir,
            name_str    => $name_str,
        },
        File::Spec->catfile( $working_dir, "chart.bat" )
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
