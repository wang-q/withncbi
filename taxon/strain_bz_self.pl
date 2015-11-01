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

use Path::Tiny;
use File::Find::Rule;

use AlignDB::Stopwatch;

use FindBin;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::Bin/../config.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

my $working_dir = ".";
my $seq_dir;    #  will do prep_fa() from this dir or use seqs store in $working_dir
my @keep;       # don't touch anything inside fasta files

my $target_id;
my @query_ids;

my $name_str = "working";

# All taxons in this project (may also contain unused taxons)
my $taxon_file = "strains_taxon_info.csv";

my $msa    = 'mafft';    # Default alignment program
my $length = 1000;

my $aligndb  = path( $Config->{run}{aligndb} )->stringify;     # alignDB path
my $egaz     = path( $Config->{run}{egaz} )->stringify;        # egaz path
my $egas     = path( $Config->{run}{egas} )->stringify;        # egas path
my $blast    = path( $Config->{run}{blast} )->stringify;       # blast path
my $circos   = path( $Config->{run}{circos} )->stringify;      # circos path
my $bat_dir  = $Config->{run}{bat};                            # Windows alignDB path

# Use name instead of taxon_id as identifier. These names should only contain
# alphanumeric value and match with sequence directory names.
# For strains not recorded in NCBI taxonomy, you should assign them fake ids.
# If this option set to be true, all $target_id, @query_ids is actually names.
my $use_name;

# RepeatMasker has been done.
my $norm;

# Don't do stat stuffs
my $nostat;

# Useless here, just keep compatibility with strain_bz.pl
my $norawphylo;

# run in parallel mode
my $parallel = $Config->{run}{parallel};

my $man  = 0;
my $help = 0;

GetOptions(
    'help'            => \$help,
    'man'             => \$man,
    'file=s'          => \$taxon_file,
    'w|working_dir=s' => \$working_dir,
    's|seq_dir=s'     => \$seq_dir,
    'keep=s'          => \@keep,
    't|target_id=s'   => \$target_id,
    'q|query_ids=s'   => \@query_ids,
    'length=i'        => \$length,
    'n|name_str=s'    => \$name_str,
    'un|use_name'     => \$use_name,
    'msa=s'           => \$msa,
    'norm'            => \$norm,
    'nostat'          => \$nostat,
    'norawphylo'      => \$norawphylo,
    'parallel=i'      => \$parallel,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Writing strains summary...");

die "$taxon_file doesn't exist\n" unless -e $taxon_file;

# prepare working dir
{
    print "Working on $name_str\n";
    $working_dir = path( $working_dir, $name_str )->absolute;
    $working_dir->mkpath;
    $working_dir = $working_dir->stringify;
    print " " x 4, "Working dir is $working_dir\n";

    path( $working_dir, 'Genomes' )->mkpath;
    path( $working_dir, 'Pairwise' )->mkpath;
    path( $working_dir, 'Processing' )->mkpath;
    path( $working_dir, 'Results' )->mkpath;
}

my @data;
{
    if ( !$use_name ) {
        @data
            = map { taxon_info( $_, $working_dir ) } ( $target_id, @query_ids );
    }
    else {
        @data = map { taxon_info_name( $_, $working_dir ) } ( $target_id, @query_ids );
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

        my $original_dir = path( $seq_dir, $id )->stringify;
        my $cur_dir = path( $working_dir, 'Genomes', $id );
        $cur_dir->mkpath;
        $cur_dir = $cur_dir->stringify;

        my @fa_files
            = File::Find::Rule->file->name( '*.fna', '*.fa', '*.fas', '*.fasta' )
            ->in($original_dir);

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

my $acc_of = {};
{
    for my $id ( $target_id, @query_ids ) {
        my $cur_dir = path( $working_dir, $id )->stringify;
        my @fa_files
            = File::Find::Rule->file->name( '*.fna', '*.fa', '*.fas', '*.fasta' )->in($cur_dir);

        $acc_of->{$id} = [];
        for my $fa_file (@fa_files) {
            my $basename = basename( $fa_file, '.fna', '.fa', '.fas', '.fasta' );
            push @{ $acc_of->{$id} }, $basename;
        }
    }
}

{
    print "Create id2name.csv\n";
    my $fh = path( $working_dir, "id2name.csv" )->openw;
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
    $tt->process( \$text, { data => \@data, }, path( $working_dir, "taxon.csv" )->stringify )
        or die Template->error;

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
taxon_id,chr,length,name,assembly
[% FOREACH item IN data -%]
[% item.taxon %],chrUn,999999999,[% item.name %]
[% END -%]
DELIMITER

if [ -f real_chr.csv ]; then
    rm real_chr.csv;
fi;

[% FOREACH item IN data -%]
if [ ! -f [% item.dir %]/chr.sizes ]; then
    faops size [% item.dir %]/*.fa > [% item.dir %]/chr.sizes;
fi;
perl -aln -F"\t" -e 'print qq{[% item.taxon %],$F[0],$F[1],[% item.name %]}' [% item.dir %]/chr.sizes >> real_chr.csv;
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
                target_id   => $target_id,
                query_ids   => \@query_ids,
            },
            path( $working_dir, $sh_name )->stringify
        ) or die Template->error;
    }

    # self_cmd.sh
    $sh_name = "3_self_cmd.sh";
    print "Create $sh_name\n";
    $text = <<'EOF';
#!/bin/bash
# strain_bz_self.pl
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]

sleep 1;

[% FOREACH id IN all_ids -%]
#----------------------------------------------------------#
# [% id %]
#----------------------------------------------------------#
if [ -d [% working_dir %]/Pairwise/[% id %]vsselfalign ]
then
    rm -fr [% working_dir %]/Pairwise/[% id %]vsselfalign
fi

#----------------------------#
# self bz
#----------------------------#
perl [% egaz %]/bz.pl \
    --is_self \
    -s set01 -C 0 --noaxt \
    -dt [% working_dir %]/Genomes/[% id %] \
    -dq [% working_dir %]/Genomes/[% id %] \
    -dl [% working_dir %]/Pairwise/[% id %]vsselfalign \
    --parallel [% parallel %]

#----------------------------#
# lpcna
#----------------------------#
perl [% egaz %]/lpcna.pl \
    -dt [% working_dir %]/Genomes/[% id %] \
    -dq [% working_dir %]/Genomes/[% id %] \
    -dl [% working_dir %]/Pairwise/[% id %]vsselfalign \
    --parallel [% parallel %]

[% END -%]

EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            egaz        => $egaz,
            name_str    => $name_str,
            all_ids     => [ $target_id, @query_ids ],
        },
        path( $working_dir, $sh_name )->stringify
    ) or die Template->error;

    # proc_cmd.sh
    $sh_name = "4_proc_cmd.sh";
    print "Create $sh_name\n";
    $text = <<'EOF';
#!/bin/bash
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]

sleep 1;

[% FOREACH id IN all_ids -%]
#----------------------------------------------------------#
# [% id %]
#----------------------------------------------------------#
if [ -d [% working_dir %]/Processing/[% id %] ]
then
    find [% working_dir %]/Processing/[% id %] -type f -not -name "circos.conf" | parallel --no-run-if-empty rm
else
    mkdir -p [% working_dir %]/Processing/[% id %]
fi

if [ ! -d [% working_dir %]/Results/[% id %] ]
then
    mkdir -p [% working_dir %]/Results/[% id %]
fi

cd [% working_dir %]/Processing/[% id %]

#----------------------------#
# quick and dirty coverage
#----------------------------#
perl [% egas %]/gather_info_axt.pl -l [% length %] -d [% working_dir %]/Pairwise/[% id %]vsselfalign --nomatch
perl [% egas %]/compare_runlist.pl -op intersect -f1 [% id %]vsselfalign.runlist.0.yml -f2 [% id %]vsselfalign.runlist.1.yml -o [% id %]vsselfalign.runlist.i.yml
perl [% egas %]/compare_runlist.pl -op xor -f1 [% id %]vsselfalign.runlist.0.yml -f2 [% id %]vsselfalign.runlist.1.yml -o [% id %]vsselfalign.runlist.x.yml
perl [% egas %]/compare_runlist.pl -op union -f1 [% id %]vsselfalign.runlist.0.yml -f2 [% id %]vsselfalign.runlist.1.yml -o [% id %]vsselfalign.runlist.u.yml

for op in 0 1 i x u
do
    perl [% egas %]/stat_runlist.pl --size [% working_dir %]/Genomes/[% id %]/chr.sizes -f [% id %]vsselfalign.runlist.$op.yml;
done

echo "group,name,length,size,coverage" > [% working_dir %]/Results/[% id %]/[% id %].[% length %].csv
for op in 0 1 i x u
do
    OP=$op perl -nl -e '/^name/ and next; print qq{$ENV{OP},$_};' [% id %]vsselfalign.runlist.$op.yml.csv;
done >> [% working_dir %]/Results/[% id %]/[% id %].[% length %].csv

for op in 0 1 i x u
do
    rm [% id %]vsselfalign.runlist.$op.yml;
    rm [% id %]vsselfalign.runlist.$op.yml.csv;
done

#----------------------------#
# genome blast
#----------------------------#
find [% working_dir %]/Genomes/[% id %] -type f -name "*.fa" \
    | sort | xargs cat \
    | perl -nl -e '/^>/ or uc; print' \
    > [% id %].genome.fasta

echo "* build genome blast db [% id %].genome.fasta"
[% blast %]/bin/formatdb -p F -o T -i [% id %].genome.fasta

perl [% egas %]/gather_seq_axt.pl -l [% length %] -d [% working_dir %]/Pairwise/[% id %]vsselfalign

echo "* blast [% id %]vsselfalign.axt.fasta"
[% blast %]/bin/blastall -p blastn -F "m D" -m 0 -b 10 -v 10 -e 1e-3 -a 8 -i [% id %]vsselfalign.axt.fasta -d [% id %].genome.fasta -o [% id %]vsselfalign.axt.blast

#----------------------------#
# paralog sequences
#----------------------------#
# XXX Omit genome locations in .axt
# XXX There are errors, especially for queries
perl [% egas %]/blastn_genome_location.pl -f [% id %]vsselfalign.axt.blast -m 0 -i 90 -c 0.95

#----------------------------#
# paralog blast
#----------------------------#
echo "* build paralog blast db [% id %]vsselfalign.gl.fasta"
[% blast %]/bin/formatdb -p F -o T -i [% id %]vsselfalign.gl.fasta

echo "* blast [% id %]vsselfalign.gl.fasta"
[% blast %]/bin/blastall -p blastn -F "m D" -m 0 -b 10 -v 10 -e 1e-3 -a 8 -i [% id %]vsselfalign.gl.fasta -d [% id %]vsselfalign.gl.fasta -o [% id %]vsselfalign.gl.blast

#----------------------------#
# merge
#----------------------------#
perl [% egas %]/blastn_paralog.pl -f [% id %]vsselfalign.gl.blast -m 0 -i 90 -c 0.9

perl [% egas %]/merge_node.pl    -v -f [% id %]vsselfalign.blast.tsv -o [% id %]vsselfalign.merge.yml -c 0.9
perl [% egas %]/paralog_graph.pl -v -f [% id %]vsselfalign.blast.tsv -m [% id %]vsselfalign.merge.yml --nonself -o [% id %]vsselfalign.merge.graph.yml
perl [% egas %]/cc.pl               -f [% id %]vsselfalign.merge.graph.yml
perl [% egas %]/proc_cc_chop.pl     -f [% id %]vsselfalign.cc.yml --size [% working_dir %]/Genomes/[% id %]/chr.sizes --msa [% msa %]
perl [% egas %]/proc_cc_stat.pl     -f [% id %]vsselfalign.cc.yml --size [% working_dir %]/Genomes/[% id %]/chr.sizes

perl [% egas %]/stat_runlist.pl --size [% working_dir %]/Genomes/[% id %]/chr.sizes -f [% id %]vsselfalign.cc.chr.runlist.yml;

#----------------------------#
# result
#----------------------------#
cp [% id %]vsselfalign.cc.yml [% working_dir %]/Results/[% id %]
mv [% id %]vsselfalign.cc.csv [% working_dir %]/Results/[% id %]
cp [% id %]vsselfalign.cc.chr.runlist.yml.csv [% working_dir %]/Results/[% id %]/[% id %]vsselfalign.chr.csv

cp [% id %]vsselfalign.cc.pairwise.fas [% working_dir %]/Results/[% id %]

#----------------------------#
# clean
#----------------------------#
find [% working_dir %]/Processing/[% id %] -type f -name "*genome.fasta*" | parallel --no-run-if-empty rm
find [% working_dir %]/Processing/[% id %] -type f -name "*gl.fasta*" | parallel --no-run-if-empty rm
find [% working_dir %]/Processing/[% id %] -type f -name "*.blast" | parallel --no-run-if-empty rm

[% END -%]
EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            aligndb     => $aligndb,
            egaz        => $egaz,
            egas        => $egas,
            blast       => $blast,
            msa         => $msa,
            name_str    => $name_str,
            all_ids     => [ $target_id, @query_ids ],
            data        => \@data,
            length      => $length,
        },
        path( $working_dir, $sh_name )->stringify
    ) or die Template->error;

    # circos.conf
    $text = <<'EOF';
<image>
dir*   = [% working_dir %]/Results/[% taxon_id %]
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
                name_str    => $name_str,
                taxon_id    => $taxon_id,
            },
            path( $working_dir, 'Processing', "${taxon_id}", "circos.conf" )->stringify
        ) or die Template->error;
    }

    # circos_cmd.sh
    $sh_name = "5_circos_cmd.sh";
    print "Create $sh_name\n";
    $text = <<'EOF';
#!/bin/bash
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]

[% FOREACH id IN all_ids -%]
#----------------------------------------------------------#
# [% id %]
#----------------------------------------------------------#

#----------------------------#
# karyotype
#----------------------------#
cd [% working_dir %]/Processing/[% id %]

# generate karyotype file
perl -anl -e '$i++; print qq{chr - $F[0] $F[0] 0 $F[1] chr$i}' [% working_dir %]/Genomes/[% id %]/chr.sizes > karyotype.[% id %].txt

#----------------------------#
# gff to highlight
#----------------------------#
# coding and other features
perl -anl -e '
    /^#/ and next;
    $F[0] =~ s/\.\d+//;
    $color = q{};
    $F[2] eq q{CDS} and $color = q{chr9};
    $F[2] eq q{ncRNA} and $color = q{dark2-8-qual-1};
    $F[2] eq q{rRNA} and $color = q{dark2-8-qual-2};
    $F[2] eq q{tRNA} and $color = q{dark2-8-qual-3};
    $F[2] eq q{tmRNA} and $color = q{dark2-8-qual-4};
    $color and ($F[4] - $F[3] > 49) and print qq{$F[0] $F[3] $F[4] fill_color=$color};
    ' \
    [% working_dir %]/Genomes/[% id %]/*.gff \
    > highlight.features.[% id %].txt

# repeats
perl -anl -e '
    /^#/ and next;
    $F[0] =~ s/\.\d+//;
    $color = q{};
    $F[2] eq q{region} and $F[8] =~ /mobile_element|Transposon/i and $color = q{chr15};
    $F[2] =~ /repeat/ and $F[8] !~ /RNA/ and $color = q{chr15};
    $color and ($F[4] - $F[3] > 49) and print qq{$F[0] $F[3] $F[4] fill_color=$color};
    ' \
    [% working_dir %]/Genomes/[% id %]/*.gff \
    > highlight.repeats.[% id %].txt

#----------------------------#
# run circos
#----------------------------#
perl [% circos %]/bin/circos -noparanoid -conf circos.conf

[% END -%]

EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            circos      => $circos,
            name_str    => $name_str,
            all_ids     => [ $target_id, @query_ids ],
            data        => \@data,
        },
        path( $working_dir, $sh_name )->stringify
    ) or die Template->error;

    # feature_cmd.sh
    $sh_name = "6_feature_cmd.sh";
    print "Create $sh_name\n";
    $text = <<'EOF';
#!/bin/bash
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]

[% FOREACH id IN all_ids -%]
#----------------------------------------------------------#
# [% id %]
#----------------------------------------------------------#

#----------------------------#
# gff to feature
#----------------------------#
cd [% working_dir %]/Processing/[% id %]

# coding
perl -anl -e '
    /^#/ and next;
    $F[0] =~ s/\.\d+//;
    $F[2] eq q{CDS} and print qq{$F[0]\t$F[3]\t$F[4]};
    ' \
    [% working_dir %]/Genomes/[% id %]/*.gff \
    > feature.coding.[% id %].bed

# repeats
perl -anl -e '
    /^#/ and next;
    $F[0] =~ s/\.\d+//;
    $F[2] eq q{region} and $F[8] =~ /mobile_element|Transposon/i and print qq{$F[0]\t$F[3]\t$F[4]};
    ' \
    [% working_dir %]/Genomes/[% id %]/*.gff \
    > feature.repeats.[% id %].bed
perl -anl -e '
    /^#/ and next;
    $F[0] =~ s/\.\d+//;
    $F[2] =~ /repeat/ and $F[8] !~ /RNA/ and print qq{$F[0]\t$F[3]\t$F[4]};
    ' \
    [% working_dir %]/Genomes/[% id %]/*.gff \
    >> feature.repeats.[% id %].bed

# others
perl -anl -e '
    /^#/ and next;
    $F[0] =~ s/\.\d+//;
    $F[2] eq q{ncRNA} and print qq{$F[0]\t$F[3]\t$F[4]};
    ' \
    [% working_dir %]/Genomes/[% id %]/*.gff \
    > feature.ncRNA.[% id %].bed
perl -anl -e '
    /^#/ and next;
    $F[0] =~ s/\.\d+//;
    $F[2] eq q{rRNA} and print qq{$F[0]\t$F[3]\t$F[4]};
    ' \
    [% working_dir %]/Genomes/[% id %]/*.gff \
    > feature.rRNA.[% id %].bed
perl -anl -e '
    /^#/ and next;
    $F[0] =~ s/\.\d+//;
    $F[2] eq q{tRNA} and print qq{$F[0]\t$F[3]\t$F[4]};
    ' \
    [% working_dir %]/Genomes/[% id %]/*.gff \
    > feature.tRNA.[% id %].bed

#----------------------------#
# merge bed and stat
#----------------------------#
for ftr in coding repeats ncRNA rRNA tRNA
do
    if [ -s feature.$ftr.[% id %].bed ]
    then
        # there are some data in .bed file
        perl [% aligndb %]/ofg/bed_op.pl --op merge_to_runlist --file feature.$ftr.[% id %].bed --name feature.$ftr.[% id %].yml;
    else
        # .bed file is empty
        # create empty runlists from chr.sizes
        perl -ane'BEGIN { print qq{---\n} }; print qq{$F[0]: "-"\n}; END {print qq{\n}};' [% working_dir %]/Genomes/[% id %]/chr.sizes > feature.$ftr.[% id %].yml;
    fi;
    perl [% egas %]/stat_runlist.pl --size [% working_dir %]/Genomes/[% id %]/chr.sizes -f feature.$ftr.[% id %].yml;
done

echo "feature,name,length,size,coverage" > [% working_dir %]/Results/[% id %]/[% id %].feature.csv
for ftr in coding repeats ncRNA rRNA tRNA
do
    FTR=$ftr perl -nl -e '/^name/ and next; print qq{$ENV{FTR},$_};' feature.$ftr.[% id %].yml.csv;
done >> [% working_dir %]/Results/[% id %]/[% id %].feature.csv

for ftr in coding repeats ncRNA rRNA tRNA
do
    perl [% egas %]/compare_runlist.pl -op intersect --mk -f1 [% id %]vsselfalign.cc.runlist.yml -f2 feature.$ftr.[% id %].yml -o [% id %]vsselfalign.cc.runlist.$ftr.yml
done

for ftr in coding repeats ncRNA rRNA tRNA
do
    perl [% egas %]/stat_runlist.pl --mk --size [% working_dir %]/Genomes/[% id %]/chr.sizes -f [% id %]vsselfalign.cc.runlist.$ftr.yml;
done

echo "feature,copy,name,length,size,coverage" > [% working_dir %]/Results/[% id %]/[% id %]vsselfalign.feature.copies.csv
for ftr in coding repeats ncRNA rRNA tRNA
do
    FTR=$ftr perl -nl -e '/^key/ and next; /\,all\,/ or next; print qq{$ENV{FTR},$_};' [% id %]vsselfalign.cc.runlist.$ftr.yml.csv;
done >> [% working_dir %]/Results/[% id %]/[% id %]vsselfalign.feature.copies.csv

for ftr in coding repeats ncRNA rRNA tRNA
do
    rm feature.$ftr.[% id %].bed;
    rm feature.$ftr.[% id %].yml;
    rm feature.$ftr.[% id %].yml.csv;
    rm [% id %]vsselfalign.cc.runlist.$ftr.yml;
    rm [% id %]vsselfalign.cc.runlist.$ftr.yml.csv;
done

[% END -%]

EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            aligndb     => $aligndb,
            egas        => $egas,
            name_str    => $name_str,
            all_ids     => [ $target_id, @query_ids ],
            data        => \@data,
        },
        path( $working_dir, $sh_name )->stringify
    ) or die Template->error;

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
    | grep -v -E "(chr_length|id2name|taxon|fake_taxon)\.csv$" \
    | grep -v -F "fake_tree.nwk" \
    > file_list.temp.txt

tar -czvf [% name_str %].tar.gz -T file_list.temp.txt

rm file_list.temp.txt

EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            egaz        => $egaz,
            target_id   => $target_id,
            name_str    => $name_str,
        },
        path( $working_dir, $sh_name )->stringify
    ) or die Template->error;

    if ( !$nostat ) {
        $sh_name = "7_pair_stat.sh";
        print "Create $sh_name\n";
        $text = <<'EOF';
#!/bin/bash
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]

[% FOREACH id IN all_ids -%]
#----------------------------------------------------------#
# [% id %]
#----------------------------------------------------------#
# gen_alignDB
perl [% aligndb %]/extra/multi_way_batch.pl \
    -d [% id %]vs[% id %] \
    -da [% working_dir %]/Results/[% id %] \
    --gff_files [% FOREACH acc IN acc_of.${id} %][% working_dir %]/Genomes/[% item.taxon %]/[% acc %].gff,[% END %] \
    --rm_gff_files [% FOREACH acc IN acc_of.${id} %][% working_dir %]/Genomes/[% item.taxon %]/[% acc %].rm.gff,[% END %] \
    --id [% working_dir %]/id2name.csv \
    -taxon [% working_dir %]/taxon.csv \
    -chr [% working_dir %]/chr_length.csv \
    -lt 1000 --parallel [% parallel %] --run 1-5,21,40

[% END -%]

#----------------------------------------------------------#
# [% name_str %]
#----------------------------------------------------------#
# init db
perl [% aligndb %]/extra/multi_way_batch.pl \
    -d [% name_str %]_paralog \
    --id [% working_dir %]/id2name.csv \
    -taxon [% working_dir %]/taxon.csv \
    -chr [% working_dir %]/chr_length.csv \
    -r 1

#----------------------------#
# gen_alignDB.pl
#----------------------------#

[% FOREACH id IN all_ids -%]
# [% id %]
# gen_alignDB to existing database
perl [% aligndb %]/extra/multi_way_batch.pl \
    -d [% name_str %]_paralog \
    --id [% working_dir %]/id2name.csv \
    -da [% working_dir %]/Results/[% id %] \
    -lt 1000 --parallel [% parallel %] --run 2

[% END -%]

#----------------------------#
# rest steps
#----------------------------#
perl [% aligndb %]/extra/two_way_batch.pl \
    -d [% name_str %]_paralog \
[% FOREACH id IN all_ids -%]
    --gff_files [% FOREACH acc IN acc_of.${id} %][% working_dir %]/[% item.taxon %]/[% acc %].gff,[% END %] \
    --rm_gff_files [% FOREACH acc IN acc_of.${id} %][% working_dir %]/[% item.taxon %]/[% acc %].rm.gff,[% END %] \
[% END -%]
    -lt 1000 --parallel [% parallel %] --batch 5 \
    --run 5,10,21,30-32,40-42,44

EOF
        $tt->process(
            \$text,
            {   stopwatch   => $stopwatch,
                parallel    => $parallel,
                working_dir => $working_dir,
                aligndb     => $aligndb,
                name_str    => $name_str,
                all_ids     => [ $target_id, @query_ids ],
                data        => \@data,
                acc_of      => $acc_of,
            },
            path( $working_dir, $sh_name )->stringify
        ) or die Template->error;

        # chart.bat
        print "Create chart.bat\n";
        $text = <<'EOF';
REM strain_bz_self.pl
REM perl [% stopwatch.cmd_line %]

REM basicstat
perl [% bat_dir %]/fig_table/collect_common_basic.pl -d .

REM common chart
if exist [% name_str %]_paralog.common.xlsx perl [% bat_dir %]/alignDB/stat/common_chart_factory.pl -i [% name_str %]_paralog.common.xlsx

REM multi chart
if exist [% name_str %]_paralog.multi.xlsx  perl [% bat_dir %]/alignDB/stat/multi_chart_factory.pl -i [% name_str %]_paralog.multi.xlsx

REM gc chart
if exist [% name_str %]_paralog.gc.xlsx     perl [% bat_dir %]/alignDB/stat/gc_chart_factory.pl --add_trend 1 -i [% name_str %]_paralog.gc.xlsx

EOF
        $tt->process(
            \$text,
            {   stopwatch   => $stopwatch,
                parallel    => $parallel,
                working_dir => $working_dir,
                bat_dir     => $bat_dir,
                name_str    => $name_str,
            },
            path( $working_dir, "chart.bat" )->stringify
        ) or die Template->error;
    }

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
        col_names => [ map { ( $_, $_ . "_id" ) } qw{strain species genus family order} ],
    };

    my $query = qq{ SELECT strain_id, strain, genus, species FROM t0 WHERE strain_id = ? };
    my $sth   = $dbh->prepare($query);
    $sth->execute($taxon_id);
    my ( $taxonomy_id, $organism_name, $genus, $species ) = $sth->fetchrow_array;
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
        dir     => path( $working_dir, 'Genomes', $taxon_id )->stringify,
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
        col_names => [ map { ( $_, $_ . "_id" ) } qw{strain species genus family order} ],
    };

    my $query = qq{ SELECT strain_id, strain, genus, species FROM t0 WHERE strain = ? };
    my $sth   = $dbh->prepare($query);
    $sth->execute($name);
    my ( $taxonomy_id, $organism_name, $genus, $species ) = $sth->fetchrow_array;
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
        dir     => path( $working_dir, 'Genomes', $name )->stringify,
    };
}

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

perl strain_bz_self.pl 
