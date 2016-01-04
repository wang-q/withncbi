#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use DBI;
use Text::CSV_XS;
use DateTime::Format::Natural;
use List::MoreUtils qw(any all uniq);
use Template;

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

strain_bz_self.pl - Full procedure for self genome alignments.

=head1 SYNOPSIS

    perl strain_bz_self.pl [options]
      Options:
        --help          -?          brief help message
        --file              STR     All taxons in this project (may also contain unused taxons)
        --working_dir   -w  STR     Default is [.]
        --seq_dir       -s  STR     Will do prep_fa() from this dir or use seqs store in $working_dir
        --keep              @STR    don't touch anything inside fasta files
        --target_id     -t  STR
        --query_ids     -q  @STR
        --length            INT     Minimal length of paralogous fragments
        --name_str      -n  STR     Default is []
        --use_name      -un         Use name instead of taxon_id as identifier. These names should only contain
                                    alphanumeric value and match with sequence directory names.
                                    For strains not recorded in NCBI taxonomy, you should assign them fake ids.
                                    If this option set to be true, all $target_id, @query_ids are actually names.        
        --msa               STR     Aligning program for refine. Default is [mafft]
        --noblast                   Don't blast against genomes
        --norm                      RepeatMasker has been done.
        --nostat                    Don't do stat stuffs
        --norawphylo                Skip rawphylo
        --parallel          INT     number of child processes

=cut

my $aligndb = path( $Config->{run}{aligndb} )->stringify;    # alignDB path
my $egaz    = path( $Config->{run}{egaz} )->stringify;       # egaz path
my $egas    = path( $Config->{run}{egas} )->stringify;       # egas path
my $blast   = path( $Config->{run}{blast} )->stringify;      # blast path
my $circos  = path( $Config->{run}{circos} )->stringify;     # circos path
my $bat_dir = $Config->{run}{bat};                           # Windows alignDB path

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'file=s'          => \( my $taxon_file  = "strains_taxon_info.csv" ),
    'working_dir|w=s' => \( my $working_dir = "." ),
    'seq_dir|s=s'     => \my $seq_dir,
    'keep=s'          => \my @keep,
    'target_id|t=s'   => \my $target_id,
    'query_ids|q=s'   => \my @query_ids,
    'length=i'        => \( my $length      = 1000 ),
    'name_str|n=s'    => \( my $name_str    = "working" ),
    'use_name|un'     => \my $use_name,
    'msa=s'           => \( my $msa         = 'mafft' ),
    'noblast'           => \my $noblast,
    'norm'            => \my $norm,
    'nostat'          => \my $nostat,
    'norawphylo'      => \my $norawphylo,
    'parallel=i'      => \( my $parallel    = $Config->{run}{parallel} ),
) or HelpMessage(1);

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
# genome sequences
#----------------------------#
echo "* Recreate genome.fa"
sleep 1;
find [% working_dir %]/Genomes/[% id %] -type f -name "*.fa" \
    | sort | xargs cat \
    | perl -nl -e '/^>/ or $_ = uc; print' \
    > genome.fa
faops size genome.fa > chr.sizes

#----------------------------#
# Correct genomic positions
#----------------------------#
echo "* Correct genomic positions"
sleep 1;

echo "* axt2fas"
fasops axt2fas [% working_dir %]/Pairwise/[% id %]vsselfalign/axtNet/*.axt.gz -l [% length %] -o stdout > axt.fas
fasops separate axt.fas -o [% working_dir %]/Processing/[% id %] --nodash -s .sep.fasta

echo "* Target positions"
perl [% egas %]/sparsemem_exact.pl -f target.sep.fasta -g genome.fa \
    --length 500 -o replace.target.tsv
fasops replace axt.fas replace.target.tsv -o axt.target.fas

echo "* Query positions"
perl [% egas %]/sparsemem_exact.pl -f query.sep.fasta -g genome.fa \
    --length 500 -o replace.query.tsv
fasops replace axt.target.fas replace.query.tsv -o axt.correct.fas

#----------------------------#
# Coverage stats
#----------------------------#
echo "* Coverage stats"
sleep 1;
fasops covers axt.correct.fas -o axt.correct.yml
runlist split axt.correct.yml -s .temp.yml
runlist compare --op union target.temp.yml query.temp.yml -o axt.union.yml
runlist stat --size chr.sizes axt.union.yml -o [% working_dir %]/Results/[% id %]/S288c.union.csv

# links by lastz-chain
fasops links axt.correct.fas -o stdout \
    | perl -nl -e 's/(target|query)\.//g; print;' \
    > links.lastz.tsv

# remove species names
fasops separate axt.correct.fas --nodash -o stdout \
    | perl -nl -e '/^>/ and s/^>(target|query)\./\>/; print;' \
    > axt.gl.fasta

[% IF noblast -%]
#----------------------------#
# Lastz paralogs
#----------------------------#
cat axt.gl.fasta > axt.all.fasta
[% ELSE -%]
#----------------------------#
# Get more paralogs
#----------------------------#
echo "* Get more paralogs"
perl [% egas %]/blastn_genome.pl -c 0.95 -f axt.gl.fasta -g genome.fa -o axt.bg.fasta
cat axt.gl.fasta axt.bg.fasta > axt.all.fasta
[% END -%]

#----------------------------#
# Link paralogs
#----------------------------#
echo "* Link paralogs"
sleep 1;
perl [% egas %]/blastn_paralog.pl -f axt.all.fasta -c 0.95 -o links.blast.tsv

#----------------------------#
# Merge paralogs
#----------------------------#
echo "* Merge paralogs"
sleep 1;

perl [% egas %]/merge_node.pl -v -c 0.95 -o [% id %].merge.yml \
[% IF noblast -%]
    -f links.lastz.tsv
[% ELSE -%]
    -f links.lastz.tsv -f links.blast.tsv
[% END -%]

perl [% egas %]/paralog_graph.pl -v -m [% id %].merge.yml --nonself -o [% id %].merge.graph.yml \
[% IF noblast -%]
    -f links.lastz.tsv
[% ELSE -%]
    -f links.lastz.tsv -f links.blast.tsv
[% END -%]

echo "* CC sequences and stats"
perl [% egas %]/cc.pl           -f [% id %].merge.graph.yml
perl [% egas %]/proc_cc_chop.pl -f [% id %].cc.yml --size chr.sizes --genome genome.fa --msa [% msa %]
perl [% egas %]/proc_cc_stat.pl -f [% id %].cc.yml --size chr.sizes

echo "* Coverage figure"
runlist stat --size chr.sizes [% id %].cc.chr.runlist.yml;
perl [% egas %]/cover_figure.pl --size chr.sizes -f [% id %].cc.chr.runlist.yml;

#----------------------------#
# result
#----------------------------#
echo "* Results"
sleep 1;

mv [% id %].cc.pairwise.fas         [% working_dir %]/Results/[% id %]
cp [% id %].cc.yml                  [% working_dir %]/Results/[% id %]
mv [% id %].cc.csv                  [% working_dir %]/Results/[% id %]
mv [% id %].cc.chr.runlist.yml.csv  [% working_dir %]/Results/[% id %]/[% id %].chr.csv
mv [% id %].cc.chr.runlist.yml      [% working_dir %]/Results/[% id %]/[% id %].paralog.yml
mv [% id %].cc.chr.runlist.png      [% working_dir %]/Results/[% id %]/[% id %].chr.png

#----------------------------#
# clean
#----------------------------#
echo "* Clean up"
sleep 1;

find [% working_dir %]/Processing/[% id %] -type f -name "*genome.fa*"   | parallel --no-run-if-empty rm
find [% working_dir %]/Processing/[% id %] -type f -name "*gl.fasta*"    | parallel --no-run-if-empty rm
find [% working_dir %]/Processing/[% id %] -type f -name "*.sep.fasta"   | parallel --no-run-if-empty rm
find [% working_dir %]/Processing/[% id %] -type f -name "axt.*"         | parallel --no-run-if-empty rm
find [% working_dir %]/Processing/[% id %] -type f -name "replace.*.tsv" | parallel --no-run-if-empty rm

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
            noblast       => $noblast,
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
file          = [% taxon_id %].cc.link4.txt
radius        = 0.88r
bezier_radius = 0.2r
color         = purple
thickness     = 2
ribbon        = yes
stroke_color  = purple
stroke_thickness = 2
</link>

<link>
file          = [% taxon_id %].cc.link3.txt
radius        = 0.88r
bezier_radius = 0.1r
color         = dgreen
thickness     = 3
ribbon        = yes
stroke_color  = dgreen
stroke_thickness = 2
</link>

<link>
file          = [% taxon_id %].cc.link2.txt
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
file = [% taxon_id %].cc.linkN.txt
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
    runlist stat --size chr.sizes feature.$ftr.[% id %].yml;
done

echo "feature,name,length,size,coverage" > [% working_dir %]/Results/[% id %]/[% id %].feature.csv
for ftr in coding repeats ncRNA rRNA tRNA
do
    FTR=$ftr perl -nl -e '/^name/ and next; print qq{$ENV{FTR},$_};' feature.$ftr.[% id %].yml.csv;
done >> [% working_dir %]/Results/[% id %]/[% id %].feature.csv

for ftr in coding repeats ncRNA rRNA tRNA
do
    runlist compare --op intersect --mk [% id %].cc.runlist.yml feature.$ftr.[% id %].yml -o [% id %].cc.runlist.$ftr.yml
done

for ftr in coding repeats ncRNA rRNA tRNA
do
    runlist stat --mk --size chr.sizes [% id %].cc.runlist.$ftr.yml;
done

echo "feature,copy,name,length,size,coverage" > [% working_dir %]/Results/[% id %]/[% id %].feature.copies.csv
for ftr in coding repeats ncRNA rRNA tRNA
do
    FTR=$ftr perl -nl -e '/^key/ and next; /\,all\,/ or next; print qq{$ENV{FTR},$_};' [% id %].cc.runlist.$ftr.yml.csv;
done >> [% working_dir %]/Results/[% id %]/[% id %].feature.copies.csv

for ftr in coding repeats ncRNA rRNA tRNA
do
    rm feature.$ftr.[% id %].bed;
    rm feature.$ftr.[% id %].yml;
    rm feature.$ftr.[% id %].yml.csv;
    rm [% id %].cc.runlist.$ftr.yml;
    rm [% id %].cc.runlist.$ftr.yml.csv;
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
