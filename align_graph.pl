#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use Bio::Seq;
use Bio::Graphics::Panel;
use Bio::Graphics::Feature;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::Ensembl;
use AlignDB::GC;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server     = $Config->{database}->{server};
my $port       = $Config->{database}->{port};
my $username   = $Config->{database}->{username};
my $password   = $Config->{database}->{password};
my $db         = $Config->{database}->{db};
my $ensembl_db = $Config->{database}->{ensembl};

# graph init values
my $CLASS = $Config->{graph}->{class};    # GD or GD::SVG
my $ftr   = 'Bio::Graphics::Feature';

my $normal_figure = $Config->{graph}->{normal_figure};
my $gc_figure     = $Config->{graph}->{gc_figure};
my $gc_wave_csv   = $Config->{graph}->{gc_wave_csv};

my $wave_window_size = $Config->{gc}->{wave_window_size};
my $wave_window_step = $Config->{gc}->{wave_window_step};
my $vicinal_size     = $Config->{gc}->{vicinal_size};
my $fall_range       = $Config->{gc}->{fall_range};

my $align_id = 1;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'server=s'   => \$server,
    'db=s'       => \$db,
    'username=s' => \$username,
    'password=s' => \$password,
    'ensembl=s'  => \$ensembl_db,
    'align_id=s' => \$align_id,
    'CLASS=s'    => \$CLASS,
    'size=i'  => \$wave_window_size,
    'step=i' => \$wave_window_step,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init objects and SQL queries
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Drawing picture for $db...");

my $gc_obj = AlignDB::GC->new(
    mysql            => "$db:$server",
    user             => $username,
    passwd           => $password,
    wave_window_size => $wave_window_size,
    wave_window_step => $wave_window_step,
    vicinal_size     => $vicinal_size,
    fall_range       => $fall_range,
);

my $ensembl = AlignDB::Ensembl->new(
    server => $server,
    db     => $ensembl_db,
    user   => $username,
    passwd => $password,
);

my $dbh = $gc_obj->dbh;

# get target and query names via AlignDB methods
my ( $target_name, $query_name ) = $gc_obj->get_names($align_id);

# get alignment info
my $align_query = qq{
    SELECT c.chr_name,
           a.align_length,
           s.chr_start,
           s.chr_end,
           t.target_seq,
           q.query_seq,
           a.comparable_runlist,
           a.indel_runlist,
           t.target_runlist,
           q.query_runlist
    FROM align a, target t, query q, sequence s, chromosome c
    WHERE a.align_id = ?
    AND a.align_id = q.align_id
    AND a.align_id = t.align_id
    AND t.seq_id = s.seq_id
    AND s.chr_id = c.chr_id
};
my $align_query_sth = $dbh->prepare($align_query);

# select all isws in this alignment
my $isw_query = q{
    SELECT i.isw_id, i.isw_start, i.isw_end, i.isw_pi,
           e.isw_feature1, e.isw_feature2,
           (i.isw_target_gc_ratio + i.isw_query_gc_ratio) / 2 isw_gc
    FROM isw i, isw_extra e, indel
    WHERE i.indel_id = indel.indel_id
    AND i.isw_id = e.isw_id
    AND indel.align_id = ?
};
my $isw_query_sth = $dbh->prepare($isw_query);

#----------------------------------------------------------#
# Get data
#----------------------------------------------------------#
print "Fetch sequences\n";
$align_query_sth->execute($align_id);

my ($chr_name,           $align_length,  $chr_start,
    $chr_end,            $target_seq,    $query_seq,
    $comparable_runlist, $indel_runlist, $target_runlist,
    $query_runlist
) = $align_query_sth->fetchrow_array;

$chr_name =~ s/chr0?//i;

# generate sets one by one
print "Generate sets\n";
my $align_set      = AlignDB::IntSpan->new("1-$align_length");
my $comparable_set = AlignDB::IntSpan->new($comparable_runlist);
my $indel_set      = AlignDB::IntSpan->new($indel_runlist);
my $target_set     = AlignDB::IntSpan->new($target_runlist);
my $query_set      = AlignDB::IntSpan->new($query_runlist);

# align_position to chr_position transforming array
# the first element [0] will be ignored
my @chr_pos;
my $gap_insert_count = 0;
for ( my $i = 1; $i <= $align_length; $i++ ) {
    my $current_base = substr( $target_seq, $i - 1, 1 );
    if ( $current_base eq '-' ) {
        $gap_insert_count++;
    }
    $chr_pos[$i] = $chr_start + $i - 1 - $gap_insert_count;
}

# make a new ensembl slice object
$ensembl->set_slice( $chr_name, $chr_start, $chr_end );

# isw
print "Prepare isw\n";
my $isw_seg        = [];
my $isw_coding_seq = [];
my $isw_repeat_seq = [];
my $isw_gc_seq     = [];
$isw_query_sth->execute($align_id);
while ( my @row4 = $isw_query_sth->fetchrow_array ) {
    my ( $isw_id, $isw_start, $isw_end, $isw_pi, $isw_coding, $isw_repeat,
        $isw_gc )
        = @row4;
    my $isw_ftr = $ftr->new(
        -start => $isw_start,
        -stop  => $isw_end,
        -type  => 'isw',
        -score => $isw_pi
    );
    push @$isw_seg, $isw_ftr;

    my $isw_coding_ftr = $ftr->new(
        -start => $isw_start,
        -stop  => $isw_end,
        -type  => 'isw',
        -score => $isw_coding
    );
    push @$isw_coding_seq, $isw_coding_ftr;

    my $isw_repeat_ftr = $ftr->new(
        -start => $isw_start,
        -stop  => $isw_end,
        -type  => 'isw',
        -score => $isw_repeat
    );
    push @$isw_repeat_seq, $isw_repeat_ftr;

    my $isw_gc_ftr = $ftr->new(
        -start => $isw_start,
        -stop  => $isw_end,
        -type  => 'isw',
        -score => $isw_gc
    );
    push @$isw_gc_seq, $isw_gc_ftr;
}

# gsw
print "Prepare gsw\n";
my @slidings = $gc_obj->gc_wave( $align_id, $comparable_set );

# write $gc_wave_csv file
if ($gc_wave_csv) {
    open my $csv_fh, ">", "$db-align-$align_id.gc.csv";
    for (@slidings) {
        $_->{set}->runlist =~ /^(\d+)/;
        print {$csv_fh} $1, ",", $_->{gc}, ",", $_->{high_low_flag}, "\n";
    }
    close $csv_fh;
}

my $sliding_seg = [];
my $crest_seg   = [];
my $trough_seg  = [];

for (@slidings) {

    if ( $_->{high_low_flag} eq 'trough' or $_->{high_low_flag} eq 'crest' ) {
        $_->{high_low_flag} =~ s/trough/T/i;
        $_->{high_low_flag} =~ s/crest/C/i;
    }

    my $sliding_ftr = $ftr->new(
        -segments => [ $_->{set}->spans ],
        -type     => 'slidings',
        -score    => $_->{gc},
    );
    push @$sliding_seg, $sliding_ftr;

    if ( $_->{high_low_flag} eq 'T' ) {
        push @$trough_seg, $sliding_ftr;
    }
    if ( $_->{high_low_flag} eq 'C' ) {
        push @$crest_seg, $sliding_ftr;
    }
}

#----------------------------------------------------------#
# Draw the normal figure
#----------------------------------------------------------#
print "Draw picture\n";

if ($normal_figure) {
    my $align_segment = $ftr->new(
        -start => 1,
        -end   => $align_length,
        -name  => $align_id,
        -type  => 'alignment'
    );
    my $panel = Bio::Graphics::Panel->new(
        -grid        => 1,
        -gridcolor   => 'lightcyan',
        -segment     => $align_segment,
        -spacing     => 15,
        -width       => 560,
        -pad_top     => 20,
        -pad_bottom  => 20,
        -pad_left    => 20,
        -pad_right   => 20,
        -key_style   => 'none',
        -image_class => $CLASS,
    );

    # text
    {
        $panel->add_track(
            $align_segment,
            -glyph => 'text_in_box',
            -text  => "DB: $db | ALIGN: $align_id | "
                . "Chr$chr_name $chr_start-$chr_end",
            -text_bgcolor => 'lightcyan',
            -height       => 10,
            -bgcolor      => 'yellow',
            -text_pad     => 4,
        );
    }

    # arrow
    {
        $panel->add_track(
            $align_segment,
            -glyph      => 'arrow',
            -double     => 1,
            -fgcolor    => 'red',
            -bump       => 0,
            -height     => 10,
            -arrowstyle => 'regular',
            -tick       => 2,
            -linewidth  => 1,
        );
    }

    # alignment
    {
        my $target_segment = $ftr->new(
            -segments => [ $target_set->spans ],
            -name     => 'target',
            -type     => 'alignment'
        );
        $panel->add_track(
            $target_segment,
            -glyph     => 'segments',
            -label     => $target_name,
            -bump      => 0,
            -height    => 10,
            -font      => 'gdSmallFont',
            -linewidth => 1,
            -bgcolor   => 'lightblue',
            -connector => 'solid',
        );
        my $query_segment = $ftr->new(
            -segments => [ $query_set->spans ],
            -name     => 'query',
            -type     => 'alignment'
        );
        $panel->add_track(
            $query_segment,
            -glyph     => 'segments',
            -label     => $query_name,
            -bump      => 0,
            -height    => 10,
            -font      => 'gdSmallFont',
            -linewidth => 1,
            -bgcolor   => 'lightblue',
            -connector => 'solid',
        );
    }

    # comparable_set
    {
        my $comparable_segment = $ftr->new(
            -segments => [ $comparable_set->spans ],
            -name     => 'comparable',
            -source   => 'confirmed',
            -type     => 'alignment'
        );
        $panel->add_track(
            $comparable_segment,
            -glyph     => 'segments',
            -label     => 'comparable',
            -bump      => 0,
            -height    => 10,
            -font      => 'gdSmallFont',
            -linewidth => 1,
            -bgcolor   => 'turquoise',
            -connector => 'solid',
        );
    }

    # indel
    {
        my $indel_segment = $ftr->new(
            -segments => [ $indel_set->spans ],
            -name     => 'indel',
            -source   => 'confirmed',
        );
        $panel->add_track(
            $indel_segment,
            -label   => 'indel',
            -bgcolor => 'yellow',
            -glyph   => 'triangle',
            -point   => 1,
            -orient  => 'N',
        );
    }

    # isw
    {

        # Pi
        my $isw_segment = $ftr->new(
            -segments => $isw_seg,
            -name     => 'isw',
            -source   => 'confirmed',
        );
        $panel->add_track(
            $isw_segment,
            -glyph      => 'xyplot',
            -label      => 'isw_pi',
            -graph_type => 'boxes',
            -bgcolor    => 'lightgreen',
            -scale      => 'left',
            -height     => 100,
        );

        # GC
        $isw_segment = $ftr->new(
            -segments => $isw_gc_seq,
            -name     => 'isw',
            -source   => 'confirmed',
        );
        $panel->add_track(
            $isw_segment,
            -glyph      => 'xyplot',
            -label      => 'isw_gc',
            -graph_type => 'boxes',
            -bgcolor    => 'green',
            -scale      => 'left',
            -height     => 100,
        );

        # coding
        $isw_segment = $ftr->new(
            -segments => $isw_coding_seq,
            -name     => 'isw_coding',
            -source   => 'confirmed',
        );
        $panel->add_track(
            $isw_segment,
            -glyph      => 'xyplot',
            -label      => 'isw_coding',
            -graph_type => 'boxes',
            -fgcolor    => 'orange',
            -bgcolor    => 'orange',
            -scale      => 'none',
            -height     => 10,
        );

        # repeat
        $isw_segment = $ftr->new(
            -segments => $isw_repeat_seq,
            -name     => 'isw_repeat',
            -source   => 'confirmed',
        );
        $panel->add_track(
            $isw_segment,
            -glyph      => 'xyplot',
            -label      => 'isw_repeat',
            -graph_type => 'boxes',
            -fgcolor    => 'pink',
            -bgcolor    => 'pink',
            -scale      => 'none',
            -height     => 10,
        );
    }

    my $gd = $panel->gd;
    my $type = ( $CLASS eq 'GD' ) ? 'png' : 'svg';
    open my $pic_fh, ">", "$db-align-$align_id.$type";
    binmode $pic_fh;
    print {$pic_fh} $gd->$type;
    close $pic_fh;
}

#----------------------------------------------------------#
# Draw the gc figure
#----------------------------------------------------------#
if ($gc_figure) {
    my $align_segment = $ftr->new(
        -start => 1,
        -end   => $align_length,
        -name  => $align_id,
        -type  => 'alignment'
    );
    my $panel = Bio::Graphics::Panel->new(
        -grid        => 1,
        -gridcolor   => 'lightcyan',
        -segment     => $align_segment,
        -spacing     => 15,
        -width       => 560,
        -pad_top     => 20,
        -pad_bottom  => 20,
        -pad_left    => 20,
        -pad_right   => 20,
        -key_style   => 'none',
        -image_class => $CLASS,
    );

    # text
    {
        $panel->add_track(
            $align_segment,
            -glyph => 'text_in_box',
            -text  => "DB: $db | ALIGN: $align_id | "
                . "Chr$chr_name $chr_start-$chr_end",
            -text_bgcolor => 'lightcyan',
            -height       => 10,
            -bgcolor      => 'yellow',
            -text_pad     => 4,
        );
    }

    # text
    {
        $panel->add_track(
            $align_segment,
            -glyph => 'text_in_box',
            -text  => "size: $wave_window_size | step: $wave_window_step | "
                . "fall range: $fall_range",
            -text_bgcolor => 'lightcyan',
            -height       => 10,
            -bgcolor      => 'yellow',
            -text_pad     => 4,
        );
    }

    # arrow
    {
        $panel->add_track(
            $align_segment,
            -glyph      => 'arrow',
            -double     => 1,
            -fgcolor    => 'red',
            -bump       => 0,
            -height     => 10,
            -arrowstyle => 'regular',
            -tick       => 2,
            -linewidth  => 1,
        );
    }

    # alignment
    {
        my $target_segment = $ftr->new(
            -segments => [ $target_set->spans ],
            -name     => 'target',
            -type     => 'alignment'
        );
        $panel->add_track(
            $target_segment,
            -glyph     => 'segments',
            -label     => $target_name,
            -bump      => 0,
            -height    => 10,
            -font      => 'gdSmallFont',
            -linewidth => 1,
            -bgcolor   => 'lightblue',
            -connector => 'solid',
        );
        my $query_segment = $ftr->new(
            -segments => [ $query_set->spans ],
            -name     => 'query',
            -type     => 'alignment'
        );
        $panel->add_track(
            $query_segment,
            -glyph     => 'segments',
            -label     => $query_name,
            -bump      => 0,
            -height    => 10,
            -font      => 'gdSmallFont',
            -linewidth => 1,
            -bgcolor   => 'lightblue',
            -connector => 'solid',
        );
    }

    # comparable_set
    {
        my $comparable_segment = $ftr->new(
            -segments => [ $comparable_set->spans ],
            -name     => 'comparable',
            -source   => 'confirmed',
            -type     => 'alignment'
        );
        $panel->add_track(
            $comparable_segment,
            -glyph     => 'segments',
            -label     => 'comparable',
            -bump      => 0,
            -height    => 10,
            -font      => 'gdSmallFont',
            -linewidth => 1,
            -bgcolor   => 'turquoise',
            -connector => 'solid',
        );
    }

    # indel
    {
        my $indel_segment = $ftr->new(
            -segments => [ $indel_set->spans ],
            -name     => 'indel',
            -source   => 'confirmed',
        );
        $panel->add_track(
            $indel_segment,
            -label   => 'indel',
            -bgcolor => 'yellow',
            -glyph   => 'triangle',
            -point   => 1,
            -orient  => 'N',
        );
    }

    # crest
    {
        my $crest_segment = $ftr->new(
            -segments => $crest_seg,
            -name     => 'crest',
            -source   => 'confirmed',
        );
        $panel->add_track(
            $crest_segment,
            -label   => 'crest',
            -bump    => 0,
            -bgcolor => 'violet',
            -glyph   => 'triangle',
            -point   => 0,
            -orient  => 'S',
        );
    }

    # sliding, GC
    {
        my $sliding_segment = $ftr->new(
            -segments => $sliding_seg,
            -name     => 'sliding',
            -source   => 'confirmed',
        );
        $panel->add_track(
            $sliding_segment,
            -glyph      => 'xyplot',
            -label      => 'gc_wave',
            -graph_type => 'line',
            -bgcolor    => 'magenta',
            -fgcolor    => 'magenta',
            -scale      => 'left',
            -height     => 120,
        );
    }

    # trough
    {
        my $trough_segment = $ftr->new(
            -segments => $trough_seg,
            -name     => 'trough',
            -source   => 'confirmed',
        );
        $panel->add_track(
            $trough_segment,
            -label   => 'trough',
            -bump    => 0,
            -bgcolor => 'violet',
            -glyph   => 'triangle',
            -point   => 0,
            -orient  => 'N',
        );
    }

    my $gd = $panel->gd;
    my $type = ( $CLASS eq 'GD' ) ? 'png' : 'svg';
    open my $pic_fh, ">", "$db-align-$align_id.gc.$type";
    binmode $pic_fh;
    print {$pic_fh} $gd->$type;
    close $pic_fh;
}

$stopwatch->end_message;
exit;

__END__


=head1 NAME

    align_graph.pl - Generate graph for one alignment in alignDB

=head1 SYNOPSIS

    align_graph.pl [options]
     Options:
       --help            brief help message
       --man             full documentation
       --server          MySQL server IP/Domain name
       --db              database name
       --username        username
       --password        password
       --ensembl         ensembl database name
       --align_id        align_id
       

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do someting
useful with the contents thereof.

=cut
