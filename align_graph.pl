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
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db       = $Config->{database}{db};

# graph init values
my $CLASS = $Config->{graph}{class};    # GD or GD::SVG
my $width = $Config->{graph}{width};

my $figure      = $Config->{graph}{figure};
my $gc_wave_csv = $Config->{graph}{gc_wave_csv};

my $wave_window_size = $Config->{gc}{wave_window_size};
my $wave_window_step = $Config->{gc}{wave_window_step};
my $vicinal_size     = $Config->{gc}{vicinal_size};
my $fall_range       = $Config->{gc}{fall_range};

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
    'align_id=s' => \$align_id,
    'CLASS=s'    => \$CLASS,
    'width=i'    => \$width,
    'size=i'     => \$wave_window_size,
    'step=i'     => \$wave_window_step,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init objects and SQL queries
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Drawing picture for $db...");

my $obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);
Moo::Role->apply_roles_to_object( $obj, qw{ AlignDB::GC } );
my %opt = (
    wave_window_size => $wave_window_size,
    wave_window_step => $wave_window_step,
    vicinal_size     => $vicinal_size,
    fall_range       => $fall_range,
);
for my $key ( sort keys %opt ) {
    $obj->$key( $opt{$key} );
}

my $dbh = $obj->dbh;

#----------------------------------------------------------#
# Get data
#----------------------------------------------------------#
my $ftr = 'Bio::Graphics::Feature';

print "Fetch info\n";

# get target and query names via AlignDB methods
my ( $target_name, $query_name ) = $obj->get_names($align_id);

my $target_info  = $obj->get_target_info($align_id);
my $chr_name     = $target_info->{chr_name};
my $chr_start    = $target_info->{chr_start};
my $chr_end      = $target_info->{chr_end};
my $align_length = $target_info->{align_length};
my $target_set   = AlignDB::IntSpan->new( $target_info->{seq_runlist} );

my $query_info = $obj->get_query_info($align_id);
my $query_set  = AlignDB::IntSpan->new( $query_info->{seq_runlist} );

my ( undef, $comparable_set, $indel_set ) = @{ $obj->get_sets($align_id) };

# align features
my ( $coding_set, $repeat_set );
{
    print "Prepare align\n";
    my $query = q{
        SELECT align_coding_runlist, align_repeats_runlist
        FROM align
        WHERE align_id = ?
    };
    my $sth = $dbh->prepare($query);

    $sth->execute($align_id);
    my ( $coding_runlist, $repeat_runlist ) = $sth->fetchrow_array;
    $coding_set = AlignDB::IntSpan->new($coding_runlist);
    $repeat_set = AlignDB::IntSpan->new($repeat_runlist);
}

# isw
my $isw_pi_seg = [];
my ( $isw_pi_max, $isw_pi_min ) = ( 0, 0 );
{
    print "Prepare isw\n";
    my $isw_query = q{
        SELECT i.isw_start, i.isw_end, i.isw_pi
        FROM isw i
        INNER JOIN indel ON i.indel_id = indel.indel_id
        WHERE indel.align_id = ?
        ORDER BY i.isw_start
    };
    my $isw_query_sth = $dbh->prepare($isw_query);

    $isw_query_sth->execute($align_id);
    while ( my @row = $isw_query_sth->fetchrow_array ) {
        my ( $isw_start, $isw_end, $isw_pi ) = @row;
        my $isw_ftr = $ftr->new(
            -start => $isw_start,
            -stop  => $isw_end,
            -score => $isw_pi,
        );
        push @$isw_pi_seg, $isw_ftr;
        $isw_pi_max = $isw_pi > $isw_pi_max ? $isw_pi : $isw_pi_max;
        $isw_pi_min = $isw_pi < $isw_pi_min ? $isw_pi : $isw_pi_min;
    }
}

# gsw
my $sliding_seg = [];
my $crest_seg   = [];
my $trough_seg  = [];
{
    print "Prepare gsw\n";
    my @slidings = $obj->gc_wave( $align_id, $comparable_set );

    # write $gc_wave_csv file
    if ($gc_wave_csv) {
        open my $csv_fh, ">", "$db-align-$align_id.gc.csv";
        for (@slidings) {
            $_->{set}->runlist =~ /^(\d+)/;
            print {$csv_fh} $1, ",", $_->{gc}, ",", $_->{high_low_flag}, "\n";
        }
        close $csv_fh;
    }

    for (@slidings) {
        if (   $_->{high_low_flag} eq 'trough'
            or $_->{high_low_flag} eq 'crest' )
        {
            $_->{high_low_flag} =~ s/trough/T/i;
            $_->{high_low_flag} =~ s/crest/C/i;
        }

        my $sliding_ftr = $ftr->new(
            -segments => [ $_->{set}->spans ],
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
}

#----------------------------------------------------------#
# Draw the figure
#----------------------------------------------------------#
if ($figure) {
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
        -width       => $width,
        -pad_top     => 20,
        -pad_bottom  => 20,
        -pad_left    => 20,
        -pad_right   => 20,
        -key_style   => 'none',
        -image_class => $CLASS,
    );

    {    # text
        $panel->add_track(
            $align_segment,
            -glyph => 'text_in_box',
            -text  => "DB: $db | Align: $align_id | "
                . "$chr_name $chr_start-$chr_end",
            -text_bgcolor => 'lightcyan',
            -height       => 10,
            -bgcolor      => 'yellow',
            -text_pad     => 4,
        );
    }

    {    # text
        $panel->add_track(
            $align_segment,
            -glyph => 'text_in_box',
            -text  => "GC wave | size: $wave_window_size | step: $wave_window_step | "
                . "fall range: $fall_range",
            -text_bgcolor => 'lightcyan',
            -height       => 10,
            -bgcolor      => 'yellow',
            -text_pad     => 4,
        );
    }

    {    # arrow
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

    {    # alignment
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

        # coding
        my $coding_segment = $ftr->new(
            -segments => [ $coding_set->spans ],
            -name     => 'coding',
            -type     => 'alignment'
        );
        $panel->add_track(
            $coding_segment,
            -glyph     => 'segments',
            -label     => 'coding',
            -bump      => 0,
            -height    => 10,
            -font      => 'gdSmallFont',
            -linewidth => 1,
            -bgcolor   => 'turquoise',
            -connector => 'dashed',
        );

        # repeat
        my $repeat_segment = $ftr->new(
            -segments => [ $repeat_set->spans ],
            -name     => 'repeat',
            -type     => 'alignment'
        );
        $panel->add_track(
            $repeat_segment,
            -glyph     => 'segments',
            -label     => 'repeat',
            -bump      => 0,
            -height    => 10,
            -font      => 'gdSmallFont',
            -linewidth => 1,
            -bgcolor   => 'green',
            -fgcolor   => 'green',
            -connector => 'dashed',
        );
    }

    {    # indel
        my $indel_segment = $ftr->new(
            -segments => [ $indel_set->spans ],
            -name     => 'indel',
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

    {    # isw
        my $isw_segment = $ftr->new(
            -segments => $isw_pi_seg,
            -name     => 'isw',
        );
        $panel->add_track(
            $isw_segment,
            -glyph      => 'xyplot',
            -label      => 'isw_pi',
            -graph_type => 'boxes',
            -fgcolor    => 'orange',
            -bgcolor    => 'orange',
            -scale      => 'left',
            -height     => 80,
            -max_score  => $isw_pi_max,
            -min_score  => $isw_pi_min,
        );
    }

    {    # crest
        my $crest_segment = $ftr->new(
            -segments => $crest_seg,
            -name     => 'crest',
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

    {    # sliding, GC
        my $sliding_segment = $ftr->new(
            -segments => $sliding_seg,
            -name     => 'sliding',
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

    {    # trough
        my $trough_segment = $ftr->new(
            -segments => $trough_seg,
            -name     => 'trough',
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
    open my $pic_fh, ">", "$db-align-$align_id.$type";
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
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --db                database name
        --username          username
        --password          password
        --align_id          align_id

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
