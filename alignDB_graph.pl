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

use Number::Format qw(:subs);

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::Ensembl;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server   = $Config->{database}->{server};
my $port     = $Config->{database}->{port};
my $username = $Config->{database}->{username};
my $password = $Config->{database}->{password};
my $db       = $Config->{database}->{db};

# graph init values
my $CLASS = $Config->{graph}->{class};    # GD or GD::SVG
my $width = $Config->{graph}->{width};

# separate from the one of database section
# default if undef
my $ensembl_db = $Config->{graph}->{ensembl};
my $graph_file = "";

# coverage on target or query
my $goal = "target";

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'server=s'   => \$server,
    'port=s'     => \$port,
    'db=s'       => \$db,
    'username=s' => \$username,
    'password=s' => \$password,
    'ensembl=s'  => \$ensembl_db,
    'output=s'   => \$graph_file,
    'class=s'    => \$CLASS,
    'width=i'    => \$width,
    'goal=s'     => \$goal,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Drawing coverage on $goal for $db...");

# alignDB object
my $alignDB_obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);

my $dbh = $alignDB_obj->dbh;

# ensembl object
my $kary_adaptor;
if ( length $ensembl_db > 0 ) {
    my $ensembl = AlignDB::Ensembl->new(
        server => $server,
        db     => $ensembl_db,
        user   => $username,
        passwd => $password,
    );
    my $db_adaptor = $ensembl->get_db_adaptor;
    $kary_adaptor = $db_adaptor->get_KaryotypeBandAdaptor;
}

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#
# get target and query names via AlignDB methods
my ( $target_name, $query_name ) = $alignDB_obj->get_names;

# select all $goal chromosomes in this database
my @chrs = @{ $alignDB_obj->get_chrs($goal) };

# select all _GOAL_ sequence in one chromosome
my $goal_query = qq{
    SELECT s.chr_start,
           s.chr_end
    FROM _GOAL_ G, sequence s
    WHERE s.chr_id = ?
    AND G.seq_id = s.seq_id
};
$goal_query =~ s/_GOAL_/$goal/;
my $goal_query_sth = $dbh->prepare($goal_query);

my @chr_infos;
my $largest = 0;

# for each chromosome
for my $row (@chrs) {
    my ( $chr_id, $chr_name, $chr_length ) = @$row;

    my $chr_set      = AlignDB::IntSpan->new("1-$chr_length");
    my $goal_chr_set = AlignDB::IntSpan->new;
    my $overlap_set  = AlignDB::IntSpan->new;

    # for each target sequence
    $goal_query_sth->execute($chr_id);
    while ( my @row2 = $goal_query_sth->fetchrow_array ) {
        my ( $chr_start, $chr_end ) = @row2;
        my $runlist = "$chr_start-$chr_end";

        my $intersect = $goal_chr_set->intersect($runlist);
        if ( $intersect->is_not_empty ) {
            print "Found overlap at $chr_name:$intersect\n";
            $overlap_set->add($intersect);
        }
        $goal_chr_set->add($runlist);
    }

    my $covered_length = $goal_chr_set->cardinality;
    my $coverage       = $covered_length / $chr_length * 100;
    $coverage = sprintf "%.2f", $coverage;

    if ( $chr_length > $largest ) {
        $largest = $chr_length;
    }

    my $chr_name2 = $chr_name;
    $chr_name2 =~ s/chr0?//;
    my $band;
    if ( defined $kary_adaptor and length $chr_name2 < 5 ) {
        $band = $kary_adaptor->fetch_all_by_chr_name($chr_name2);
    }

    my $info = {
        chr_name     => $chr_name,
        chr_id       => $chr_id,
        chr_length   => $chr_length,
        chr_set      => $chr_set,
        goal_chr_set => $goal_chr_set,
        overlap_set  => $overlap_set,
        coverage     => $coverage,
        band         => $band,
    };

    push @chr_infos, $info;
}

#----------------------------#
# write a chr yaml
#----------------------------#
my $full_chr_set        = {};
my $full_overlap_set    = {};
my $full_nonoverlap_set = {};
for my $chr_info (@chr_infos) {
    my $chr_name     = $chr_info->{chr_name};
    my $goal_chr_set = $chr_info->{goal_chr_set};
    my $overlap_set  = $chr_info->{overlap_set};

    $full_chr_set->{$chr_name}        = $goal_chr_set;
    $full_overlap_set->{$chr_name}    = $overlap_set;
    $full_nonoverlap_set->{$chr_name} = $goal_chr_set->diff($overlap_set);
}
DumpFile( "$db.$goal.yml",            $full_chr_set );
DumpFile( "$db.$goal.overlap.yml",    $full_overlap_set );
DumpFile( "$db.$goal.nonoverlap.yml", $full_nonoverlap_set );

#----------------------------#
# draw pic
#----------------------------#
print "Draw picture...\n";
my $ftr = 'Bio::Graphics::Feature';
{
    my $largest_chr = $ftr->new(
        -start => 1,
        -end   => $largest,
        -name  => $db,
        -type  => 'alignment'
    );
    my $panel = Bio::Graphics::Panel->new(
        -grid        => 1,
        -gridcolor   => 'lightcyan',
        -segment     => $largest_chr,
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
        my $title
            = "DB: $db"
            . " | Target: $target_name"
            . " | Query: $query_name"
            . " | On $goal";
        $panel->add_track(
            $largest_chr,
            -glyph        => 'text_in_box',
            -text         => $title,
            -text_bgcolor => 'lightcyan',
            -height       => 10,
            -bgcolor      => 'yellow',
            -text_pad     => 4,
            -font         => 'gdMediumBoldFont',
        );
    }

    {    # arrow
        $panel->add_track(
            $largest_chr,
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

    foreach my $chr_info (@chr_infos) {
        my $chr_name     = $chr_info->{chr_name};
        my $chr_set      = $chr_info->{chr_set};
        my $chr_length   = $chr_info->{chr_length};
        my $goal_chr_set = $chr_info->{goal_chr_set};
        my $coverage     = $chr_info->{coverage};
        my $band         = $chr_info->{band};

        my $chr_segment = $ftr->new(
            -segments => [ $chr_set->spans ],
            -name     => $chr_name,
            -type     => 'alignment'
        );

        {    # text
            my $title
                = "$chr_name: "
                . format_bytes($chr_length) . " bp"
                . " | Coverage: $coverage%";
            $panel->add_track(
                $chr_segment,
                -glyph        => 'text_in_box',
                -text         => $title,
                -text_bgcolor => 'lightcyan',
                -height       => 10,
                -bgcolor      => 'yellow',
                -text_pad     => 4,
            );
        }

        {    # alignment
            $panel->add_track(
                $chr_segment,
                -glyph     => 'segments',
                -label     => $chr_name,
                -bump      => 0,
                -height    => 10,
                -font      => 'gdSmallFont',
                -linewidth => 1,
                -bgcolor   => 'blue',
                -connector => 'solid',
            );

            my $target_segment = $ftr->new(
                -segments => [ $goal_chr_set->spans ],
                -name     => 'coverage',
                -type     => 'alignment'
            );
            $panel->add_track(
                $target_segment,
                -glyph     => 'segments',
                -label     => 'coverage',
                -bump      => 0,
                -height    => 10,
                -font      => 'gdSmallFont',
                -linewidth => 1,
                -bgcolor   => 'lightgreen',
                -fgcolor   => 'lightgreen',
                -connector => 'solid',
            );
        }

        # bands
        if ( defined $band ) {
            my $band_seg   = {};
            my $band_color = {
                acen    => 'green',
                gneg    => 'cyan',
                gpos100 => 'gray',
                gpos75  => 'gray',
                gpos50  => 'gray',
                gpos25  => 'gray',
                gvar    => 'purple',
                stalk   => 'green',
            };
            foreach (@$band) {
                my $source   = $_->stain;
                my $band_ftr = Bio::Graphics::Feature->new(
                    -start => $_->start,
                    -stop  => $_->end,
                    -name  => $_->name,
                    -type  => 'band',

                    #-source => $_->stain,
                );
                if ( exists $band_seg->{$source} ) {
                    push @{ $band_seg->{$source} }, $band_ftr;

                }
                else {
                    $band_seg->{$source} = [$band_ftr];
                }
            }

            foreach my $source ( sort keys %$band_seg ) {
                my $band_segment = $ftr->new(
                    -segments => $band_seg->{$source},
                    -name     => $source,
                );

                $panel->add_track(
                    $band_segment,
                    -glyph   => 'segments',
                    -label   => $source,
                    -fgcolor => $band_color->{$source},
                    -bgcolor => 'white',
                    -height  => 10,
                );
            }
        }

        # line
        $panel->add_track(
            $largest_chr,
            -glyph     => 'line',
            -bump      => 0,
            -height    => 20,
            -linewidth => 1,
            -bgcolor   => 'turquoise',
            -fgcolor   => 'turquoise',
        );
    }

    my $gd = $panel->gd;
    my $type = ( $CLASS eq 'GD' ) ? 'png' : 'svg';
    $graph_file = ">$db.$goal.$type" unless $graph_file;
    open PICFH, ">$graph_file";
    binmode PICFH;
    print PICFH $gd->$type;
    close PICFH;
}

$stopwatch->end_message;
exit;

__END__


=head1 NAME

    alignDB_graph.pl - Generate graph for chromosome coverage in alignDB

=head1 SYNOPSIS

    alignDB_graph.pl [options]
        Options:
            --help              brief help message
            --man               full documentation
            --server            MySQL server IP/Domain name
            --db                database name
            --username          username
            --password          password
            --class             GD or GD::SVG
            --width             width + 40 = figure width

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
