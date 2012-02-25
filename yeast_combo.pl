#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use List::MoreUtils qw(zip);
use Math::Combinatorics;
use File::Copy::Recursive qw(rmove);
use Template;

use AlignDB::Stopwatch;

use FindBin;

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

# running options
my $k        = 2;                                   # which combination
my $phase    = 0;                                   # diffirent working phases
my $base_dir = '/home/wangq/Date/Alignment/yeast6/';

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    'k=i'          => \$k,
    'p|phase=s'    => \$phase,
    'b|base_dir=s' => \$base_dir,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
my $tt = Template->new;

my @names    = qw{S288C RM11 YJM789 DBVPG6765 SK1 Y55};
my @orders   = ( 0 .. 5 );
my %order_of = zip @names, @orders;

my $cmd_1 = qq{
    perl $FindBin::Bin/../multi/fasta_malignDB.pl \\
    --dir=[% dir %] \\
    --db=[% db %] \\
    --parallel=4
};
my $cmd_2 = qq{
    perl $FindBin::Bin/../multi/update_multi.pl \\
    --db=[% db %] \\
    --ensembl=yeast_50
};
my $cmd_3 = qq{
    perl $FindBin::Bin/../stat/mvar_stat_factory.pl \\
    --db=[% db %] \\
    --output=[% mvar_xls %]
};
my $cmd_4 = qq{
    perl $FindBin::Bin/../stat/multi_stat_factory.pl \\
    --db=[% db %] \\
    --output=[% multi_xls %]
};
my $cmd_5 = qq{
    perl $FindBin::Bin/../stat/indel_content.pl \\
    --db=[% db %] \\
    --min_length=1 \\
    --max_length=50 \\
    --output=[% content_xls %]
};

my $dispatch = {};
my $common   = {
    1 => $cmd_1,
    2 => $cmd_2,
    3 => $cmd_3,
    4 => $cmd_4,
    5 => $cmd_5,
};

for ( 2 .. 6 ) {
    $dispatch->{$_} = $common;
}

$dispatch->{2}{db} = q{[% combo.0 %]vs[% combo.1 %]refSpar};
$dispatch->{2}{0} = qq{
    perl ref_outgroup.pl \\
    --first_db=[% combo.0 %]vs[% combo.1 %] \\
    --second_db=[% combo.0 %]vsSpar \\
    --goal_db=[% db %] \\
    --target=0target --query=0query --outgroup=1query \\
    --no_insert=1 --trimmed_fasta=1 \\
    --chr_id_runlist=1-1000
};

$dispatch->{3}{db} = q{[% combo.0 %]vs[% combo.1 %]vs[% combo.2 %]refSpar};
$dispatch->{3}{0} = qq{
    perl ref_outgroup.pl \\
    --first_db=[% combo.0 %]vs[% combo.1 %] \\
    --second_db=[% combo.0 %]vs[% combo.2 %] \\
    --third_db=[% combo.0 %]vsSpar \\
    --goal_db=[% db %] \\
    --third=1query --outgroup=2query \\
    --no_insert=1 --trimmed_fasta=1 \\
    --chr_id_runlist=1-1000
};

$dispatch->{4}{db}
    = q{[% combo.0 %]vs[% combo.1 %]vs[% combo.2 %]vs[% combo.3 %]refSpar};
$dispatch->{4}{0} = qq{
    perl ref_outgroup.pl \\
    --first_db=[% combo.0 %]vs[% combo.1 %] \\
    --second_db=[% combo.0 %]vs[% combo.2 %] \\
    --third_db=[% combo.0 %]vsSpar \\
    --forth_db=[% combo.0 %]vs[% combo.3 %] \\
    --goal_db=[% db %] \\
    --third=1query --outgroup=2query \\
    --no_insert=1 --trimmed_fasta=1 \\
    --chr_id_runlist=1-1000
};

$dispatch->{5}{db}
    = q{[% combo.0 %]vs[% combo.1 %]vs[% combo.2 %]vs[% combo.3 %]vs[% combo.4 %]refSpar};
$dispatch->{5}{0} = qq{
    perl ref_outgroup.pl \\
    --first_db=[% combo.0 %]vs[% combo.1 %] \\
    --second_db=[% combo.0 %]vs[% combo.2 %] \\
    --third_db=[% combo.0 %]vsSpar \\
    --forth_db=[% combo.0 %]vs[% combo.3 %] \\
    --fifth_db=[% combo.0 %]vs[% combo.4 %] \\
    --goal_db=[% db %] \\
    --third=1query --outgroup=2query \\
    --no_insert=1 --trimmed_fasta=1 \\
    --chr_id_runlist=1-1000
};

$dispatch->{6}{db}
    = q{[% combo.0 %]vs[% combo.1 %]vs[% combo.2 %]vs[% combo.3 %]vs[% combo.4 %]vs[% combo.5 %]refSpar};
$dispatch->{6}{0} = qq{
    perl ref_outgroup.pl \\
    --first_db=[% combo.0 %]vs[% combo.1 %] \\
    --second_db=[% combo.0 %]vs[% combo.2 %] \\
    --third_db=[% combo.0 %]vsSpar \\
    --forth_db=[% combo.0 %]vs[% combo.3 %] \\
    --fifth_db=[% combo.0 %]vs[% combo.4 %] \\
    --sixth_db=[% combo.0 %]vs[% combo.5 %] \\
    --goal_db=[% db %] \\
    --third=1query --outgroup=2query \\
    --no_insert=1 --trimmed_fasta=1 \\
    --chr_id_runlist=1-1000
};

#----------------------------#
# Use the AlignDB::Multi way to generate multi-sequence fasta files
#----------------------------#
my $combinat = Math::Combinatorics->new(
    count => $k,
    data  => [@names],
);

$stopwatch->block_message(
    "Combinations of $k from: " . join( " ", @names ) );

while ( my @combo = $combinat->next_combination ) {
    @combo = sort { $order_of{$a} <=> $order_of{$b} } @combo;
    print "# ", join( ',', @combo ), "\n";

    my ( $goal_db, $goal_dir, $mvar_xls, $multi_xls, $content_xls, $cmd );

    # use the dispatch template to generate $goal_db
    $tt->process( \$dispatch->{$k}{db}, { combo => \@combo }, \$goal_db )
        or die Template->error;

    $goal_dir    = "$base_dir$k/$goal_db";
    $mvar_xls    = "$goal_db.mvar.xls";
    $multi_xls   = "$goal_db.multi.xls";
    $content_xls = "$goal_db.indel.length1-50.xls";

    # use the dispatch template to generate $cmd
    $tt->process(
        \$dispatch->{$k}{$phase},
        {   combo       => \@combo,
            db          => $goal_db,
            dir         => $goal_dir,
            mvar_xls    => $mvar_xls,
            multi_xls   => $multi_xls,
            content_xls => $content_xls,
        },
        \$cmd
    ) or die Template->error;

    print $cmd, "\n";
    system $cmd;

    if ( $phase == 0 ) {
        rmove( $goal_db, $goal_dir );
    }
    elsif ( $phase == 3 ) {
        rmove( $mvar_xls, "$k/mvar/$goal_db.mvar.xls" );
    }
    elsif ( $phase == 4 ) {
        rmove( $multi_xls, "$k/multi/$goal_db.multi.xls" );
    }
    elsif ( $phase == 5 ) {
        rmove( $content_xls, "$k/content/$goal_db.indel.length1-50.xls" );
    }
}
print "\n";

$stopwatch->end_message;
exit;

__END__

perl yeast_combo.pl -k=6 -p=0
