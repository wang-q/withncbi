#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
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

# Database init values
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db       = $Config->{database}{db};

my $goal_db;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'server=s'   => \$server,
    'port=i'     => \$port,
    'db=s'       => \$db,
    'username=s' => \$username,
    'password=s' => \$password,
    'goal=s'     => \$goal_db,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

unless ($goal_db) {
    $goal_db = $db . "_dup";
}

#----------------------------------------------------------#
# call mysql
#----------------------------------------------------------#
$stopwatch->start_message("Operation start...");

my $str       = " -h$server -P$port -u$username -p$password ";
my $drop      = "mysql $str -e \"DROP DATABASE IF EXISTS $goal_db;\"";
my $create    = "mysql $str -e \"CREATE DATABASE $goal_db;\"";
my $duplicate = "mysqldump $str $db | mysql $str $goal_db";

print "#drop\n$drop\n";
system($drop);
print "#create\n$create\n";
system($create );
print "#duplicate\n$duplicate\n";
system($duplicate);

$stopwatch->end_message;

## store program running meta info to database
#END {
#    $obj->add_meta_stopwatch($stopwatch);
#}

__END__


=head1 NAME

    init_alignDB.pl - Initiate alignDB

=head1 SYNOPSIS

    init_alignDB.pl [options]
      Options:
        --help          brief help message
        --man           full documentation
        --server        MySQL server IP/Domain name
        --port          MySQL server port
        --db            database name
        --username      username
        --password      password
        --init_sql      init sql filename
      
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
