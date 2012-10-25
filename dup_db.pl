#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

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

# Database init values
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db       = $Config->{database}{db};

my $goal_db;
my $file_dump;

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
    'file=s'     => \$file_dump,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

unless ($goal_db) {
    $goal_db = $db . "_dup";
}
unless ($file_dump) {
    $file_dump = $db . ".dump.sql";
}

#----------------------------------------------------------#
# call mysql
#----------------------------------------------------------#
$stopwatch->start_message("Operation start...");

my $str       = " -h$server -P$port -u$username -p$password ";
my $drop      = "mysql $str -e \"DROP DATABASE IF EXISTS $goal_db;\"";
my $create    = "mysql $str -e \"CREATE DATABASE $goal_db;\"";
my $dump      = "mysqldump $str $db > $file_dump";
my $duplicate = "mysql $str $goal_db < $file_dump";

print "#drop\n$drop\n\n";
system($drop);
print "#create\n$create\n\n";
system($create );
print "#dump\n$dump\n\n";
system($dump);
print "#duplicate\n$duplicate\n\n";
system($duplicate);

$stopwatch->end_message;

__END__


=head1 NAME

    dup_db.pl - Duplicate a database

=head1 SYNOPSIS

    dup_db.pl [options]
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
