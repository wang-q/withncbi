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
    'help|?'       => \$help,
    'man'          => \$man,
    's|server=s'   => \$server,
    'P|port=i'     => \$port,
    'd|db=s'       => \$db,
    'u|username=s' => \$username,
    'p|password=s' => \$password,
    'f|file=s'     => \$file_dump,
    'g|goal=s'     => \$goal_db,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# call mysql
#----------------------------------------------------------#
$stopwatch->start_message("Operation start...");

my $str = " -h$server -P$port -u$username -p$password ";

# dump db
if ($file_dump) {
    if ( !-e $file_dump ) {
        my $dump = "mysqldump $str $db > $file_dump";
        print "#dump\n$dump\n\n";
        system($dump);
    }
}
else {
    $file_dump = $db . ".dump.sql";
    my $dump = "mysqldump $str $db > $file_dump";
    print "#dump\n$dump\n\n";
    system($dump);
}

# load dump
if ($goal_db) {
    my $drop      = "mysql $str -e \"DROP DATABASE IF EXISTS $goal_db;\"";
    my $create    = "mysql $str -e \"CREATE DATABASE $goal_db;\"";
    my $duplicate = "mysql $str $goal_db < $file_dump";

    print "#drop\n$drop\n\n";
    system($drop);
    print "#create\n$create\n\n";
    system($create );
    print "#duplicate\n$duplicate\n\n";
    system($duplicate);
}

$stopwatch->end_message;

__END__

=head1 NAME

    dup_db.pl - Duplicate a database

=head1 SYNOPSIS

    dup_db.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --port              MySQL server port
        --db                database name
        --username          username
        --password          password
        --file              dump file name
        --goal              dup db name

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
