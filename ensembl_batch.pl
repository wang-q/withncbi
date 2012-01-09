#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Find::Rule;
use File::Spec;
use Path::Class;
use File::Basename;
use String::Compare;
use Set::Light;

use Template;

use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";

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

# write_axt parameter
my $yml_file = "$FindBin::Bin/fungi.yml";

my $run;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    'server=s'     => \$server,
    'port=i'       => \$port,
    'username=s'   => \$username,
    'password=s'   => \$password,
    'y|yml_file=s' => \$yml_file,
    'r|run'        => \$run,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Ensembl batching for $yml_file");

my $dispatch = LoadFile($yml_file);

my $mysql_dir = dir( $dispatch->{dir}{mysql} );
my $fasta_dir = dir( $dispatch->{dir}{fasta} );

my $species_ref = $dispatch->{species};

my $version     = $dispatch->{meta}{version};
my $initrc_file = $dispatch->{meta}{initrc_file};

#----------------------------------------------------------#
# Write .axt files from alignDB
#----------------------------------------------------------#
my $initrc_str;
for my $sp ( sort keys %{$species_ref} ) {
    print "$sp\n";

    my (@aliases);
    {
        my ( $genus, $species, $strain ) = split " ", $sp;

        if ($species) {
            my $str = lc( substr( $genus, 0, 1 ) . substr( $species, 0, 3 ) );
            push @aliases, $str;
            push @aliases, $genus . '_' . $species;
            push @aliases, substr( $genus, 0, 1 ) . $species;
            push @aliases, substr( $genus, 0, 1 ) . '_' . $species;
        }
        else {
            my $str = lc $genus;
            push @aliases, $str;
            push @aliases, $genus;
        }
    }

    #print Dump { dirs => \@dirs, db_base => $db_base, aliases => \@aliases};

    # build databases
    for my $group (qw{core funcgen otherfeatures variation compara}) {
        if ( $species_ref->{$sp}{$group} ) {
            my $ensembl_db = build_db( $sp, $group, $mysql_dir, $run );
            if ( $species_ref->{$sp}{initrc} ) {
                my $aliases_ref = [@aliases];
                push @{$aliases_ref}, $aliases[0] . "_$group" . "_$version";
                push @{$aliases_ref}, $aliases[0] . "_$version"
                    if $group eq 'core';
                if ( ref $species_ref->{$sp}{$group} eq 'HASH'
                    and exists $species_ref->{$sp}{$group}{aliases} )
                {
                    push @{$aliases_ref},
                        @{ $species_ref->{$sp}{$group}{aliases} };
                }
                $initrc_str
                    .= build_initrc( $sp, $group, $ensembl_db, $aliases_ref );
            }
        }
    }
}

{
    open my $out_fh, ">", $initrc_file;
    print {$out_fh} $initrc_str;
    close $out_fh;
}

$stopwatch->end_message;

sub build_db {
    my $sp        = shift;
    my $group     = shift;
    my $mysql_dir = shift;
    my $run       = shift;

    my %score_of;
    my ( $genus, $species, $strain ) = split " ", $sp;
    my $set = Set::Light->new(qw{ core funcgen otherfeatures variation });
    my $str;
    if ( $set->has($group) ) {
        $str = lc join "_", grep {$_} ( $genus, $species, $strain, $group );
    }
    else {
        $str = lc join "_", ( 'ensembl', $group, $genus  );
    }

    while ( my $child = $mysql_dir->next ) {
        next unless -d $child;
        my $dir      = $child->stringify;
        my $base_dir = basename($dir);
        my $score    = String::Compare::char_by_char( $base_dir, $str );
        $score_of{$dir} = $score;
    }

    my ($ensembl_dir)
        = sort { $score_of{$b} <=> $score_of{$a} } keys %score_of;
    my $ensembl_db = basename($ensembl_dir);
    print "Find $ensembl_db for $str\n";

    my $cmd
        = "perl $FindBin::Bin/build_ensembl.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " --initdb"
        . " --db=$ensembl_db"
        . " --ensembl=$ensembl_dir";

    $stopwatch->block_message("Build ensembl for $ensembl_db");
    $stopwatch->block_message($cmd);
    system $cmd if $run;
    $stopwatch->block_message("Finish build");

    return $ensembl_db;
}

sub build_initrc {
    my $sp          = shift;
    my $group       = shift;
    my $db          = shift;
    my $aliases_ref = shift;

    my $text = <<'EOF';
{    # [% sp %]
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host    => $host,
        -user    => $user,
        -pass    => $pass,
        -port    => $port,
        -species => '[% sp %]',
        -group   => '[% group %]',
        -dbname  => '[% db %]',
    );

    my @aliases = (
[% FOREACH item IN aliases -%]
        '[% item %]',
[% END -%]
    );

    Bio::EnsEMBL::Utils::ConfigRegistry->add_alias(
        -species => '[% sp %]',
        -alias   => \@aliases,
    );
}

EOF

    my $tt = Template->new;
    my $output;
    $tt->process(
        \$text,
        {   sp      => $sp,
            group   => $group,
            db      => $db,
            aliases => $aliases_ref,
        },
        \$output
    ) or die Template->error;

    return $output;
}

__END__

=head1 NAME

    write_axt_slice.pl - extract alignment slices from alignDB

=head1 SYNOPSIS

    write_axt_slice.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --db                database name
        --username          username
        --password          password
        -y, --yaml_dir      dir of yaml

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

=cut

