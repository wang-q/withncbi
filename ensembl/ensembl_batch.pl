#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use Path::Tiny;
use String::Compare;
use Set::Light;
use Template;

use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

=head1 NAME

ensembl_batch.pl - Ensembl batching

=head1 SYNOPSIS

    perl ensembl_batch.pl [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port      -P  INT     MySQL server port
        --username  -u  STR     username
        --password  -p  STR     password
        --file      -i  STR     YAML file

    perl ensembl_batch.pl -i ensembl_82.yml

=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'server|s=s'   => \( my $server   = $Config->{database}{server} ),
    'port|P=i'     => \( my $port     = $Config->{database}{port} ),
    'username|u=s' => \( my $username = $Config->{database}{username} ),
    'password|p=s' => \( my $password = $Config->{database}{password} ),
    'file|i=s'     => \( my $yml_file = "$FindBin::Bin/ensembl_82.yml" ),
) or HelpMessage(1);

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
$stopwatch->start_message("Ensembl batching for $yml_file");

my $dispatch = LoadFile($yml_file);

my $mysql_dir = path( $dispatch->{dir}{mysql} )->absolute->stringify;
my $fasta_dir = path( $dispatch->{dir}{fasta} )->absolute->stringify;

my $species_ref = $dispatch->{species};

my $version     = $dispatch->{meta}{version};
my $initrc_file = $dispatch->{meta}{initrc_file};
my $sh_file     = $dispatch->{meta}{sh_file};

#----------------------------------------------------------#
# Write .axt files from alignDB
#----------------------------------------------------------#
my $initrc_str = q{
use strict;
use warnings;
use autodie;

use Bio::EnsEMBL::Utils::ConfigRegistry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

my $host = 'localhost';
my $port = 3306;
my $user = 'alignDB';
my $pass = 'alignDB';

};
my $sh_str;
for my $sp ( sort keys %{$species_ref} ) {
    print "$sp\n";
    if ( $species_ref->{$sp}{skip} ) {
        print " " x 4, "Skip $sp\n";
        next;
    }

    my (@aliases);
    {
        my ( $genus, $species, $strain ) = split " ", $sp;

        if ($species) {
            my $str = lc( substr( $genus, 0, 1 ) . substr( $species, 0, 3 ) );
            push @aliases, $str;
            push @aliases, ucfirst $str;
            push @aliases, $genus . '_' . $species;
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
            my ( $ensembl_db, $cmd ) = build_cmd( $sp, $group, $mysql_dir );
            $sh_str .= $cmd;
            if ( $species_ref->{$sp}{initrc} ) {
                # Don't affect other groups' aliases
                my $aliases_ref = [@aliases];
                push @{$aliases_ref}, $aliases[0] . "_$group" . "_$version";
                push @{$aliases_ref}, $aliases[0] . "_$version"
                    if $group eq 'core';
                push @{$aliases_ref}, $ensembl_db;

                if ( ref $species_ref->{$sp}{$group} eq 'HASH'
                    and exists $species_ref->{$sp}{$group}{aliases} )
                {
                    push @{$aliases_ref}, @{ $species_ref->{$sp}{$group}{aliases} };
                }
                $initrc_str .= build_initrc( $sp, $group, $ensembl_db, $aliases_ref );
            }
        }
    }
}

$initrc_str .= q{
1;
};

#----------------------------#
# Write outfiles
#----------------------------#
{
    path($initrc_file)->spew($initrc_str);
    path($sh_file)->spew($sh_str);
}

$stopwatch->end_message;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#
sub build_cmd {
    my $sp        = shift;
    my $group     = shift;
    my $mysql_dir = shift;

    my %score_of;
    my ( $genus, $species, $strain ) = split " ", $sp;
    my $set = Set::Light->new(qw{ core funcgen otherfeatures variation });
    my $str;
    if ( $set->has($group) ) {
        $str = lc join "_", grep {$_} ( $genus, $species, $strain, $group );
    }
    else {
        $str = lc join "_", ( 'ensembl', $group, $genus );
    }

    my $iter = path($mysql_dir)->iterator( { recurse => 0, } );
    while ( my $child = $iter->() ) {
        next unless $child->is_dir;
        my $base_dir = $child->basename;
        my $score = String::Compare::char_by_char( $base_dir, $str );
        $score_of{ $child->stringify } = $score;
    }

    my ($ensembl_dir)
        = sort { $score_of{$b} <=> $score_of{$a} } keys %score_of;
    my $ensembl_db = path($ensembl_dir)->basename;
    print " " x 4, "Find $ensembl_db for $str\n";
    if ( index( $ensembl_db, $str ) != 0 ) {
        print " " x 4, "Be careful for [$str]\n";
    }

    my $cmd
        = "perl $FindBin::Bin/build_ensembl.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " --initdb"
        . " --db $ensembl_db"
        . " --ensembl $ensembl_dir" . "\n\n";

    return ( $ensembl_db, $cmd );
}

sub build_initrc {
    my $sp          = shift;
    my $group       = shift;
    my $db          = shift;
    my $aliases_ref = shift;

    my $text = <<'EOF';
# [% sp %]
{
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
