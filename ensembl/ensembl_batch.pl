#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw();
use Config::Tiny;
use FindBin;
use YAML::Syck qw();

use Path::Tiny qw();
use String::Compare qw();
use Set::Scalar;
use Template;

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

ensembl_batch.pl - Ensembl batching

=head1 SYNOPSIS

    perl ensembl_batch.pl [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port          INT     MySQL server port
        --username  -u  STR     username
        --password  -p  STR     password
        --file      -i  STR     YAML file

    perl ensembl_batch.pl -i ensembl_82.yml

=cut

Getopt::Long::GetOptions(
    'help|?'       => sub { Getopt::Long::HelpMessage(0) },
    'server|s=s'   => \( my $server = $Config->{database}{server} ),
    'port=i'       => \( my $port = $Config->{database}{port} ),
    'username|u=s' => \( my $username = $Config->{database}{username} ),
    'password|p=s' => \( my $password = $Config->{database}{password} ),
    'file|i=s'     => \( my $yml_file = "$FindBin::Bin/ensembl_98.yml" ),
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
$stopwatch->start_message("Ensembl batching for $yml_file");

my $dispatch = YAML::Syck::LoadFile($yml_file);

my $mysql_dir = Path::Tiny::path( $dispatch->{dir}{mysql} )->absolute->stringify;
my $fasta_dir = Path::Tiny::path( $dispatch->{dir}{fasta} )->absolute->stringify;
my $gff3_dir  = Path::Tiny::path( $dispatch->{dir}{gff3} )->absolute->stringify;
my $dest_dir  = Path::Tiny::path( $dispatch->{dir}{dest} )->absolute->stringify;

my $species_ref = $dispatch->{species};

my $version = $dispatch->{meta}{version};

#----------------------------------------------------------#
# Write .axt files from alignDB
#----------------------------------------------------------#
my $initrc_str = qq{
use strict;
use warnings;
use autodie;

use Bio::EnsEMBL::Utils::ConfigRegistry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

my \$host = '$server';
my \$port = $port;
my \$user = '$username';
my \$pass = '$password';

};

my $build_str = q{#!/usr/bin/env bash

};
my $fasta_str = q{#!/usr/bin/env bash

};
my $anno_str = q{#!/usr/bin/env bash

};

for my $sp ( sort keys %{$species_ref} ) {
    print "$sp\n";
    if ( $species_ref->{$sp}{skip} ) {
        print " " x 4, "Skip $sp\n";
        next;
    }
    my $info = $species_ref->{$sp};

    my (@aliases);
    my $core_dbname;
    {
        my ( $genus, $species, undef ) = split " ", $sp;

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
        if ( $info->{$group} ) {
            my ( $ensembl_base, $cmd ) = build_cmd( $sp, $group );
            $core_dbname = $ensembl_base if $group eq "core";
            $build_str .= $cmd;
            if ( $info->{initrc} ) {

                # Don't affect other groups' aliases
                my $aliases_ref = [@aliases];
                push @{$aliases_ref}, $aliases[0] . "_$group" . "_$version";
                push @{$aliases_ref}, $aliases[0] . "_$version"
                    if $group eq 'core';
                push @{$aliases_ref}, $ensembl_base;

                if ( ref $info->{$group} eq 'HASH' and exists $info->{$group}{aliases} ) {
                    push @{$aliases_ref}, @{ $info->{$group}{aliases} };
                }
                $initrc_str .= build_initrc( $sp, $group, $ensembl_base, $aliases_ref );
            }
        }
    }

    # build fasta
    if ( $info->{fasta} ) {
        my $alias   = $aliases[1];
        my $pattern = "*dna_sm.toplevel*";
        my $append;
        if ( ref $info->{fasta} eq 'HASH' ) {
            if ( exists $info->{fasta}{alias} ) {
                $alias = $info->{fasta}{alias};
            }
            if ( exists $info->{fasta}{pattern} ) {
                $pattern = $info->{fasta}{pattern};
            }
            if ( exists $info->{fasta}{append} ) {
                $append = $info->{fasta}{append} . "\n";
            }
        }
        $fasta_str .= build_fasta( $sp, $alias, $pattern, $append );

        # build anno
        $anno_str .= build_anno( $sp, $alias );
    }
}

$initrc_str .= q{
1;
};

#----------------------------#
# Write outfiles
#----------------------------#
{
    Path::Tiny::path( $dispatch->{meta}{initrc_file} )->spew($initrc_str);
    Path::Tiny::path( $dispatch->{meta}{build_file} )->spew($build_str);
    Path::Tiny::path( $dispatch->{meta}{fasta_file} )->spew($fasta_str);
    Path::Tiny::path( $dispatch->{meta}{anno_file} )->spew($anno_str);
}

$stopwatch->end_message;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#
sub build_fasta {
    my $sp      = shift;
    my $alias   = shift;
    my $pattern = shift;
    my $append  = shift;

    my %score_of;
    my ( $genus, $species, $strain ) = split " ", $sp;
    my $str = lc join "_", grep {$_} ( $genus, $species, $strain );

    my $iter = Path::Tiny::path($fasta_dir)->iterator( { recurse => 0, } );
    while ( my Path::Tiny $child = $iter->() ) {
        next unless $child->is_dir;
        my $base_dir = $child->basename;
        my $score = String::Compare::char_by_char( $base_dir, $str );
        $score_of{ $child->stringify } = $score;
    }

    my ($ensembl_dir)
        = sort { $score_of{$b} <=> $score_of{$a} } keys %score_of;
    my $ensembl_base = Path::Tiny::path($ensembl_dir)->basename;
    printf " " x 4 . "fasta: [%s] for [%s]\n", $ensembl_base, $str;
    if ( index( $ensembl_base, $str ) != 0 ) {
        print " " x 4, "*** Be careful for [$str]\n";
    }

    if ($append) {
        $append = join "\n", map { " " x 4 . $_ } split( "\n", $append );
    }

    my $text = <<'EOF';
# [% sp %]
if [ ! -d [% dest %]/[% alias %] ]; then
    echo "==> [% sp %]"

    mkdir -p [% dest %]/[% alias %]
    cd [% dest %]/[% alias %]

    find [% dir %]/dna/ -name "[% pattern %]" |
        xargs gzip -d -c > toplevel.fa

    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

[% append %]

else
    echo "==> [% dest %]/[% alias %] exists"
fi

EOF

    my $tt = Template->new;
    my $output;
    $tt->process(
        \$text,
        {   sp      => $sp,
            dir     => $ensembl_dir,
            dest    => $dest_dir,
            alias   => $alias,
            pattern => $pattern,
            append  => $append,
        },
        \$output
    ) or die Template->error;

    return $output;
}

sub build_anno {
    my $sp      = shift;
    my $alias   = shift;

    my %score_of;
    my ( $genus, $species, $strain ) = split " ", $sp;
    my $str = lc join "_", grep {$_} ( $genus, $species, $strain );

    my $iter = Path::Tiny::path($gff3_dir)->iterator( { recurse => 0, } );
    while ( my Path::Tiny $child = $iter->() ) {
        next unless $child->is_dir;
        my $base_dir = $child->basename;
        my $score = String::Compare::char_by_char( $base_dir, $str );
        $score_of{ $child->stringify } = $score;
    }

    my ($ensembl_dir)
        = sort { $score_of{$b} <=> $score_of{$a} } keys %score_of;
    my $ensembl_base = Path::Tiny::path($ensembl_dir)->basename;
    printf " " x 4 . "gff3: [%s] for [%s]\n", $ensembl_base, $str;
    if ( index( $ensembl_base, $str ) != 0 ) {
        print " " x 4, "*** Be careful for [$str]\n";
    }

    my $text = <<'EOF';
# [% sp %]
if [ -d [% dest %]/[% alias %] ]; then
    echo "==> [% sp %]"

    if [ -f [% dest %]/[% alias %]/anno.yml ]; then
        echo "==> [% dest %]/[% alias %]/anno.yml exists"
    else
        cd [% dest %]/[% alias %]

        find [% dir %]/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        spanr gff chr.gff --tag CDS -o cds.yml

        faops masked *.fa |
            spanr cover stdin -o repeat.yml

        spanr merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> [% dest %]/[% alias %] does not exist"
fi

EOF

    my $tt = Template->new;
    my $output;
    $tt->process(
        \$text,
        {   sp      => $sp,
            dir     => $ensembl_dir,
            dest    => $dest_dir,
            alias   => $alias,
        },
        \$output
    ) or die Template->error;

    return $output;
}

sub build_cmd {
    my $sp    = shift;
    my $group = shift;

    my %score_of;
    my ( $genus, $species, $strain ) = split " ", $sp;
    my $set = Set::Scalar->new(qw{ core funcgen otherfeatures variation });
    my $str;
    if ( $set->has($group) ) {
        $str = lc join "_", grep {$_} ( $genus, $species, $strain, $group );
    }
    else {
        $str = lc join "_", ( 'ensembl', $group, $genus );
    }

    my $iter = Path::Tiny::path($mysql_dir)->iterator( { recurse => 0, } );
    while ( my Path::Tiny $child = $iter->() ) {
        next unless $child->is_dir;
        my $base_dir = $child->basename;
        my $score = String::Compare::char_by_char( $base_dir, $str );
        $score_of{ $child->stringify } = $score;
    }

    my ($ensembl_dir)
        = sort { $score_of{$b} <=> $score_of{$a} } keys %score_of;
    my $ensembl_base = Path::Tiny::path($ensembl_dir)->basename;
    printf " " x 4 . "db: [%s] for [%s]\n", $ensembl_base, $str;
    if ( index( $ensembl_base, $str ) != 0 ) {
        print " " x 4, "*** Be careful for [$str]\n";
    }

    my $cmd
        = "# $sp\n"
        . "perl $FindBin::Bin/build_ensembl.pl" . " \\\n"
        . " " x 4
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password" . " \\\n"
        . " " x 4
        . " --initdb"
        . " --db $ensembl_base" . " \\\n"
        . " " x 4
        . " --ensembl $ensembl_dir" . "\n\n";

    return ( $ensembl_base, $cmd );
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
