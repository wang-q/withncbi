#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw();
use YAML::Syck qw();
use Path::Tiny qw();

use Template;
use Text::CSV_XS;

use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new();

=head1 NAME

assembly_prep.pl - prepare ASSEMBLY materials

=head1 SYNOPSIS

    perl assembly_prep.pl [options]
      Options:
        --help, -?              brief help message

        --file, -f      STR     tab seperated file containing wgs prefix and name
        --outdir, -o    STR     output dir
        --csvonly               only generate the csv file

    perl assembly_prep.pl -f trichoderma.assembly.tsv -o ASSEMBLY

    #name   ftp_path    organism    assembly_level

    Three files will be generated:

        trichoderma.assembly.rsync.sh
        rsync.tsv
        trichoderma.assembly.collect.sh

    The latter one will create:

        trichoderma.assembly.collect.csv

=cut

Getopt::Long::GetOptions(
    'help|?'     => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s'   => \my $infile,
    'outdir|o=s' => \( my $outdir = "." ),
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Prepare NCBI ASSEMBLY");

#----------------------------#
# Read
#----------------------------#
$stopwatch->block_message("Load $infile.");
my $basename = Path::Tiny::path($infile)->basename( ".txt", ".tab", ".tsv" );

my $ftp_of = {};
my @orig_orders;
{
    my @lines = Path::Tiny::path($infile)->lines;
    for my $line (@lines) {
        chomp $line;
        $line =~ /^#/ and next;
        my ( $name, $ftp ) = split /\t/, $line;
        $ftp or next;
        $ftp_of->{$name} = $ftp;
        push @orig_orders, $name;
    }
}

Path::Tiny::path($outdir)->mkpath();

#----------------------------#
# rsync
#----------------------------#
{
    $stopwatch->block_message("Generate .rsync.sh");

    my $file_url = Path::Tiny::path( $outdir, "rsync.tsv" );
    $file_url->remove if $file_url->is_file;

    for my $key (@orig_orders) {
        my $ftp = $ftp_of->{$key};

        # ftp   - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/167/675/GCA_000167675.2_v2.0
        # rsync - ftp.ncbi.nlm.nih.gov::genomes/all/GCA/000/167/675/GCA_000167675.2_v2.0

        my $rsync = $ftp;
        $rsync =~ s/(ftp|https?):\/\/ftp.ncbi.nlm.nih.gov\//ftp.ncbi.nlm.nih.gov::/;
        if ( $rsync eq $ftp ) {
            die "Check the ftp url: [$key] $ftp\n";
        }

        $file_url->append( sprintf( "%s\t%s\n", $key, $rsync ) );
    }

    my $file_rsync = Path::Tiny::path( $outdir, "$basename.rsync.sh" );
    $file_rsync->remove if $file_rsync->is_file;

    $file_rsync->append(
        <<'EOF'
#!/bin/bash

BASE_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd ${BASE_DIR}

signaled () {
    echo Interrupted
    exit 1
}
trap signaled TERM QUIT INT

cat rsync.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
        echo >&2
        echo >&2 "==> {1}"
        mkdir -p {1}
        rsync -avP {2}/ {1}/
    '

EOF
    );

}

#----------------------------#
# collect
#----------------------------#
{
    $stopwatch->block_message("Generate .collect.sh");

    my $file_collect = Path::Tiny::path( $outdir, "$basename.collect.sh" )->stringify;

    my @columns = (
        "Organism name",
        "Taxid",
        "Assembly name",
        "Infraspecific name",
        "BioSample",
        "BioProject",
        "Submitter",
        "Date",
        "Assembly type",
        "Release type",
        "Assembly level",
        "Genome representation",
        "WGS project",
        "Assembly method",
        "Genome coverage",
        "Sequencing technology",
        "RefSeq category",
        "RefSeq assembly accession",
        "GenBank assembly accession",
    );

    my $template = <<'EOF';
#!/bin/bash

BASE_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd ${BASE_DIR}

signaled () {
    echo Interrupted
    exit 1
}
trap signaled TERM QUIT INT

echo "name,[% header %]" \
    > [% basename %].collect.csv

[% FOREACH n IN names -%]
echo >&2
echo >&2 "==> [% n %]"
find [% n %] -type f -name "*_assembly_report.txt" |
    xargs cat |
    perl -nl -e '
        BEGIN { our %stat = (); }

        m{^#\s+} or next;
        s/^#\s+//;
        @O = split /\:\s*/;
        scalar @O == 2 or next;
        $O[0] =~ s/\s*$//g;
        $O[0] =~ s/\W/_/g;
        $O[1] =~ /([\w =.-]+)/ or next;
        $stat{$O[0]} = $1;

        END {
            my @c;
            for my $key ( qw{ [% columns %] } ) {
                if (exists $stat{$key}) {
                    push @c, $stat{$key};
                }
                else {
                    push @c, q{};
                }
            }
            print join(q{,}, q{[% n %]}, @c);
        }
    ' \
    >> [% basename %].collect.csv

[% END -%]

EOF
    my $tt = Template->new;
    $tt->process(
        \$template,
        {   basename => $basename,
            names    => \@orig_orders,
            header   => join( ",", map { s/\s+/_/g; $_ } @columns ),
            columns  => join( " ", map { s/\s+/_/g; $_ } @columns ),
        },
        $file_collect
    ) or die Template->error;
}

$stopwatch->end_message;

exit;

__END__
