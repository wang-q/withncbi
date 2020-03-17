#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw();
use YAML::Syck qw();
use Path::Tiny qw();

use Template;
use Text::CSV_XS;
use WWW::Mechanize;
use HTML::TableExtract;
use Bio::DB::Taxonomy;

use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new();

=head1 NAME

wgs_prep.pl - prepare WGS materials

=head1 SYNOPSIS

    perl wgs_prep.pl [options]
      Options:
        --help      -?          brief help message

        --file, -f      STR     tab seperated file containing wgs prefix and name
        --outdir, -o    STR     output dir
        --fix                   sometimes WGS records miss assigning strain id
        --nofix         @STR    Skip some strains
        --csvonly               only generate the csv file

    perl wgs_prep.pl -f trichoderma.wgs.tsv -o WGS

    #name	prefix	organism	contigs

    Three files will be generated.

        trichoderma.wgs.csv
        trichoderma.wgs.rsync.sh
        trichoderma.wgs.data.yml

=cut

# for unrecorded strains, give them arbitrary ids
my $arbitrary = 100_000_000;

Getopt::Long::GetOptions(
    'help|?'     => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s'   => \my $infile,
    'outdir|o=s' => \( my $outdir = "." ),
    'fix'        => \my $fix_strain,
    'nofix=s'    => \my @nofix,
    'csvonly'    => \my $csvonly,
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Prepare NCBI WGS");

my $taxon_db = Bio::DB::Taxonomy->new( -source => 'entrez', );

#----------------------------#
# Read
#----------------------------#
$stopwatch->block_message("Load $infile.");
my $basename = Path::Tiny::path($infile)->basename( ".txt", ".tab", ".tsv" );

my $wgsid_of = {};
my @orig_orders;
{
    my @lines = Path::Tiny::path($infile)->lines;
    for my $line (@lines) {
        chomp $line;
        $line =~ /^#/ and next;
        my ( $name, $prefix ) = split /\t/, $line;
        $prefix or next;
        $wgsid_of->{$name} = $prefix;
        push @orig_orders, $name;
    }
}

#----------------------------#
# scraper
#----------------------------#
$stopwatch->block_message("Scrapping NCBI WGS...");
my $master = {};
{
    for my $key (@orig_orders) {
        print "$key\n";
        my $prefix = $wgsid_of->{$key};
        my $info   = wgs_worker($prefix);
        $info->{name} = $key;
        $master->{$key} = $info;
    }

    print "\n", "=" x 20, "\n";
    print "Finish scrapping\n";
}

#----------------------------#
# csv and rsync
#----------------------------#
$stopwatch->block_message("Generate .csv for info and .rsync.sh for downloading ");
{
    mkdir $outdir unless -d $outdir;

    my $csv = Text::CSV_XS->new( { binary => 1 } )
        or die "Cannot use CSV: " . Text::CSV_XS->error_diag;
    $csv->eol("\n");

    my $file_csv = Path::Tiny::path( $outdir, "$basename.csv" )->stringify;

    open my $csv_fh, ">", $file_csv;

    my @columns = (
        'prefix',                'taxon_id',     'name',         'Organism',
        'Biosource',             'BioProject',   'Keywords',     'Genome_Coverage',
        'Sequencing_Technology', '#_of_Contigs', 'Total_length', 'Assembly_Method',
        'Assembly_Name',         'Update_date',  'pubmed',
    );

    $csv->print( $csv_fh, \@columns );

    for my $key (@orig_orders) {

        # Don't use hashref here, because I want use hash slices.
        my %info = %{ $master->{$key} };

        if ( !$csvonly and $fix_strain ) {
            if ( grep { $_ eq $key } @nofix ) {
                print "Skip $info{name} as you don't want fix it\n";
            }
            elsif ( $info{Organism} =~ /$info{Biosource}/ ) {
                print "Don't need fixing for $info{name}\n";
            }
            else {
                # Sometimes the uploader didn't create a new strain, assign its own id and marked
                # species name as strain name.
                # So try looking up this strain in taxonomy dumps
                print "Fix strain taxon info for $info{name}\n";
                $info{Organism} = $info{Organism} . " " . $info{Biosource};

                my $node;
                eval { $node = $taxon_db->get_taxon( -name => $info{Organism} ); };
                if ( $@ or !$node ) {
                    print " " x 4, "Can't find taxon for $info{name}\n";
                    $arbitrary++;
                    print " " x 4, "Give it arbitrary id as $arbitrary\n";
                    $info{original_id} = $info{taxon_id};
                    $info{taxon_id}    = $arbitrary;
                }
                else {
                    $info{taxon_id} = $node->id;
                }
            }
        }
        $csv->print( $csv_fh, [ @info{@columns} ] );
        $master->{$key} = \%info;

    }
    close $csv_fh;

    print "\n", "=" x 20, "\n";
    print ".csv generated.\n";

}

#----------------------------#
# .rsync.sh and .aria2.sh
#----------------------------#
if ( !defined $csvonly ) {
    # rsync
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

EOF
    );

    for my $key (@orig_orders) {
        my %info = %{ $master->{$key} };

        my ($url) = grep {/\.fsa_nt\.gz/} @{ $info{download} };
        $url =~ s/\/[\w.]+fsa_nt\.gz$//;
        $url =~ /wgs_aux\/(.+)$/;
        my $rsync = $1;
        $rsync = "ftp.ncbi.nlm.nih.gov::sra/wgs_aux/" . $rsync;

        # ftp   - ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/JY/NM/JYNM02/JYNM02.1.fsa_nt.gz
        # rsync - ftp.ncbi.nlm.nih.gov::sra/wgs_aux/JY/NM/JYNM02/JYNM02.1.fsa_nt.gz
        # https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/BC/GH/BCGH01/BCGH01.1.fsa_nt.gz

        $file_rsync->append(
            <<"EOF"
echo >&2
echo >&2 "==> $key"
mkdir -p $key
rsync -avP $rsync/ $key/

EOF
        );
    }

    # aria2 for missing files
    my $file_aria2 = Path::Tiny::path( $outdir, "$basename.aria2.sh" );
    $file_aria2->remove if $file_rsync->is_file;

    $file_aria2->append(
        <<'EOF'
#!/bin/bash

BASE_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd ${BASE_DIR}

signaled () {
    echo Interrupted
    exit 1
}
trap signaled TERM QUIT INT

EOF
    );

    for my $key (@orig_orders) {
        my %info = %{ $master->{$key} };

        my ($url) = grep {/\.fsa_nt\.gz/} @{ $info{download} };
        $url =~ /\/([\w.]+fsa_nt\.gz)$/;
        my $filename = $1;

        $file_aria2->append(
            <<"EOF"
echo >&2
echo >&2 "==> $key"
if [ ! -e $key/$filename ]; then
    aria2c -UWget -x 6 -s 3 -c $url -d $key
fi

EOF
        );
    }

}

#----------------------------#
# @data yaml
#----------------------------#
if ( !defined $csvonly ) {
    $stopwatch->block_message("Generate .yml");

    my $file_data = Path::Tiny::path( $outdir, "$basename.data.yml" )->stringify;

    my $text = <<'EOF';
---
data:
[% FOREACH name IN names -%]
  - taxon: [% master.$name.taxon_id %]
    name: [% master.$name.name %]
    sciname: [% master.$name.Organism %]
    prefix: [% master.$name.prefix %]
    coverage: [% master.$name.Genome_Coverage %] [% master.$name.Sequencing_Technology %]
[% IF master.$name.original_id -%]
    original_id: [% master.$name.original_id %]
[% END -%]
[% END -%]

EOF
    my $tt = Template->new;
    $tt->process( \$text, { names => [@orig_orders], master => $master, }, $file_data )
        or die Template->error;

    print ".data.yml generated.\n";

    YAML::Syck::DumpFile( Path::Tiny::path( $outdir, "$basename.master.yml" )->stringify, $master );
    print ".master.yml generated.\n";
}

$stopwatch->end_message;

exit;

sub wgs_worker {
    my $term = shift;

    my $mech = WWW::Mechanize->new;
    $mech->stack_depth(0);    # no history to save memory

    # local shadowsocks proxy
    if ( $ENV{SSPROXY} ) {
        $mech->proxy( [ 'http', 'https' ], 'socks://127.0.0.1:1080' );
    }

    my $url_part = "http://www.ncbi.nlm.nih.gov/Traces/wgs/";
    my $url      = $url_part . '?val=' . $term;
    warn " " x 4 . $url . "\n";

    my $info = { prefix => $term, };
    $mech->get($url);

    {    # extract from tables
        my $page    = $mech->content;
        my @tables  = qw{ master-table structured-comments };
        my @columns = (
            '#_of_Contigs',    'Total_length',
            'Update_date',     'BioProject',
            'Keywords',        'Organism',
            'Assembly_Method', 'Assembly_Name',
            'Genome_Coverage', 'Sequencing_Technology',
            'Biosource',
        );

        for my $table (@tables) {
            print " " x 4 . "Extract from table ", $table, "\n";
            my $te = HTML::TableExtract->new( attribs => { class => $table, }, );
            $te->parse($page);

            for my $ts ( $te->table_states ) {
                for my $row ( $ts->rows ) {
                    for my $cell (@$row) {
                        if ($cell) {
                            $cell =~ s/[,:]//g;
                            $cell =~ s/^\s+//g;
                            $cell =~ s/\s+$//g;
                            $cell =~ s/\s+/ /g;
                        }
                    }
                    next unless $row->[0];
                    $row->[0] =~ s/\s+/_/g;
                    next unless grep { $row->[0] eq $_ } @columns;

                    $row->[1] =~ s/\s+.\s+show.+lineage.+$//g;
                    if ( $row->[0] eq 'Biosource' ) {
                        my ($biosource_strain)
                            = grep {/strain = /} grep {defined} split /\//,
                            $row->[1];

                        #print $row->[1], "\n";
                        $biosource_strain =~ s/strain = //;

                        #print $biosource_strain, "\n";
                        $info->{ $row->[0] } = $biosource_strain;
                    }
                    else {
                        $info->{ $row->[0] } = $row->[1];
                    }
                }
            }
        }
    }

    {    # taxon id
        my @links = $mech->find_all_links( url_regex => => qr{wwwtax}, );
        if ( @links and $links[0]->url =~ /\?id=(\d+)/ ) {
            $info->{taxon_id} = $1;
        }
    }

    {    # pubmed id
        my @links = $mech->find_all_links( url_regex => => qr{\/pubmed\/}, );
        if ( @links and $links[0]->url =~ /\/pubmed\/(\d+)/ ) {
            $info->{pubmed} = $1;
        }
    }

    {    # downloads
        my @links = $mech->find_all_links(
            text_regex => qr{$term},
            url_regex  => qr{\.gz$},
        );
        $info->{download} = [ map { $_->url } @links ];
    }

    return $info;
}

__END__
