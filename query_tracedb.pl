#!/usr/bin/perl -w
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use threads;
use threads::shared;
use LWP::UserAgent;
use HTTP::Request::Common 'POST';
use File::Spec::Functions qw(catfile);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $outdir    = undef;
my $parallel  = 5;
my $species   = 'ACYRTHOSIPHON PISUM';    # pea aphid
my $page_size = 1000;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'outdir=s'   => \$outdir,
    'parallel=i' => \$parallel,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

unless ( defined $outdir ) {
    $outdir = $species;
    $outdir =~ s/\s/_/g;
    print "Write to directory $outdir\n";
    mkdir $outdir unless -d $outdir;
}

my $http_proxy = 'http://127.0.0.1:8088/';

$ENV{'LANG'}   = 'C';
$ENV{'LC_ALL'} = 'C';

#----------------------------------------------------------#
# Download
#----------------------------------------------------------#
$|++;

my $count = get_count($species);
print "$species has $count traces\n";

my $total_page = int( $count / $page_size );
my @pages;
share(@pages);
@pages = ( 0 .. $total_page );
print "Page_size is $page_size; there are total $total_page pages\n";

my $download_tgz = sub {
    my $thread_id   = shift;
    my $thread_name = "thread $thread_id: ";

GET: while (1) {
        my $page_number;
        {
            lock @pages;
            last GET if @pages == 0;
            $page_number = shift @pages;
        }
        print $thread_name, "getting page: $page_number\n";

        my $filename = catfile( $outdir, "page$page_number.tgz" );
        my $is_success
            = get_data( $species, $page_size, $page_number, $filename,
            $thread_name );
        if ($is_success) {
            print $thread_name, "Done.\n";
        }
        else {
            lock @pages;
            unshift @pages, $page_number;
            print $thread_name, "Broken file downloaded. Retry...\n";
            redo GET;
        }
    }

    return;
};

my @workers;
for ( 1 .. $parallel ) {
    $workers[$_] = threads->new( $download_tgz, $_ );
}
for ( 1 .. $parallel ) {
    $workers[$_]->join;
}

exit;

#----------------------------------------------------------#
# Subroutine
#----------------------------------------------------------#

sub get_count {
    my $species = shift;

    my $query      = "query count species_code='$species'";
    my $req_result = trace_request($query);
    $req_result =~ m{(\d+)};

    return $1;
}

sub get_data {
    my $species     = shift;
    my $page_size   = shift;
    my $page_number = shift;
    my $filename    = shift;
    my $thread_name = shift;

    $thread_name ||= ' ' x 4;

    my $bin_query
        = "query page_size $page_size page_number $page_number binary species_code='$species'";
    my $bin_result = trace_request($bin_query);
    return unless $bin_result;
    print $thread_name, "Got bin page $page_number\n";

    my $tgz_query  = "retrieve_tgz fasta 0b" . $bin_result;
    my $tgz_result = trace_request($tgz_query);
    return unless $tgz_result;
    print $thread_name, "Got tgz page $page_number\n";

    open my $fh, '>', $filename
        or die "Cannot open file: $!";
    binmode $fh;
    print {$fh} $tgz_result;
    close $fh;
    print $thread_name, "Write file page $page_number\n";

    return 1;
}

sub trace_request {
    my $query = shift;

    my $req = POST 'http://trace.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=raw',
        [ query => $query ];

    my $ua = LWP::UserAgent->new;
    $ua->proxy( [ 'http', 'ftp' ], $http_proxy );

    my $res = $ua->request($req);

    #print Dump($res), "\n";
    if ( $res->is_success ) {
        return $res->content;    # or whatever
    }
    else {
        warn "Couldn't connect to TRACE server\n", $res->status_line, "\n";
        return;
    }
}

__END__

Reference:
[1] http://www.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=show&f=doc&m=obtain&s=stips
[2] ftp://ftp.ncbi.nih.gov/pub/TraceDB/misc/query_tracedb

 How to download large data sets

The number of records which can be obtained on a single request is limited.
Currently this number is set to 40,000. In order to download more records, you
would need to place several requests accordingly. Although it is generally
possible to download all needed data with a browser, the best approach to do
this job is to use our Perl script query_tracedb. After copying this script,
don't forget to make it executable.

All records in the archive are assigned a unique identifier - TI, and
therefore, first, you would need to obtain all identifiers which comply to
your query. Using these identifiers you can then retrieve the actual data that
you need. Let's see how this works on a real example (please note that this
page is static, and all the numbers shown in the example may not reflect the
current status of the archive):

   1. The first step is to count all available records:
      query_tracedb "query count species_code='AEDES AEGYPTI'"
      122116
   2. A simple calculation shows that to retrieve all records we will need to
   make at least 4 requests, so let's obtain the identifiers. Please note that
   the identifiers are in network (BIG ENDIAN) format:
      query_tracedb "query page_size 40000 page_number 0 binary species_code='AEDES AEGYPTI'" > page1.bin
      query_tracedb "query page_size 40000 page_number 1 binary species_code='AEDES AEGYPTI'" > page2.bin
      ...
      query_tracedb "query page_size 40000 page_number 3 binary species_code='AEDES AEGYPTI'" > page4.bin
   3. You can now retrieve the data in the submission form (tarball). Pointer
   "0b" shows that following data are in binary format.
      (echo -n "retrieve_tgz all 0b"; cat page1.bin) | query_tracedb > data1.tgz
      ...
      (echo -n "retrieve_tgz all 0b"; cat page4.bin) | query_tracedb > data4.tgz
      The above will retrieve all files from the archive: fasta, quality scores, chromatograms in scf format, mate_pairs, and ancillary files.
   4. *Note: steps 2 and 3 can be done at the same time:
      (echo -n "retrieve_tgz all 0b"; query_tracedb "query page_size 40000 page_number 0 binary species_code='AEDES AEGYPTI'") | query_tracedb > data1.tgz

For more information please apply 'query_tracedb help' for available data
formats, and 'query_tracedb usage' for usage examples.

If you need to save only TI numbers for future reference, you might want to
obtain them in text form:
query_tracedb "query page_size 40000 page_number 0 text species_code='AEDES AEGYPTI'" > page1.txt
