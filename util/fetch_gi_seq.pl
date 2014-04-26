#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use WWW::Mechanize;

use AlignDB::Run;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $gi_file          = '';
my $length_threshold = 10;
my $outfile          = undef;
my $failed_file      = undef;
my $parallel         = 5;
my $batch_number     = 100;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'length=i'   => \$length_threshold,
    'infile=s'   => \$gi_file,
    'outfile=s'  => \$outfile,
    'failed=s'   => \$failed_file,
    'parallel=i' => \$parallel,
    'batch=i'    => \$batch_number,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

unless ( defined $outfile ) {
    $outfile = $gi_file;
    $outfile =~ s/\.\w+?$/\.fasta/;
    print "Write to $outfile\n";
}
unless ( defined $failed_file ) {
    $failed_file = $gi_file;
    $failed_file =~ s/\.\w+?$/\.failed/;
}

#----------------------------------------------------------#
# Download
#----------------------------------------------------------#
$|++;

my @jobs;
{    # read gi file
    open my $infh, "<", $gi_file;
    my $file_content = do { local $/; <$infh> };
    close $infh;
    my @gids = $file_content =~ /gi\:(\d+)/g;

    while ( scalar @gids ) {
        my @batching = splice @gids, 0, $batch_number;
        push @jobs, [@batching];
    }
}

#----------------------------------------------------------#
# Worker
#----------------------------------------------------------#
my $worker = sub {
    my $job = shift;
    my $opt = shift;

    my @gids  = @$job;
    my $retry = $opt->{retry};

    my $mech = WWW::Mechanize->new;
    $mech->stack_depth(0);    # no history to save memory

    my $url_part_1 = "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?";
    my $url_part_2 = "dopt=fasta&dispmax=5&sendto=t&db=nucleotide";
    my $url_part_3 = "&qty=1&c_start=1&list_uids=";
    my $url_part_4 = "&from=begin&to=end";

    my $error = 0;
    my $content;
GET: while ( scalar @gids ) {
        my $gi = shift @gids;

        my $url = join '',
            ( $url_part_1, $url_part_2, $url_part_3, $gi, $url_part_4 );

        print "Getting gi: [$gi]\n";
        $mech->get($url);

        my $html_content = $mech->content;
        if ( $html_content =~ /^>.+\n\n$/s ) {
            $error = 0;
            if ( length $html_content > $length_threshold ) {
                $content .= $html_content;
                print "Download [", length $html_content, "] bytes.\n";
            }
            else {
                print "Seq too short. Jump to next\n";
            }
        }
        else {
            $error++;
            unshift @gids, $gi;
            if ( $error < $retry ) {
                print "Broken file downloaded. Retry...\n";
                redo GET;
            }
            else {
                print "Connection lost.\n";
                last GET;
            }
        }
    }

    open my $outfh, ">>", $outfile
        or die "Cannot append to file: $!";
    print {$outfh} $content;
    close $outfh;
    print "Write seqs to $outfile\n";

    if ( @gids > 0 ) {
        open my $failed_fh, ">>", $failed_file
            or die "Cannot append to file: $!";
        print {$failed_fh} join '', map {"gi:$_\n"} @gids;
        close $failed_fh;
        print "Write failed gids to $failed_file\n";
    }

    return;
};

#----------------------------------------------------------#
# Run
#----------------------------------------------------------#
my $run = AlignDB::Run->new(
    parallel => $parallel,
    jobs     => \@jobs,
    code     => $worker,
    opt      => { retry => $batch_number * 5 },
);
$run->run;

exit;

__END__

=head1 NAME

    fetch_gi_seq.pl - Fetch fasta file of gi's from NCBI

=head1 SYNOPSIS

    fetch_gi_seq.pl [options]
        Options:
            --help              brief help message
            --man               full documentation
            --length            length threshold
            --infile            gi filename
            --outfile            output dir of fasta files
            --length            length threshold
            --parallel          run in parallel mode

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
