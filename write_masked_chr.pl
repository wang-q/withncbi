#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML::XS qw(Dump Load DumpFile LoadFile);

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;

use File::Find::Rule;
use File::Basename;
use String::Compare;
use List::Util qw(reduce);
use List::MoreUtils qw(zip);

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB::Ensembl;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server     = $Config->{database}{server};
my $port       = $Config->{database}{port};
my $username   = $Config->{database}{username};
my $password   = $Config->{database}{password};
my $ensembl_db = $Config->{database}{ensembl};

my $feature = "repeat";
my $yaml_file;
my $dir_fa;
my $linelen = 60;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'        => \$help,
    'man'           => \$man,
    'server=s'      => \$server,
    'port=i'        => \$port,
    'username=s'    => \$username,
    'password=s'    => \$password,
    'ensembl=s'     => \$ensembl_db,
    'feature=s'     => \$feature,
    'y|yaml_file=s' => \$yaml_file,
    'dir=s'         => \$dir_fa,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Write masked chr from $ensembl_db...");

#----------------------------#
# Get runlist of $feature
#----------------------------#
my $ftr_of = {};

if ($yaml_file) {
    $ftr_of = LoadFile($yaml_file);
}
else {

    # ensembl handler
    my $ensembl = AlignDB::Ensembl->new(
        server => $server,
        db     => $ensembl_db,
        user   => $username,
        passwd => $password,
    );

    # ensembl handler
    my $db_adaptor = $ensembl->db_adaptor;

    my $slice_adaptor = $db_adaptor->get_SliceAdaptor;
    my @slices        = @{ $slice_adaptor->fetch_all('chromosome') };
    while ( my $chr_slice = shift @slices ) {
        my $chr   = $chr_slice->seq_region_name;
        my $start = $chr_slice->start;
        my $end   = $chr_slice->end;

        print "$chr:$start:$end\n";
        my $slice_set = AlignDB::IntSpan->new("$start-$end");    # repeat set

        eval { $ensembl->set_slice( $chr, $start, $end ); };
        if ($@) {
            warn "Can't get annotation\n";
            next;
        }

        my $slice       = $ensembl->slice;
        my $ftr_chr_set = $slice->{"_$feature\_set"};

        $ftr_of->{$chr} = $ftr_chr_set->runlist;
    }

    DumpFile( "${ensembl_db}_${feature}.yml", $ftr_of );
}

#----------------------------#
# Soft mask
#----------------------------#
if ($dir_fa) {
    my @files = sort File::Find::Rule->file->name('*.fa')->in($dir_fa);

    my @chrs = map { basename $_ , ".fa" } @files; # strip dir and suffix
    while (1) {
        my $lcss = lcss(@chrs);
        last unless $lcss;
        print "LCSS [$lcss]\n";
        my $rx = quotemeta $lcss;
        $chrs[$_] =~ s/$rx// for 0 .. $#chrs;
    }
    my $file_of = { zip( @chrs, @files ) };

    for my $file_chr ( sort @chrs ) {

        # use the most similar chr name
        my ($ftr_chr) = map { $_->[0] }
            sort { $b->[1] <=> $a->[1] }
            map { [ $_, compare( $_, $file_chr ) ] } keys %{$ftr_of};

        printf "Write masked seq for ftr_chr %8s file_chr %8s\n", $ftr_chr,
            $file_chr;

        my $ftr_set = AlignDB::IntSpan->new( $ftr_of->{$ftr_chr} );

        my ( $seq_of, $seq_names ) = read_fasta( $file_of->{$file_chr} );
        my $seq = $seq_of->{ $seq_names->[0] };

        my @sets = $ftr_set->sets;
        for my $set (@sets) {
            my $offset = $set->min - 1;
            my $length = $set->size;

            my $str = substr $seq, $offset, $length;
            $str = lc $str;
            substr $seq, $offset, $length, $str;
        }

        open my $out_fh, '>', $file_of->{$file_chr} . ".masked";
        print {$out_fh} ">chr$ftr_chr\n";
        print {$out_fh} substr( $seq, 0, $linelen, '' ) . "\n" while ($seq);
        close $out_fh;
    }
}

$stopwatch->end_message;

# comes from
# http://stackoverflow.com/questions/499967/how-do-i-determine-the-longest-similar-portion-of-several-strings
sub lcss {
    return '' unless @_;
    return $_[0] if @_ == 1;
    my $i          = 0;
    my $first      = shift;
    my $min_length = length($first);
    for (@_) {
        $min_length = length($_) if length($_) < $min_length;
    }
INDEX: for my $ch ( split //, $first ) {
        last INDEX unless $i < $min_length;
        for my $string (@_) {
            last INDEX if substr( $string, $i, 1 ) ne $ch;
        }
    }
    continue { $i++ }
    return substr $first, 0, $i;
}

__END__

=head1 NAME

    write_masked_chr.pl - just like RepeatMasker does, but use ensembl
                          annotations.
                          And it change fasta headers for you.

=head1 SYNOPSIS

    write_masked_chr.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --username          username
        --password          password
        --ensembl           ensembl database name
        --feature           mask which feature, default is "repeat"
        -y, --yaml_file     use a exists yaml file as annotation
        --dir               .fa dir
    
    > perl write_masked_chr.pl -e arabidopsis_58 
    > perl write_masked_chr.pl --dir e:\data\alignment\arabidopsis\ath_58\ -y e:\wq\Scripts\alignDB\util\arabidopsis_58_repeat.yml
    
    $ perl write_masked_chr.pl -e nip_65
    $ perl write_masked_chr.pl --dir /home/wangq/data/alignment/rice/nip_65 -y nip_65_repeat.yml

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

