#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use Template;

use FindBin;
use lib "$FindBin::Bin/../lib";
use MyUtil qw(replace_home);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::Bin/../config.ini");

my $seq_dir;
my $taxon_file;

my $withncbi = replace_home( $Config->{run}{withncbi} );    # withncbi path

my $parallel = 4;

my $man  = 0;
my $help = 0;

GetOptions(
    'help'         => \$help,
    'man'          => \$man,
    'seq_dir=s'    => \$seq_dir,
    'taxon_file=s' => \$taxon_file,
    'parallel=i'   => \$parallel,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# start
#----------------------------------------------------------#
my $tt = Template->new;

my $text = <<'EOF';
# [% name %]
perl [% withncbi %]/taxon/strain_bz.pl \
    --file [% taxon_file %] \
[% IF seq_dir -%]
    --seq_dir [% seq_dir %] \
[% END -%]
    --use_name \
    --name [% name %] \
    --parallel [% parallel %] \
[% IF o -%]
    -o [% o %] \
[% END -%]
[% FOREACH q IN qs -%]
    -q [% q %] \
[% END -%]
    -t [% t %]

EOF

while (<>) {
    chomp;
    /^#/ and next;

    my @fields = split /\t/;
    next unless scalar @fields >= 3;

    $tt->process(
        \$text,
        {   seq_dir    => $seq_dir,
            taxon_file => $taxon_file,
            withncbi   => $withncbi,
            parallel   => $parallel,
            name       => $fields[0],
            t          => $fields[1],
            qs         => [ split /,/, $fields[2] ],
            o          => $fields[3],
        },
        \*STDOUT
    ) or die Template->error;
}

exit;

__END__

=head1 NAME

cmd_template.pl - Simple template for strain_bz.pl

=head1 SYNOPSIS

    cat <file> | perl cmd_template.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --seq_dir           ~/data/organelle/plastid_genomes
        --taxon_file        ~/data/organelle/plastid_genomes/plastid_ncbi.csv
        --parallel          number of child processes

=cut
