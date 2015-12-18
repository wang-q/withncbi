#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use Template;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::RealBin/../config.ini");

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

my $withncbi = path( $Config->{run}{withncbi} )->stringify;    # withncbi path

GetOptions(
    'help|?'       => sub { HelpMessage(0) },
    'seq_dir=s'    => \my $seq_dir,
    'taxon_file=s' => \my $taxon_file,
    'parallel=i' => \( my $parallel = 4 ),
) or HelpMessage(1);

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
    -q [% o %] \
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
