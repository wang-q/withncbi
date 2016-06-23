#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use Template;
use Path::Tiny;

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

pop_prep.pl - prepare pop: file, rm and info.

=head1 SYNOPSIS

    perl pop_prep.pl [options]
      Options:
        --help      -?          brief help message
        --file      -i  STR     input yaml
        --parallel      INT     number of child processes

=cut

my $withncbi = path( $Config->{run}{withncbi} )->stringify;
my $egaz     = path( $Config->{run}{egaz} )->stringify;

GetOptions(
    'help|?'     => sub { Getopt::Long::HelpMessage(0) },
    'file|i=s'   => \my $file_yaml,
    'parallel=i' => \my $parallel,
) or Getopt::Long::HelpMessage(1);

die "Need a YAML file" unless $file_yaml;

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
$stopwatch->start_message("Prepare pop");

my $yml = LoadFile($file_yaml);

# these three are needed
my $group_name = $yml->{group_name};
my $base_dir   = path( $yml->{base_dir} )->stringify;
my $data_dir   = path( $yml->{data_dir} )->stringify;

# Can load from .ini file
my $split_about
    = exists $yml->{split_about}
    ? $yml->{split_about}
    : $Config->{pop}{split_about};
my $min_contig
    = exists $yml->{min_contig}
    ? $yml->{min_contig}
    : $Config->{pop}{min_contig};
my $per_seq_min_contig
    = exists $yml->{per_seq_min_contig}
    ? $yml->{per_seq_min_contig}
    : $Config->{pop}{per_seq_min_contig};

if ( !$parallel ) {
    $parallel
        = exists $yml->{parallel} ? $yml->{parallel} : $Config->{run}{parallel};
}

my $rm_species;
if ( exists $yml->{rm_species} ) {
    $rm_species = $yml->{rm_species};
}

my $phylo_tree;
if ( exists $yml->{phylo_tree} ) {
    $phylo_tree = $yml->{phylo_tree};
}

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#
print "Create dir\n";
my @data = @{ $yml->{data} };
for my $item (@data) {
    if ( exists $item->{skip} ) {
        printf " " x 4 . $item->{name} . " SKIP! %s\n", $item->{skip};
    }
    else {
        printf " " x 4 . $item->{name} . "\n";
        my $dir = path( $data_dir, 'Genomes', $item->{name} );
        $dir->mkpath;
        $item->{dir} = $dir->stringify;
    }
}

{
    my $tt = Template->new;
    my $text;
    my $sh_name;

    # 01_file.sh
    $sh_name = "01_file.sh";
    print "Create $sh_name\n";
    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------------------------------------#
# unzip, filter and split
#----------------------------------------------------------#

[% FOREACH item IN data -%]
#----------------------------#
# [% item.name %]
#----------------------------#
echo [% item.name %]
cd [% item.dir %]

[% IF item.skip -%]
echo '    SKIP! [% item.skip %]'
[% ELSIF item.downloaded -%]
echo '    Downloaded files.'
[% IF item.pre_dir -%]
find [% item.pre_dir %] -type f | parallel --no-run-if-empty cp {} .
[% END -%]
[% ELSE -%]
echo '    Unzip, filter and split.'
gzip -d -c --force [% item.fasta %] > toplevel.fa
perl -p -i -e '/\>gi\|/ and s/\>gi\|(\d+).*/\>gi_$1/' toplevel.fa
[% IF item.per_seq -%]
faops count toplevel.fa | perl -aln -e 'next if $F[0] eq 'total'; print $F[0] if $F[1] > [% per_seq_min_contig * 10 %]; print $F[0] if $F[1] > [% per_seq_min_contig %]  and $F[6]/$F[1] < 0.05' | uniq > listFile
faops some toplevel.fa listFile toplevel.filtered.fa
faops split-name toplevel.filtered.fa .
[% ELSE -%]
faops count toplevel.fa | perl -aln -e 'next if $F[0] eq 'total'; print $F[0] if $F[1] > [% min_contig * 10 %]; print $F[0] if $F[1] > [% min_contig %]  and $F[6]/$F[1] < 0.05' | uniq > listFile
faops some toplevel.fa listFile toplevel.filtered.fa
faops split-about toplevel.filtered.fa [% split_about %] .
[% END -%]
rm toplevel.fa toplevel.filtered.fa listFile
[% END -%]

[% END -%]

EOF

    $tt->process(
        \$text,
        {   data               => \@data,
            data_dir           => $data_dir,
            split_about        => $split_about,
            min_contig         => $min_contig,
            per_seq_min_contig => $per_seq_min_contig,
        },
        path( $data_dir, $sh_name )->stringify
    ) or die Template->error;

    # 02_rm.sh
    $sh_name = "02_rm.sh";
    print "Create $sh_name\n";
    $text = <<'EOF';
#!/bin/bash

#----------------------------------------------------------#
# RepeatMasker
#----------------------------------------------------------#
cd [% data_dir %]
echo Doing RepeatMasker

[% FOREACH item IN data -%]
#----------------------------#
# [% item.name %]
#----------------------------#
echo [% item.name %]

[% IF item.skip -%]
echo '    SKIP! [% item.skip %]'
[% ELSE -%]
cd [% item.dir %]

for f in `find [% item.dir%] -name "*.fa"` ; do
    rename 's/fa$/fasta/' $f ;
done

for f in `find [% item.dir%] -name "*.fasta"` ; do
    RepeatMasker $f [% IF rm_species %]-species [% rm_species %][% END %] -xsmall --parallel [% parallel %] ;
done

for f in `find [% item.dir%] -name "*.fasta.out"` ; do
    rmOutToGFF3.pl $f > `dirname $f`/`basename $f .fasta.out`.rm.gff;
done
[% END -%]

[% END -%]

#----------------------------------------------------------#
# Clean RepeatMasker
#----------------------------------------------------------#
cd [% data_dir %]
echo Cleaning RepeatMasker

[% FOREACH item IN data -%]
#----------------------------#
# [% item.name %]
#----------------------------#
echo [% item.name %]

[% IF item.skip -%]
echo '    SKIP! [% item.skip %]'
[% ELSE -%]
cd [% item.dir %]
for f in `find [% item.dir%] -name "*.fasta"` ; do
    if [ -f $f.masked ];
    then
        rename 's/fasta.masked$/fa/' $f.masked;
        find [% item.dir%] -type f -name "`basename $f`*" | xargs rm;
    else
        rename 's/fasta$/fa/' $f;
        echo `date` "RepeatMasker on $f failed.\n" >> RepeatMasker.log
        find [% item.dir%] -type f -name "`basename $f`*" | xargs rm;
    fi;
done;
[% END -%]

[% END -%]

EOF

    $tt->process(
        \$text,
        {   data       => \@data,
            data_dir   => $data_dir,
            rm_species => $rm_species,
            parallel   => $parallel,
        },
        path( $data_dir, $sh_name )->stringify
    ) or die Template->error;

    # 03_strain_info.sh
    $sh_name = "03_strain_info.sh";
    print "Create $sh_name\n";
    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------------------------------------#
# generate taxon file
#----------------------------------------------------------#
perl [% withncbi %]/taxon/strain_info.pl \
[% FOREACH item IN data -%]
[% IF ! item.skip -%]
    --id   [% item.taxon %] \
    --name [% item.taxon %]=[% item.name %] \
[% IF item.original_id -%]
    --species [% item.taxon %]=[% item.original_id %] \
[% END -%]
[% END -%]
[% END -%]
    --file [% data_dir %]/[% group_name %].taxon.csv

EOF

    $tt->process(
        \$text,
        {   data       => \@data,
            group_name => $group_name,
            data_dir   => $data_dir,
            base_dir   => $base_dir,
            withncbi   => $withncbi,
            parallel   => $parallel,
        },
        path( $data_dir, $sh_name )->stringify
    ) or die Template->error;

    # plan_ALL.sh
    $sh_name = "plan_ALL.sh";
    print "Create $sh_name\n";
    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------------------------------------#
# alignment plan of all genomes
#----------------------------------------------------------#
# Don't copy sequences (RepeatMasker done)
# This plan includes all genomes.
[% IF phylo_tree -%]
# Use [% phylo_tree %] as guide tree for other plans.
[% ELSE -%]
# Use the generated phylogenetic tree in this step as guide tree for other plans.
[% END -%]

# plan_ALL
cd [% data_dir %]
perl [% egaz %]/multi_batch.pl \
    --file [% data_dir %]/[% group_name %].taxon.csv \
    -w     [% base_dir %] \
    --name [% group_name %] \
    --multi_name plan_ALL \
    --use_name \
    --parallel [% parallel %] \
    --norm \
[% IF phylo_tree -%]
    --phylo_tree [% phylo_tree %] \
[% END -%]
[% FOREACH item IN data -%]
[% IF loop.index != 0 -%]
[% IF ! item.skip -%]
    -q [% item.name %] \
[% END -%]
[% END -%]
[% END -%]
    -t [% data.0.name %]

EOF

    $tt->process(
        \$text,
        {   data       => \@data,
            group_name => $group_name,
            data_dir   => $data_dir,
            base_dir   => $base_dir,
            phylo_tree => $phylo_tree,
            withncbi   => $withncbi,
            egaz       => $egaz,
            parallel   => $parallel,
        },
        path( $data_dir, $sh_name )->stringify
    ) or die Template->error;

    if ( exists $yml->{plans} ) {
        print "Create .sh for each plans\n";
        my @plans = @{ $yml->{plans} };
        for my $plan (@plans) {
            my $plan_name = $plan->{name};

            $sh_name = "plan_$plan_name.sh";
            print " " x 4, "$sh_name\n";
            $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------------------------------------#
# alignment plan for [% plan_name %]
#----------------------------------------------------------#
perl [% egaz %]/multi_batch.pl \
    --file [% data_dir %]/[% group_name %].taxon.csv \
    -w     [% base_dir %] \
    --name [% group_name %] \
    --multi_name [% plan_name %] \
    --use_name \
    --parallel [% parallel %] \
    --norm \
[% IF phylo_tree -%]
    --phylo_tree [% phylo_tree %] \
[% END -%]
[% IF plan.o -%]
    -o [% plan.o %] \
[% END -%]
[% FOREACH q IN plan.qs -%]
    -q [% q %] \
[% END -%]
    -t [% plan.t %]

EOF

            $tt->process(
                \$text,
                {   group_name => $group_name,
                    data_dir   => $data_dir,
                    base_dir   => $base_dir,
                    withncbi   => $withncbi,
                    egaz       => $egaz,
                    parallel   => $parallel,
                    plan_name  => $plan_name,
                    plan       => $plan,
                },
                path( $data_dir, $sh_name )->stringify
            ) or die Template->error;
        }
    }
}

$stopwatch->end_message;

exit;

__END__
