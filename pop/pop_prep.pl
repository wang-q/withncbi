#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use Template;
use File::Basename;
use File::Find::Rule;
use File::Spec;
use String::Compare;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use MyUtil qw(replace_home);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::Bin/../config.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

my $withncbi = replace_home( $Config->{run}{withncbi} );    # withncbi path

my $file_yaml;

my $parallel;

my $man  = 0;
my $help = 0;

GetOptions(
    'help'       => \$help,
    'man'        => \$man,
    'i|file=s'   => \$file_yaml,
    'parallel=i' => \$parallel,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
$stopwatch->start_message("Prepare pop");
die "Need a YAML file" unless $file_yaml;

my $yml = LoadFile($file_yaml);

# these three are needed
my $group_name = $yml->{group_name};
my $base_dir   = replace_home( $yml->{base_dir} );
my $data_dir   = replace_home( $yml->{data_dir} );

# Can load from .ini file
my $split_about
    = exists $yml->{split_about}
    ? $yml->{split_about}
    : $Config->{pop}{split_about};
my $min_contig
    = exists $yml->{min_contig}
    ? $yml->{min_contig}
    : $Config->{pop}{min_contig};

if ( !$parallel ) {
    $parallel
        = exists $yml->{parallel} ? $yml->{parallel} : $Config->{run}{parallel};
}

my $rm_species;
if ( exists $yml->{rm_species} ) {
    $rm_species = $yml->{rm_species};
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
        my $dir = File::Spec->catdir( $data_dir, $item->{name} );
        mkdir $dir if !-e $dir;
        $item->{dir} = $dir;
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

#----------------------------#
# unzip, filter and split
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %] 
echo [% item.name %]

[% IF item.skip -%]
echo '    SKIP! [% item.skip %]'
[% ELSIF item.downloaded -%]
echo '    Downloaded files. [% item.skip %]'
[% ELSE -%]
cd [% item.dir %]
gzip -d -c [% item.fasta %] > toplevel.fa
perl -p -i -e '/>/ and s/\>gi\|(\d+).*/\>gi_$1/' toplevel.fa
faops count toplevel.fa | perl -aln -e 'next if $F[0] eq 'total'; print $F[0] if $F[1] > [% min_contig * 10 %]; print $F[0] if $F[1] > [% min_contig %]  and $F[6]/$F[1] < 0.05' | uniq > listFile
faops some toplevel.fa listFile toplevel.filtered.fa
[% IF item.per_seq -%]
faops split-name toplevel.filtered.fa .
[% ELSE -%]
faops split-about toplevel.filtered.fa [% split_about %] .
[% END -%]
rm toplevel.fa toplevel.filtered.fa listFile

rename 's/fa$/fasta/' *.fa;
[% END -%]

[% END -%]

EOF

    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            split_about => $split_about,
            min_contig  => $min_contig,
        },
        File::Spec->catfile( $data_dir, $sh_name )
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
# [% item.name %] [% item.coverage %]
#----------------------------#
echo [% item.name %]

[% IF item.skip -%]
echo '    SKIP! [% item.skip %]'
[% ELSE -%]
cd [% item.dir %]
RepeatMasker [% item.dir %]/*.fasta [% IF rm_species %]-species [% rm_species %][% END %] -xsmall --parallel [% parallel %]
[% END -%]

[% END -%]

#----------------------------------------------------------#
# Clean RepeatMasker
#----------------------------------------------------------#
cd [% data_dir %]
echo Cleaning RepeatMasker

[% FOREACH item IN data -%]
#----------------------------#
# [% item.name %] [% item.coverage %]
#----------------------------#
echo [% item.name %]

[% IF item.skip -%]
echo '    SKIP! [% item.skip %]'
[% ELSE -%]
cd [% item.dir %]
for i in *.fasta;
do
    if [ -f $i.masked ];
    then
        rename 's/fasta.masked$/fa/' $i.masked;
        find [% item.dir %] -type f -name "`basename $i`*" | xargs rm 
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
        File::Spec->catfile( $data_dir, $sh_name )
    ) or die Template->error;

    # 03_strain_info.sh
    $sh_name = "03_strain_info.sh";
    print "Create $sh_name\n";
    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# generate taxon file
#----------------------------#
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
    --file [% data_dir %]/[% group_name %].csv

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
        File::Spec->catfile( $data_dir, $sh_name )
    ) or die Template->error;

    # 04_plan_ALL.sh
    $sh_name = "04_plan_ALL.sh";
    print "Create $sh_name\n";
    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# alignment plan of all genomes 
#----------------------------#
# Don't copy sequences (RepeatMasker done)
# This plan includes all genomes, use the generated phylogenetic tree as guide
#   tree for other plans

# plan_ALL
cd [% data_dir %]
perl [% withncbi %]/taxon/strain_bz.pl \
    --file [% data_dir %]/[% group_name %].csv \
    -w     [% base_dir %] \
    --name [% group_name %] \
    --multi_name plan_ALL \
    --use_name \
    --parallel [% parallel %] \
    --norm \
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
            withncbi   => $withncbi,
            parallel   => $parallel,
        },
        File::Spec->catfile( $data_dir, $sh_name )
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

#----------------------------#
# alignment plan for [% plan_name %]
#----------------------------#
perl [% withncbi %]/taxon/strain_bz.pl \
    --file [% data_dir %]/[% group_name %].csv \
    -w     [% base_dir %] \
    --name [% group_name %] \
    --multi_name [% plan_name %] \
    --use_name \
    --parallel [% parallel %] \
    --norm \
    --phylo_tree [% data_dir %]/plan_ALL_phylo/plan_ALL.nwk \
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
                    parallel   => $parallel,
                    plan_name  => $plan_name,
                    plan       => $plan,
                },
                File::Spec->catfile( $data_dir, $sh_name )
            ) or die Template->error;
        }
    }
}

$stopwatch->end_message;

exit;

__END__

=head1 NAME

pop_prep.pl - prepare pop

=head1 SYNOPSIS

    perl pop_prep.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        -i, --file          input yaml
        --parallel          number of child processes

=cut
