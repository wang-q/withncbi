#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Template;
use File::Basename;
use File::Find::Rule;
use File::Spec;
use String::Compare;
use YAML qw(Dump Load DumpFile LoadFile);

use FindBin;
use lib "$FindBin::Bin/../lib";
use MyUtil qw(replace_home);

my $conf_file = "$FindBin::Bin/trichoderma_data.yml";
my $yml       = LoadFile($conf_file);

my $group_name = $yml->{group_name};
my $base_dir   = replace_home( $yml->{base_dir} );
my $data_dir   = replace_home( $yml->{data_dir} );
my $pl_dir     = replace_home( $yml->{pl_dir} );
my $parallel   = $yml->{parallel} || 4;

# NCBI WGS
my $fasta_dir = replace_home( $yml->{fasta_dir} );

my @data = @{ $yml->{data} };

my @subdirs_fasta = File::Find::Rule->file->name('*.fsa_nt.gz')->in($fasta_dir);
for my $item (@data) {
    printf "Matching [%s]\n", $item->{name};
    if ( exists $item->{skip} ) {
        printf " " x 4 . "SKIP! %s\n", $item->{skip};
        next;
    }

    # match the most similar name
    my ($fasta) = map { $_->[0] }
        sort { $b->[1] <=> $a->[1] }
        map { [ $_, compare( basename($_), $item->{prefix} ) ] } @subdirs_fasta;
    $item->{fasta} = $fasta;

    # prepare working dir
    my $dir = File::Spec->catdir( $data_dir, $item->{name} );
    mkdir $dir if !-e $dir;
    $item->{dir} = $dir;
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
[% ELSE -%]
cd [% item.dir %]
gzip -d -c [% item.fasta %] > toplevel.fa
perl -p -i -e '/>/ and s/\>gi\|(\d+).*/\>gi_$1/' toplevel.fa
faops count toplevel.fa | perl -aln -e 'next if $F[0] eq 'total'; print $F[0] if $F[1] > 100000; print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05' | uniq > listFile
faops some toplevel.fa listFile toplevel.filtered.fa
[% IF item.per_seq -%]
faops split-name toplevel.filtered.fa .
[% ELSE -%]
faops split-about toplevel.filtered.fa 10000000 .
[% END -%]
rm toplevel.fa toplevel.filtered.fa listFile

rename 's/fa$/fasta/' *.fa;
[% END -%]

[% END -%]

EOF

    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
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
echo 'SKIP! [% item.skip %]'
[% ELSE -%]
cd [% item.dir %]
RepeatMasker [% item.dir %]/*.fasta -species Fungi -xsmall --parallel [% parallel %]
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
echo 'SKIP! [% item.skip %]'
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
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
            parallel => $parallel,
        },
        File::Spec->catfile( $data_dir, $sh_name )
    ) or die Template->error;

    # 03_prepare.sh
    $sh_name = "03_prepare.sh";
    print "Create $sh_name\n";
    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# generate taxon file
#----------------------------#
perl [% pl_dir %]/withncbi/taxon/strain_info.pl \
[% FOREACH item IN data -%]
    --id   [% item.taxon %] \
    --name [% item.taxon %]=[% item.name %] \
[% IF item.original_id -%]
    --species [% item.taxon %]=[% item.original_id %] \
[% END -%]
[% END -%]
    --file [% data_dir %]/[% group_name %].csv

#----------------------------#
# multi genome alignment plan
#----------------------------#
# don't copy sequences (RepeatMasker done)
# Execute the following lines by copy & paste.

# Trichoderma_10way
cd [% data_dir %]
perl [% pl_dir %]/withncbi/taxon/strain_bz.pl \
    --file [% data_dir %]/[% group_name %].csv \
    -w     [% base_dir %] \
    --name [% group_name %] \
    --multi_name Trichoderma_10way \
    --use_name \
    --parallel [% parallel %] \
    --norm \
    -t Tatr_IMI_2206040 \
    -q Tatr_XS215 \
    -q Thar_B05 \
    -q Thar_T6776 \
    -q Tlon_SMF2 \
    -q Tpar \
    -q Tree_QM6a \
    -q Tree_RUT_C_30 \
    -q Tvir_FT_333 \
    -q Tvir_Gv29_8

EOF

    $tt->process(
        \$text,
        {   data       => \@data,
            group_name => $group_name,
            data_dir   => $data_dir,
            base_dir   => $base_dir,
            pl_dir     => $pl_dir,
            parallel   => $parallel,
        },
        File::Spec->catfile( $data_dir, $sh_name )
    ) or die Template->error;

}

__END__

# create pop/trichoderma.tsv manually, be careful with tabs and spaces.
# http://www.ncbi.nlm.nih.gov/Traces/wgs/?page=1&term=trichoderma

mkdir -p ~/data/alignment/trichoderma
cd ~/data/alignment/trichoderma

perl ~/Scripts/withncbi/util/wgs_prep.pl \
    -f ~/Scripts/withncbi/pop/trichoderma.tsv \
    --fix \
    -o WGS \
    -a 

aria2c -x 6 -s 3 -c -i WGS/trichoderma.url.txt

find WGS -name "*.gz" | xargs gzip -t 

# rsync --progress -av wangq@139.162.23.84:/home/wangq/data/alignment/trichoderma/ ~/data/alignment/trichoderma

# Add some contents to WGS/trichoderma.data.yml, get pop/trichoderma_data.yml

perl ~/Scripts/withncbi/pop/trichoderma.pl
sh 01_file.sh
sh 02_rm.sh

# execute 03_prepare.sh by copy & paste  

# for each multi_name, execute the following bash file
sh 1_real_chr.sh
sh 3_pair_cmd.sh
sh 4_rawphylo.sh
sh 5_multi_cmd.sh
sh 6_var_list.sh
sh 7_multi_db_only.sh
sh 9_pack_it_up.sh
