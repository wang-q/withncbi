#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use File::Basename;

use AlignDB::IntSpan;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $op = "merge_to_runlist";    # or bed_diff, bed_intersect

my @files;
my $filename;
my $remove_chr;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'   => \$help,
    'man'      => \$man,
    'o|op=s'   => \$op,
    'f|file=s' => \@files,
    'n|name=s' => \$filename,
    'r|remove' => \$remove_chr,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#

if ( $op eq "merge_to_runlist" ) {
    if ( !$filename ) {
        $filename = "$op.yml";
    }

    my $chr_data_set = {};
    for my $file (@files) {
        print "Target File [$file]\n";
        load_bed( $file, $chr_data_set, $remove_chr );
    }

    print "Transform AlignDB::IntSpan objects to runlists\n";
    for my $key ( keys %{$chr_data_set} ) {
        $chr_data_set->{$key} = $chr_data_set->{$key}->runlist;
    }

    print "Start dumping $filename\n";
    DumpFile( $filename, $chr_data_set );
}
elsif ( $op eq "bed_diff" or $op eq "bed_intersect" ) {
    if ( !$filename ) {
        die "Need a yaml or bed file for --name\n";
    }

    my $chr_data_set = {};
    if ( $filename =~ /\.ya?ml$/ ) {
        print "Op File [$filename]\n";
        $chr_data_set = load_yaml($filename);
    }
    elsif ( $filename =~ /\.bed$/ ) {
        print "Op File [$filename]\n";
        load_bed( $filename, $chr_data_set, $remove_chr );
    }
    else {
        die "Not a yaml or bed file\n";
    }

    for my $file (@files) {
        print "Target File [$file]\n";
        write_filtered_bed( $file, $op, $chr_data_set, $remove_chr );
    }
}

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#
sub read_bed {
    my $file       = shift;
    my $remove_chr = shift;

    my @data;
    open my $fh, '<', $file;
    while ( my $string = <$fh> ) {
        next unless defined $string;
        chomp $string;
        my ( $chr, $start, $end )
            = ( split /\t/, $string )[ 0, 1, 2 ];
        next unless $chr =~ /^\w+$/;
        $chr =~ s/chr0?//i if $remove_chr;
        next unless $start =~ /^\d+$/;
        next unless $end =~ /^\d+$/;
        if ( $start > $end ) {
            ( $start, $end ) = ( $end, $start );
        }
        next if $end - $start < 10;
        my $set = AlignDB::IntSpan->new("$start-$end");
        push @data, { chr => $chr, set => $set, };
    }
    close $fh;

    print " " x 4, scalar @data, " records\n";
    return \@data;
}

sub write_bed {
    my $file = shift;
    my $data = shift;

    open my $fh, '>', $file;
    for my $item ( @{$data} ) {
        print {$fh} $item->{chr}, "\t", $item->{set}->min, "\t",
            $item->{set}->max,
            "\n";
    }
    close $fh;

    return;
}

sub load_bed {
    my $filename     = shift;
    my $chr_data_set = shift;    # need be a hashref
    my $remove_chr   = shift;

    my $data = read_bed( $filename, $remove_chr );
    for my $item ( @{$data} ) {
        my $chr = $item->{chr};
        my $set = $item->{set};
        if ( !exists $chr_data_set->{$chr} ) {
            $chr_data_set->{$chr} = AlignDB::IntSpan->new;
        }
        $chr_data_set->{$chr}->merge($set);
    }

    return;
}

sub load_yaml {
    my $filename = shift;

    my $chr_data_set = {};
    $chr_data_set = LoadFile($filename);
    for my $key ( sort keys %{$chr_data_set} ) {

        # compatible with Object dump or runlists
        $chr_data_set->{$key}
            = AlignDB::IntSpan->new( $chr_data_set->{$key} );
    }

    return $chr_data_set;
}

sub write_filtered_bed {
    my $filename     = shift;
    my $op           = shift;
    my $chr_data_set = shift;
    my $remove_chr   = shift;

    my $data = read_bed( $filename, $remove_chr );

    my @filtered;
    for my $item ( @{$data} ) {
        my $chr = $item->{chr};
        my $set = $item->{set};

        if ( $op eq "bed_diff" ) {
            if ( $chr_data_set->{$chr}->intersect($set)->is_empty ) {
                push @filtered, $item;
            }
        }
        elsif ( $op eq "bed_intersect" ) {
            if ( $chr_data_set->{$chr}->intersect($set)->is_not_empty ) {
                push @filtered, $item;
            }
        }
    }

    my $basename = basename($filename);
    my $outfile  = "$basename.$op";
    print "Write [$outfile]\n";
    write_bed( $outfile, \@filtered );

    return;
}
