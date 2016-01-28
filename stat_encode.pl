#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Basename;
use File::Find::Rule;
use List::MoreUtils qw(any all uniq zip);
use Set::Scalar;

use DBI;
use Text::CSV_XS;

use AlignDB::Stopwatch;

use FindBin;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

# executable file location
#my $dir_result = "/home/wangq/data/encode/process";
my $dir_result    = "/home/wangq/data/encode/process";
my $term_info     = "/home/wangq/data/encode/cv.ra";
my $raw_stat_file = "encode_raw_stat.csv";
my $term_yml_file = "encode_term.yml";

# filter by filename regex
my $regex;

# filter by item counts
my $count;

# filter by meta info
my ( @dataType, @cell, @lab, @antibody );

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'regex=s'    => \$regex,
    'count=i'    => \$count,
    'dataType=s' => \@dataType,
    'cell=s'     => \@cell,
    'lab=s'      => \@lab,
    'antibody=s' => \@antibody,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# bed info
#----------------------------------------------------------#
if ( !-e $term_yml_file ) {
    print "Gather term info\n";
    open my $in_fh, '<', $term_info;
    my $info_of = {};

    my $content = '';
    while ( my $line = <$in_fh> ) {
        next if $line =~ /^#/;
        if ( $line =~ /^\s+$/ and $content =~ /\S/ ) {
            my @lines = grep {/\S/} split /\n/, $content;
            $content = '';

            my $info = {};
            for (@lines) {
                /^(\w+) (.+)$/;    # term L1210
                next unless $1;
                $info->{$1} = $2;
            }

            next unless $info->{term};
            if ( $info->{deprecated} ) {
                print "Deprecated!\n";
                print Dump \@lines;
                print "\n";
            }

            $info_of->{ $info->{term} } = $info;
        }
        else {
            $content .= $line;
        }
    }
    close $in_fh;

    DumpFile( $term_yml_file, $info_of );
}

my @heads
    = qw{dataType cell cell_tag antibody antibody_tag view itemCount average_size filename};
my $term_info_of = LoadFile($term_yml_file);
if ( !-e $raw_stat_file ) {
    print "Gather raw stat\n";

    my @files_yml
        = sort File::Find::Rule->file->name( '*.yml', '*.yaml' )
        ->in($dir_result);
    printf "\n----Total .yml Files: %4s----\n\n", scalar @files_yml;

    my @ymls;
    for my $file (@files_yml) {
        my $yml = LoadFile($file);
        next unless $yml->{dataType};
        if ( $yml->{dataType} eq "ChIP-seq" ) {
            $yml->{dataType} = "ChipSeq";
        }
        if ( $yml->{cell} ) {
            $yml->{cell_tag} = $term_info_of->{ $yml->{cell} }->{tag};
            if ( !$yml->{cell_tag} ) {
                $yml->{cell_tag} = $yml->{cell};
                $yml->{cell_tag} =~ s/_\(.+\)$//;
            }
            if ( $term_info_of->{ $yml->{cell} }->{deprecated} ) {
                $term_info_of->{ $yml->{cell} }->{deprecated} =~ /same as (.+)/;
                my $cell = $1;
                if ( $term_info_of->{$cell}->{tag} ) {
                    $yml->{cell_tag} = $term_info_of->{$cell}->{tag};
                }
                else {
                    printf "Skip deprecated cell %s\n", $yml->{cell};
                    next;
                }
                printf "Set cell line %s to %s\n", $yml->{cell},
                    $yml->{cell_tag};
            }
        }
        if ( $yml->{antibody} ) {
            $yml->{antibody_tag} = $term_info_of->{ $yml->{antibody} }->{tag};
            if ( !$yml->{antibody_tag} ) {
                $yml->{antibody_tag} = $yml->{antibody};
                $yml->{antibody_tag} =~ s/_\(.+\)$//;
            }
            if ( $term_info_of->{ $yml->{antibody} }->{deprecated} ) {
                printf "Skip deprecated antibody %s\n", $yml->{antibody};
                next
                    ; # mostly uneffective antibody or data withdrawn, except TCF4
            }
        }
        if ( $yml->{dataType} eq "ChipSeq" ) {
            if ( $yml->{antibody_tag} =~ /^H\d/ ) {
                $yml->{dataType} = "Histone";
            }
            else {
                $yml->{dataType} = "TFBS";
            }
        }
        push @ymls, $yml;
    }
    printf "Final ymls: %d\n", scalar @ymls;

    my $csv = Text::CSV_XS->new( { binary => 1 } );
    $csv->eol("\n");

    my @rows;
    for my $yml (@ymls) {
        my @row;
        for my $head (@heads) {
            push @row, $yml->{$head};
        }
        push @rows, [@row];
    }

    open my $fh, ">", $raw_stat_file;
    $csv->print( $fh, $_ ) for @rows;
    close $fh;
}

{
    print "Term stat\n";

    my $dbh = DBI->connect("DBI:CSV:");

    #----------------------------#
    # load tab sep. txt files
    #----------------------------#
    $dbh->{csv_tables}->{t0} = {
        eol            => "\n",
        sep_char       => ",",
        file           => $raw_stat_file,
        skip_first_row => 0,
        col_names      => [@heads],
    };

    {
        my $csv = Text::CSV_XS->new( { binary => 1 } );
        $csv->eol("\n");

        my $query = qq{
            SELECT 
                t0.dataType,
                SUM(t0.itemCount),
                AVG(t0.average_size),
                count(*)
            FROM   t0
            WHERE 1 = 1
            AND t0.itemCount >= 1000
            GROUP BY t0.dataType
        };
        my $sth = $dbh->prepare($query);
        $sth->execute;
        open my $fh, ">", "encode_dataType.csv";
        $csv->print( $fh, [qw{dataType sum_itemCount average_size count}] );
        while ( my @row = $sth->fetchrow_array ) {
            $csv->print( $fh, [@row] );
        }
        close $fh;
    }

    {
        my $csv = Text::CSV_XS->new( { binary => 1 } );
        $csv->eol("\n");

        my $query = qq{
            SELECT 
                t0.cell_tag,
                t0.dataType
            FROM   t0
            WHERE 1 = 1
            AND t0.itemCount >= 1000
            AND t0.dataType IN ('DnaseSeq', 'FaireSeq', 'Histone', 'RepliChip', 'RepliSeq')
            ORDER BY t0.cell_tag, t0.dataType
        };
        my $sth = $dbh->prepare($query);
        $sth->execute;
        open my $fh, ">", "encode_cell_dataType.csv";
        $csv->print( $fh, [qw{cell_tag dataType}] );
        while ( my @row = $sth->fetchrow_array ) {
            $csv->print( $fh, [@row] );
        }
        close $fh;
    }

    {
        my $csv = Text::CSV_XS->new( { binary => 1 } );
        $csv->eol("\n");

        my $query = qq{
            SELECT 
                t0.cell_tag,
                SUM(t0.itemCount),
                AVG(t0.average_size),
                count(*)
            FROM   t0
            WHERE 1 = 1
            AND t0.itemCount >= 1000
            GROUP BY t0.cell_tag
        };
        my $sth = $dbh->prepare($query);
        $sth->execute;
        open my $fh, ">", "encode_cell_tag.csv";
        $csv->print( $fh, [qw{cell_tag sum_itemCount average_size count}] );
        while ( my @row = $sth->fetchrow_array ) {
            $csv->print( $fh, [@row] );
        }
        close $fh;
    }

    {
        my $csv = Text::CSV_XS->new( { binary => 1 } );
        $csv->eol("\n");

        my $query = qq{
            SELECT 
                t0.antibody_tag,
                SUM(t0.itemCount),
                AVG(t0.average_size),
                count(*)
            FROM   t0
            WHERE 1 = 1
            AND t0.antibody_tag IS NOT NULL
            AND t0.itemCount >= 1000
            GROUP BY t0.antibody_tag
        };
        my $sth = $dbh->prepare($query);
        $sth->execute;
        open my $fh, ">", "encode_antibody_tag.csv";
        $csv->print( $fh,
            [qw{antibody_tag sum_itemCount average_size count used}] );
        while ( my @row = $sth->fetchrow_array ) {
            $csv->print( $fh, [@row] );
        }
        close $fh;
    }

    {
        my $csv = Text::CSV_XS->new( { binary => 1 } );
        $csv->eol("\n");

        my $query = qq{
            SELECT 
                t0.cell_tag,
                SUM(t0.itemCount),
                AVG(t0.average_size),
                count(*)
            FROM   t0
            WHERE 1 = 1
            AND t0.dataType = 'TFBS'
            AND t0.itemCount >= 1000
            GROUP BY t0.cell_tag
        };
        my $sth = $dbh->prepare($query);
        $sth->execute;
        open my $fh, ">", "encode_tfbs_cell.csv";
        $csv->print( $fh, [qw{cell_tag sum_itemCount average_size count}] );
        while ( my @row = $sth->fetchrow_array ) {
            $csv->print( $fh, [@row] );
        }
        close $fh;
    }

    #{
    #    my $csv = Text::CSV_XS->new( { binary => 1 } );
    #    $csv->eol("\n");
    #
    #    my ( @cells, @anitbodies );
    #    {
    #        my $query = qq{
    #            SELECT 
    #                t0.cell_tag,
    #                count(*) 
    #            FROM   t0
    #            WHERE 1 = 1
    #            AND t0.dataType = 'TFBS'
    #            AND t0.itemCount >= 1000
    #            GROUP BY t0.cell_tag
    #        };
    #        my $sth = $dbh->prepare($query);
    #        $sth->execute;
    #        while ( my @row = $sth->fetchrow_array ) {
    #            if ( $row[1] >= 3 ) {
    #                push @cells, $row[0];
    #            }
    #        }
    #    }
    #    {
    #        my $query = qq{
    #            SELECT 
    #                t0.antibody_tag,
    #                count(*)
    #            FROM   t0
    #            WHERE 1 = 1
    #            AND t0.dataType = 'TFBS'
    #            AND t0.itemCount >= 1000
    #            GROUP BY t0.antibody_tag
    #        };
    #        my $sth = $dbh->prepare($query);
    #        $sth->execute;
    #        while ( my @row = $sth->fetchrow_array ) {
    #            if ( $row[1] >= 3 ) {
    #                push @anitbodies, $row[0];
    #            }
    #        }
    #    }
    #
    #    {
    #        my $in_cells
    #            = " AND cell_tag IN ("
    #            . join( ",", map { $dbh->quote($_) } @cells ) . ")\n";
    #        my $in_antibodies
    #            = " AND antibody_tag IN ("
    #            . join( ",", map { $dbh->quote($_) } @anitbodies ) . ")\n";
    #        my $query = qq{
    #            SELECT *
    #            FROM   t0
    #            WHERE 1 = 1
    #            AND t0.dataType = 'TFBS'
    #            AND t0.itemCount >= 1000
    #        } . $in_antibodies . $in_cells;# . $in_antibodies;
    #        my $sth = $dbh->prepare($query);
    #        $sth->execute;
    #        open my $fh, ">", "encode_tfbs_cell_antibody.csv";
    #
    #       #$csv->print( $fh, [qw{cell_tag sum_itemCount average_size count}] );
    #        while ( my @row = $sth->fetchrow_array ) {
    #            $csv->print( $fh, [@row] );
    #        }
    #        close $fh;
    #    }
    #}

    {
        my @used = qw{
            ATF3 BATF BCL11A BCL3 BCLAF1 BDP1 BHLHE40 BRCA1 BRF1 BRF2 CCNT2
            CEBPB CHD2 CTBP2 CTCF CTCFL E2F1 E2F4 E2F6 EBF1 EGR1 ELF1 ELK4 EP300
            ESR1 ESRRA ETS1 FOS FOSL1 FOSL2 FOXA1 FOXA2 GABPA GATA1 GATA2 GATA3
            GTF2B GTF2F1 GTF3C2 HDAC2 HDAC8 HMGN3 HNF4A HNF4G HSF1 IRF1 IRF3
            IRF4 JUN JUNB JUND MAFF MAFK MAX MEF2A MEF2C MXI1 MYC NANOG NFE2
            NFKB1 NFYA NFYB NR2C2 NR3C1 NRF1 PAX5 PBX3 POLR2A POLR3A POLR3G
            POU2F2 POU5F1 PPARGC1A PRDM1 RAD21 RDBP REST RFX5 RXRA SETDB1 SIN3A
            SIRT6 SIX5 SMARCA4 SMARCB1 SMARCC1 SMARCC2 SMC3 SP1 SP2 SPI1 SREBF1
            SRF STAT1 STAT2 STAT3 SUZ12 TAF1 TAF7 TAL1 TBP TCF12 TCF7L2 TFAP2A
            TFAP2C THAP1 TRIM28 USF1 USF2 WRNIP1 YY1 ZBTB33 ZBTB7A ZEB1 ZNF143
            ZNF263 ZNF274 ZZZ3
        };
        my $used_set = Set::Scalar->new(@used);
        my $csv = Text::CSV_XS->new( { binary => 1 } );
        $csv->eol("\n");

        my $query = qq{
            SELECT 
                t0.antibody_tag,
                SUM(t0.itemCount),
                AVG(t0.average_size),
                count(*)
            FROM   t0
            WHERE 1 = 1
            AND t0.dataType = 'TFBS'
            AND t0.itemCount >= 1000
            GROUP BY t0.antibody_tag
        };
        my $sth = $dbh->prepare($query);
        $sth->execute;
        open my $fh, ">", "encode_tfbs_antibody.csv";
        $csv->print( $fh,
            [qw{antibody_tag sum_itemCount average_size count used}] );
        while ( my @row = $sth->fetchrow_array ) {
            if ( $used_set->has( $row[0] ) ) {
                $csv->print( $fh, [ @row, "used" ] );
                $used_set->delete( $row[0] );
            }
            else {
                $csv->print( $fh, [@row] );
            }
        }
        close $fh;
        print "Unused " . $used_set->members . "\n";
        print "$_ " for sort $used_set->members;
        print "\n";
    }

    {
        my $csv = Text::CSV_XS->new( { binary => 1 } );
        $csv->eol("\n");

        my $query = qq{
            SELECT 
                t0.antibody_tag,
                SUM(t0.itemCount),
                AVG(t0.average_size),
                count(*)
            FROM   t0
            WHERE 1 = 1
            AND t0.dataType = 'Histone'
            GROUP BY t0.antibody_tag
        };
        my $sth = $dbh->prepare($query);
        $sth->execute;
        open my $fh, ">", "encode_histone_antibody.csv";
        $csv->print( $fh,
            [qw{antibody_tag sum_itemCount average_size count }] );
        while ( my @row = $sth->fetchrow_array ) {
            $csv->print( $fh, [@row] );
        }
        close $fh;
    }

    {
        my $csv = Text::CSV_XS->new( { binary => 1 } );
        $csv->eol("\n");

        my $query = qq{
            SELECT 
                t0.cell_tag,
                SUM(t0.itemCount),
                AVG(t0.average_size),
                count(*)
            FROM   t0
            WHERE 1 = 1
            AND t0.dataType = 'Histone'
            GROUP BY t0.cell_tag
        };
        my $sth = $dbh->prepare($query);
        $sth->execute;
        open my $fh, ">", "encode_histone_cell.csv";
        $csv->print( $fh, [qw{cell_tag sum_itemCount average_size count}] );
        while ( my @row = $sth->fetchrow_array ) {
            $csv->print( $fh, [@row] );
        }
        close $fh;
    }
}

__END__

=head1 SYNOPSIS

    stat_encode.pl [options]
      Options:
        --help            brief help message
        --man             full documentation
        --file            csv file name
        --rep             times of bootstrap simulations

=cut
