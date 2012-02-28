#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use WWW::Mechanize;
use Regexp::Common qw(balanced);
use List::MoreUtils qw(uniq zip);
use HTML::TableExtract;

my $srs_worker = sub {
    my $term = shift;

    my $mech = WWW::Mechanize->new;
    $mech->stack_depth(0);    # no history to save memory

    my $url_part1 = "http://www.ncbi.nlm.nih.gov/biosample/?term=";
    my $url_part2 = "&from=begin&to=end&dispmax=200";
    my $url       = $url_part1 . $term . $url_part2;
    print $url, "\n";

    my @srx;

    $mech->get($url);

    # this link exists in both summary and detailed pages
    $mech->follow_link(
        text_regex => qr{SRS\d+},
        url_regex  => => qr{sample},
    );

    {
        my @links = $mech->find_all_links(
            text_regex => qr{SRX\d+},
            url_regex  => qr{report},
        );

        printf "OK, get %d SRX\n", scalar @links;
        @srx = map { $_->text } @links;
    }

    return \@srx;
};

my $srx_worker = sub {
    my $term = shift;

    my $mech = WWW::Mechanize->new;
    $mech->stack_depth(0);    # no history to save memory

    my $url_part = "http://www.ncbi.nlm.nih.gov/sra?term=";
    my $url      = $url_part . $term;
    print $url, "\n";

    my $info = {
        sample   => "",
        library  => "",
        platform => "",
        layout   => "",
    };

    $mech->get($url);

    my $page = $mech->content;
    my $te   = HTML::TableExtract->new;
    $te->parse($page);

    my %srr_info;
    for my $ts ( $te->table_states ) {
        for my $row ( $ts->rows ) {
            for my $cell (@$row) {
                $cell =~ s/,//g;
                $cell =~ s/\s+//g;
            }
            next unless $row->[0] =~ /\d+/;
            $srr_info{ $row->[1] } = { spot => $row->[2], base => $row->[3] };
        }
    }
    $info->{srr_info} = \%srr_info;

    {
        my @links = $mech->find_all_links(
            text      => "sra",
            url_regex => => qr{ftp},
        );
        $info->{ftp_base} = $links[0]->url;
        ( $info->{srx} ) = reverse grep {$_} split /\//, $info->{ftp_base};
    }

    {
        my ( @srr, @downloads );
        my @links = $mech->find_all_links( text_regex => qr{SRR}, );
        printf "OK, get %d SRR\n", scalar @links;

        @srr       = map { $_->text } @links;
        @downloads = map { $info->{ftp_base} . "/$_/$_.sra" } @srr;

        $info->{srr}       = \@srr;
        $info->{downloads} = \@downloads;
    }

    {
        my @links = $mech->find_all_links(
            text      => "Study",
            url_regex => => qr{study},
        );
        ( $info->{srp} ) = reverse grep {$_} split /\=/, $links[0]->url;
    }

    {
        my @links = $mech->find_all_links(
            text_regex => qr{SRS},
            url_regex  => => qr{sample},
        );
        $info->{srs} = $links[0]->text;
    }

    {
        my $content = $mech->content;
        $content =~ s/^.+Accession\://s;
        $content =~ s/Download reads.+$//s;
        $content =~ s/$RE{balanced}{-parens=>'<>'}/ /g;
        $content =~ s/$RE{balanced}{-parens=>'()'}/\n/g;
        $content =~ s/ +/ /g;
        $content =~ s/\n+/\n/g;
        $content =~ s/\s{2,}/\n/g;
        my @lines = grep {$_} split /\n/, $content;

        while (@lines) {
            my $line = shift @lines;
            if ( $line =~ /(sample|library|platform)\:\s+(.+)$/i ) {
                $info->{ lc $1 } = $2;
            }
            if ( $line =~ /(layout)\:\s+(\w+)/i ) {
                $info->{ lc $1 } = $2;
            }
        }

        #DumpFile(
        #    "$term.yml",
        #    {   content => $content,
        #    }
        #);
    }

    return $info;
};

my @ids = (
    391, 517, 375, 321, 40,  852, 176, 57,  380, 208, 38,  138, 727, 443,
    757, 897, 738, 332, 181, 406, 177, 320, 737, 392, 357, 492, 377, 502,
    508, 381, 426, 370, 373, 892, 405, 857, 359, 837, 461, 908, 42,  142,
    555, 383, 513, 887, 855, 441, 235, 75,  531, 707, 882, 439, 907, 714,
    790, 822, 799, 440, 639, 427, 705, 894, 83,  765, 732, 812, 237, 712,
    805, 890, 21,  365, 49,  399, 59,  810, 358, 776, 352, 787, 509, 318,
    796, 535, 761, 804, 26,  491, 818, 85,  195, 93,  45,  356, 859, 820,
    239, 73,  136, 28,  280, 256, 563, 161, 703, 149, 642, 832, 69,  730,
    911, 91,  105, 783, 88,  41,  217, 374, 884, 309, 589, 367, 595, 310,
    591, 721, 716, 233, 646, 808, 101, 861, 228, 227, 313, 325, 158, 801,
    786, 879, 379, 208, 371, 386, 338, 301, 517, 427, 802, 437, 315, 129,
    774, 307, 765, 229, 391, 109, 350, 639, 313, 379, 707, 409, 852, 324,
    153, 362, 714, 317, 555, 365, 730, 712, 786, 859, 375, 303, 705, 799,
    774, 820, 287, 335, 358, 360, 399, 732, 357, 486, 437, 304, 362, 380,
    306,
);
@ids = uniq(@ids);

my $master = {};
for my $id (@ids) {
    my $name = "DGRP-$id";
    print "$name\n";
    my @srx = @{ $srs_worker->($name) };
    print "@srx\n";

    my $sample = {};
    for (@srx) {
        $sample->{$_} = $srx_worker->($_);
    }
    $master->{$name} = $sample;
    print "\n";
}

DumpFile( "DGRP.yml", $master );
