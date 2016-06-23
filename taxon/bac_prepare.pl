#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use File::Find::Rule;
use Path::Tiny;
use Text::Table;
use List::MoreUtils qw(uniq);
use Archive::Extract;

use DBI;
use Template;

use Bio::Taxon;
use Bio::DB::Taxonomy;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/../lib";
use MyUtil qw(find_ancestor);

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

bac_prepare.pl

=head1 SYNOPSIS

    perl bac_prepare.pl --base_dir d:/bacteria/bacteria_101015 --parent 562

=cut

# paths
my $td_dir  = path( $Config->{path}{td} )->stringify;   # taxdmp
my $nb_dir  = path( $Config->{path}{nb} )->stringify;   # NCBI genomes bac
my $nbd_dir = path( $Config->{path}{nbd} )->stringify;  # NCBI genomes bac draft
my $ngbd_dir
    = path( $Config->{path}{ngbd} )->stringify; # NCBI genbank genomes bac draft

my $egaz = path( $Config->{run}{egaz} )->stringify;    # egaz path

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'server|s=s'    => \( my $server      = $Config->{database}{server} ),
    'port=i'        => \( my $port        = $Config->{database}{port} ),
    'db|d=s'        => \( my $db_name     = $Config->{database}{db} ),
    'username|u=s'  => \( my $username    = $Config->{database}{username} ),
    'password|p=s'  => \( my $password    = $Config->{database}{password} ),
    'seq_dir=s'     => \( my $seq_dir     = "~/data/bacteria/bac_seq_dir" ),
    'working_dir=s' => \( my $working_dir = "." ),
    'parent|p=s' => \( my $parent = "562,585054" ),   # E.coli and E. fergusonii
    'target|t=i' => \my $target,
    'outgroup|o=i' => \my $outgroup,
    'exclude|e=s'  => \( my $exclude = '0' ),
    'name_str|n=s' => \
        my $name_str
    , # use custom name_str, working dir and goal db name. mysql restrict db name length 64
    'strains' => \my $require_strains,    # skip taxonomy_id == species_id
    'is_self' => \my $is_self,            # is self alignment (paralog)
    'length=i' => \( my $paralog_length = 1000 ),    # paralog length
    'get_seq'  => \
        my $get_seq,    # download sequences via get_seq.pl if not existing
    'scaffold' => \my $scaffold,    # including scaffolds and contigs
    'parallel=i' => \( my $parallel = $Config->{run}{parallel} ),
) or Getopt::Long::HelpMessage(1);

$seq_dir = path($seq_dir)->stringify;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Preparing whole group...");

# Database handler
my $dsn
    = "dbi:mysql:"
    . "database="
    . $db_name
    . ";host="
    . $server
    . ";port="
    . $port;
my DBI $dbh = DBI->connect( $dsn, $username, $password )
    or die $DBI::errstr;

my $id_str;
{
    $stopwatch->block_message("Load taxdmp and expand --parent");

    my $taxon_db = Bio::DB::Taxonomy->new(
        -source    => 'flatfile',
        -directory => $td_dir,
        -nodesfile => "$td_dir/nodes.dmp",
        -namesfile => "$td_dir/names.dmp",
    );
    my @parents = split /,/, $parent;

    my $sub_id_set = AlignDB::IntSpan->new;
    for my $p_id (@parents) {
        $sub_id_set->add($p_id);
        my $p_ids = $taxon_db->get_taxon( -taxonid => $p_id );

        my @taxa = $taxon_db->get_all_Descendents($p_ids);
        for my $taxon (@taxa) {
            $sub_id_set->add( $taxon->id );
        }
    }

    my $db_id_set = AlignDB::IntSpan->new;
    {
        my $query
            = $scaffold
            ? q{ SELECT taxonomy_id FROM gr WHERE 1 = 1  }
            : q{ SELECT taxonomy_id FROM gr WHERE status NOT IN ('Contig', 'Scaffold') };
        $query .= "AND taxonomy_id != species_id" if $require_strains;
        my DBI $sth = $dbh->prepare($query);
        $sth->execute;
        while ( my ($id) = $sth->fetchrow_array ) {
            $db_id_set->add($id);
        }
    }

    my $id_set = $sub_id_set->intersect($db_id_set);
    $id_set->remove( split /,/, $exclude );
    $id_str = '(' . ( join ",", $id_set->as_array ) . ')';

    die "Wrong id_str $id_str\n" unless $id_str =~ /\d+/;
}

{
    $stopwatch->block_message("making working dir");

    if ( !$name_str ) {
        my $query = qq{
            SELECT DISTINCT species
            FROM gr
            WHERE taxonomy_id IN $id_str
        };
        my DBI $sth = $dbh->prepare($query);
        $sth->execute;

        while ( my ($name) = $sth->fetchrow_array ) {
            $name_str .= "_$name";
        }
        $name_str =~ s/\W/_/g;
        $name_str =~ s/^_+//g;
        $name_str =~ s/\s+/_/g;
    }

    print "Working on $name_str\n";
    $working_dir = path( $working_dir, $name_str )->absolute;
    if ( !-d $working_dir ) {
        $working_dir->mkpath;
    }
    $working_dir = $working_dir->stringify;
    print "Working dir is $working_dir\n";
}

my @query_ids;
{
    $stopwatch->block_message("find all strains' taxon ids");

    # select all strains in this species
    my $query = qq{
        SELECT taxonomy_id, organism_name, released_date, status, code
        FROM gr
        WHERE taxonomy_id IN $id_str
        ORDER BY released_date, status, code
    };
    my DBI $sth = $dbh->prepare($query);
    $sth->execute;

    # header line
    my @strains;
    my $table = Text::Table->new( @{ $sth->{NAME} } );
    while ( my @row = $sth->fetchrow_array ) {
        push @strains, [@row];
        $table->load( [@row] );
    }

    my $fh = path( $working_dir, "table.txt" )->openw;
    print {$fh} $table, "\n";
    print $table, "\n";

    {
        my $message = "There are " . scalar @strains . " strains\n";
        print {$fh} $message;
        print $message;
    }

    if ($target) {
        my ($exist) = grep { $_->[0] == $target } @strains;
        if ( defined $exist ) {
            my $message
                = "Use [$exist->[1] ($exist->[0])] as target, as you wish.\n";
            print {$fh} $message;
            print $message;
        }
        else {
            print "Taxon [$target] doesn't exist, please check.\n";
            exit;
        }
    }
    else {
        $target = $strains[0]->[0];
        my $message
            = "Use [$strains[0]->[1] ($strains[0]->[0])] as target, the oldest strain on NCBI.\n";
        print {$fh} $message;
        print $message;
    }

    @query_ids = map { $_->[0] == $target ? () : $_->[0] } @strains;

    if ($outgroup) {
        my ($exist) = grep { $_ == $outgroup } @query_ids;
        if ( defined $exist ) {
            my $message = "Use [$exist] as outgroup, as you wish.\n";
            print {$fh} $message;
            print $message;
        }
        else {
            print "Taxon [$outgroup] doesn't exist, please check.\n";
        }
    }

    print "\n\n";
    print {$fh} "perl " . $stopwatch->cmd_line, "\n";

    close $fh;
}

my %id_missing_file;
my @ids_missing;
{
    $stopwatch->block_message("build fasta files");

    # read all filenames, then grep
    print "Reading file list\n";
    my @fna_files = File::Find::Rule->file->name('*.fna')->in($nb_dir);
    my @gff_files = File::Find::Rule->file->name('*.gff')->in($nb_dir);
    my ( @scaff_files, @contig_files );
    if ($scaffold) {
        @scaff_files
            = File::Find::Rule->file->name('*.scaffold.fna.tgz')->in($nbd_dir);
        @contig_files
            = File::Find::Rule->file->name('*.contig.fna.tgz')->in($ngbd_dir);
    }

    print "Rewrite seqs for every strains\n";
ID: for my $taxon_id ( $target, @query_ids ) {
        print "taxon_id $taxon_id\n";
        my $id_dir = path( $seq_dir, $taxon_id );
        if ( !-e $id_dir ) {
            $id_dir->mkpath;
        }

        my @accs;    # complete accessions
        {
            my DBI $sth = $dbh->prepare(
                q{
                SELECT chr FROM gr WHERE taxonomy_id = ?
                }
            );
            $sth->execute($taxon_id);
            my ($acc) = $sth->fetchrow_array;
            push @accs,
                ( map { s/\.\d+$//; $_ } grep {defined} ( split /\|/, $acc ) );
        }

      # for NZ_CM***, NZ_CP*** accessions, the following codes will find nothing
        for my $acc ( grep {defined} @accs ) {
            my ($fna_file) = grep {/$acc/} @fna_files;
            if ( !defined $fna_file ) {
                if ($get_seq) {
                    print
                        "Download from entrez. id: [$taxon_id]\tseq: [$acc]\n";
                    if ( -e "$id_dir/$acc.gb" ) {
                        print "Sequence [$id_dir/$acc.gb] exists, next\n";
                        next;
                    }

                    # download
                    system "perl $FindBin::Bin/get_seq.pl $acc $id_dir";
                }
                else {
                    warn ".fna file for [$acc] of [$taxon_id] doesn't exist.\n";
                    $id_missing_file{$taxon_id}++;
                    push @ids_missing, $taxon_id;
                    next ID;
                }
            }
            else {
                path($fna_file)->copy($id_dir);

                my ($gff_file) = grep {/$acc/} @gff_files;
                path($gff_file)->copy($id_dir);
            }
        }

        # prep_scaff() will find the scaffolds for NZ_CM***, NZ_CP***
        if ($scaffold) {
            my ($wgs) = get_taxon_wgs($taxon_id);

            next unless $wgs;

            $wgs =~ s/\d+$//;
            my $rc = prep_scaff( \@scaff_files, "NZ_$wgs", $id_dir );
            if ($rc) {
                warn " " x 4, $rc;
                print " " x 4, "Try contig files\n";
                my $rc2 = prep_scaff( \@contig_files, $wgs, $id_dir );
                warn " " x 4, $rc2 if $rc2;
            }
        }
    }
}

# report missing
if (@ids_missing) {
    @query_ids = grep { !$id_missing_file{$_} } @query_ids;

    my $fh = path( $working_dir, "table.txt" )->opena;
    print {$fh} "Can't find files for the following ID:\n@ids_missing\n";
    close $fh;
}

{
    my $tt = Template->new;
    my $text;

    # prepare.sh
    print "Create prepare.sh\n";
    $text = <<'EOF';
#!/bin/bash

cd [% working_dir %]

[% IF ! redo -%]
perl [% findbin %]/strain_info.pl \
    --file [% working_dir %]/info.csv \
[% FOREACH id IN query_ids -%]
    --id [% id %] \
[% END -%]
    --id [% target %]
[% END -%]

[% IF ! is_self -%]
perl [% egaz %]/multi_batch.pl \
    --file [% working_dir %]/info.csv \
    -w [% working_dir %]/.. \
[% IF ! redo -%]
    --seq_dir [% seq_dir %] \
[% END -%]
    --name [% name_str %] \
    --parallel [% parallel %] \
[% IF outgroup -%]
    -o [% outgroup %] \
[% END -%]
[% FOREACH id IN query_ids -%]
    -q [% id %] \
[% END -%]
    -t [% target %]

[% ELSE -%]
perl [% egaz %]/self_batch.pl \
    --file [% working_dir %]/info.csv \
    -w [% working_dir %]/.. \
[% IF ! redo -%]
    --seq_dir [% seq_dir %] \
[% END -%]
    --name [% name_str %] \
    --parallel [% parallel %] \
    --length [% length %] \
[% FOREACH id IN query_ids -%]
    -q [% id %] \
[% END -%]
    -t [% target %]

[% END -%]

EOF
    $tt->process(
        \$text,
        {   findbin     => $FindBin::RealBin,
            egaz        => $egaz,
            working_dir => $working_dir,
            seq_dir     => $seq_dir,
            name_str    => $name_str,
            parallel    => $parallel,
            target      => $target,
            query_ids   => \@query_ids,
            outgroup    => $outgroup,
            is_self     => $is_self,
            length      => $paralog_length,
        },
        path( $working_dir, "prepare.sh" )->stringify
    ) or die Template->error;

    print "Create redo_prepare.sh\n";
    $tt->process(
        \$text,
        {   findbin     => $FindBin::RealBin,
            egaz        => $egaz,
            working_dir => $working_dir,
            seq_dir     => $seq_dir,
            name_str    => $name_str,
            parallel    => $parallel,
            target      => $target,
            query_ids   => \@query_ids,
            outgroup    => $outgroup,
            is_self     => $is_self,
            length      => $paralog_length,
            redo        => 1,                # If RepeatMasker has been executed
                                             # don't pass $seq_dir
                                             # don't gather taxon info
        },
        path( $working_dir, "redo_prepare.sh" )->stringify
    ) or die Template->error;
}

$stopwatch->end_message;
exit;

sub get_taxon_wgs {
    my $taxon_id = shift;

    my DBI $sth = $dbh->prepare(
        q{
        SELECT wgs FROM gr WHERE taxonomy_id = ?
        }
    );
    $sth->execute($taxon_id);
    my ($wgs) = $sth->fetchrow_array;

    return $wgs;
}

sub prep_scaff {
    my $all_files = shift;
    my $wgs       = shift;
    my $dir       = shift;

    my ($wgs_file) = grep {/$wgs/} @{$all_files};
    if ( !$wgs_file ) {
        return "Can't find fasta file for $wgs\n";
    }

    my $ae = Archive::Extract->new( archive => $wgs_file );
    my $ok = $ae->extract( to => $dir );

    if ( !$ok ) {
        return $ae->error;
    }

    my (@files) = map { path($_)->absolute($dir)->stringify } @{ $ae->files };

    for my $file (@files) {
        unless ( -e $file ) {
            return "$file not exists!\n";
        }

        # file size less than 1k
        if ( ( stat($file) )[7] < 1024 ) {
            next;
        }

        my $basename = path($file)->basename(".fna");
        my $fa_file = path( $dir, "$basename.fa" );
        path($file)->copy($fa_file);
    }

    unlink $_ for @files;

    return;
}

__END__
