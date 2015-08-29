#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

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

use FindBin;
use lib "$FindBin::Bin/../lib";
use MyUtil qw(replace_home find_ancestor);

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

# running options
my $seq_dir     = "~/data/bacteria/bac_seq_dir";
my $working_dir = ".";

my $parent_id = "562,585054";    # E.coli and E. fergusonii
my $target_id;
my $outgroup_id;
my $exclude_ids = '0';

# is self alignment (paralog)
my $is_self;

# paralog length
my $paralog_length = 1000;

# use custom name_str
# working dir and goal db name
# mysql restrict db name length 64
my $name_str;

# Database init values
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db_name  = $Config->{database}{db};

# download sequences via get_seq.pl if not existing
my $get_seq;

# including scaffolds and contigs
my $scaffold;

# paths
my $td_dir  = replace_home( $Config->{path}{td} );     # taxdmp
my $nb_dir  = replace_home( $Config->{path}{nb} );     # NCBI genomes bac
my $nbd_dir = replace_home( $Config->{path}{nbd} );    # NCBI genomes bac draft
my $ngbd_dir
    = replace_home( $Config->{path}{ngbd} );           # NCBI genbank genomes bac draft

# run in parallel mode
my $parallel = $Config->{run}{parallel};

my $man  = 0;
my $help = 0;

GetOptions(
    'help'           => \$help,
    'man'            => \$man,
    's|server=s'     => \$server,
    'P|port=i'       => \$port,
    'u|username=s'   => \$username,
    'p|password=s'   => \$password,
    'd|db=s'         => \$db_name,
    'seq_dir=s'      => \$seq_dir,
    'working_dir=s'  => \$working_dir,
    'p|parent_id=s'  => \$parent_id,
    't|target_id=i'  => \$target_id,
    'o|r|outgroup=i' => \$outgroup_id,
    'e|exclude=s'    => \$exclude_ids,
    'n|name_str=s'   => \$name_str,
    'is_self'        => \$is_self,
    'length=i'       => \$paralog_length,
    'get_seq'        => \$get_seq,
    'scaffold'       => \$scaffold,
    'parallel=i'     => \$parallel,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

$seq_dir = replace_home($seq_dir);

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Preparing whole group...");

my $dbh = DBI->connect( "dbi:mysql:$db_name:$server", $username, $password );

my $id_str;
{
    $stopwatch->block_message("Load taxdmp and expand --parent_id");

    my $taxon_db = Bio::DB::Taxonomy->new(
        -source    => 'flatfile',
        -directory => $td_dir,
        -nodesfile => "$td_dir/nodes.dmp",
        -namesfile => "$td_dir/names.dmp",
    );
    my @parent_ids = split /,/, $parent_id;

    my $sub_id_set = AlignDB::IntSpan->new;
    for my $p_id (@parent_ids) {
        $sub_id_set->add($p_id);
        my $parent = $taxon_db->get_taxon( -taxonid => $p_id );

        my @taxa = $taxon_db->get_all_Descendents($parent);
        for my $taxon (@taxa) {
            $sub_id_set->add( $taxon->id );
        }
    }

    my $db_id_set = AlignDB::IntSpan->new;
    {
        my $query
            = $scaffold
            ? q{ SELECT taxonomy_id FROM gr WHERE 1 = 1 }
            : q{ SELECT taxonomy_id FROM gr WHERE status NOT IN ('Contig', 'Scaffold') };
        my $sth = $dbh->prepare($query);
        $sth->execute;
        while ( my ($id) = $sth->fetchrow_array ) {
            $db_id_set->add($id);
        }
    }

    my $id_set = $sub_id_set->intersect($db_id_set);
    $id_set->remove( split /,/, $exclude_ids );
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
        my $sth = $dbh->prepare($query);
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
    my $sth = $dbh->prepare($query);
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

    if ($target_id) {
        my ($exist) = grep { $_->[0] == $target_id } @strains;
        if ( defined $exist ) {
            my $message = "Use [$exist->[1]] as target, as you wish.\n";
            print {$fh} $message;
            print $message;
        }
        else {
            print "Taxon $target_id doesn't exist, please check.\n";
            exit;
        }
    }
    else {
        $target_id = $strains[0]->[0];
        my $message = "Use [$strains[0]->[1]] as target, the oldest strain on NCBI.\n";
        print {$fh} $message;
        print $message;
    }

    @query_ids = map { $_->[0] == $target_id ? () : $_->[0] } @strains;

    if ($outgroup_id) {
        my ($exist) = grep { $_ == $outgroup_id } @query_ids;
        if ( defined $exist ) {
            my $message = "Use [$exist] as reference, as you wish.\n";
            print {$fh} $message;
            print $message;

            # make $outgroup_id first
            @query_ids = map { $_ == $outgroup_id ? () : $_ } @query_ids;
            unshift @query_ids, $outgroup_id;
        }
        else {
            print "Taxon $outgroup_id doesn't exist, please check.\n";
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
        @scaff_files  = File::Find::Rule->file->name('*.scaffold.fna.tgz')->in($nbd_dir);
        @contig_files = File::Find::Rule->file->name('*.contig.fna.tgz')->in($ngbd_dir);
    }

    print "Rewrite seqs for every strains\n";
ID: for my $taxon_id ( $target_id, @query_ids ) {
        print "taxon_id $taxon_id\n";
        my $id_dir = path( $seq_dir, $taxon_id );
        if ( !-e $id_dir ) {
            $id_dir->mkpath;
        }

        my @accs;    # complete accessions

        my $query = qq{ SELECT chr FROM gr WHERE taxonomy_id = ? };
        my $sth   = $dbh->prepare($query);
        $sth->execute($taxon_id);
        my ($acc) = $sth->fetchrow_array;
        push @accs, ( map { s/\.\d+$//; $_ } grep {defined} ( split /,/, $acc ) );

        # for NZ_CM*** accessions, the following prep_fa() will find nothing
        # AND is $scaffold, prep_scaff() will find the scaffolds
        for my $acc ( grep {defined} @accs ) {
            my ($fna_file) = grep {/$acc/} @fna_files;
            if ( !defined $fna_file ) {
                if ($get_seq) {
                    print "Download from entrez. id: [$taxon_id]\tseq: [$acc]\n";
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

        if ($scaffold) {
            my ($wgs) = get_taxon_wgs( $dbh, $taxon_id );

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
    --id [% target_id %]
[% END -%]

[% IF ! is_self -%]
perl [% findbin %]/strain_bz.pl \
    --file [% working_dir %]/info.csv \
    -w [% working_dir %]/.. \
[% IF ! redo -%]
    --seq_dir [% seq_dir %] \
[% END -%]
    --name [% name_str %] \
[% FOREACH id IN query_ids -%]
    -q [% id %] \
[% END -%]
    -t [% target_id %]

[% ELSE -%]
perl [% findbin %]/strain_bz_self.pl \
    --file [% working_dir %]/info.csv \
    -w [% working_dir %]/.. \
[% IF ! redo -%]
    --seq_dir [% seq_dir %] \
[% END -%]
    --name [% name_str %] \
    --length [% length %] \
[% FOREACH id IN query_ids -%]
    -q [% id %] \
[% END -%]
    -t [% target_id %]

[% END -%]

EOF
    $tt->process(
        \$text,
        {   findbin     => $FindBin::Bin,
            working_dir => $working_dir,
            seq_dir     => $seq_dir,
            name_str    => $name_str,
            target_id   => $target_id,
            query_ids   => \@query_ids,
            is_self     => $is_self,
            length      => $paralog_length,
        },
        path( $working_dir, "prepare.sh" )->stringify
    ) or die Template->error;

    print "Create redo_prepare.sh\n";
    $tt->process(
        \$text,
        {   findbin     => $FindBin::Bin,
            working_dir => $working_dir,
            seq_dir     => $seq_dir,
            name_str    => $name_str,
            target_id   => $target_id,
            query_ids   => \@query_ids,
            is_self     => $is_self,
            length      => $paralog_length,
            redo        => 1,                 # If RepeatMasker has been executed
                                              # don't pass $seq_dir
                                              # don't gather taxon info
        },
        path( $working_dir, "redo_prepare.sh" )->stringify
    ) or die Template->error;
}

$stopwatch->end_message;
exit;

sub prep_fa {
    my $all_files = shift;
    my $acc       = shift;
    my $dir       = shift;

    my ($fna_file) = grep {/$acc/} @{$all_files};
    if ( !$fna_file ) {
        return "Can't find fasta file for $acc\n";
    }

    open my $in_fh, '<', $fna_file;
    open my $out_fh, '>', path( $dir, "$acc.fa" )->openw;
    while (<$in_fh>) {
        if (/>/) {
            print {$out_fh} ">$acc\n";
        }
        else {
            print {$out_fh} $_;
        }
    }
    close $out_fh;
    close $in_fh;

    return;
}

sub get_taxon_wgs {
    my $dbh      = shift;
    my $taxon_id = shift;

    my $query = qq{ SELECT wgs FROM gr WHERE taxonomy_id = ? };
    my $sth   = $dbh->prepare($query);
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

perl bac_prepare.pl --base_dir d:/bacteria/bacteria_101015 --parent 562
perl d:/wq/Scripts/tool/replace.pl -d d:/wq/Scripts/alignDB/bac -p "cmd.bat" -f /home/wangq -r d:/wq
