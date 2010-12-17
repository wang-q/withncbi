#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use Bio::SeqIO;
use Bio::Location::Simple;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Utils::Exception qw(throw);

use FindBin;
use AlignDB::Stopwatch;

use Readonly;
Readonly my $MAX_SEQUENCE_LEN => 2e7;
Readonly my $GENE_LOGIC_NAMES => {
    CDS  => 'genbank_cds',
    tRNA => 'genbank_trna',
    rRNA => 'genbank_rrna'
};
Readonly my $REPEAT_LOGIC_NAME    => 'genbank_repeat';
Readonly my $COORD_SYSTEM_VERSION => 'wangq1';

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server   = $Config->{database}->{server};
my $port     = $Config->{database}->{port};
my $username = $Config->{database}->{username};
my $password = $Config->{database}->{password};

my $do_init_db = 0;
my $ensembl_db = "try_embl";
my $infile     = undef;
my $format     = "genbank";
my $gene_type  = 'known';
my $to_fasta   = 0;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'server=s'   => \$server,
    'port=i'     => \$port,
    'username=s' => \$username,
    'password=s' => \$password,
    'infile|i=s' => \$infile,
    'ensembl=s'  => \$ensembl_db,
    'genetype=s' => \$gene_type,
    'init_db=s'  => \$do_init_db,
    'to_fasta=s'  => \$to_fasta,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

throw("--ensembl argument required") if ( !defined($ensembl_db) );
throw("--infile argument required")  if ( !defined($infile) );

#----------------------------------------------------------#
# init database
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Build ensembl $ensembl_db...");

if ($do_init_db) {
    my $cmd    = "mysql -h$server -P$port -u$username -p$password ";
    my $drop   = "-e \"DROP DATABASE IF EXISTS $ensembl_db;\"";
    my $create = "-e \"CREATE DATABASE $ensembl_db;\"";
    my $init   = "$FindBin::Bin/../ensembl.sql";

    print "#drop\n" . "$cmd $drop\n";
    system("$cmd $drop");
    print "#create\n" . "$cmd $create\n";
    system("$cmd $create");
    print "#init\n" . "$cmd $ensembl_db < $init\n";
    system("$cmd $ensembl_db < $init");
}

#----------------------------------------------------------#
# init objects
#----------------------------------------------------------#
# Get Ensembl DB adaptor
my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host   => $server,
    -user   => $username,
    -pass   => $password,
    -port   => $port,
    -dbname => $ensembl_db,
);

# Create new analysis object indicating a curated origin
my $gene_analyses;
for my $which ( keys %$GENE_LOGIC_NAMES ) {
    $gene_analyses->{$which} = Bio::EnsEMBL::Analysis->new(
        -logic_name => $GENE_LOGIC_NAMES->{$which} );
}

my $repeat_analysis
    = Bio::EnsEMBL::Analysis->new( -logic_name => $REPEAT_LOGIC_NAME );

# Create sequence object and load it
my $seqin = Bio::SeqIO->new(
    -file   => $infile,
    -format => $format,
);

# Main loading loop
while ( my $seq = $seqin->next_seq ) {

    # Split sequence and return a slice
    my $slice = load_sequence( $dba, $seq );

    foreach my $f ( $seq->get_SeqFeatures ) {
        if ( $f->primary_tag =~ /(CDS)|([rt]RNA)/i ) {
            store_gene( $dba, $f, $slice, $to_fasta );
        }

        if ( $f->primary_tag =~ /repeat/ ) {

            # store_repeat($f, $slice);
        }
    }
}

$stopwatch->end_message;
exit;

##################################################
# Usage      : my $slice = load_sequence( $dba, $seq );
# Purpose    : Loads the sequence for this organism into the database,
#              creating Coordinate system and seq_region as necessary
# Returns    : Bio::EnsEMBL::Slice
# Parameters : Bio::EnsEMBL::DBAdaptor $db
#              Bio::SeqI $seq
# Throws     : no exceptions
# Comments   : none
# See Also   : n/a
sub load_sequence {
    my $dba = shift;
    my $seq = shift;

    my ( $csname, $seq_region ) = get_seq_region_data($seq);

    my $csa = $dba->get_CoordSystemAdaptor;
    my $cs  = $csa->fetch_by_name($csname);

    if ( !$cs ) {
        $cs = Bio::EnsEMBL::CoordSystem->new(
            -NAME           => $csname,
            -SEQUENCE_LEVEL => 1,
            -VERSION        => $COORD_SYSTEM_VERSION,
            -RANK           => 1,
            -DEFAULT        => 1
        );
        $csa->store($cs);
    }

    my $slice = Bio::EnsEMBL::Slice->new(
        -SEQ_REGION_NAME   => $seq_region,
        -COORD_SYSTEM      => $cs,
        -START             => 1,
        -END               => $seq->length,
        -SEQ_REGION_LENGTH => $seq->length
    );

    my $slice_adaptor = $dba->get_SliceAdaptor;

    my $sequence = $seq->seq;
    $slice_adaptor->store( $slice, \$sequence );

    return $slice;
}

##################################################
# Usage      : my ($csname, $seq_region_name) = get_seq_region_data($seq);
# Purpose    : Gets the coordinate system name (e.g. 'chromosome', 'plasmid')
#              and the name of the seq_region
# Returns    : pair of strings
# Parameters : Bio::SeqI
# Throws     : no exceptions
# Comments   : none
# See Also   : n/a
sub get_seq_region_data {
    my $seq  = shift;
    my $type = 'chromosome';
    my $name = "00";

    foreach my $f ( $seq->get_SeqFeatures ) {
        if ( $f->primary_tag eq 'source' ) {
            if ( $f->start == 1 and $f->end == $seq->length ) {

                if ( $f->has_tag('plasmid') ) {
                    $type = 'plasmid';
                    my @vals = $f->get_tag_values('plasmid');
                    $name = $vals[0];

                    #strip off plasmid prefix if present
                    $name =~ s/plasmid\s+//i;
                }
                elsif ( $f->has_tag('chromosome') ) {
                    $type = 'chromosome';
                    my @vals = $f->get_tag_values('chromosome');
                    $name = $vals[0];

                    #strip off chromosome prefix if present
                    $name =~ s/chr(omosome)\s+//i;
                }
                last;
            }
        }
    }

    return ( $type, $name );
}

##################################################
# Usage      : store_gene( $dba, $f, $slice )
# Purpose    : Stores a gene feature (both protein- and RNA-coding)
# Returns    : void
# Parameters : Bio::SeqFeatureI, Bio::EnsEMBL::DBSQL::SliceAdaptor
# Throws     : Failed loading $gene_stable_id
# Comments   : none
# See Also   : n/a
sub store_gene {
    my ( $db, $f, $slice, $to_fasta ) = @_;

    my $gene_stable_id = get_gene_stable_id($f);

    my $gene = Bio::EnsEMBL::Gene->new;

    $gene->analysis( $gene_analyses->{ $f->primary_tag } );
    $gene->biotype($gene_type);
    $gene->stable_id($gene_stable_id);
    $gene->version(1);
    $gene->slice($slice);

    #$gene->description(undef) if !$gene->description;
    #$gene->status(undef) if !$gene->status;

    print STDERR sprintf( "Found CDS with ID %s\n", $gene_stable_id );

    my $tcount = 0;
    my $transcript_id = sprintf( "%s.%d", $gene_stable_id, $tcount++ );

    my $transcript = Bio::EnsEMBL::Transcript->new;
    $transcript->stable_id($transcript_id);
    $transcript->version(1);
    $transcript->slice($slice);

    $gene->add_Transcript($transcript);

    # Add the exons
    my @exons  = create_exons($f);
    my $ecount = 0;
    foreach my $exon (@exons) {
        $exon->stable_id( sprintf( "%s.%d", $gene_stable_id, $ecount++ ) );
        $exon->version(1);
        $exon->slice($slice);
        $transcript->add_Exon($exon);
    }

    if ( $f->primary_tag =~ /RNA/ ) {
        foreach my $exon (@exons) {
            $exon->phase(-1);
            $exon->end_phase(-1);
        }
    }

    # Handle protein CDS features
    if ( $f->primary_tag eq 'CDS' ) {

        # Based on /codon_start EMBL qualifier
        my $frame = get_initial_frame($f);

        @exons = @{ $transcript->get_all_Exons };

        # This code assumes no UTRs
        my $translation = Bio::EnsEMBL::Translation->new;
        my $rcount      = 0;
        $translation->stable_id(
            sprintf( "%s.%d", $gene_stable_id, $rcount++ ) );
        $translation->version(1);
        $translation->start_Exon( $exons[0] );
        $translation->start( 1 + $frame );
        $translation->end_Exon( $exons[$#exons] );
        $translation->end( $translation->end_Exon->length );

        set_exon_phases( $translation, @exons );

        foreach my $exon (@exons) {
            print STDERR
                sprintf(
                "Added exon start: %d end: %d strand: %d phase: %d\n",
                $exon->start, $exon->end, $exon->strand, $exon->phase );
        }

        $transcript->translation($translation);

        my $mrna_seq = Bio::Seq->new(
            -seq      => $transcript->translateable_seq,
            -moltype  => "dna",
            -alphabet => 'dna',
            -id       => $translation->stable_id
        );

        # Translate args: stop char, unknown aa char, frame,
        # table, full CDS, throw
        my $aa_seq = $transcript->translate->seq;

        if ( $aa_seq =~ /\*/ ) {
            print STDERR
                sprintf( "Failed translating %s after phase setting\n",
                $translation->stable_id );
        }
    }

    eval { $db->get_GeneAdaptor->store($gene); };
    throw( sprintf( "Failed loading %s\n%s\n", $gene_stable_id, $@ ) )
        if ($@);
    
    if ($to_fasta) {
        my ($trans) = @{ $gene->get_all_Transcripts };
        my $tag = $f->primary_tag;
        my $seq;
        if ($tag =~ /CDS/i) {
            $seq = $trans->translateable_seq;
        }
        else {
            $seq = $trans->spliced_seq;
        }
        
        open my $out_fh, '>>', "$tag.fasta";
        print {$out_fh} $seq, 'N' x 10, "\n";
        close $out_fh;
    }
}

##################################################
# Usage      : my $stable_id = get_gene_stable_id($f);
# Purpose    : Tries to get a sensible stable identifier from a bioperl
#              feature created from an embl flat file.
# Returns    : string
# Parameters : Bio::SeqFeatureI $f
# Throws     : throw if cannot determine stable identifier
# Comments   : none
# See Also   : n/a
sub get_gene_stable_id {
    my $f = shift;

    my $stable_id;

    if ( $f->has_tag('protein_id') ) {
        my @vals = $f->get_tag_values('protein_id');
        ($stable_id) = split( /\s+/, $vals[0] );
    }
    if ( !$stable_id && $f->has_tag('locus_tag') ) {
        my @vals = $f->get_tag_values('locus_tag');
        ($stable_id) = split( /\s+/, $vals[0] );
    }

    if ( !$stable_id && $f->has_tag('gene') ) {
        my @vals = $f->get_tag_values('gene');
        ($stable_id) = split( /\s+/, $vals[0] );
    }

    if ( !$stable_id && $f->has_tag('product') ) {
        my @vals = $f->get_tag_values('product');
        ($stable_id) = split( /\s+/, $vals[0] );
    }

    throw("Could not determine gene identifier\n") if !$stable_id;

    return $stable_id;
}

##################################################
# Usage      : my @exons  = create_exons($f);
# Purpose    : Returns a list of exons created from the location of
#              feature.
# Returns    : List of Bio::EnsEMBL::Exons
# Parameters : Bio::SeqFeatureI
# Throws     : no exceptions
# Comments   : none
# See Also   : n/a
sub create_exons {
    my $f = shift;

    my @exons;

    foreach my $loc ( $f->location->each_Location ) {
        push(
            @exons,
            Bio::EnsEMBL::Exon->new(
                -start  => $loc->start,
                -end    => $loc->end,
                -strand => $loc->strand
            )
        );

        print STDERR sprintf( "Creating exon at %d..%d on strand %d\n",
            $loc->start, $loc->end, $loc->strand );
    }

    if ( $f->has_tag('pseudo') ) {
        foreach my $exon (@exons) {
            $exon->phase(-1);
            $exon->end_phase(-1);
        }
    }
    else {
        foreach my $exon (@exons) {
            $exon->end_phase(0);
        }
    }

    return @exons;
}

##################################################
# Usage      : my $frame = get_initial_frame($f);
# Purpose    : Returns the frame specified by the codon_start tag of
#              feature.
# Returns    : int
# Parameters : Bio::SeqFeatureI
# Throws     : no exceptions
# Comments   : none
# See Also   : n/a
sub get_initial_frame {
    my $f = shift;

    my $frame = 0;

    if ( $f->has_tag('codon_start') ) {
        my @vals = $f->get_tag_values('codon_start');
        $frame = $vals[0] - 1;
    }

    return $frame;
}

##################################################
# Usage      : set_exon_phases( $translation, @exons );
# Purpose    : Sets the start and end phases of exons.
# Returns    : void
# Parameters : Bio::Otter::AnnotatedTranscript, Bio::Ensembl::Exons
# Throws     : no exceptions
# Comments   : none
# See Also   : n/a
sub set_exon_phases {
    my $translation = shift;
    my @exons       = @_;

    my $found_start = 0;
    my $found_end   = 0;
    my $phase       = 0;

    foreach my $exon (@exons) {

        # Internal and end exons
        if ( $found_start && !$found_end ) {
            $exon->phase($phase);
            $exon->end_phase( ( $exon->length + $exon->phase ) % 3 );
            $phase = $exon->end_phase;
        }

        if ( $translation->start_Exon == $exon ) {
            $exon->phase($phase);
            $exon->end_phase(
                (   ( $exon->length - $translation->start + 1 ) + $exon->phase
                ) % 3
            );
            $phase       = $exon->end_phase;
            $found_start = 1;
        }

        if ( $translation->end_Exon == $exon ) {
            $found_end = 1;
        }
    }
}

