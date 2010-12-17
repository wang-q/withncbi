#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use Data::Dumper;

use Bio::SeqIO;
use Bio::Location::Simple;
use Bio::Annotation::Collection;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Attribute;

use Bio::EnsEMBL::Utils::Exception qw(throw);

###
### Script globals.
###

my $MAX_SEQUENCE_LEN = 2e7;

#Hard-coded should be command line options:
my $CHROM_PREFIX     = "R";
my $GENE_LOGIC_NAMES = {
    CDS  => 'submittergeneannotation',
    tRNA => 'genbank_trna',
    rRNA => 'genbank_rrna'
};
my $REPEAT_LOGIC_NAME    = 'genbank_repeat';
my $COORD_SYSTEM_VERSION = 'TIGR3';

#Organelles -- some not rice
my %skip_accession = map { $_, 1 } qw( AB042240 AB076665 AB076666
    D32052 NC_001320 X86563);

###
### Command line options
###
my $gene_type = 'known';

my $ensembl_species = $ENV{ENSEMBL_SPECIES};

#Argument Processing
{    
    my $help = 0;
    my $man  = 0;
    GetOptions(
        "help|?"     => \$help,
        "man"        => \$man,
        "species=s"  => \$ensembl_species,
        'genetype:s' => \$gene_type
    ) or pod2usage(2);
    pod2usage( -verbose => 2 ) if $man;
    pod2usage(1) if $help;
}

throw("file argument required") unless @ARGV;

###
### Get Ensembl DB adaptor
###

$ENV{'ENSEMBL_SPECIES'} = $ensembl_species;

my $reg = "Bio::EnsEMBL::Registry";    #Use this to get adaptors
$reg->load_all("$ENV{GrameneEnsemblDir}/conf/init_rw")
    or die "load_all failed";

my $slice_adaptor = $reg->get_adaptor( $ensembl_species, 'core', 'Slice' )
    or die "can't get Slice adaptor for $ensembl_species";
my $attribute_adaptor
    = $reg->get_adaptor( $ensembl_species, 'core', 'Attribute' )
    or die "can't get Attribute adaptor for $ensembl_species";

my $dba = $slice_adaptor->db;          #DBAdaptor
my $dbc = $dba->dbc;                   #DBConnection
warn "user " . $dbc->username . ", db " . $dbc->dbname . "\n";

my $dbh = $dbc->db_handle;    #you can use this as a DBI database handle
warn "\$dbh is " . ref($dbh) . "\n";

my $assembly_sth = $dbc->prepare(
    "insert into assembly
      (asm_seq_region_id ,cmp_seq_region_id ,asm_start ,asm_end ,cmp_start ,cmp_end ,ori)
      values  (?, ?,?,?,?,?,?)"
);

###
### Create new analysis objects indicating a curated origin
###
my $repeat_analysis
    = Bio::EnsEMBL::Analysis->new( -logic_name => $REPEAT_LOGIC_NAME );
my $gene_analyses;
for my $what ( keys %$GENE_LOGIC_NAMES ) {
    $gene_analyses->{$what} = Bio::EnsEMBL::Analysis->new(
        -logic_name => $GENE_LOGIC_NAMES->{$what} );
}

###
# Get our meta right
###
$dbh->do(
    "insert into meta (meta_key,meta_value)
	    values  ('assembly.mapping'
	             ,'chromosome:$COORD_SYSTEM_VERSION|clone')"
);

###
# More globals!
###
my %gene_stable_id;    #Cavalierly assuming nothing like ours in db
my $default_gene_stable_id = "GRMG1000000000";
my $fasta_out;

###
### Create sequence object and load it
###
for my $emblfile (@ARGV) {
    warn "\n==== $emblfile ====";
    my $seqi = Bio::SeqIO->new(
        -file   => $emblfile,
        -format => 'genbank'
    );
    ###
    # Fasta output
    ###
    $fasta_out = Bio::SeqIO->new( -file => ">$emblfile.fasta",
        '-format' => 'fasta' );
    ###
    ### Main loading loop
    ###
    while ( my $seq = $seqi->next_seq ) {
        next if $skip_accession{ $seq->accession };
        warn $seq->accession;

        # Split sequence and return a slice
        my $slice = load_sequence( $dba, $seq );

        foreach my $f ( $seq->get_SeqFeatures ) {
            if (   $f->primary_tag =~ /CDS/
                || $f->primary_tag =~ /[rt]RNA/ )
            {
                store_gene( $dba, $f, $slice );
            }

            if ( $f->primary_tag =~ /repeat/ ) {

                # store_repeat($f, $slice);
            }
        }
    }

    #close io
    undef $seqi;
    undef $fasta_out;
}

=head2 store_gene

 Title    : store_gene
 Function : Stores a gene feature (both protein- and RNA-coding)
 Returns  : void
 Argument : Bio::SeqFeatureI, Bio::EnsEMBL::DBSQL::SliceAdaptor

=cut

sub store_gene {
    my ( $db, $f, $slice ) = @_;

    my $gene_stable_id = get_gene_stable_id($f);

    my $gene = Bio::EnsEMBL::Gene->new();

    $gene->analysis( $gene_analyses->{ $f->primary_tag } );
    $gene->type($gene_type);
    $gene->stable_id($gene_stable_id);
    $gene->version(1);
    $gene->slice($slice);

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
        eval {  #because have overlapping exons, e.g. in  CAE04768 on AL662978
            $transcript->add_Exon($exon);
        };
        warn "add_Exon: $@" if $@;
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
        my $aa_seq = $transcript->translate()->seq();

        print STDERR sprintf( "Translation is: %s\n", $aa_seq );

        if ( $aa_seq =~ /\*/ ) {
            print STDERR
                sprintf( "Failed translating %s after phase setting\n",
                $translation->stable_id );
        }
    }

    eval { $db->get_GeneAdaptor->store($gene); };
    throw( sprintf( "Failed loading %s\n%s\n", $gene_stable_id, $@ ) )
        if ($@);
}

=head2 get_gene_stable_id

  Arg [1]    : Bio::SeqFeatureI $f
  Example    : my $stable_id = get_gene_stable_id($f);
  Description: Tries to get a sensible stable identifier from a bioperl feature
               created from an embl flat file.
  Returntype : string
  Exceptions : throw if cannot determine stable identifier
  Caller     :

=cut

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

    $stable_id ||= $default_gene_stable_id++;

    while ( $gene_stable_id{$stable_id}++ ) {
        my $old = $stable_id++;
        $stable_id .= ".1" if $stable_id eq $old;
    }

    return $stable_id;

}

=head2 create_exons

 Title    : create_exons
 Function : Returns a list of exons created from the location of
            feature.
 Returns  : List of Bio::EnsEMBL::Exons
 Argument : Bio::SeqFeatureI

=cut

sub create_exons {
    my $f = shift;

    my @exons;

    foreach my $loc ( $f->location->each_Location ) {
        warn "Skipping tiny exon at " . $loc->start . ".." . $loc->end
            and next
            if abs( $loc->start - $loc->end ) <= 2;
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

=head2 get_initial_frame

 Title    : get_initial_frame
 Function : Returns the frame specified by the codon_start tag of
            feature.
 Returns  : int
 Argument : Bio::SeqFeatureI

=cut

sub get_initial_frame {
    my $f = shift;

    my $frame = 0;

    if ( $f->has_tag('codon_start') ) {
        my @vals = $f->get_tag_values('codon_start');
        $frame = $vals[0] - 1;
    }

    return $frame;
}

=head2 set_exon_phases

 Title    : set_exon_phases
 Function : Sets the start and end phases of exons.
 Returns  : void
 Argument : Bio::Otter::AnnotatedTranscript, Bio::Ensembl::Exons

=cut

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

=head2 get_seq_region_data

  Arg [1]    : Bio::SeqI
  Example    : my ($csname, $seq_region_name) = get_seq_region_data($seq);
  Description: Gets the coordinate system name (e.g. 'chromosome', 'plasmid')
               and the name of the seq_region
  Returntype : pair of strings
  Exceptions : none
  Caller     : load_sequence

=cut

sub get_seq_region_data {
    my $seq  = shift;
    my $type = 'chromosome';
    my $name = "00";

    foreach my $f ( $seq->get_SeqFeatures ) {
        if ( $f->primary_tag eq 'source' ) {
            if ( $f->start == 1 and $f->end == $seq->length ) {
                if ( $f->has_tag('chromosome') ) {
                    $type = 'chromosome';
                    my @vals = $f->get_tag_values('chromosome');
                    $name = $vals[0];
                    $name =~ s/chr(omosome)\s+//i
                        ;    #strip off chromosome prefix if present
                    last;
                }

                if ( $f->has_tag('plasmid') ) {
                    $type = 'plasmid';
                    my @vals = $f->get_tag_values('plasmid');
                    $name = $vals[0];
                    $name =~ s/plasmid\s+//i
                        ;    #strip off plasmid prefix if present
                    last;
                }
            }
        }
    }
    warn "seq_region_data $name $type\n";

    return ( $type, $name );
}

=head2 load_sequence

  Arg [1]    : Bio::EnsEMBL::DBAdaptor $db
  Arg [2]    : Bio::SeqI $seq
  Example    : load_sequence
  Description: Loads the sequence for this organism into the database, creating
               Coordinate system and seq_region as necessary
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub load_sequence {
    my $db  = shift;
    my $seq = shift;

    my $accession = $seq->accession,

        my ( $nomchrom, $clone );
    if ( my $source = source_feature($seq) ) {
        ($nomchrom) = $source->get_tag_values('chromosome')
            if $source->has_tag('chromosome');
        if ( $source->has_tag('clone') ) {
            ($clone) = $source->get_tag_values('clone');
            $clone =~ s/\s*\bBAC\b\s*//i;
            $clone =~ s/\s*\bPLASMID\b\s*//i;
            $clone =~ s/\s*\bCLONE\b\s*//i;
            $clone =~ s/ (PCR) (PRODUCT)/_$1_$2/i
                ;    #bizarre, but OSJNBa0010P23_PCR_product in TIGR rice assy
        }
    }
    unless ( defined $nomchrom ) {
        $nomchrom = $1 if $seq->description =~ /\bCHROMOSOME\s+([^ ,;]+),/i;
    }
    unless ($clone) {
        warn "look for clone in " . $seq->description;
        $clone = $1 if $seq->description =~ /\bCLONE\s+([^ ,;]+),/i;
        $clone =~ s/ (PCR) (PRODUCT)/_$1_$2/i
            ;        #bizarre, but OSJNBa0010P23_PCR_product in TIGR rice assy
    }
    $clone ||= $seq->accession;    #For TIGR Rice Assembly BX842241, BX957224
    my ( $csname, $seq_region ) = get_seq_region_data($seq);

    my $csa = $db->get_CoordSystemAdaptor();

    my $topcs = $csa->fetch_by_name($csname);
    if ( !$topcs ) {
        $topcs = Bio::EnsEMBL::CoordSystem->new(
            -NAME           => $csname,
            -SEQUENCE_LEVEL => 0,
            -VERSION        => $COORD_SYSTEM_VERSION,
            -RANK           => 1,
            -DEFAULT        => 1
        );
        $csa->store($topcs);
    }
    my $topslice = Bio::EnsEMBL::Slice->new(
        -SEQ_REGION_NAME   => "$CHROM_PREFIX${nomchrom}_$accession",
        -COORD_SYSTEM      => $topcs,
        -START             => 1,
        -END               => $seq->length(),
        -SEQ_REGION_LENGTH => $seq->length()
    );
    $slice_adaptor->store($topslice);

    my $cs = $csa->fetch_by_name('clone');
    if ( !$cs ) {
        $cs = Bio::EnsEMBL::CoordSystem->new(
            -NAME           => 'clone',
            -SEQUENCE_LEVEL => 1,
            -RANK           => 2,
            -DEFAULT        => 1
        );
        $csa->store($cs);
    }

    my $slice = Bio::EnsEMBL::Slice->new(
        -SEQ_REGION_NAME   => $accession,
        -COORD_SYSTEM      => $cs,
        -START             => 1,
        -END               => $seq->length(),
        -SEQ_REGION_LENGTH => $seq->length()
    );

    my $sequence = uc $seq->seq();
    $sequence =~ tr/ACGTN/N/c;
    $seq->seq($sequence);           #for fasta file
    $seq->id( $seq->accession );    #for fasta file
    $fasta_out->write_seq($seq);
    $slice_adaptor->store( $slice, \$sequence );

    write_assembly( $topslice, $slice );

    my @attributes;

    #bacpac: clone version phase site gi chrom length
    #maybe add later.

    push @attributes,
        Bio::EnsEMBL::Attribute->new(
        -CODE        => 'nomchrom',
        -NAME        => 'nominal chromosome',
        -DESCRIPTION => 'chromosome according to source of clone',
        -VALUE       => $nomchrom
        ) if defined $nomchrom;

    push @attributes,
        Bio::EnsEMBL::Attribute->new(
        -CODE        => 'version',
        -NAME        => 'GenBank version number',
        -DESCRIPTION => 'GenBank version number',
        -VALUE       => $seq->version
        ) if $seq->version;    # $seq->seq_version seems to be the same

    push @attributes,
        Bio::EnsEMBL::Attribute->new(
        -CODE        => 'gi',
        -NAME        => 'GenBank gi number',
        -DESCRIPTION => 'GenBank gi number',
        -VALUE       => $seq->primary_id
        ) if $seq->primary_id;

    push @attributes,
        Bio::EnsEMBL::Attribute->new(
        -CODE        => 'clone',
        -NAME        => 'clone name',
        -DESCRIPTION => 'clone name',
        -VALUE       => $clone
        ) if $clone;

    warn join( ",", get_PubMed($seq) ) . " --pm--\n";
    push @attributes, map {
        Bio::EnsEMBL::Attribute->new(
            -CODE => 'pubmed',
            -NAME => 'pubmed',
            -DESCRIPTION =>
                'pubmed id of a reference in clone genbank record',
            -VALUE => $_
            )
    } get_PubMed($seq);

    {
        my $phase = HTG_Phase($seq);
        push @attributes,
            Bio::EnsEMBL::Attribute->new(
            -CODE        => 'HTGphase',
            -NAME        => 'HTG phase',
            -DESCRIPTION => 'High Throughput Genome phase',
            -VALUE       => $phase
            ) if defined $phase;
    }

    {
        my $site = getSite($seq);
        warn "site ", $site || 'none', "\n";
        push @attributes,
            Bio::EnsEMBL::Attribute->new(
            -CODE        => 'seqsite',
            -NAME        => 'Sequencing Site',
            -DESCRIPTION => 'Sequencing Site',
            -VALUE       => $site
            ) if defined $site;
    }

    warn join( "; ", map { $_->code . "=" . $_->value } @attributes ) . "\n";

    $attribute_adaptor->store_on_Slice( $slice, \@attributes );

    return $slice;
}

sub source_feature
{    #given sequence,returns pointer to hash of source attributes
    my ($seq) = @_;
    for my $f ( $seq->get_SeqFeatures ) {
        return $f if $f->primary_tag eq 'source';
    }
    return undef;
}

sub HTG_Phase {    # High Throughput Genomics Phase
                   # may be 0, 1, 2 or 3 (level of quality)
                   # should default to 0 for non HTG clones

    my ($seq) = @_;

    return 3 if $seq->division eq 'PLN';

    for my $kw ( $seq->annotation->get_Annotations('keyword') ) {
        return $1 if $kw =~ /PHASE([0-3])/;
    }
    return 0;
}

sub getSite {

    my ($seq) = @_;
    my $site;

    foreach my $ref ( $seq->annotation()->get_Annotations('reference') ) {
        my $journal = $ref->location;

        warn "reference location $journal\n";

        # these should be made more specific and expanded
        #   for handling species other than rice:
        if ( $journal =~ /Japan/ ) {
            $site = "RGP(Japan)";

        }
        elsif ( $journal =~ /CHINA/ || $journal =~ /China/ ) {
            $site = "NCGR(China)";

        }
        elsif ( $journal =~ /Korea/ ) {
            $site = "KRGRP(Korea)";

        }
        elsif ( $journal =~ /Taiwan/ ) {
            $site = "ASPGC(Taiwan)";

        }
        elsif ( $journal =~ /Thailand/ ) {
            $site = "BIOTEC(Thailand)";

        }
        elsif ( $journal =~ /FRANCE/ || $journal =~ /France/ ) {
            $site = "Genoscope(France)";

        }
        elsif ( $journal =~ /MD/ ) {
            $site = "TIGR(USA)";

        }
        elsif ( $journal =~ /NY/ ) {
            $site = "CSHL(USA)";

        }
        elsif ( $journal =~ /SC/ ) {
            $site = "CUGI(USA)";
        }
        elsif ( $journal =~ /Arizona\s+Genomics\s+Institute/i ) {
            $site = "AGI(USA)";

        }
        elsif ($journal =~ /Plant\s+Genome\s+Initiative\s+at\s+Rutgers/
            || $journal =~ /\bWaksman\s*Institute.*Rutgers/ )
        {
            $site = "PGIR(USA)";

        }
        elsif ( $journal =~ /WI/ ) {
            $site = "GCW(USA)";

        }
        elsif ( $journal =~ /MO/ ) {
            $site = "GSC(USA)";

        }
        elsif ($journal =~ /\bIndian\s+Initiative\b/i
            || $journal =~ /\bIIRGS\b/
            || $journal =~ /University of Delhi South Campus.*India/
            || $journal
            =~ /Indian Agricultural Research Institute, LBS Centre/ )
        {
            $site = "IIRGS(India)";
        }
        elsif ( $journal =~ /\bBotany\b.*University of Georgia\b/i ) {
            $site = "LGBUGA(USA)";
        }
        print STDERR "j $journal -> $site ", $seq->accession, "\n" if $site;
        last if $site;
    }
    return $site;
}

sub get_PubMed {
    my ($seq) = @_;
    my $ac = $seq->annotation();

    my @ref = $ac->get_Annotations('reference');
    warn scalar(@ref) . " references\n";

    #    foreach my $ref ( @ref ) {
    #        warn join(" | ",ref($ref),$ref->can('pubmed'),$ref->pubmed,"\n")
    #	      .Dumper($ref)."\n";
    #    }
    return grep {$_} map {
        $_->pubmed

            #|| $_->medline ?? does this make sense?
    } @ref;
}

sub write_assembly {
    my ( $asmslice, $cmpslice ) = @_;
    $assembly_sth->execute(
        $asmslice->get_seq_region_id,
        $cmpslice->get_seq_region_id,
        $asmslice->start, $asmslice->end, $cmpslice->start, $cmpslice->end, 1
    );
}
