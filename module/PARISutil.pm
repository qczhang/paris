## --------------- package PARIS ---------------
package PARISutil;

use strict;
use warnings;

use Carp;
use File::Basename;
use Data::Dumper;
#use Set::IntervalTree;

use base 'Exporter';
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK );

our $VERSION     = '0.01';
our @EXPORT      = ();
our @EXPORT_OK   = qw( readBED readIcSHAPE loadGenome readGTF_ensembl_new getExonID get5primeLen getCDSLen get3primeLen getBioType parseCigar reverseComplement localAlignment );

my $_debug = 0;
my %enviroment = (
    maxChr	      => 999999999999,
    matchEnergy       => { 'AA' => -999, 'AC' => -999, 'AG' => -999, 'AT' => 1,
                           'CA' => -999, 'CC' => -999, 'CG' => 1, 'CT' => -999,
                           'GA' => -999, 'GC' => 1, 'GG' => -999, 'GT' => 1,
                           'TA' => 1, 'TC' => -999, 'TG' => 1, 'TT' => -999 },
    gapPenalty        => -1
);


sub readBED 
{
    my %parameters = @_;
}

sub loadGenome
{
    my $fasta = shift;
    my %parameters = @_;

    my %chr_seq = ();
    my $lineCount = 0;
    open ( FA, $fasta ) or die ( "Error in reading fasta file $fasta!\n" );
    print "read fasta file $fasta...\n\tTime:", `date`;
    if ( $parameters{faidx} ) {

    }
    elsif ( $parameters{simple} ) {
        while ( my $line = <FA> ) {
            $lineCount++;
            print "line: $lineCount\n\t", `date` if ( $lineCount % 1000000 == 0 );
            chomp $line;
            my $id = substr ( $line, 1 );
            $id =~ s/^(\s+)//g;

            $line = <FA>;
            chomp $line;
            $chr_seq{$id} = $line;
        }
    }
    else {
        my $id = "";
        my $seq = "";
        my $line = "";
        while ( $line = <FA> ) {
            $lineCount++;
            print "line: $lineCount\n\t", `date` if ( $lineCount % 1000000 == 0 );
            chomp $line;
            if ( $line =~ /^>/ ) {
                if ( $seq ) {
                    $chr_seq{$id} = $seq;
                    $seq = "";
                }
                $id = substr ( $line, 1 );
                ( $id ) = split ( /\s/, $id );
            }
            else {
                $seq .= $line;
            }
        }
        $chr_seq{$id} = $seq if ( $seq );
    }
    close FA;

    print "$lineCount lines read from file $fasta.\n\tTime:", `date`, "\n";
    return \%chr_seq;
}

sub readIcSHAPE
{
    my $icSHAPE = shift;
    my %parameters = @_;

    my %trans_icSHAPE = ();
    my $lineCount = 0;
    open ( SH, $icSHAPE ) or die ( "Error in reading icSHAPE file $icSHAPE!\n" );
    print "read icSHAPE file $icSHAPE...\n";
    while ( my $line = <SH> ) {
        $lineCount++;
        chomp $line;

        my ( $id, $length, $rpkm, @scores ) = split ( /\t/, $line );
        $trans_icSHAPE{$id} = \@scores;
    }
    close SH;

    return \%trans_icSHAPE;
}


sub readGTF_ensembl_new {
    my $gtfFile = shift;
    my %parameters = @_;

    my %chr_size = (); my %gene_info = (); my %transcript_info = (); my %exon_info = ();

    my $lineCount = 0;
    open ( GTF, $gtfFile ) or die ( "Error in reading GTF file $gtfFile!\n" );
    print "read genomic annotations from file: $gtfFile\n\t", `date`;
    my $geneID = "";  my $transcriptID = "";  my $exonID = "";  my $exonNum = 0;
    my $geneType = "";  my $geneName = "";  my $transcriptType = "";  my $transcriptName = "";
    while ( my $line = <GTF> ) {
        next if ( $line =~ /^#/ ); chomp $line;
        $lineCount++; print "  line: $lineCount\n" if ( ( $parameters{verbose} ) and ( $lineCount % 100000 == 0 ) );

        # 1     nonsense_mediated_decay     stop_codon      4782680     4782682     .       -       0       
        # gene_id "ENSMUSG00000033845"; transcript_id "ENSMUST00000045689"; exon_number "2"; 
        # gene_name "Mrpl15"; gene_biotype "protein_coding"; transcript_name "Mrpl15-003";
        # exon_id ...
        ## this is the statistics of feature field of mouse.74.gtf (in total 1242614 lines):
        #   432247 CDS
        #   628052 exon
        #    41170 start_codon
        #    41145 stop_codon
        ## this is the statistics of attribute field of mouse.74.gtf (in total 1242614 lines):
        #  1142614 gene_id
        #  1142614 gene_name
        #  1142614 gene_biotype
        #  1142614 transcript_id
        #  1142614 transcript_name
        #  1142614 exon_number
        #   628052 exon_id
        #   432247 protein_id
	$geneID = "";  $transcriptID = "";  $exonID = "";  $exonNum = 0;  $geneType = "";  $geneName = "";  $transcriptType = "";  $transcriptName = "";

        my ( $chr, $class, $feature, $start, $end, $score, $strand, $frame, $attribute ) = split ( /\t/, $line );
	if ( not defined $chr_size{$chr} ) { $chr_size{$chr} = $end; }
	else { $chr_size{$chr} = $end if ( $end > $chr_size{$chr} ); }

        my @data = split ( /; /, $attribute );
        foreach my $field ( @data ) {
            my $index = index ( $field, " " );
            next if ( $index <= 0 );

            my $type = substr ( $field, 0, $index ); 
            if ( $type eq "gene_id" ) { $geneID = substr ( $field, $index+2, -1 ); }
            elsif ( $type eq "gene_type" ) { $geneType = substr ( $field, $index+2, -1 ); }
            elsif ( $type eq "gene_name" ) { $geneName = substr ( $field, $index+2, -1 ); }
            elsif ( $type eq "transcript_id" ) { $transcriptID = substr ( $field, $index+2, -1 ); }
            elsif ( $type eq "transcript_type" ) { $transcriptType = substr ( $field, $index+2, -1 ); }
            elsif ( $type eq "transcript_name" ) { $transcriptName = substr ( $field, $index+2, -1 ); }
            elsif ( $type eq "exon_id" ) { $exonID = substr ( $field, $index+2, -1 ); }
            elsif ( $type eq "exon_number" ) { ( $exonNum ) = ( substr ( $field, $index+1 ) =~ /(\d+)/ ); }
        }
        if ( ( $feature eq "gene" ) and ( $geneID ) ) {                   # not defined in ensembl
            if ( defined $gene_info{$geneID} ) { print STDERR "Warnning! skipping line $lineCount of repeated geneID: $geneID\n\t$line\n"; next; }
            $gene_info{$geneID}{chr} = $chr; $gene_info{$geneID}{strand} = $strand;
            $gene_info{$geneID}{start} = $start; $gene_info{$geneID}{end} = $end;
            $gene_info{$geneID}{bioType} = $geneType; $gene_info{$geneID}{geneName} = $geneName;
        }
        if ( ( $feature eq "transcript" ) and ( $transcriptID ) ) {       # not defined in ensembl
            if ( not $geneID ) { print STDERR "Warnning! skipping line $lineCount of no geneID:\n\t$line\n"; next; }
            if ( defined $transcript_info{$transcriptID} ) 
                { print STDERR "Warnning! skipping line $lineCount of repeated transcriptID: $transcriptID\n\t$line\n"; next; }
            $transcript_info{$transcriptID}{gene} = $geneID;
            $transcript_info{$transcriptID}{start} = $start; $transcript_info{$transcriptID}{end} = $end;
            $transcript_info{$transcriptID}{bioType} = $transcriptType;
            $transcript_info{$transcriptID}{transcriptName} = $transcriptName;

            if ( not defined $gene_info{$geneID} ) {
                $gene_info{$geneID}{chr} = $chr; $gene_info{$geneID}{strand} = $strand;
                $gene_info{$geneID}{start} = $start; $gene_info{$geneID}{end} = $end;
            }
            else {
                $gene_info{$geneID}{start} = $start if ( $start < $gene_info{$geneID}{start} ); 
                $gene_info{$geneID}{end} = $end if ( $end > $gene_info{$geneID}{end} );
            }
            push @{$gene_info{$geneID}{transcript}}, $transcriptID;
       }
        if ( ( $feature eq "start_codon" ) and ( $transcriptID ) ) {
            if ( not $geneID ) { print STDERR "Warnning! skipping line $lineCount of no geneID:\n\t$line\n"; next; }
            if ( defined $transcript_info{$transcriptID}{startCodon} ) 
                { print STDERR "Warnning! skipping line $lineCount of repeated start_codon of transcriptID: $transcriptID\n\t$line\n"; next; }
            $transcript_info{$transcriptID}{startCodon} = $start;
        }
        if ( ( $feature eq "stop_codon" ) and ( $transcriptID ) ) {
            if ( not $geneID ) { print STDERR "Warnning! skipping line $lineCount of no geneID:\n\t$line\n"; next; }
            if ( defined $transcript_info{$transcriptID}{stopCodon} ) 
                { print STDERR "Warnning! skipping line $lineCount of repeated stop_codon of transcriptID: $transcriptID\n\t$line\n"; next; }
            $transcript_info{$transcriptID}{stopCodon} = $start;
        }
        if ( $feature eq "exon" ) {                   
            if ( not $geneID ) { print STDERR "Warnning! skipping line $lineCount of no geneID:\n\t$line\n"; next; }
            if ( not $transcriptID ) { print STDERR "Warnning! skipping line $lineCount of no transcriptID:\n\t$line\n"; next; }
            if ( not $exonID ) { print STDERR "Warnning! skipping line $lineCount of no exonID:\n\t$line\n"; next; }

            if ( defined $exon_info{$exonID} ) {
                if ( ( $exon_info{$exonID}{start} != $start ) or ( $exon_info{$exonID}{end} != $end ) ) 
                    { print STDERR "Error! line $lineCount of inconsistent exonID $exonID:\n\t$line\n"; next; }
            }
            if ( defined $transcript_info{$transcriptID} ) {
                if ( not defined $gene_info{$geneID} ) 
                    { print STDERR "Warnning! in consistent transcript annotation of transcript $transcriptID in $line\n"; next; }
                if ( $transcript_info{$transcriptID}{gene} ne $geneID ) 
                    { print STDERR "Warnning! in consistent transcript annotation of transcript $transcriptID in $line\n"; next; }
                if ( defined $transcript_info{$transcriptID}{exon}{$exonNum} ) 
                    { print STDERR "Warnning! skipping line $lineCount of repeated exon_number $exonNum:\n\t$line\n"; next; }
            }
            if ( defined $gene_info{$geneID} ) {
                if ( ( $gene_info{$geneID}{chr} ne $chr ) or ( $gene_info{$geneID}{strand} ne $strand ) )
                    { print STDERR "Warnning! in consistent location annotation of gene $geneID in $line\n"; next; }
            }

            if ( not defined $gene_info{$geneID} ) {
                $gene_info{$geneID}{chr} = $chr; $gene_info{$geneID}{strand} = $strand;
                $gene_info{$geneID}{start} = $start; $gene_info{$geneID}{end} = $end;
            }
            else {
                $gene_info{$geneID}{start} = $start if ( $start < $gene_info{$geneID}{start} ); 
                $gene_info{$geneID}{end} = $end if ( $end > $gene_info{$geneID}{end} );
            }

            if ( not defined $transcript_info{$transcriptID} ) {
                $transcript_info{$transcriptID}{gene} = $geneID;
                $transcript_info{$transcriptID}{start} = $start; $transcript_info{$transcriptID}{end} = $end;
                $transcript_info{$transcriptID}{length} = $end - $start + 1;
                $transcript_info{$transcriptID}{exonNum} = 1;
                push @{$gene_info{$geneID}{transcript}}, $transcriptID;
            }
            else {
                $transcript_info{$transcriptID}{start} = $start if ( $start < $transcript_info{$transcriptID}{start} ); 
                $transcript_info{$transcriptID}{end} = $end if ( $end > $transcript_info{$transcriptID}{end} );
                $transcript_info{$transcriptID}{length} += ($end - $start + 1);
                $transcript_info{$transcriptID}{exonNum}++;
            }

            $transcript_info{$transcriptID}{exon}{$exonNum} = $exonID;
            $exon_info{$exonID}{start} = $start; $exon_info{$exonID}{end} = $end;
        }
    }
    close GTF;

    return { 
        chr_size            => \%chr_size,
        gene_info           => \%gene_info, 
        transcript_info     => \%transcript_info, 
        exon_info           => \%exon_info 
    };
}

sub getExonID
{
    my $ref_annotation = shift;
    my $transID = shift;
    my $start = shift;
    my $end = shift;

    my @exonID = ();
    my @exonLen = ();
    my $geneID = $ref_annotation->{transcript_info}{$transID}{gene};
    if ( $ref_annotation->{gene_info}{$geneID}{strand} eq "+" ) {
        for ( my $idxExon = 1; $idxExon <= $ref_annotation->{transcript_info}{$transID}{exonNum}; $idxExon++ ) {
            if ( $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{end} < $start ) {
                next;
            }
            elsif ( $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{start} > $end ) {
                last;
            }
            else {
                push @exonID, $idxExon;
		my $len = $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{end} - $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{start} + 1;
                push @exonLen, $len;
            }
        }
    }
    elsif ( $ref_annotation->{gene_info}{$geneID}{strand} eq "-" ) {
        for ( my $idxExon = 1; $idxExon <= $ref_annotation->{transcript_info}{$transID}{exonNum}; $idxExon++ ) {
            if ( $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{start} > $end ) {
                next;
            }
            elsif ( $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{end} < $start ) {
                last;
            }
            else {
                push @exonID, $idxExon;
		my $len = $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{end} - $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{start} + 1;
                push @exonLen, $len;
            }
        }
    }

    return ( \@exonID, \@exonLen );
}

sub getCDSLen
{
    my $ref_annotation = shift;
    my $transID = shift;

    return $ref_annotation->{transcript_info}{$transID}{length} - get5primeLen ( $ref_annotation, $transID ) - get3primeLen ( $ref_annotation, $transID );
}

sub get5primeLen 
{
    my $ref_annotation = shift;
    my $transID = shift;

    my $fivePrimeLen = 0;
    if ( not defined $ref_annotation->{transcript_info}{$transID}{startCodon} )  { return $fivePrimeLen };

    my $geneID = $ref_annotation->{transcript_info}{$transID}{gene};
    if ( $ref_annotation->{gene_info}{$geneID}{strand} eq "+" ) {
        my $TSS = $ref_annotation->{transcript_info}{$transID}{startCodon};
        for ( my $idxExon = 1; $idxExon <= $ref_annotation->{transcript_info}{$transID}{exonNum}; $idxExon++ ) {
            if ( $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{end} >= $TSS ) {
		$fivePrimeLen += $TSS - $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{start}; 
                last;
            }
            else {
		$fivePrimeLen += $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{end} - $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{start} + 1; 
            }
        }
    }
    elsif ( $ref_annotation->{gene_info}{$geneID}{strand} eq "-" ) {
        my $TSS = $ref_annotation->{transcript_info}{$transID}{startCodon} + 2;
        for ( my $idxExon = 1; $idxExon <= $ref_annotation->{transcript_info}{$transID}{exonNum}; $idxExon++ ) {
            if ( $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{start} < $TSS ) {
		$fivePrimeLen += $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{end} - $TSS; 
                last;
            }
            else {
		$fivePrimeLen += $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{end} - $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{start} + 1; 
            }
        }
    }

    return $fivePrimeLen;
}

sub get3primeLen 
{
    my $ref_annotation = shift;
    my $transID = shift;

    my $threePrimeLen = 0;
    if ( not defined $ref_annotation->{transcript_info}{$transID}{stopCodon} )  { return $threePrimeLen };

    my $geneID = $ref_annotation->{transcript_info}{$transID}{gene};
    if ( $ref_annotation->{gene_info}{$geneID}{strand} eq "-" ) {
        my $TES = $ref_annotation->{transcript_info}{$transID}{stopCodon};
        for ( my $idxExon = $ref_annotation->{transcript_info}{$transID}{exonNum}; $idxExon >= 1; $idxExon-- ) {
            if ( $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{end} >= $TES ) {
		$threePrimeLen += $TES - $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{start}; 
                last;
            }
            else {
		$threePrimeLen += $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{end} - $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{start} + 1; 
            }
        }
    }
    elsif ( $ref_annotation->{gene_info}{$geneID}{strand} eq "+" ) {
        my $TES = $ref_annotation->{transcript_info}{$transID}{stopCodon} + 2;
        for ( my $idxExon = $ref_annotation->{transcript_info}{$transID}{exonNum}; $idxExon >= 1; $idxExon-- ) {
            if ( $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{start} <= $TES ) {
		$threePrimeLen += $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{end} - $TES; 
                last;
            }
            else {
		$threePrimeLen += $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{end} - $ref_annotation->{exon_info}{$ref_annotation->{transcript_info}{$transID}{exon}{$idxExon}}{start} + 1; 
            }
        }
    }

    return $threePrimeLen;
}


sub getBioType 
{
    my $ref_annotation = shift;
    my $featureID = shift;
    my %parameters = @_;

    my $featureType = ( defined $parameters{featureType} ) ? $parameters{featureType} : "transcript_info";
    my $biotypeString = "NULL";
    if ( defined $ref_annotation->{$featureType}{$featureID} ) {
        if ( defined $ref_annotation->{$featureType}{$featureID}{bioType} ) {
            $biotypeString = $ref_annotation->{$featureType}{$featureID}{bioType};
        }

        $biotypeString =~ s/^_//; $biotypeString =~ s/_/ /g;
        $biotypeString =~ s/PSEUDO/pseudogene/;
        $biotypeString =~ s/^gene$/unclassified_gene/i; $biotypeString =~ s/unclassified gene$/unclassified_gene/i; $biotypeString =~ s/ gene$//;
    }

    return $biotypeString;
}

sub loadGTF_interval
{
    my $gtfFile = shift;
    my %parameters = @_;

    my %gene_info = (); my %transcript_info = (); my %exon_info = ();
    my %gene_interval = (); my %intergenic_interval = (); my %transcript_interval = (); my %exon_interval = (); my %intron_interval = ();

    my $lineCount = 0;
    open ( GTF, $gtfFile ) or die ( "Error in reading GTF file $gtfFile!\n" );
    print STDERR "read GTF file $gtfFile...\n\tTime:", `date`;
    while ( my $line = <GTF> ) {
        next if ( $line =~ /^#/ );
        $lineCount++;
        print "line: $lineCount\n\t", `date` if ( $lineCount % 100000 == 0 );
        last if ( $_debug and ( $lineCount > 48 ) ); 

        my ( $chr, $source, $class, $start, $end, $reserve, $strand, $reserve2, $info ) = split ( /\t/, $line );
        if ( not defined $gene_interval{$chr} ) { 
            $gene_interval{$chr}{"+"} = Set::IntervalTree->new;
            $gene_interval{$chr}{"-"} = Set::IntervalTree->new;
            $transcript_interval{$chr}{"+"} = Set::IntervalTree->new;
            $transcript_interval{$chr}{"-"} = Set::IntervalTree->new;
            $exon_interval{$chr}{"+"} = Set::IntervalTree->new;
            $exon_interval{$chr}{"-"} = Set::IntervalTree->new;
        }

        my $geneID = "";  my $transcriptID = "";  my $exonID = "";  my $exonNum = 0;
        my @data = split ( /; /, $info );
        foreach my $field ( @data ) {
            my $index = index ( $field, " " );
            next if ( $index <= 0 );
            my $type = substr ( $field, 0, $index ); 
            if ( $type eq "gene_id" ) { $geneID = substr ( $field, $index+2, -1 ); }
            elsif ( $type eq "transcript_id" ) { $transcriptID = substr ( $field, $index+2, -1 ); }
            elsif ( $type eq "exon_id" ) { $exonID = substr ( $field, $index+2, -1 ); }
            elsif ( $type eq "exon_number" ) { $exonNum = substr ( $field, $index+1 ); }
        }

        if ( ( $class eq "gene" ) and ( $geneID ) ) {
            if ( defined $gene_info{$geneID} ) {
                print STDERR "Warnning! skipping line $lineCount of repeated geneID: $geneID\n\t$line\n";
                next;
            }
            $gene_info{$geneID}{chr} = $chr;
            $gene_info{$geneID}{strand} = $strand;
            $gene_info{$geneID}{start} = $start;
            $gene_info{$geneID}{end} = $end;
            $gene_interval{$chr}{$strand}->insert ( { 
                    id => $geneID, 
                    start => $start,
                    end => $end }, $start, $end+1 );
        }
        if ( ( $class eq "transcript" ) and ( $transcriptID ) ) {
            if ( not $geneID ) {
                print STDERR "Warnning! skipping line $lineCount of no geneID:\n\t$line\n";
                next;
            }
            if ( defined $transcript_info{$transcriptID} ) {
                print STDERR "Warnning! skipping line $lineCount of repeated transcriptID: $transcriptID\n\t$line\n";
                next;
            }
            push @{$gene_info{$geneID}{transcript}}, $transcriptID;
            $transcript_info{$transcriptID}{gene} = $geneID;
            $transcript_info{$transcriptID}{start} = $start;
            $transcript_info{$transcriptID}{end} = $end;
            $transcript_interval{$chr}{$strand}->insert ( { 
                    id => $transcriptID, 
                    gene => $geneID,
                    start => $start,
                    end => $end }, $start, $end+1 );
        }
        if ( ( $class eq "exon" ) and ( $exonID ) ) {
            if ( not $transcriptID ) {
                print STDERR "Warnning! skipping line $lineCount of no transcriptID:\n\t$line\n";
                next;
            }
            elsif ( defined $transcript_info{$transcriptID}{exon}{$exonNum} ) {
                print STDERR "Warnning! skipping line $lineCount of repeated exon_number $exonNum:\n\t$line\n";
                next;
            }

            if ( defined $exon_info{$exonID} ) {
                $transcript_info{$transcriptID}{exon}{$exonNum} = $exonID;
                if ( ( $exon_info{$exonID}{start} != $start ) or ( $exon_info{$exonID}{end} != $end ) ) {
                    print STDERR "Warnning! line $lineCount of repeated and inconsistent exonID $exonID:\n\t$line\n";
                    $exonID .= "." . $lineCount;
                    print STDERR "\tthe exon is renamed to $exonID.\n";
                }
            }

            $transcript_info{$transcriptID}{exon}{$exonNum} = $exonID;
            $exon_info{$exonID}{start} = $start;
            $exon_info{$exonID}{end} = $end;
            $exon_interval{$chr}{$strand}->insert ( { 
                    id => $exonID, 
                    transcript => $transcriptID,
                    start => $start,
                    end => $end }, $start, $end+1 );
        }
    }
    close GTF;

    getIntergenic ( \%intergenic_interval, \%gene_info );
    getIntrons ( \%intron_interval, \%gene_info, \%transcript_info, \%exon_info );

    if ( $_debug ) {
        print STDERR "gene info:\n"; print STDERR Dumper \%gene_info;
        print STDERR "gene interval:\n"; printAnnotation ( \%gene_interval );
        print STDERR "transcript info:\n"; print STDERR Dumper \%transcript_info;
        print STDERR "transcript interval:\n"; printAnnotation ( \%transcript_interval );
        print STDERR "exon info:\n"; print STDERR Dumper \%exon_info;
        print STDERR "exon interval:\n"; printAnnotation ( \%exon_interval );
        print STDERR "intron interval:\n"; printAnnotation ( \%intron_interval );
    }

    return { 
        gene_info           => \%gene_info, 
        transcript_info     => \%transcript_info, 
        exon_info           => \%exon_info, 

        gene_interval       => \%gene_interval, 
        intergenic_interval => \%intergenic_interval, 
        transcript_interval => \%transcript_interval, 
        exon_interval       => \%exon_interval, 
        intron_interval     => \%intron_interval 
    };
}


sub getIntergenic
{
    1;
}

sub getIntrons
{
    my $ref_intron_interval = shift;
    my ( $ref_gene_info, $ref_transcript_info, $ref_exon_info ) = @_;

    foreach my $geneID ( keys %{$ref_gene_info} ) {
        my $chr = $ref_gene_info->{$geneID}{chr};
        my $strand = $ref_gene_info->{$geneID}{strand};
        if ( not defined $ref_intron_interval->{$chr}{$strand} ) { $ref_intron_interval->{$chr}{$strand} = Set::IntervalTree->new; }

        foreach my $transcriptID ( @{$ref_gene_info->{$geneID}{transcript}} ) {
            my $tStart = $ref_transcript_info->{$transcriptID}{start};
            my $tEnd = $ref_transcript_info->{$transcriptID}{end};

            my @exons = ();
            if ( $strand eq "+" ) { @exons = sort { $a <=> $b } ( keys %{$ref_transcript_info->{$transcriptID}{exon}} ); }
            elsif ( $strand eq "-" ) { @exons = sort { $b <=> $a } ( keys %{$ref_transcript_info->{$transcriptID}{exon}} ); }

            $ref_transcript_info->{$transcriptID}{exonInError} = "";
            my $exonID = $ref_transcript_info->{$transcriptID}{exon}{$exons[0]};
            if ( $ref_exon_info->{$exonID}{start} != $tStart ) { 
                print STDERR "ERROR in first exon for transcript: $transcriptID.\n"; 
                $ref_transcript_info->{$transcriptID}{exonInError} .= "firstExonError;";
            }
            $exonID = $ref_transcript_info->{$transcriptID}{exon}{$exons[-1]};
            if ( $ref_exon_info->{$exonID}{end} != $tEnd ) { 
                print STDERR "ERROR in last exon for transcript: $transcriptID.\n"; 
                $ref_transcript_info->{$transcriptID}{exonInError} .= "lastExonError;";
            }

            my $orderCheck = 0; my $eStart = 0; my $eEnd = 0;
            foreach my $exonNum ( @exons ) {
                my $exonID = $ref_transcript_info->{$transcriptID}{exon}{$exonNum};
                if ( $eStart > $ref_exon_info->{$exonID}{start} ) {
                    print STDERR "ERROR in exon order for transcript: $transcriptID.\n";  
                    $orderCheck = 1;
                    last;
                }
                else { $eStart = $ref_exon_info->{$exonID}{start}; }
                if ( $eEnd > $ref_exon_info->{$exonID}{end} ) {
                    print STDERR "ERROR in exon order for transcript: $transcriptID.\n";  
                    $orderCheck = 1;
                    last;
                }
                else { $eEnd = $ref_exon_info->{$exonID}{end}; }
            }
            if ( $orderCheck ) { $ref_transcript_info->{$transcriptID}{exonInError} .= "exonOrderError;"; }
            if ( $ref_transcript_info->{$transcriptID}{exonInError} ) { 
                $ref_transcript_info->{$transcriptID}{exonInError} =~ s/;$//; 
                print STDERR "\t...exon in error, introns will not be calculated.\n";  
                next; 
            }

            my $intronCount = 0; my $iStart = 0; my $iEnd = 0;
            foreach my $exonNum ( @exons ) {
                my $exonID = $ref_transcript_info->{$transcriptID}{exon}{$exonNum};
                $eStart = $ref_exon_info->{$exonID}{start};
                $eEnd = $ref_exon_info->{$exonID}{end};

                if ( $eStart > $tStart ) { $iEnd = $eStart-1; }
                if ( $iEnd ) {
                    $intronCount++;
                    $ref_intron_interval->{$chr}{$strand}->insert ( { 
                            id => $intronCount, 
                            transcript => $transcriptID,
                            start => $iStart,
                            end => $iEnd }, $iStart, $iEnd+1 );
                }
                if ( $eEnd < $tEnd ) { $iStart = $eEnd+1; }
            }
        }
    }

    1;
}


sub printArray
{
    my $ref_array = shift;

    for ( my $idx = 0; $idx < scalar ( @{$ref_array} ); $idx++ ) {
        if ( not defined $ref_array->[$idx] ) { print STDERR "\tnotDefined"; }
        elsif ( $ref_array->[$idx] eq "" ) { print STDERR "\tblank"; }
        else { print STDERR "\t", $ref_array->[$idx]; }
    }

    1;
}

sub reverseComplement
{
    my $inseq = shift;
    my $outseq = reverse ( $inseq );

    $outseq =~ tr/AGTCagtc/TCAGtcag/;
    return $outseq;
}

sub parseCigar
{
    my $cigar = shift;
    my %parameters = @_;

    my @match = split ( /[0-9]+/, $cigar );             # CIGAR: \*|([0-9]+[MIDNSHPX=])+
    shift @match;
    my @matchSize = split ( /[MIDNSHPX=]/, $cigar );

    if ( $_debug ) {
        print STDERR "CIGAR\t", $cigar, "\nmatchOper";
        printArray ( \@match );
        print STDERR "\nmatchSize";
        printArray ( \@matchSize );
        print STDERR "\n";
    }

    if ( $parameters{getLargestM} ) {
        my $largestM = 0;
        for ( my $idx = 0; $idx < scalar ( @match ); $idx++ ) {
            if ( $match[$idx] eq "M" ) {
                if ( $matchSize[$idx] > $largestM ) {
                    $largestM = $matchSize[$idx];
                }
            }
        }

        return $largestM;
    }
    if ( $parameters{getLeadingS} ) {
        if ( $match[0] ne "S" ) {
            print STDERR "Warning! unexpected CIGAR string: $cigar!\n";
            return 0;
        }
        else { return $matchSize[0]; }
    }
    elsif ( $parameters{getMatchLen} ) {
    }
    else {
        return ( \@match, \@matchSize );
    }

    1;
}

sub localAlignment
{
    my ($seq1, $seq2) = @_;

    # # initialization
    my @matrix;
    $matrix[0][0]{score}   = 0;
    $matrix[0][0]{pointer} = "none";
    for(my $j = 1; $j <= length($seq1); $j++) {
        $matrix[0][$j]{score}   = 0;
        $matrix[0][$j]{pointer} = "none";
    }
    for (my $i = 1; $i <= length($seq2); $i++) {
        $matrix[$i][0]{score}   = 0;
        $matrix[$i][0]{pointer} = "none";
    }

    # fill
    my $maxScore = 0; my $maxI = 0; my $maxJ = 0;
    for(my $i = 1; $i <= length($seq2); $i++) {
        for(my $j = 1; $j <= length($seq1); $j++) {
            my ($diagonalScore, $leftScore, $upScore);

            # calculate match score
            my $pair = substr($seq1, $j-1, 1) . substr($seq2, $i-1, 1);
            $diagonalScore = $matrix[$i-1][$j-1]{score}+$enviroment{matchEnergy}{$pair};

            # calculate gap scores
	    $upScore   = $matrix[$i-1][$j]{score}+$enviroment{gapPenalty};
            $leftScore = $matrix[$i][$j-1]{score}+$enviroment{gapPenalty};

            if ($diagonalScore <= 0 and $upScore <= 0 and $leftScore <= 0) {
                $matrix[$i][$j]{score}   = 0;
                $matrix[$i][$j]{pointer} = "none";
                next; # terminate this iteration of the loop
            }

            # choose best score
	    if ($diagonalScore >= $upScore) {
                if ($diagonalScore >= $leftScore) {
                    $matrix[$i][$j]{score}   = $diagonalScore;
                    $matrix[$i][$j]{pointer} = "diagonal";
                }
                else {
                    $matrix[$i][$j]{score}   = $leftScore;
                    $matrix[$i][$j]{pointer} = "left";
                }
            } 
            else {
                if ($upScore >= $leftScore) {
                    $matrix[$i][$j]{score}   = $upScore;
                    $matrix[$i][$j]{pointer} = "up";
                }
                else {
                    $matrix[$i][$j]{score}   = $leftScore;
                    $matrix[$i][$j]{pointer} = "left";
                }
            }

            # set maximum score
	    if ($matrix[$i][$j]{score} > $maxScore) {
                $maxI     = $i;
                $maxJ     = $j;
                $maxScore = $matrix[$i][$j]{score};
            }
        }
    }

    # trace-back
    my $align1 = ""; my $align2 = "";
    my $i = $maxI; my $j = $maxJ;
    while (1) {
        last if $matrix[$i][$j]{pointer} eq "none";

        if ($matrix[$i][$j]{pointer} eq "diagonal") {
            $align1 .= substr($seq1, $j-1, 1);
            $align2 .= substr($seq2, $i-1, 1);
            $i--; $j--;
        }
        elsif ($matrix[$i][$j]{pointer} eq "left") {
            $align1 .= substr($seq1, $j-1, 1);
            $align2 .= "-";
            $j--;
        }
        elsif ($matrix[$i][$j]{pointer} eq "up") {
            $align1 .= "-";
            $align2 .= substr($seq2, $i-1, 1);
            $i--;
        }   
    }

    my $align = ( reverse $align1 ) . ":" . $align2;
    if ( $_debug ) {
        print STDERR "local alignment input seq1: ", $seq1, ", input seq2: ", reverse ( $seq2 ), "\n";
        print STDERR "alignment: ", $align, ", max score: ", $maxScore, "\n";
    }

    return ( $maxScore, $align, $j+1, $maxJ, $i+1, $maxI);
}

sub printAnnotation 
{
    my $ref_annotation = shift;
    foreach my $chr ( keys %{$ref_annotation} ) {
        foreach my $strand ( keys %{$ref_annotation->{$chr}} ) {
        print STDERR "\t$chr\t$strand\n";
            my $all = $ref_annotation->{$chr}{$strand}->fetch ( 1, $enviroment{maxChr} );
            print STDERR Dumper $all;
        }
    }

    1;
}
