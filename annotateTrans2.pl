#! /usr/bin/perl
#
use strict;
use warnings;

use Data::Dumper;
use Time::HiRes;
use Memory::Usage;
use Set::IntervalTree;
#
use lib "/Share/home/zhangqf/cliff/icshape/module";
use icSHAPEutil qw( &readGTF_ensembl_new );

##--------------------------------------------------
#
my $_debug = 0;

my %enviroment = (
    annotationFile => "/Share/home/zhangqf/cliff/paris/data/ref/gencode.v21.chr_patch_hapl_scaff.annotation.gtf"
);

my %global = (
    bw => 100000,
    filter => "onlyInter"
);

#
##--------------------------------------------------
#
my $readGroup = shift;
my $supportSAM = shift;
my $output = shift;

&main ( $readGroup, $supportSAM, $output );
sub main
{
    my $readGroupFile = shift;
    my $supportSamFile = shift;
    my $outputFile = shift;

    my %readGroup = ();
    my %supportReads = ();

    my $ref_annotation = readGTF_ensembl_new ( $enviroment{annotationFile} );
    my $ref_bin = binize ( $ref_annotation->{gene_info}, $ref_annotation->{chr_size}, bw => $global{bw} );

    my %supportReads = ();
    my $memoryUsage = Memory::Usage->new(); $memoryUsage->record('starting work');

    my $validGroup = loadReadGroup ( $readGroupFile, \%readGroup, \%supportReads );
    #my $totalGroup = loadReadGroup ( $readGroupFile, $outputFile, \%supportReads, $ref_annotation, $ref_bin, bw => $global{bw}, onlyInter => 1 );
    if ( defined $supportSamFile ) {
        my $validReads = loadSupportSam ( $supportSamFile, \%supportReads  );
    }

    my $filteredSamFile = $supportSamFile; $filteredSamFile =~ s/.sam$/.filtered.sam/;
    printRG ( \%readGroup, $outputFile, minSupport => $global{support} );
    printSupportSam ( \%supportReads, $filteredSamFile );

    $memoryUsage->record('final memory usage'); $memoryUsage->dump();

    1;
}

## ----------------------------------
sub loadReadGroup
{
    my $readGroupFile = shift;
    my $outputFile = shift;
    my $ref_supportReads = shift;
    my $ref_annotation = shift;
    my $ref_bin = shift;
    my %parameters = @_;

    my $readGroupCount = 0;
    my %readGroup_interval = ();
    open ( RG, $readGroupFile ) or die "Cannot open $readGroupFile for reading!\n";
    open ( OUT, ">$outputFile" ) or die "Cannot open $outputFile for writing!\n";
    OUTER: while ( my $line = <RG> ) {
        next if ( $line =~ /^#/ );
	if ( $line =~ /^Group/ ) {
	    $readGroupCount++;
	    my ( $chr1, $strand1, $start1, $end1, $chr2, $strand2, $start2, $end2, $support ) = ( $line =~ /position (.+)\(([+-])\):(\d+)-(\d+)\|(.+)\(([+-])\):(\d+)-(\d+), support (.+)\./ );
	    print OUT $line;

	    my @overlapRegion = ();
	    &convert ( \@overlapRegion, $chr1, $strand1, $start1, $end1, $ref_annotation, $ref_bin, bw => $parameters{bw} );
    	    my ( $transcript1, $range1 ) = &printPosition ( \@overlapRegion );
	    @overlapRegion = ();
	    &convert ( \@overlapRegion, $chr2, $strand2, $start2, $end2, $ref_annotation, $ref_bin, bw => $parameters{bw} );
    	    my ( $transcript2, $range2 ) = &printPosition ( \@overlapRegion );

	    my $transcriptLine = "";
	    if ( $parameters{onlyIntra} ) {
		if ( $transcript1 and $transcript2 and ( $transcript1 eq $transcript2 ) ) {
		    $transcriptLine = $line . $transcript1 . $range1 . "\n" . $transcript2 . $range2 . "\n";
		}
	    }
	    elsif ( $parameters{onlyInter} ) {
		if ( $transcript1 and $transcript2 and ( $transcript1 ne $transcript2 ) ) {
		    $transcriptLine = $line . $transcript1 . $range1 . "\n" . $transcript2 . $range2 . "\n";
		}
	    } 
	    elsif ( $parameters{requireBothAnnotation} ) {
		if ( $transcript1 and $transcript2 ) {
		    $transcriptLine = $line . $transcript1 . $range1 . "\n" . $transcript2 . $range2 . "\n";
		}
	    } 
	    elsif ( $parameters{requireAnnotation} ) {
		if ( $transcript1 ) {
		    $transcriptLine = $line . $transcript1 . $range1 . "\n";
		}
		if ( $transcript2 ) {
		    $transcriptLine = $line . $transcript2 . $range2 . "\n";
		}
	    }
	    else {
		$transcriptLine = $line . "no transcript annotated\n";
	    }

	    print OUT $transcriptLine if ( $transcriptLine );

	    $line = <RG>;
	    print OUT $line;
	    INNER: while ( $line =<RG> ) {
	       print OUT $line if ( $transcriptLine );

	       my @support = split ( /\s+/, $line );
	       if ( $support[0] ) { $ref_supportReads->{$support[0]} = 1; }
	       last INNER if ( not ( $line =~ /\S/ ) );
       	    }
	}
    }
    close RG;
    print "in total $readGroupCount read groups read from $readGroupFile.\n";

    return $readGroupCount;
}

## ------------------------------------
sub binize
{
    my $ref_featurePos = shift;
    my $ref_chr_size = shift;
    my %parameters = @_;

    print STDERR "Binize genome to speed up searching.\n\t", `date`;
    my $bw = $parameters{bw};
    my %bin = ();
    foreach my $chr ( keys %{$ref_chr_size} ) {
        my $count = int ( $ref_chr_size->{$chr} / $bw ) + 1;
        for ( my $idx = 0; $idx < $count; $idx++ ) { $bin{$chr}{"+"}[$idx] = ();  $bin{$chr}{"-"}[$idx] = ();  }
    }

    foreach my $featureID ( keys %{$ref_featurePos} ) {
        my $chr = $ref_featurePos->{$featureID}{chr};
        my $strand = $ref_featurePos->{$featureID}{strand};
        my $start = int ( $ref_featurePos->{$featureID}{start} / $bw );
        my $end = int ( $ref_featurePos->{$featureID}{end} / $bw );
        for ( my $idx = $start; $idx <= $end; $idx++ ) { push @{$bin{$chr}{$strand}[$idx]}, $featureID; }
    }

    return \%bin;
}

sub convert
{
    my $ref_overlapRegion = shift;
    my $chr = shift;  my $strand = shift;  my $start = shift;  my $end = shift;
    my $ref_annotation = shift;
    my $ref_bin = shift;
    my %parameters = @_;

    my $bw = $parameters{bw};
    my $overlapCount = 0;

        for ( my $idxBin = int ( $start / $bw ); $idxBin <= int ( $end / $bw ); $idxBin++ ) {
            ## get genes in the bin
            foreach my $gene ( @{$ref_bin->{$chr}{$strand}[$idxBin]} ) {
                # get transcripts for each gene
                next if ( ( $end < $ref_annotation->{gene_info}{$gene}{start} ) or ( $start > $ref_annotation->{gene_info}{$gene}{end} ) );
                foreach my $transID ( @{$ref_annotation->{gene_info}{$gene}{transcript}} ) {
                    next if ( ( $end < $ref_annotation->{transcript_info}{$transID}{start} ) or ( $start > $ref_annotation->{transcript_info}{$transID}{end} ) );
                    $overlapCount += &overlapTrans ( $ref_overlapRegion, $ref_annotation, $transID, $strand, $start, $end );
                }
            }
        }

    return $overlapCount;
}

sub overlapTrans
{
    my $ref_overlapRegion = shift;
    my $ref_annotation = shift; my $transID = shift;
    my $strand = shift; my $absStart = shift; my $absEnd = shift; 

    my $transcript = $ref_annotation->{transcript_info}{$transID};
    my $numExon = scalar ( keys %{$transcript->{exon}} );

    my $overlapCount = 0;
    my $relExonStart = 0;  my $startPosiInExon = 0; my $endPosiInExon = 0; my $absStartInExon = 0; my $absEndInExon = 0; 
    for ( my $idxExon = 1; $idxExon <= $numExon; $idxExon++ ) {
        my $idxExonWithStrand = ( $strand eq "+" ) ? $idxExon : ( $numExon - $idxExon + 1 );
        my $exonID = $ref_annotation->{transcript_info}{$transID}{exon}{$idxExon};
        my $exonStart = $ref_annotation->{exon_info}{$exonID}{start};
        my $exonEnd = $ref_annotation->{exon_info}{$exonID}{end};
        my $exonLength = $exonEnd - $exonStart + 1;

        if ( ( $exonStart <= $absEnd ) && ( $exonEnd >= $absStart ) ) {
            ## overlapped
            if ( $strand eq "+" ) {
                $startPosiInExon = $absStart - $exonStart + 1;
                $startPosiInExon = 1 if ( $startPosiInExon < 1 );
                $endPosiInExon = $absEnd - $exonStart + 1;
                $endPosiInExon = $exonLength if ( $endPosiInExon > $exonLength );
            }
            elsif ( $strand eq "-" ) {
                $startPosiInExon = $exonEnd - $absEnd + 1;
                $startPosiInExon = 1 if ( $startPosiInExon < 1 );
                $endPosiInExon = $exonEnd - $absStart + 1;
                $endPosiInExon = $exonLength if ( $endPosiInExon > $exonLength );
            }

            my $relStart += $relExonStart + $startPosiInExon;
            my $relEnd += $relExonStart + $endPosiInExon;

            $absStartInExon = ( $absStart >= $exonStart ) ? $absStart : $exonStart; 
            $absEndInExon = ( $absEnd >= $exonEnd ) ? $exonEnd : $absEnd; 

            my $overlapString = join ( "\t", $absStartInExon, $absEndInExon, $transID, $relStart, $relEnd );
            push @{$ref_overlapRegion}, $overlapString;
            $overlapCount++;
        }

        $relExonStart += $exonLength;
    }

    return $overlapCount;
}

sub printPosition
{
    my $ref_overlapRegion = shift;

    my $outline = "";
    foreach my $overlap ( @{$ref_overlapRegion} ) 
        { $outline .= $overlap . "\n"; }

    return $outline;
}
