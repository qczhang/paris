#! /usr/bin/perl
#
use strict;
use warnings;

use Data::Dumper;
use Time::HiRes;
use Memory::Usage;
use Set::IntervalTree;
use Getopt::Std;

use lib "module";
use PARISutil qw( &loadGTF &loadGenome &parseCigar &reverseComplement &localAlignment );
#
##--------------------------------------------------
#
my $_debug = 0;

my %environment = (
    genomeFile        => "/Share/home/zhangqf/cliff/paris/data/ref/grch38.fna",
    annotationFile    => "/Share/home/zhangqf/cliff/paris/data/ref/gencode.v21.chr_patch_hapl_scaff.annotation.gtf",
    transcriptomeFile => "/Share/home/zhangqf/cliff/paris/data/ref/transcriptome.fa",

    maxChr            => 999999999999,
    intronFlanking    => 3,
    minOverhang	      => 20,
    minSupport        => 1
);

my %global = (
    genomeSeq               => {},
    annotation              => {},

    dsPairCount             => 0,
    dsPair                  => [],
    dsPairCluster           => {},

    dspInterval             => {},
    dspIntervalClusterCount => 0,
    dspIntervalCluster      => {},
    dspIntervalClusterRead  => {}
);

#
##--------------------------------------------------
#

use vars qw ($opt_h $opt_V $opt_D $opt_i $opt_j $opt_s $opt_o $opt_g $opt_a $opt_t $opt_l $opt_p );
&getopts('hVDi:j:s:o:g:a:t:l:p:');

my $usage = <<_EOH_;
## --------------------------------------
call base pair groups from PARIS sequencing
Command:
$0 -i input_sam_file -j chiastic_junction_file -s chiastic_support_sam_file -o output_read_group
# what it is:
 -i     input sam file of spliced alignment
 -j     chiastic junction file
 -s	chiastic junction support alignment file
 -o	output base pair groups
# more options:
 -g     genome file
 -a     annotation file
 -t     transcriptome file

 -l	minimum overhang length for valid mapping
 -p	minimum number of supporting reads 
_EOH_
;

&main();


sub main
{
    my %parameters = &init();

    my $samFileList = $parameters{samFiles};
    my $chiasticFileList = $parameters{chiasticFiles};
    my $chiasticSamFileList = $parameters{chiasticSamFiles};
    my $outputFile = $parameters{output};

    my $memoryUsage = Memory::Usage->new(); $memoryUsage->record('starting work');
    $global{annotation} = loadGTF ( $parameters{annotationFile} );
    $memoryUsage->record('Memory after GTF loading'); $memoryUsage->dump();
    $global{genomeSeq} = loadGenome ( $parameters{genomeFile} );
    $memoryUsage->record('Memory after genome loading'); $memoryUsage->dump();

    my $allSupportSam = "";
    if ( $samFileList ne "NULL" ) {
        $allSupportSam = $samFileList;
        my @samFiles = split ( /:/, $samFileList );
        foreach my $samFile ( @samFiles ) { genPairClusterFromSamFile ( $samFile ); }
    }
    $memoryUsage->record('Memory after sam file process'); $memoryUsage->dump();
    if ( $chiasticFileList ne "NULL" ) {
        if ( $allSupportSam ) { $allSupportSam = $samFileList . ":" . $chiasticSamFileList;  }
        else {  $allSupportSam = $chiasticSamFileList;  }
        my @chiasticFiles = split ( /:/, $chiasticFileList );
        foreach my $chiasticFile ( @chiasticFiles ) { genPairClusterFromJunctionFile ( $chiasticFile ); }
    }
    $memoryUsage->record('Memory after chiastic junction file process'); $memoryUsage->dump();

    sortCluster ( minSupport => $parameters{minSupport} );
    $memoryUsage->record('Memory after cluster sorting'); $memoryUsage->dump();
    printCluster ( $outputFile, supportSam => 1, inputSam => $allSupportSam );
    $memoryUsage->record('final memory usage'); $memoryUsage->dump();

    1;
}

## ----------------------------------
sub init {
    my %parameters = ();

    die $usage if ( $opt_h || ( not $opt_i ) || ( ( not $opt_j ) && ( not $opt_s ) ) || ( not $opt_o ) );
    $opt_V = 0 if ( not defined $opt_V );
    $opt_D = 0 if ( not defined $opt_D );

    if ( defined $opt_i ) { $parameters{samFiles} = $opt_i; }
    if ( defined $opt_j ) { $parameters{chiasticFiles} = $opt_j; }
    if ( defined $opt_s ) { $parameters{chiasticSamFiles} = $opt_s; }
    if ( defined $opt_o ) { $parameters{output} = $opt_o; }

    if ( defined $opt_g ) { $parameters{genomeFile} = $opt_g; }
    if ( defined $opt_a ) { $parameters{annotationFile} = $opt_a; }
    if ( defined $opt_t ) { $parameters{transcriptomeFile} = $opt_t; }

    if ( defined $opt_l ) { $environment{minOverhang} = $opt_l; }
    if ( defined $opt_p ) { $parameters{minSupport} = $opt_p; }

    return ( %parameters );
}

sub genPairClusterFromSamFile
{
    my $samFile = shift;

    my $lineCount = 0;
    my $validCount = 0;
    open ( SAM, $samFile ) or die ( "Error in reading sam file $samFile!\n" );
    print "read sam file $samFile...\n\tTime: ", `date`;
    while ( my $line = <SAM> ) {
        next if ( $line =~ /^#/ );
        if ( $line =~ /^@/ ) { }
        else {
            $lineCount++;
            if ( $lineCount % 100000 == 0 ) {
                print "line: $lineCount\n";
                print "\tvalid line: $validCount\n\t", `date`;
            }

            #last if ( $lineCount > 10 );
            $validCount += genPairClusterFromSamLine ( $line );
        }
    }
    close SAM;
    print "$lineCount lines read from sam file $samFile.\n";
    print "among which $validCount lines generate supports for base pairing.\n\tTime: ", `date`, "\n";

    return ( $lineCount, $validCount );
}

sub genPairClusterFromSamLine
{
    my $line = shift;

    my @data = split ( /\t/, $line );
    if ( $data[5] !~ /[ND]/ ) {
        print STDERR "Skip read that does not contain cleavage!\n\t$line\n";
        return 0;
    }

    my $strand = ( $data[1] & 16 ) ? "-" : "+";
    my ( $alignment, $pair1s, $pair1e, $pair2s, $pair2e ) = getSamPair ( $data[2], $strand, $data[3], $data[5], minOverhang => $environment{minOverhang} );
    if ( not $alignment ) {
        print STDERR "\t$line\n" if ( $opt_V );
        return 0;
    }

    $global{dsPair}[$global{dsPairCount}][2] = $alignment;
    $global{dsPair}[$global{dsPairCount}][3] = $data[2];  
    $global{dsPair}[$global{dsPairCount}][4] = $strand;
    $global{dsPair}[$global{dsPairCount}][5] = $pair1s; 
    $global{dsPair}[$global{dsPairCount}][6] = $pair1e;
    $global{dsPair}[$global{dsPairCount}][7] = $data[2];  
    $global{dsPair}[$global{dsPairCount}][8] = $strand;
    $global{dsPair}[$global{dsPairCount}][9] = $pair2s; 
    $global{dsPair}[$global{dsPairCount}][10] = $pair2e;
    $global{dsPair}[$global{dsPairCount}][11] = $data[0];
    $global{dsPair}[$global{dsPairCount}][12] = $data[5];
    $global{dsPair}[$global{dsPairCount}][13] = $data[3];

    my ( $contained, $id1, $id2 ) = newPair ( );
    if ( $contained ) {
        $global{dsPair}[$global{dsPairCount}][0] = $id1;
        $global{dsPair}[$global{dsPairCount}][1] = $id2;
    }
    else {
        $global{dsPair}[$global{dsPairCount}][0] = $global{dspIntervalClusterCount};
        $global{dsPair}[$global{dsPairCount}][1] = $global{dspIntervalClusterCount}+1;
        $global{dspIntervalClusterCount} += 2;
    }
    $global{dsPairCount}++;

    1;
}

sub genPairClusterFromJunctionFile
{
    my $junctionFile = shift;

    my $lineCount = 0;
    my $validCount = 0;
    open ( JUNC, $junctionFile ) or die ( "Error in reading junction file $junctionFile!\n" );
    print "read junction file $junctionFile...\n\tTime: ", `date`;
    while ( my $line = <JUNC> ) {
        next if ( $line =~ /^#/ );
        $lineCount++;
        if ( $lineCount % 100000 == 0 ) {
            print "line: $lineCount\n";
            print "\tvalid line: $validCount\n\t", `date`;
        }

        print "line: $lineCount\n\t", `date` if ( $lineCount % 100000 == 0 );
        chomp $line;
        $validCount += genPairClusterFromOneJunction ( $line );
    }
    close JUNC;
    print "$lineCount lines read from junction file $junctionFile.\n";
    print "among which $validCount lines generate supports for base pairing.\n\tTime: ", `date`, "\n";

    return $lineCount;
}

sub genPairClusterFromOneJunction
{
    my $line = shift;

    my @data = split ( /\t/, $line );
    my $cigar = "";  my $isChiastic = 0;
    if ( ( $data[0] ne $data[3] ) or ( $data[2] ne $data[5] ) ) {
        $isChiastic = 2;
        $cigar = $data[11] . "|" . $data[13];
    }
    else {
        if ( $data[2] eq "+" )  {
            if ( $data[1] > $data[4] ) { $isChiastic = 1; $cigar = getNewCigar ( 1, "+", $data[10], $data[12], $data[11], $data[13] ); }
            else { $cigar = getNewCigar ( 0, "+", $data[10], $data[12], $data[11], $data[13] ); }
        }
        elsif ( $data[2] eq "-" ) {
            if ( $data[1] < $data[4] ) { $isChiastic = 1; $cigar = getNewCigar ( 1, "-", $data[10], $data[12], $data[11], $data[13] ); }
            else { $cigar = getNewCigar ( 0, "-", $data[10], $data[12], $data[11], $data[13] ); }
        }
        else { print STDERR "\t$line\n"; return 0; }
    }
    if ( not $cigar ) {  print STDERR "\tSkip line of inapproprieate alignment: $line\n";  return 0;  }

    my ( $alignment, $pair1s, $pair1e, $pair2s, $pair2e ) = 
	getJuncPair ( $data[0], $data[2], $data[1], $data[10], $data[11], $data[3], $data[5], $data[4], $data[12], $data[13], minOverhang => $environment{minOverhang} );
    if ( not $alignment ) { print STDERR "\t$line\n"; return 0; }

    $global{dsPair}[$global{dsPairCount}][2] = $alignment;
    $global{dsPair}[$global{dsPairCount}][3] = $data[0];  
    $global{dsPair}[$global{dsPairCount}][4] = $data[2];
    $global{dsPair}[$global{dsPairCount}][5] = $pair1s; 
    $global{dsPair}[$global{dsPairCount}][6] = $pair1e;
    $global{dsPair}[$global{dsPairCount}][7] = $data[3];  
    $global{dsPair}[$global{dsPairCount}][8] = $data[5];
    $global{dsPair}[$global{dsPairCount}][9] = $pair2s; 
    $global{dsPair}[$global{dsPairCount}][10] = $pair2e;
    $global{dsPair}[$global{dsPairCount}][11] = $data[9];
    $global{dsPair}[$global{dsPairCount}][12] = $cigar . ":" . $isChiastic;
    if ( $isChiastic < 2 ) { $global{dsPair}[$global{dsPairCount}][13] = ( $data[10] < $data[12] ) ? $data[10] : $data[12]; }
    else {
        $global{dsPair}[$global{dsPairCount}][13] = $data[10];
        $global{dsPair}[$global{dsPairCount}][14] = $data[12];
    }

    my ( $contained, $id1, $id2 ) = newPair ( );
    if ( $contained ) {
        $global{dsPair}[$global{dsPairCount}][0] = $id1;
        $global{dsPair}[$global{dsPairCount}][1] = $id2;
    }
    else {
        $global{dsPair}[$global{dsPairCount}][0] = $global{dspIntervalClusterCount};
        $global{dsPair}[$global{dsPairCount}][1] = $global{dspIntervalClusterCount}+1;
        $global{dspIntervalClusterCount} += 2;
    }
    $global{dsPairCount}++;

    1;
}

sub getNewCigar
{
    my $isChiastic = shift;
    my $strand = shift;
    my ( $start1, $start2, $frag1Cigar, $frag2Cigar ) = @_;

    my ( $ref_match1, $ref_matchSize1 ) = parseCigar ( $frag1Cigar );
    my ( $ref_match2, $ref_matchSize2 ) = parseCigar ( $frag2Cigar );

    my $cigar = "";
    my $newReadFrag1 = "";  my $newReadFrag2 = "";
    if ( not $isChiastic ) {
        if ( $strand eq "+" ) {
            my $lenN = $start2 - $start1;
            for ( my $idx = 0; $idx < scalar ( @{$ref_match1} - 1 ); $idx++ ) { 
                $cigar .= $ref_matchSize1->[$idx] . $ref_match1->[$idx];  
                if ( ( $ref_match1->[$idx] eq "M" ) or ( $ref_match1->[$idx] eq "=" ) or ( $ref_match1->[$idx] eq "X" ) or ( $ref_match1->[$idx] eq "D" ) or ( $ref_match1->[$idx] eq "N" ) ) { 
                    $lenN -= $ref_matchSize1->[$idx];
                }
            }
	    #if ( $ref_match2->[0] eq "S" ) { 
	    #    $lenN -= $ref_matchSize2->[0];
	    #}
            if ( $lenN < 0 ) {  $cigar = 0;  }
            else {
                $cigar .= $lenN . "N";
                for ( my $idx = 1; $idx < scalar ( @{$ref_match2} ); $idx++ ) { $cigar .= $ref_matchSize2->[$idx] . $ref_match2->[$idx]; }
            }
        }
        elsif ( $strand eq "-" ) {
            my $lenN = $start1 - $start2;
            for ( my $idx = 0; $idx < scalar ( @{$ref_match2} - 1 ); $idx++ ) { 
                $cigar .= $ref_matchSize2->[$idx] . $ref_match2->[$idx]; 
                if ( ( $ref_match2->[$idx] eq "M" ) or ( $ref_match2->[$idx] eq "=" ) or ( $ref_match2->[$idx] eq "X" ) or ( $ref_match2->[$idx] eq "D" ) or ( $ref_match2->[$idx] eq "N" ) ) { 
                    $lenN -= $ref_matchSize2->[$idx];
                }
            }
	    #if ( $ref_match1->[0] eq "S" ) { 
	    #    $lenN -= $ref_matchSize1->[0];
	    #}
            if ( $lenN < 0 ) {  $cigar = 0;  }
            else {
                $cigar .= $lenN . "N";
                for ( my $idx = 1; $idx < scalar ( @{$ref_match1} ); $idx++ ) { $cigar .= $ref_matchSize1->[$idx] . $ref_match1->[$idx]; }
            }
        }
    }
    elsif ( $isChiastic == 1 ) {
        if ( $strand eq "+" ) {
            my $lenN = $start1 - $start2;
            for ( my $idx = 1; $idx < scalar ( @{$ref_match2} ); $idx++ ) { 
                if ( $ref_match2->[$idx] eq "S" )  { $cigar .= $ref_matchSize2->[$idx] . "M"; }
		else  { $cigar .= $ref_matchSize2->[$idx] . $ref_match2->[$idx]; }

                if ( ( $ref_match2->[$idx] eq "M" ) or ( $ref_match2->[$idx] eq "=" ) or ( $ref_match2->[$idx] eq "X" ) or ( $ref_match2->[$idx] eq "D" ) or ( $ref_match2->[$idx] eq "N" ) or ( $ref_match2->[$idx] eq "S" ) ) { 
                    $lenN -= $ref_matchSize2->[$idx];
                }
            }
	    if ( $ref_match1->[0] eq "S" )  { $lenN -= $ref_matchSize1->[0]; }
            if ( $lenN < 0 ) {  $cigar = 0;  }
            else {
                $cigar .= $lenN . "N";
                for ( my $idx = 0; $idx < scalar ( @{$ref_match1} - 1 ); $idx++ ) { 
                    if ( $ref_match1->[$idx] eq "S" )  { $cigar .= $ref_matchSize1->[$idx] . "M"; }
		    else  { $cigar .= $ref_matchSize1->[$idx] . $ref_match1->[$idx]; }
		}
            }
        }
        if ( $strand eq "-" ) {
            my $lenN = $start2 - $start1;
            for ( my $idx = 1; $idx < scalar ( @{$ref_match1} ); $idx++ ) { 
                if ( $ref_match1->[$idx] eq "S" )  { $cigar .= $ref_matchSize1->[$idx] . "M"; }
		else  { $cigar .= $ref_matchSize1->[$idx] . $ref_match1->[$idx]; }

                if ( ( $ref_match1->[$idx] eq "M" ) or ( $ref_match1->[$idx] eq "=" ) or ( $ref_match1->[$idx] eq "X" ) or ( $ref_match1->[$idx] eq "D" ) or ( $ref_match1->[$idx] eq "N" ) or ( $ref_match1->[$idx] eq "S" ) ) { 
                    $lenN -= $ref_matchSize1->[$idx];
                }
            }
	    if ( $ref_match2->[0] eq "S" )  { $lenN -= $ref_matchSize2->[0]; }
            if ( $lenN < 0 ) {  $cigar = 0;  }
            else {
                $cigar .= $lenN . "N";
                for ( my $idx = 0; $idx < scalar ( @{$ref_match2} - 1 ); $idx++ ) { 
                    if ( $ref_match2->[$idx] eq "S" )  { $cigar .= $ref_matchSize2->[$idx] . "M"; }
		    else  { $cigar .= $ref_matchSize2->[$idx] . $ref_match2->[$idx]; }
		}
            }
        }
    }

    return $cigar; 
}

sub getSamPair
{
    my $chr = shift;  my $strand = shift;  my $pos = shift;  my $cigar = shift;
    my %parameters = @_;

    my ( $ref_match, $ref_matchSize ) = parseCigar ( $cigar );
    my ( $frag1, $frag2, $frag1Pos, $frag2Pos ) = getBothFragment ( $chr, $strand, $pos, $ref_match, $ref_matchSize );
    if ( ( not $frag1 ) or ( not $frag2 ) ) {
        print STDERR "Skip read that aligns to an intron!\n";
        return 0;
    }
    elsif ( ( defined $parameters{minOverhang} ) and ( ( length ( $frag1 ) < $parameters{minOverhang} ) or ( length ( $frag2 ) < $parameters{minOverhang} ) ) ) {
        print STDERR "Skip read that does not contain enough overhang in either end!\n" if ( $opt_V );
        return 0;
    }

    $frag2 = reverse ( $frag2 );
    my ( $maxScore, $alignment, $intv1s, $intv1e, $intv2s, $intv2e ) = localAlignment ( $frag1, $frag2 );
    if ( not $maxScore ) {
        print STDERR "Skip read that cannot be paired!\n";
        return 0; 
    }

    print STDERR "relative intervals: ", $intv1s, ",", $intv1e, ",", $intv2s, ",", $intv2e, "\n" if ( $_debug );
    if ( $strand eq "+" ) {
        $intv1s += $frag1Pos-1;
        $intv1e += $frag1Pos;
        my $tmpPos = length ( $frag2 )-$intv2s+$frag2Pos;
        $intv2s = length ( $frag2 )-$intv2e+$frag2Pos;
        $intv2e = $tmpPos+1;
    }
    else {
        my $tmpPos = length ( $frag1 )-$intv1s+$frag1Pos;
        $intv1s = length ( $frag1 )-$intv1e+$frag1Pos;
        $intv1e = $tmpPos+1;
        $intv2s += $frag2Pos-1;
        $intv2e += $frag2Pos;
    }
    print STDERR "alignment: ", $alignment, ", absolute intervals: ", $intv1s, ",", $intv1e, ",", $intv2s, ",", $intv2e, "\n" if ( $_debug );

    return ( $alignment, $intv1s, $intv1e, $intv2s, $intv2e );
}


sub getJuncPair
{
    my $chr1 = shift;  my $strand1 = shift;  my $donor = shift;  my $pos1 = shift;  my $cigar1 = shift;
    my $chr2 = shift;  my $strand2 = shift;  my $acceptor = shift;  my $pos2 = shift;  my $cigar2 = shift;
    my %parameters = @_;

    my $isIntron = 0;
    if ( ( $chr1 eq $chr2 ) and ( $strand1 eq $strand2 ) ) {
        if ( ( $strand1 eq "+" ) and ( $donor < $acceptor ) ) { $isIntron = checkJuncIntron ( $chr1, $strand1, $donor, $acceptor ); }
        elsif ( ( $strand1 eq "-" ) and ( $donor > $acceptor ) ) { $isIntron = checkJuncIntron ( $chr1, $strand1, $acceptor, $donor ); }
    }
    if ( $isIntron ) {
        print STDERR "Skip read that aligns to an intron!\n";
        return 0;
    }

    my $frag1 = getOneFragment ( $chr1, $strand1, $pos1, $cigar1 );
    my $frag2 = getOneFragment ( $chr2, $strand2, $pos2, $cigar2 );
    if ( ( defined $parameters{minOverhang} ) and ( ( length ( $frag1 ) < $parameters{minOverhang} ) or ( length ( $frag2 ) < $parameters{minOverhang} ) ) ) {
        print STDERR "Skip read that does not contain enough overhang in either end!\n";
        return 0;
    }

    $frag2 = reverse ( $frag2 );
    my ( $maxScore, $alignment, $intv1s, $intv1e, $intv2s, $intv2e ) = localAlignment ( $frag1, $frag2 );
    if ( not $maxScore ) {
        print STDERR "Skip junction that cannot be paired!\n";
        return 0; 
    }

    print STDERR "relative intervals: ", $intv1s, ",", $intv1e, ",", $intv2s, ",", $intv2e, "\n" if ( $_debug );
    if ( $strand1 eq "+" ) {
        $intv1s += $pos1-1;
        $intv1e += $pos1;
        my $tmpPos = length ( $frag2 )-$intv2s+$pos2;
        $intv2s = length ( $frag2 )-$intv2e+$pos2;
        $intv2e = $tmpPos+1;
    }
    else {
        my $tmpPos = length ( $frag1 )-$intv1s+$pos1;
        $intv1s = length ( $frag1 )-$intv1e+$pos1;
        $intv1e = $tmpPos+1;
        $intv2s += $pos2-1;
        $intv2e += $pos2;
    }
    print STDERR "alignment: ", $alignment, ", absolute intervals: ", $intv1s, ",", $intv1e, ",", $intv2s, ",", $intv2e, "\n" if ( $_debug );

    return ( $alignment, $intv1s, $intv1e, $intv2s, $intv2e );
}

sub checkJuncIntron 
{
    my ( $chr, $strand, $pos1, $pos2 ) = @_;

    my $overlap = [];
    if ( ( defined $global{annotation}{intron_interval}{$chr} ) and ( defined $global{annotation}{intron_interval}{$chr}{$strand} ) ) {
        $overlap = $global{annotation}{intron_interval}{$chr}{$strand}->fetch ( $pos1-$environment{intronFlanking}, $pos2+$environment{intronFlanking} );
    }
    foreach my $interval ( @{$overlap} ) {
        if ( abs ( ( $interval->{start}-$pos1 ) < $environment{intronFlanking} ) and ( abs ( $interval->{end}-$pos2 ) < $environment{intronFlanking} ) ) {
            return 1;
        }
    }

    return 0;
}

sub newPair
{
    my $chr1 = $global{dsPair}[$global{dsPairCount}][3];
    my $strand1 = $global{dsPair}[$global{dsPairCount}][4];
    my $pair1s = $global{dsPair}[$global{dsPairCount}][5];
    my $pair1e = $global{dsPair}[$global{dsPairCount}][6];
    my $chr2 = $global{dsPair}[$global{dsPairCount}][7];
    my $strand2 = $global{dsPair}[$global{dsPairCount}][8];
    my $pair2s = $global{dsPair}[$global{dsPairCount}][9];
    my $pair2e = $global{dsPair}[$global{dsPairCount}][10];

    if ( not defined $global{dspInterval}{$chr1} ) { 
        $global{dspInterval}{$chr1}{"+"} = Set::IntervalTree->new;
        $global{dspInterval}{$chr1}{"-"} = Set::IntervalTree->new;
    }
    my $pair1Overlap = $global{dspInterval}{$chr1}{$strand1}->fetch ( $pair1s, $pair1e );
    if ( not defined $global{dspInterval}{$chr2} ) { 
        $global{dspInterval}{$chr2}{"+"} = Set::IntervalTree->new;
        $global{dspInterval}{$chr2}{"-"} = Set::IntervalTree->new;
    }
    my $pair2Overlap = $global{dspInterval}{$chr2}{$strand2}->fetch ( $pair2s, $pair2e );

    my $contained = 0;  my $id1 = "";  my $id2 = "";
    my $support = 0;  my $start1 = 0; my $end1 = 0; my $start2 =0; my $end2 = 0;
    if ( ( scalar @{$pair1Overlap} ) and  ( scalar @{$pair2Overlap} ) )  { 
        ( $contained, $id1, $id2, $support, $start1, $end1, $start2, $end2 ) = checkContain ( $pair1Overlap, $pair2Overlap, $pair1s, $pair1e, $pair2s, $pair2e ); 
    }

    if ( $contained ) {
        if ( $_debug ) {
            print STDERR join ( "\t", $pair1s, $pair1e, $pair2s, $pair2e, "\n" );
            print STDERR join ( "\t", $contained, $id1, $id2, "\n" );
            print STDERR Dumper $pair1Overlap;
            print STDERR Dumper $pair2Overlap;
        }

        $global{dspInterval}{$chr1}{$strand1}->remove ( $start1, $end1, sub { my $item = shift; return ( $item->{id} == $id1 ); });
        $global{dspInterval}{$chr2}{$strand2}->remove ( $start2, $end2, sub { my $item = shift; return ( $item->{id} == $id2 ); });
        $global{dspInterval}{$chr1}{$strand1}->insert ( { 
                id => $id1, 
                pair => $id2,
                support => $support+1,
                start => $start1,
                end => $end1 }, $start1, $end1 );
        $global{dspInterval}{$chr2}{$strand2}->insert ( { 
                id => $id2,
                pair => $id1,
                support => $support+1,
                start => $start2,
                end => $end2 }, $start2, $end2 );
    }
    else {
        my $overlapped = 0; my %toDelete1 = (); my %toDelete2 = ();
        my $minS1 = 0; my $maxE1 = 0; my $minS2 = 0; my $maxE2 = 0;
        if ( ( scalar @{$pair1Overlap} ) and  ( scalar @{$pair2Overlap} ) )  { 
            ( $overlapped, $minS1, $maxE1, $minS2, $maxE2 ) = checkOverlap ( \%toDelete1, \%toDelete2, $pair1Overlap, $pair2Overlap, $pair1s, $pair1e, $pair2s, $pair2e ); 
        }
        foreach my $id ( keys %toDelete1 ) { $global{dspIntervalCluster}{$id} = $global{dspIntervalClusterCount}; }
        foreach my $id ( keys %toDelete2 ) { $global{dspIntervalCluster}{$id} = $global{dspIntervalClusterCount}+1; }

        if ( $overlapped > 1 ) {
            $global{dspInterval}{$chr1}{$strand1}->remove ( $minS1, $maxE1, sub { my $item = shift; return ( defined $toDelete1{$item->{id}} ); });
            $global{dspInterval}{$chr2}{$strand2}->remove ( $minS2, $maxE2, sub { my $item = shift; return ( defined $toDelete2{$item->{id}} ); });

            $global{dspInterval}{$chr1}{$strand1}->insert ( { 
                    id => $global{dspIntervalClusterCount}, 
                    pair => $global{dspIntervalClusterCount}+1, 
                    support => $overlapped,
                    start => $minS1,
                    end => $maxE1 }, $minS1, $maxE1 );
            $global{dspInterval}{$chr2}{$strand2}->insert ( { 
                    id => $global{dspIntervalClusterCount}+1, 
                    pair => $global{dspIntervalClusterCount}, 
                    support => $overlapped,
                    start => $minS2,
                    end => $maxE2 }, $minS2, $maxE2 );
        }
        else {
            $global{dspInterval}{$chr1}{$strand1}->insert ( { 
                    id => $global{dspIntervalClusterCount}, 
                    pair => $global{dspIntervalClusterCount}+1, 
                    support => 1,
                    start => $pair1s,
                    end => $pair1e }, $pair1s, $pair1e );
            $global{dspInterval}{$chr2}{$strand2}->insert ( { 
                    id => $global{dspIntervalClusterCount}+1, 
                    pair => $global{dspIntervalClusterCount}, 
                    support => 1,
                    start => $pair2s,
                    end => $pair2e }, $pair2s, $pair2e );
        }

        $global{dspIntervalCluster}{$global{dspIntervalClusterCount}} = $global{dspIntervalClusterCount};
        $global{dspIntervalCluster}{$global{dspIntervalClusterCount}+1} = $global{dspIntervalClusterCount}+1;
    }

    return ( $contained, $id1, $id2 );
}

sub checkContain 
{
    ## NEED CHECK: note that intervals are [,) 
    my $intervalList1 = shift;
    my $intervalList2 = shift;
    my ( $minS1, $maxE1, $minS2, $maxE2 ) = @_;

    my %list1 = ();
    for ( my $idx = 0; $idx < scalar @{$intervalList1}; $idx++ ) {
        $list1{$intervalList1->[$idx]{id}} = $idx;
    }
    for ( my $idx = 0; $idx < scalar @{$intervalList2}; $idx++ ) {
        if ( defined $list1{$intervalList2->[$idx]{pair}} ) {
            my $list1Idx = $list1{$intervalList2->[$idx]{pair}};
            my $list2Idx = $idx;

            if ( $minS1 <= $minS2 ) {
                if ( $intervalList1->[$list1Idx]{start} <= $intervalList2->[$list2Idx]{start} ) {
                    if ( ( $intervalList1->[$list1Idx]{start} <= $minS1 )
                            and ( $intervalList1->[$list1Idx]{end} >= $maxE1 )
                            and ( $intervalList2->[$list2Idx]{start} <= $minS2 )
                            and ( $intervalList2->[$list2Idx]{end} >= $maxE2 ) ) {
                        return ( 1, $intervalList1->[$list1Idx]{id}, $intervalList2->[$list2Idx]{id}, $intervalList1->[$list1Idx]{support}, $intervalList1->[$list1Idx]{start}, $intervalList1->[$list1Idx]{end}, $intervalList2->[$list2Idx]{start}, $intervalList2->[$list2Idx]{end} );
                    }
                }
                else {
                    if ( $intervalList1->[$list1Idx]{start} > $intervalList2->[$list2Idx]{start} ) {
                        if ( ( $intervalList1->[$list1Idx]{start} <= $minS1 )
                                and ( $intervalList1->[$list1Idx]{end} >= $maxE1 )
                                and ( $intervalList2->[$list2Idx]{start} <= $minS2 )
                                and ( $intervalList2->[$list2Idx]{end} >= $maxE2 ) ) {
                            return ( 1, $intervalList1->[$list1Idx]{id}, $intervalList2->[$list2Idx]{id}, $intervalList1->[$list1Idx]{support}, $intervalList1->[$list1Idx]{start}, $intervalList1->[$list1Idx]{end}, $intervalList2->[$list2Idx]{start}, $intervalList2->[$list2Idx]{end} );
                        }
                    }
                }
            }
        }
    }

    return 0;
}

sub checkOverlap 
{
    ## NEED CHECK: note that intervals are [,) 
    my $ref_toDelete1 = shift;
    my $ref_toDelete2 = shift;
    my $intervalList1 = shift;
    my $intervalList2 = shift;
    my ( $minS1, $maxE1, $minS2, $maxE2 ) = @_;

    my $overlapped = 1;
    my %list1 = ();
    for ( my $idx = 0; $idx < scalar @{$intervalList1}; $idx++ ) {
        $list1{$intervalList1->[$idx]{id}} = $idx;
    }
    for ( my $idx = 0; $idx < scalar @{$intervalList2}; $idx++ ) {
        if ( defined $list1{$intervalList2->[$idx]{pair}} ) {
            my $list1Idx = $list1{$intervalList2->[$idx]{pair}};
            my $list2Idx = $idx;

            if ( $minS1 <= $minS2 ) {
                if ( $intervalList1->[$list1Idx]{start} <= $intervalList2->[$list2Idx]{start} ) {
                    $overlapped += $intervalList1->[$list1Idx]{support};
                    $ref_toDelete1->{$intervalList1->[$list1Idx]{id}} = 1;
                    $ref_toDelete2->{$intervalList2->[$list2Idx]{id}} = 1;

                    if ( $intervalList1->[$list1Idx]{start} < $minS1 ) { $minS1 = $intervalList1->[$list1Idx]{start}; }
                    if ( $intervalList1->[$list1Idx]{end} > $maxE1 ) { $maxE1 = $intervalList1->[$list1Idx]{end}; }
                    if ( $intervalList2->[$list2Idx]{start} < $minS2 ) { $minS2 = $intervalList2->[$list2Idx]{start}; }
                    if ( $intervalList2->[$list2Idx]{end} > $maxE2 ) { $maxE2 = $intervalList2->[$list2Idx]{end}; }
                }
            }
            else {
                if ( $intervalList1->[$list1Idx]{start} > $intervalList2->[$list2Idx]{start} ) {
                    $overlapped += $intervalList1->[$list1Idx]{support};
                    $ref_toDelete1->{$intervalList1->[$list1Idx]{id}} = 1;
                    $ref_toDelete2->{$intervalList2->[$list2Idx]{id}} = 1;

                    if ( $intervalList1->[$list1Idx]{start} < $minS1 ) { $minS1 = $intervalList1->[$list1Idx]{start}; }
                    if ( $intervalList1->[$list1Idx]{end} > $maxE1 ) { $maxE1 = $intervalList1->[$list1Idx]{end}; }
                    if ( $intervalList2->[$list2Idx]{start} < $minS2 ) { $minS2 = $intervalList2->[$list2Idx]{start}; }
                    if ( $intervalList2->[$list2Idx]{end} > $maxE2 ) { $maxE2 = $intervalList2->[$list2Idx]{end}; }
                }
            }
        }
    }

    return  ( $overlapped, $minS1, $maxE1, $minS2, $maxE2 );
}

sub getOneFragment
{
    my ( $chr, $strand, $pos, $cigar ) = @_;

    my $largestMatch = parseCigar ( $cigar,  getLargestM => 1 );
    my $frag = substr ( $global{genomeSeq}->{$chr}, $pos-1, $largestMatch );
    $frag = reverseComplement ( $frag ) if ( $strand eq "-" );

    print STDERR join ( "\t", $frag, $pos, "\n" ) if ( $_debug );
    return $frag;
}

sub getBothFragment
{
    my ( $chr, $strand, $pos, $ref_match, $ref_matchSize ) = @_;

    my $frag1 = ""; my $frag2 = "";
    my $maxNonMatch = maxND ( $chr, $strand, $pos, $ref_match, $ref_matchSize );
    return ( $frag1, $frag2 ) if ( not $maxNonMatch );

    my $pos1 = $pos; 
    for ( my $idx = 0; $idx <= $maxNonMatch; $idx++ ) {
        if ( ( $ref_match->[$idx] eq "S" ) or ( $ref_match->[$idx] eq "H" ) or ( $ref_match->[$idx] eq "I" ) or ( $ref_match->[$idx] eq "P" ) )  { }
        elsif ( ( $ref_match->[$idx] eq "M" ) or ( $ref_match->[$idx] eq "=" ) or ( $ref_match->[$idx] eq "X" ) ) { 
            $frag1 .= substr ( $global{genomeSeq}->{$chr}, $pos-1, $ref_matchSize->[$idx] );
            $pos += $ref_matchSize->[$idx];
        }
        elsif ( $ref_match->[$idx] eq "D" ) {
            $frag1 .= substr ( $global{genomeSeq}->{$chr}, $pos-1, $ref_matchSize->[$idx] );
            $pos += $ref_matchSize->[$idx];
        }
        elsif ( ( $ref_match->[$idx] eq "N" ) ) { 
            $pos += $ref_matchSize->[$idx];
        };
    }

    my $pos2 = $pos;
    for ( my $idx = $maxNonMatch+1; $idx < scalar ( @{$ref_matchSize} ); $idx++ ) {
        if ( ( $ref_match->[$idx] eq "S" ) or ( $ref_match->[$idx] eq "H" ) or ( $ref_match->[$idx] eq "I" ) or ( $ref_match->[$idx] eq "P" ) )  { }
        elsif ( ( $ref_match->[$idx] eq "M" ) or ( $ref_match->[$idx] eq "=" ) or ( $ref_match->[$idx] eq "X" ) ) { 
            $frag2 .= substr ( $global{genomeSeq}->{$chr}, $pos-1, $ref_matchSize->[$idx] );
            $pos += $ref_matchSize->[$idx];
        }
        elsif ( $ref_match->[$idx] eq "D" ) {
            $frag2 .= substr ( $global{genomeSeq}->{$chr}, $pos-1, $ref_matchSize->[$idx] );
            $pos += $ref_matchSize->[$idx];
        }
        elsif ( ( $ref_match->[$idx] eq "N" ) ) { 
            $pos += $ref_matchSize->[$idx];
        }
    }

    if ( $strand eq "-" ) {
        $frag1 = reverseComplement ( $frag1 );
        $frag2 = reverseComplement ( $frag2 );
    }

    print STDERR join ( "\t", $frag1, $frag2, $pos1, $pos2, "\n" ) if ( $_debug );
    return ( $frag1, $frag2, $pos1, $pos2 );
}

sub maxND
{
    my ( $chr, $strand, $pos, $ref_match, $ref_matchSize ) = @_;

    my $totalLen = eval ( join ( "+", @{$ref_matchSize} ) );
    my $overlap = [];
    if ( ( defined $global{annotation}{intron_interval}{$chr} ) and ( defined $global{annotation}{intron_interval}{$chr}{$strand} ) ) {
        $overlap = $global{annotation}{intron_interval}{$chr}{$strand}->fetch ( $pos-$environment{intronFlanking}, $pos+$totalLen+$environment{intronFlanking} );
    }

    my $maxMatch = 0;  my $tmpMax = 0; 
    my $oldPos = $pos;
    for ( my $idx = 0; $idx < scalar ( @{$ref_matchSize} ); $idx++ ) {
        if ( ( $ref_match->[$idx] eq "S" ) or ( $ref_match->[$idx] eq "H" ) or ( $ref_match->[$idx] eq "I" ) or ( $ref_match->[$idx] eq "P" ) )  { }
        elsif ( ( $ref_match->[$idx] eq "M" ) or ( $ref_match->[$idx] eq "=" ) or ( $ref_match->[$idx] eq "X" ) ) { 
            $pos += $ref_matchSize->[$idx];
        }
        elsif ( ( $ref_match->[$idx] eq "D" ) or ( $ref_match->[$idx] eq "N" ) ) { 
            $oldPos = $pos;
            $pos += $ref_matchSize->[$idx];

            ## not intron
            my $isIntron = 0;
            foreach my $interval ( @{$overlap} ) {
                if ( abs ( ( $interval->{start}-$oldPos ) < $environment{intronFlanking} ) and ( abs ( $interval->{end}-$pos ) < $environment{intronFlanking} ) ) {
                    $isIntron = 1;
                    last;
                }
            }
            if ( ( not $isIntron ) and ( $ref_matchSize->[$idx] > $tmpMax ) ) {
                $tmpMax = $ref_matchSize->[$idx];
                $maxMatch = $idx;
            }
        }
    }

    return $maxMatch;
}

sub sortCluster
{
    my %parameters = @_;

    my %removedCluster = ();
    foreach my $chr ( keys %{$global{dspInterval}} ) {
        foreach my $strand ( keys %{$global{dspInterval}{$chr}} ) {
            my $removed = $global{dspInterval}{$chr}{$strand}->remove ( 1, $environment{maxChr}, 
                sub { my $item = shift; return ( $item->{support} < $parameters{minSupport} ); });
            for ( my $idx = 0; $idx < scalar ( @{$removed} ); $idx++ ) { $removedCluster{$removed->[$idx]{id}} = 1; }
        }
    }

    for ( my $idx = 0; $idx < scalar ( @{$global{dsPair}} ); $idx++ ) {
        next if ( ( defined $removedCluster{$global{dsPair}[$idx][0]} ) or ( defined $removedCluster{$global{dsPair}[$idx][1]} ) );

        my $id1 = $global{dsPair}[$idx][0];
        while ( $global{dspIntervalCluster}{$id1} ne $id1 ) {  $id1 = $global{dspIntervalCluster}{$id1};  }
        my $id2 = $global{dsPair}[$idx][1];
        while ( $global{dspIntervalCluster}{$id2} ne $id2 ) {  $id2 = $global{dspIntervalCluster}{$id2};  }
        if ( $id1 < $id2 ) {
            my $id = $id1 . ":" . $id2;
            if ( defined $global{dsPairCluster}{$id} ) { $global{dsPairCluster}{$id}{set} .= ":" . $idx; }
            else { $global{dsPairCluster}{$id}{set} = $idx; }
        }
        else {
            my $id = $id2 . ":" . $id1;
            if ( defined $global{dsPairCluster}{$id} ) { $global{dsPairCluster}{$id}{set} .= ":-" . $idx; }
            else { $global{dsPairCluster}{$id}{set} = "-" . $idx; }
        }
    }

    foreach my $chr ( keys %{$global{dspInterval}} ) {
        foreach my $strand ( keys %{$global{dspInterval}{$chr}} ) {
            my $all = $global{dspInterval}{$chr}{$strand}->fetch ( 1, $environment{maxChr} );
            for ( my $idx = 0; $idx < scalar ( @{$all} ); $idx++ ) {
                if ( $all->[$idx]{id} < $all->[$idx]{pair} ) {
                    my $id = $all->[$idx]{id} . ":" . $all->[$idx]{pair};
                    next if ( not defined $global{dsPairCluster}{$id}{set} );
                    $global{dsPairCluster}{$id}{chr1} = $chr;
                    $global{dsPairCluster}{$id}{strand1} = $strand;
                    $global{dsPairCluster}{$id}{start1} = $all->[$idx]{start};
                    $global{dsPairCluster}{$id}{end1} = $all->[$idx]{end}-1;
                }
                else {
                    my $id = $all->[$idx]{pair} . ":" . $all->[$idx]{id};
                    next if ( not defined $global{dsPairCluster}{$id}{set} );
                    $global{dsPairCluster}{$id}{chr2} = $chr;
                    $global{dsPairCluster}{$id}{strand2} = $strand;
                    $global{dsPairCluster}{$id}{start2} = $all->[$idx]{start};
                    $global{dsPairCluster}{$id}{end2} = $all->[$idx]{end}-1;
                }
            }
    print Dumper $all;
        }
    }

    1;
}

sub printCluster
{
    my $outputFile = shift;
    my %parameters = @_;

    open ( OUT, ">$outputFile" ) or die ( "Error in opening file $outputFile for writing!\n" );
    print "Output results to file $outputFile...\n\tTime:", `date`;
    my $clusterCount = 0;
    foreach my $cluster ( sort { $global{dsPairCluster}{$a}{chr1} cmp $global{dsPairCluster}{$b}{chr1} or
                                 $global{dsPairCluster}{$a}{strand1} cmp $global{dsPairCluster}{$b}{strand1} or
                                 $global{dsPairCluster}{$a}{start1} <=> $global{dsPairCluster}{$b}{start1} or
                                 $global{dsPairCluster}{$a}{end1} <=> $global{dsPairCluster}{$b}{end1} or
                                 $global{dsPairCluster}{$a}{chr2} cmp $global{dsPairCluster}{$b}{chr2} or
                                 $global{dsPairCluster}{$a}{strand2} cmp $global{dsPairCluster}{$b}{strand2} or
                                 $global{dsPairCluster}{$a}{start2} <=> $global{dsPairCluster}{$b}{start2} or
                                 $global{dsPairCluster}{$a}{end2} <=> $global{dsPairCluster}{$b}{end2} } keys %{$global{dsPairCluster}} ) {
        $clusterCount++;
        if ( $_debug ) { print $cluster; print Dumper $global{dsPairCluster}{$cluster}; }

        my @pairSet = split ( /:/, $global{dsPairCluster}{$cluster}{set} );
        print OUT "Group $clusterCount == position "; 
        print OUT $global{dsPairCluster}{$cluster}{chr1}, "(", $global{dsPairCluster}{$cluster}{strand1}, "):";
        print OUT $global{dsPairCluster}{$cluster}{start1}, "-", $global{dsPairCluster}{$cluster}{end1};
        print OUT "|", $global{dsPairCluster}{$cluster}{chr2}, "(", $global{dsPairCluster}{$cluster}{strand2}, "):";
        print OUT $global{dsPairCluster}{$cluster}{start2}, "-", $global{dsPairCluster}{$cluster}{end2};
        print OUT ", support ", scalar(@pairSet), ".\n----\n"; 
        for ( my $idx = 0; $idx < scalar ( @pairSet ); $idx++ ) {
            my $fragNeed2Switch = 0;
            my $pairIdx = $pairSet[$idx];
            if ( $pairSet[$idx] < 0 ) {
                $pairIdx = -$pairSet[$idx];
                $fragNeed2Switch = 1;
            }
            my $ref_toPrint = printStem ( $global{dsPair}[$pairIdx], $global{dsPairCluster}{$cluster}, $fragNeed2Switch );   
            #print OUT $global{dsPair}[$pairSet[$idx]][11], "\t", $$ref_toPrint, "\n";
            ## DOUBLE CHECK this
            print OUT $global{dsPair}[$pairIdx][11], "\t", $$ref_toPrint, "\n";

            my $pos = $global{dsPair}[$pairIdx][3] . $global{dsPair}[$pairIdx][4] . $global{dsPair}[$pairIdx][13];  # one read can have multiple positions
            if ( not defined $global{dspIntervalClusterRead}{$global{dsPair}[$pairIdx][11]}{$pos} ) { 
                $global{dspIntervalClusterRead}{$global{dsPair}[$pairIdx][11]}{$pos}{dg} = "$clusterCount"; 
                $global{dspIntervalClusterRead}{$global{dsPair}[$pairIdx][11]}{$pos}{cigar} = $global{dsPair}[$pairIdx][12]; 
                if ( $global{dsPair}[$pairIdx][12] =~ /:/ ) { 
                    my ( $cigar, $isChiastic ) = split ( /:/, $global{dspIntervalClusterRead}{$global{dsPair}[$pairIdx][11]}{$pos}{cigar} ); 
                    if ( $isChiastic == 2 ) {
                        my $pos = $global{dsPair}[$pairIdx][3] . $global{dsPair}[$pairIdx][4] . $global{dsPair}[$pairIdx][14];  # one read can have multiple positions
                        $global{dspIntervalClusterRead}{$global{dsPair}[$pairIdx][11]}{$pos}{cigar} = $global{dsPair}[$pairIdx][12]; 
                    }
                }
            }
            else { $global{dspIntervalClusterRead}{$global{dsPair}[$pairIdx][11]}{$pos}{dg} .= ",$clusterCount"; }
        }
        print OUT "\n";
    }
    print "$clusterCount clusteres output to file $outputFile.\n\tTime: ", `date`, "\n";

    if ( ( defined $parameters{supportSam} ) and $parameters{supportSam} ) {
        my $samFile = $outputFile . ".support.sam";
        outputSam ( $samFile, $parameters{inputSam} );
    }

    1;
}

sub outputSam
{
    my $outputFile = shift;
    my $inputSamFileList = shift;
    my %parameters = @_;

    open ( OUT, ">$outputFile" ) or die ( "Error in openning file $outputFile to output supporting reads!\n" );
    print "Output supporting reads to $outputFile...\n\tTime: ", `date`;
    my $lineCount = 0;
    my $validCount = 0;
    my @samFiles = split ( /:/, $inputSamFileList );
    foreach my $samFile ( @samFiles ) { 
        open ( SAM, $samFile ) or die ( "Error in reading sam file $samFile!\n" );
        print "check for supporting reads from sam file $samFile...\n\tTime: ", `date`;
        while ( my $line = <SAM> ) {
            next if ( $line =~ /^#/ );
            if ( $line =~ /^@/ ) { }
            else {
                $lineCount++;
                if ( $lineCount % 1000000 == 0 ) { print "line: $lineCount\n"; print "\tvalid line: $validCount\n\t", `date`; }

                my @data = split ( /\t/, $line );
                my $strand = ( $data[1] & 16 ) ? "-" : "+";
                my $pos = $data[2] . $strand . $data[3];
                if ( defined $global{dspIntervalClusterRead}{$data[0]}{$pos} ) { 
                    $validCount++; 
                    my $cigar = ""; my $isChiastic = 0;
                    if ( $global{dspIntervalClusterRead}{$data[0]}{$pos}{cigar} =~ /:/ ) { ( $cigar, $isChiastic ) = split ( /:/, $global{dspIntervalClusterRead}{$data[0]}{$pos}{cigar} ); }
                    else { $cigar = $global{dspIntervalClusterRead}{$data[0]}{$pos}{cigar}; }
                    if ( $isChiastic ) {
                        $data[1] = 0;
                        $data[9] = reverseRead ( $data[5], $data[9] );
                        $data[10] = reverseRead ( $data[5], $data[10] );
                    }
                    $data[5] = $cigar;
                    chomp $data[-1];
                    print OUT join ( "\t", @data, "DG:i:$global{dspIntervalClusterRead}{$data[0]}{$pos}{dg}", "XG:i:$isChiastic" ), "\n"; 
                }
            }
        }
        close SAM;
        print "in total $lineCount lines read from sam files.\n";
        print "among which $validCount lines generate supports for valid base pairing.\n\tTime: ", `date`, "\n";

    }

    return ( $lineCount, $validCount );
}


sub reverseRead 
{
    my $cigar = shift;
    my $seq = shift;

    my $reverseSeq = "";
    my $leadingS = parseCigar ( $cigar,  getLeadingS => 1 );
    $reverseSeq = substr ( $seq, $leadingS );
    $reverseSeq .= substr ( $seq, 0, $leadingS );

    return $reverseSeq;
}

sub printStem
{
    my $ref_aPair = shift;
    my $ref_aCluster = shift;
    my $needSwitch = shift;

    my $string = "";
    my $frag = $ref_aPair->[2]; $frag =~ s/-//g; $frag =~ s/://g;

    if ( not $needSwitch ) {
        for ( my $idx = $ref_aCluster->{start1}; $idx < $ref_aPair->[5]; $idx++ ) { $string .= "."; }
        for ( my $idx = 0; $idx < ( $ref_aPair->[6]-$ref_aPair->[5] ); $idx++ ) {
            $string .= substr ( $frag, $idx, 1 );
        }
        for ( my $idx = $ref_aPair->[6]; $idx < $ref_aCluster->{end1}; $idx++ ) { $string .= "."; }

        $string .= "    gap    ";
        for ( my $idx = $ref_aCluster->{start2}; $idx < $ref_aPair->[9]; $idx++ ) { $string .= "."; }
        for ( my $idx = ( $ref_aPair->[6]-$ref_aPair->[5] ); $idx < ( $ref_aPair->[6]+$ref_aPair->[10]-$ref_aPair->[5]-$ref_aPair->[9] ); $idx++ ) {
            $string .= substr ( $frag, $idx, 1 );
        }
        for ( my $idx = $ref_aPair->[10]; $idx < $ref_aCluster->{end2}; $idx++ ) { $string .= "."; }
    }
    else {
        for ( my $idx = $ref_aCluster->{start1}; $idx < $ref_aPair->[9]; $idx++ ) { $string .= "."; }
        for ( my $idx = 0; $idx < ( $ref_aPair->[10]-$ref_aPair->[9] ); $idx++ ) {
            $string .= substr ( $frag, $idx, 1 );
        }
        for ( my $idx = $ref_aPair->[10]; $idx < $ref_aCluster->{end1}; $idx++ ) { $string .= "."; }

        $string .= "    gap    ";
        for ( my $idx = $ref_aCluster->{start2}; $idx < $ref_aPair->[5]; $idx++ ) { $string .= "."; }
        for ( my $idx = ( $ref_aPair->[10]-$ref_aPair->[9] ); $idx < ( $ref_aPair->[10]+$ref_aPair->[6]-$ref_aPair->[9]-$ref_aPair->[5] ); $idx++ ) {
            $string .= substr ( $frag, $idx, 1 );
        }
        for ( my $idx = $ref_aPair->[6]; $idx < $ref_aCluster->{end2}; $idx++ ) { $string .= "."; }
    }

    return \$string;
}
