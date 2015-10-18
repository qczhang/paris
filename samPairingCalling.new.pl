#! /usr/bin/perl
#
use strict;
use warnings;

use Data::Dumper;
use Getopt::Std;

use lib "module";
#use PARISutil qw( &loadGTF &loadGenome &parseCigar &reverseComplement &localAlignment );
use PARISutil qw( &loadGenome &parseCigar &reverseComplement &localAlignment );
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

use vars qw ($opt_h $opt_V $opt_D $opt_i $opt_j $opt_s $opt_o $opt_g $opt_a $opt_t $opt_l $opt_p $opt_z $opt_c );
&getopts('hVDi:j:s:o:g:a:t:l:p:z:c:');

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
 -z	genome size file
 -a     annotation file
 -t     transcriptome file

 -l	minimum overhang length for valid mapping
 -p	minimum number of supporting reads 
 -c 	scoring method (harmonic or geometric)
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

#    $global{annotation} = loadGTF ( $parameters{annotationFile} );
    $global{genomeSeq} = loadGenome ( $parameters{genomeFile} );

    my $allSupportSam = "";
    my $duplexGroupBed = "tmp.$$.duplexGroup.bed";
    my $readClusterBed = "tmp.$$.readCluster.bed";
    if ( $samFileList ne "NULL" ) {
        $allSupportSam = $samFileList;
        my @samFiles = split ( /:/, $samFileList );
        foreach my $samFile ( @samFiles ) { genPairClusterFromSamFile ( $samFile, duplexGroup => $duplexGroupBed, readCluster => $readClusterBed, removeRedundancy => 1 ); }
    }
    if ( $chiasticFileList ne "NULL" ) {
        if ( $allSupportSam ) { $allSupportSam = $samFileList . ":" . $chiasticSamFileList;  }
        else {  $allSupportSam = $chiasticSamFileList;  }
        my @chiasticFiles = split ( /:/, $chiasticFileList );
        foreach my $chiasticFile ( @chiasticFiles ) { genPairClusterFromJunctionFile ( $chiasticFile, duplexGroup => $duplexGroupBed, readCluster => $readClusterBed ); }
    }

    my %read_tag = ();
    genDuplexGroup ( \%read_tag, $duplexGroupBed );
    nonOverlappingTag ( \%read_tag, $readClusterBed );
    sortCluster ( minSupport => $parameters{minSupport}, outputBed => 1, inputSam => $allSupportSam, genomeSizeFile => $parameters{genomeSizeFile} );
    printCluster ( $outputFile, supportSam => 1, inputSam => $allSupportSam, method => $parameters{scoringMethod} );

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
    if ( defined $opt_z ) { $parameters{genomeSizeFile} = $opt_z; }
    if ( defined $opt_a ) { $parameters{annotationFile} = $opt_a; }
    if ( defined $opt_t ) { $parameters{transcriptomeFile} = $opt_t; }

    if ( defined $opt_l ) { $environment{minOverhang} = $opt_l; }
    if ( defined $opt_p ) { $parameters{minSupport} = $opt_p; }
    if ( defined $opt_c ) { $parameters{scoringMethod} = $opt_c; }


    return ( %parameters );
}

sub genPairClusterFromSamFile
{
    my $samFile = shift;
    my %parameters = @_;

    if ( $parameters{removeRedundancy} ) {
        print STDERR "Remove redundancy in input SAM file $samFile.\t", `date`;
        my $sortedSamFile = "tmp.$$.sorted.sam";
        my $uniqSamFile = "tmp.$$.uniq.sam";
        print STDERR `grep "^@" $samFile > $sortedSamFile`;
        print STDERR `grep -v "^@" $samFile | sort -k3,3 -k4,4n -k10,10 >> $sortedSamFile`;
        uniqSam ( $sortedSamFile, $uniqSamFile );
        $samFile = $uniqSamFile;
    }

    my $lineCount = 0;
    my $validCount = 0;
    open ( SAM, $samFile ) or die ( "Error in reading sam file $samFile!\n" );
    open ( DG, ">>$parameters{duplexGroup}" ) or die ( "Error in opening $parameters{duplexGroup} for output read duplex groups!\n" );
    open ( RC, ">>$parameters{readCluster}" ) or die ( "Error in opening $parameters{readCluster} for output read clusters!\n" );
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
            my ( $duplexStemLine, $duplexIntervalLine ) = genPairClusterFromSamLine ( $line );
            if ( $duplexStemLine ) {
                print DG $duplexStemLine;
                print RC $duplexIntervalLine;
                $validCount++;
            }
        }
    }
    close SAM;
    close RC;
    print "$lineCount lines read from sam file $samFile.\n";
    print "among which $validCount lines generate supports for base pairing.\n\tTime: ", `date`, "\n";

    return ( $lineCount, $validCount );
}

sub uniqSam 
{
    my $sortedSamFile = shift;
    my $uniqSamFile = shift;

    open ( SAM, $sortedSamFile ) or die ( "Error in reading sam file $sortedSamFile!\n" );
    open ( OUT, ">$uniqSamFile" ) or die ( "Error in opening $uniqSamFile for output uniq SAM!\n" );
    print "read sam file $sortedSamFile...\n\tTime: ", `date`;
    my $seqName = "";  my $pos = 0; my $seq = "";
    my $lineCount = 0; my $uniqCount = 0;
    while ( my $line = <SAM> ) {
        next if ( $line =~ /^#/ );
        if ( $line =~ /^@/ ) { print OUT $line; }
        else {
            $lineCount++;
            if ( $lineCount % 100000 == 0 ) { print "line: $lineCount\n", `date`; }

            #last if ( $lineCount > 10 );
            my @data = split ( /\t/, $line );
            if ( ( $data[9] ne $seq ) or ( $data[3] != $pos ) or ( $data[2] ne $seqName ) ) {
                print OUT $line;
                $uniqCount++;
            }
            $seqName = $data[2];
            $pos = $data[3];
            $seq = $data[9];
        }
    }
    close SAM;
    close OUT;

    print "$lineCount lines read from sam file $sortedSamFile.\n";
    print "among which $uniqCount lines are unique.\n\tTime: ", `date`, "\n";

    1;
}

sub genPairClusterFromSamLine
{
    my $line = shift;

    my @data = split ( /\t/, $line );
    if ( $data[5] !~ /[ND]/ ) {
        print STDERR "Skip read that does not contain cleavage!\n\t$line\n" if ( $opt_V );
        return 0;
    }

    my $strand = ( $data[1] & 16 ) ? "-" : "+";
    my ( $alignment, $pair1s, $pair1e, $pair2s, $pair2e ) = getSamPair ( $data[2], $strand, $data[3], $data[5], minOverhang => $environment{minOverhang} );
    if ( not $alignment ) {
        print STDERR "No valid duplex: \t$line\n" if ( $opt_V );
        return 0;
    }
    if ( ( $pair1s == $pair2s ) or ( $pair1e == $pair2e ) 
            or ( ( $pair1s < $pair2e ) and ( $pair1e > $pair2s ) ) )  {
        print STDERR "inapproprieate duplex: \t$line\n" if ( $opt_V );
        return 0;
    }

    return 0 if ( ( $pair1s == $pair2s ) or ( $pair1e == $pair2e ) );
    my $stemBed = join ( "\t", $data[2], $pair1s, $pair1e, $data[0], "1", $strand ) . "\n";
    $stemBed .= join ( "\t", $data[2], $pair2s, $pair2e, $data[0], "2", $strand ) . "\n";
    my $intervalBed = join ( "\t", $data[2], $pair1s, $pair2e, $data[0], ".", $strand ) . "\n";
    
    return ( $stemBed, $intervalBed );

    1;
}

sub genPairClusterFromJunctionFile
{
    my $junctionFile = shift;
    my %parameters = @_;

    if ( $parameters{removeRedundancy} ) {
        print STDERR "Remove redundancy in input Junction file $junctionFile.\t", `date`;
        my $sortedJunctionFile = "tmp.$$.sorted.junction";
        my $uniqJunctionFile = "tmp.$$.uniq.junction";
        print STDERR `sort -k1,1 -k4,4 -k3,3 -k6,6 -k11,11n -k13,13n -k12,12 -k14,14 $junctionFile > $sortedJunctionFile`;
        uniqJunction ( $sortedJunctionFile, $uniqJunctionFile );
        $junctionFile = $uniqJunctionFile;
    }

    my $lineCount = 0;
    my $validCount = 0;
    open ( JUNC, $junctionFile ) or die ( "Error in reading junction file $junctionFile!\n" );
    open ( DG, ">>$parameters{duplexGroup}" ) or die ( "Error in opening $parameters{duplexGroup} for output read duplex groups!\n" );
    open ( RC, ">>$parameters{readCluster}" ) or die ( "Error in opening $parameters{readCluster} for output read clusters!\n" );
    print "read junction file $junctionFile...\n\tTime: ", `date`;
    while ( my $line = <JUNC> ) {
        next if ( $line =~ /^#/ );
        $lineCount++;
        if ( $lineCount % 100000 == 0 ) {
            print "line: $lineCount\n";
            print "\tvalid line: $validCount\n\t", `date`;
        }

        print "line: $lineCount\n\t", `date` if ( $lineCount % 100000 == 0 );
        my ( $duplexStemLine, $duplexIntervalLine ) = genPairClusterFromOneJunction ( $line );
        if ( $duplexStemLine ) {
            print DG $duplexStemLine;
            print RC $duplexIntervalLine;
            $validCount++;
        }
    }
    close JUNC;
    close RC;
    print "$lineCount lines read from junction file $junctionFile.\n";
    print "among which $validCount lines generate supports for base pairing.\n\tTime: ", `date`, "\n";

    return $lineCount;
}

sub uniqJunction 
{
    my $sortedJunctionFile = shift;
    my $uniqJunctionFile = shift;

    open ( JUNC, $sortedJunctionFile ) or die ( "Error in reading junction file $sortedJunctionFile!\n" );
    open ( OUT, ">$uniqJunctionFile" ) or die ( "Error in opening $uniqJunctionFile for output uniq junction file!\n" );
    print "read junction file $sortedJunctionFile...\n\tTime: ", `date`;
    my $seqName1 = "";  my $strand1 = "";  my $pos1 = 0;  my $cigar1 = "";
    my $seqName2 = "";  my $strand2 = "";  my $pos2 = 0;  my $cigar2 = "";
    my $lineCount = 0; my $uniqCount = 0;
    while ( my $line = <JUNC> ) {
        next if ( $line =~ /^#/ );
        if ( $line =~ /^@/ ) { print OUT $line; }
        else {
            $lineCount++;
            if ( $lineCount % 100000 == 0 ) { print "line: $lineCount\n", `date`; }

            #last if ( $lineCount > 10 );
            my @data = split ( /\t/, $line );
            if ( ( $data[10] != $pos1 ) or ( $data[12] != $pos2 ) 
                    or ( $data[11] ne $cigar1 ) or ( $data[13] ne $cigar2 ) 
                    or ( $data[0] ne $seqName1 ) or ( $data[3] ne $seqName2 ) 
                    or ( $data[2] ne $strand1 ) or ( $data[5] ne $strand2 ) ) {
                print OUT $line;
                $uniqCount++;
            }
            $seqName1 = $data[0]; $seqName2 = $data[3];
            $strand1 = $data[2]; $strand2 = $data[5];
            $pos1 = $data[10]; $pos2 = $data[12];
            $cigar2 = $data[13]; $cigar1 = $data[11];
        }
    }
    close JUNC;
    close OUT;

    print "$lineCount lines read from junction file $sortedJunctionFile.\n";
    print "among which $uniqCount lines are unique.\n\tTime: ", `date`, "\n";

    1;
}


sub genPairClusterFromOneJunction
{
    my $line = shift;

    chomp $line;
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
    if ( not $cigar ) {  print STDERR "\tSkip line of inapproprieate alignment: $line\n" if ( $opt_V );  return 0;  }

    my ( $alignment, $pair1s, $pair1e, $pair2s, $pair2e ) = getJuncPair ( $data[0], $data[2], $data[1], $data[10], $data[11], $data[3], $data[5], $data[4], $data[12], $data[13], minOverhang => $environment{minOverhang} );
    if ( not $alignment ) {
        print STDERR "No valid duplex: \t$line\n" if ( $opt_V );
        return 0;
    }
    if ( ( $pair1s == $pair2s ) or ( $pair1e == $pair2e ) 
            or ( ( $pair1s < $pair2e ) and ( $pair1e > $pair2s ) ) ) {
        print STDERR "inapproprieate duplex: \t$line\n" if ( $opt_V );
        return 0;
    }

    my $stemBed = join ( "\t", $data[0], $pair1s, $pair1e, $data[9], "1", $data[2] ) . "\n";
    $stemBed .= join ( "\t", $data[3], $pair2s, $pair2e, $data[9], "2", $data[5] ) . "\n";
    my $intervalBed = $stemBed;
    if ( ( $data[0] eq $data[3] ) and ( $data[2] eq $data[5] ) ) { 
        if ( $pair1s < $pair2s ) { $intervalBed = join ( "\t", $data[0], $pair1s, $pair2e, $data[9], ".", $data[2] ) . "\n"; }
        else { $intervalBed = join ( "\t", $data[0], $pair2s, $pair1e, $data[9], ".", $data[2] ) . "\n"; }
    }

    return ( $stemBed, $intervalBed );
    
    1;
}

sub genDuplexGroup 
{
    my $ref_read_tag = shift;
    my $duplexGroupBedFile = shift;

    my $sortedDuplexGroupBedFile = $duplexGroupBedFile . ".sorted";
    my $uniqDuplexGroupBedFile = $duplexGroupBedFile . ".uniq";
    my $duplexGroupFile = $duplexGroupBedFile . ".intersect";
    my $duplexConnectFile = $duplexGroupBedFile . ".connect";

    print STDERR `sort -k1,1 -k6,6 -k2,2n -k3,3n $duplexGroupBedFile -o $sortedDuplexGroupBedFile`;
    uniqBed ( $sortedDuplexGroupBedFile, $uniqDuplexGroupBedFile, sorted => 1);
    print STDERR `bedtools intersect -a $uniqDuplexGroupBedFile -b $uniqDuplexGroupBedFile -wa -wb -s -r -f 0.50 > $duplexGroupFile`;

    processIntersect ( $duplexGroupFile, $duplexConnectFile, minOverlap => 10 );
    ## generate proper tags for reads in $ref_read_tag 
}

sub processIntersect
{
    my $duplexGroupFile = shift;
    my $duplexConnectFile = shift;
    my %parameters = @_;

    print "Process intersect file $duplexGroupFile to generate connect file $duplexConnectFile\n";
    my %read_connect = ();
    open ( IN, $duplexGroupFile ) or die "Cannot open $duplexGroupFile for reading!\n";
    print "  Open intersect file $duplexGroupFile for reading\n";
    my $tmpConnect = "tmp.$$.connect";
    my $lineCount = 0;

    my %content = ();
    while ( my $line = <IN> ) {
        chomp $line;
        $lineCount++;
        if ( $lineCount % 10000 == 0 ) { 
            print "line: $lineCount\t", `date`; 

            for ( my $idx = 0; $idx < 10; $idx++ ) {
                my $file = "$tmpConnect.$idx";
                open ( OUT, ">>$file" );
                print OUT $content{$idx};
                close OUT;
                $content{$idx} = "";
            }
        }

        my @data = split ( /\t/, $line );
        my @reads1 = split ( /;/, $data[3] );
        my @reads2 = split ( /;/, $data[9] );
        foreach my $read1 ( @reads1 ) { 
            foreach my $read2 ( @reads2 ) { 
                next if ( $read1 eq $read2 );
                my $fileTag = substr ( $read1, -1 );
                $content{$fileTag} .= $read1 . "\t" . $read2 . "\n";
            } 
        }
    }
    close IN;

    for ( my $idx = 0; $idx < 10; $idx++ ) {
        next if ( not -e "$tmpConnect.$idx" );
        my $tmpConnectSorted = "tmp.$$.connect.sorted.$idx";
        print STDERR `sort -i $tmpConnect.$idx -o $tmpConnectSorted`;

        open ( IN, $tmpConnectSorted ) or die "Cannot open $tmpConnectSorted for reading!\n";
        open ( TCC, ">>$duplexConnectFile" ) or die "Cannot open $duplexConnectFile for writing!\n";
        print "  Open connect file $duplexConnectFile for writing\n";
        my $connectCount = 0;
        my $lastLine = "";
        while ( my $line = <IN> ) {
            if ( $line eq $lastLine ) {
                $connectCount++;
                $lastLine = $line;
            }
            else {
                if ( $lastLine ) {
                    chomp $lastLine;
                    print TCC $lastLine, "\t", $connectCount, "\n";
                    $lastLine = $line;
                    $connectCount = 1;
                }
            }
        }
        print TCC $lastLine, "\t", $connectCount, "\n";
        close IN;
        close TCC;
    }
}

sub nonOverlappingTag 
{
    my $ref_read_tag = shift;
    my $readClusterBedFile = shift;

    my $sortedReadClusterBedFile = $readClusterBedFile . ".sorted";
    my $uniqReadClusterBedFile = $readClusterBedFile . ".uniq";
    my $readClusterFile = $readClusterBedFile . ".cluster";

    print STDERR `sort -k1,1 -k6,6 -k2,2n -k3,3n $readClusterBedFile -o $sortedReadClusterBedFile`;
    uniqBed ( $sortedReadClusterBedFile, $uniqReadClusterBedFile, sorted => 1 );
    print STDERR `bedtools cluster -i $uniqReadClusterBedFile -s > $readClusterFile`;

    ## generate proper tags for reads in $ref_read_tag 
}

sub uniqBed 
{
    my $inputBed = shift;
    my $uniqBed = shift;
    my %parameters = @_;

    my $bedPos = "";  my $tag = "";  my $bedStrand = "";  my $label = "";  my $count = 0;
    open ( IN, $inputBed ) or die "Cannot open $inputBed for reading!\n";
    open ( OUT, ">$uniqBed" ) or die "Cannot open $uniqBed for writing!\n";
    if ( $parameters{sorted} ) {
        while ( my $line = <IN> ) {
            chomp $line;
            my @data = split ( /\t/, $line );
            my $tmpPos = join ( "\t", $data[0], $data[1], $data[2] );
            if ( ( $tmpPos ne $bedPos ) or ( $data[5] ne $bedStrand ) ) {
                if ( $bedPos ) { print OUT join ( "\t", $bedPos, $tag, $count, $bedStrand ), "\n"; }
                $bedPos = $tmpPos;
                $bedStrand = $data[5];
                $tag = $data[3];
                $label = $data[4];
                $count = 1;
            }
            else { $tag .= ";" . $data[3]; $count++;  }
        }
        if ( $bedPos ) { print OUT join ( "\t", $bedPos, $tag, $count, $bedStrand ), "\n"; }
    }
    else {
        print STDERR "not implemented!\n";
        return 0;
    }
    close IN;
    close OUT;
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
        print STDERR "Skip read that aligns to an intron!\n" if ( $opt_V );
        return 0;
    }
    elsif ( ( defined $parameters{minOverhang} ) and ( ( length ( $frag1 ) < $parameters{minOverhang} ) or ( length ( $frag2 ) < $parameters{minOverhang} ) ) ) {
        print STDERR "Skip read that does not contain enough overhang in either end!\n" if ( $opt_V );
        return 0;
    }

    if ( defined $parameters{useLongestPairing} ) {
        ## the following is to use longest pairing region to replace the total fragment. 
        $frag2 = reverse ( $frag2 );
        my ( $maxScore, $alignment, $intv1s, $intv1e, $intv2s, $intv2e ) = localAlignment ( $frag1, $frag2 );
        if ( not $maxScore ) {
            print STDERR "Skip read that cannot be paired!\n" if ( $opt_V );
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
    else {
        return ( 1, $frag1Pos, $frag1Pos+length($frag1), $frag2Pos, $frag2Pos+length($frag2) );
    }
}


sub getJuncPair
{
    my $chr1 = shift;  my $strand1 = shift;  my $donor = shift;  my $pos1 = shift;  my $cigar1 = shift;
    my $chr2 = shift;  my $strand2 = shift;  my $acceptor = shift;  my $pos2 = shift;  my $cigar2 = shift;
    my %parameters = @_;

    my $isIntron = 0;
    if ( ( $chr1 eq $chr2 ) and ( $strand1 eq $strand2 ) ) {
#        if ( ( $strand1 eq "+" ) and ( $donor < $acceptor ) ) { $isIntron = checkJuncIntron ( $chr1, $strand1, $donor, $acceptor ); }
#        elsif ( ( $strand1 eq "-" ) and ( $donor > $acceptor ) ) { $isIntron = checkJuncIntron ( $chr1, $strand1, $acceptor, $donor ); }
    }
    if ( $isIntron ) {
        print STDERR "Skip read that aligns to an intron!\n" if ( $opt_V );
        return 0;
    }

    my $frag1 = getOneFragment ( $chr1, $strand1, $pos1, $cigar1 );
    my $frag2 = getOneFragment ( $chr2, $strand2, $pos2, $cigar2 );
    if ( ( defined $parameters{minOverhang} ) and ( ( length ( $frag1 ) < $parameters{minOverhang} ) or ( length ( $frag2 ) < $parameters{minOverhang} ) ) ) {
        print STDERR "Skip read that does not contain enough overhang in either end!\n" if ( $opt_V );
        return 0;
    }

    if ( defined $parameters{useLongestPairing} ) {
        ## the following is to use longest pairing region to replace the total fragment. 
        $frag2 = reverse ( $frag2 );
        my ( $maxScore, $alignment, $intv1s, $intv1e, $intv2s, $intv2e ) = localAlignment ( $frag1, $frag2 );
        if ( not $maxScore ) {
            print STDERR "Skip junction that cannot be paired!\n" if ( $opt_V );
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
    else {
        return ( 1, $pos1, $pos1+length($frag1), $pos2, $pos2+length($frag2) );
    }
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

    my $tmpClusterBed = "tmp.$$.cluster.bed";
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

            print Dumper $all if ( $_debug );
        }
    }

    if ( defined $parameters{outputBed} ) {
        my $posBed = "tmp.$$.cluster.pos.bed";
        my $negBed = "tmp.$$.cluster.neg.bed";
        open ( POS, ">$posBed" ); 
        open ( NEG, ">$negBed" ); 
        foreach my $cluster ( keys %{$global{dsPairCluster}} ) {
            if ( $global{dsPairCluster}{$cluster}{strand1} eq "+" ) { print POS join ( "\t", $global{dsPairCluster}{$cluster}{chr1}, $global{dsPairCluster}{$cluster}{start1}, $global{dsPairCluster}{$cluster}{end1}, $global{dsPairCluster}{$cluster}{strand1}, $cluster ), "\n"; }
            else { print NEG join ( "\t", $global{dsPairCluster}{$cluster}{chr1}, $global{dsPairCluster}{$cluster}{start1}, $global{dsPairCluster}{$cluster}{end1}, $global{dsPairCluster}{$cluster}{strand1}, $cluster ), "\n"; }

            if ( $global{dsPairCluster}{$cluster}{strand2} eq "+" ) { print POS join ( "\t", $global{dsPairCluster}{$cluster}{chr2}, $global{dsPairCluster}{$cluster}{start2}, $global{dsPairCluster}{$cluster}{end2}, $global{dsPairCluster}{$cluster}{strand2}, $cluster ), "\n"; }
            else { print NEG join ( "\t", $global{dsPairCluster}{$cluster}{chr2}, $global{dsPairCluster}{$cluster}{start2}, $global{dsPairCluster}{$cluster}{end2}, $global{dsPairCluster}{$cluster}{strand2}, $cluster ), "\n"; }
        }
        close POS;
        close NEG;

        mergeSam ( "tmp.$$", $parameters{inputSam}, outputBam => 1 );
        print STDERR `bedtools genomecov -ibam "tmp.bam" -g $parameters{genomeSizeFile} -bg -strand + > "tmp.$$.pos.genomeCov"`; 
        print STDERR `bedtools genomecov -ibam "tmp.bam" -g $parameters{genomeSizeFile} -bg -strand - > "tmp.$$.neg.genomeCov"`; 
        print STDERR `bedtools intersect -wb -a "tmp.$$.pos.genomeCov" -b $posBed > "tmp.$$.pos.intCov"`; 
        addGenomeScore ( "tmp.$$.pos.intCov" );
        print STDERR `bedtools intersect -wb -a "tmp.$$.neg.genomeCov" -b $negBed > "tmp.$$.neg.intCov"`; 
        addGenomeScore ( "tmp.$$.neg.intCov" );
    }

    1;
}

sub addGenomeScore 
{
    my $intCovFile = shift;
    open ( ICF, $intCovFile );
    while ( my $line = <ICF> ) {
        chomp $line;
        my @data = split ( /\t/, $line );
        if ( not defined $global{dsPairCluster}{$data[8]}{$data[4]}{$data[7]}{$data[5]} ) {
            $global{dsPairCluster}{$data[8]}{$data[4]}{$data[7]}{$data[5]} = $data[3];
        }
        elsif ( $global{dsPairCluster}{$data[8]}{$data[4]}{$data[7]}{$data[5]} < $data[3] ) {
            $global{dsPairCluster}{$data[8]}{$data[4]}{$data[7]}{$data[5]} = $data[3];
        }
    }
    close ICF;
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

        my $support = scalar(@pairSet);
        my $left = $global{dsPairCluster}{$cluster}{$global{dsPairCluster}{$cluster}{chr1}}{$global{dsPairCluster}{$cluster}{strand1}}{$global{dsPairCluster}{$cluster}{start1}};
        my $right = $global{dsPairCluster}{$cluster}{$global{dsPairCluster}{$cluster}{chr2}}{$global{dsPairCluster}{$cluster}{strand2}}{$global{dsPairCluster}{$cluster}{start2}};
        my $score = supportScore ( $support, $left, $right, $parameters{method} );
        print OUT ", support $support, left $left, right $right, score $score.\n----\n";

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

sub supportScore
{
  my ( $support, $left, $right, $method ) = @_;

  my $score = -1;
  if ( $method eq "geometric" ) { $score = ( $support + 1 ) / sqrt ( ( $left + 1 ) * ( $right + 1 ) ); }
  elsif ( $method eq "harmonic" ) { $score = ( $support + 1 ) * ( 1/ ($left+1) + 1/ ($right+1) ) / 2; }
  else { $score = 0; }
  $score = 1 if ( $score > 1 );

  return $score;
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

sub mergeSam
{
    my $outputFile = shift;
    my $inputSamFileList = shift;
    my %parameters = @_;

    my $samOutputFile = ( $parameters{outputBam}) ? $outputFile . ".sam" : $outputFile;
    open ( OUT, ">$samOutputFile" ) or die ( "Error in openning file $outputFile to output merged SAM files!\n" );
    print "Output merged SAM files to $samOutputFile...\n\tTime: ", `date`;
    my $lineCount = 0;
    my $firstSam = 1;
    my @samFiles = split ( /:/, $inputSamFileList );
    foreach my $samFile ( @samFiles ) { 
        open ( SAM, $samFile ) or die ( "Error in reading sam file $samFile!\n" );
        print "check for supporting reads from sam file $samFile...\n\tTime: ", `date`;
        while ( my $line = <SAM> ) {
            next if ( $line =~ /^#/ );
            if ( $line =~ /^@/ ) { print OUT $line if ( $firstSam ); }
            else {
                $lineCount++;
                if ( $lineCount % 1000000 == 0 ) { print "line: $lineCount\n\t", `date`; }

                print OUT $line;
            }
        }
        close SAM;
        print "in total $lineCount lines read from sam files.\n";

        $firstSam = 0;
    }
    close OUT;
    if ( $parameters{outputBam} ) { print STDERR `samtools view -bS $samOutputFile | samtools sort - $outputFile`; }

    return ( $lineCount );
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
