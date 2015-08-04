#! /usr/bin/perl
#
use strict;
use warnings;

use Data::Dumper;
use Time::HiRes;
use Memory::Usage;
use Set::IntervalTree;
#
##--------------------------------------------------
#
my $_debug = 0;

my %enviroment = (
);

my %global = (
);

#
##--------------------------------------------------
#
my $inputSam = shift;
my $readGroup = shift;
my $output = shift;

&main ( $inputSam, $readGroup, $output );
sub main
{
    my $samFileList = shift;
    my $readGroupFile = shift;
    my $outputFile = shift;
    my %parameters = @_;

    my $memoryUsage = Memory::Usage->new(); $memoryUsage->record('starting work');
    my %read_group = ();
    my %read_start = ();
    my %read_end = ();
    my $totalGroup = loadReadGroup ( $readGroupFile, \%read_group, \%read_start, \%read_end );
    $memoryUsage->record('Memory after load read group file'); $memoryUsage->dump();

    my $allSupportSam = "";
    if ( $samFileList ne "NULL" ) { outputSam ( $outputFile, $samFileList, \%read_group, \%read_start, \%read_end ); }
    $memoryUsage->record('final memory usage'); $memoryUsage->dump();

    1;
}

## ----------------------------------
sub loadReadGroup
{
    my $readGroupFile = shift;
    my $ref_read_group = shift;
    my $ref_read_start = shift;
    my $ref_read_end = shift;

    my $readGroupCount = 0;
    my %readGroup_interval = ();
    open ( RG, $readGroupFile ) or die "Cannot open $readGroupFile for reading!\n";
    OUTER: while ( my $line = <RG> ) {
        next if ( $line =~ /^#/ );
        if ( $line =~ /^Group/ ) {
            $readGroupCount++;
            my ( $chr1, $strand1, $start1, $end1, $chr2, $strand2, $start2, $end2, $support ) = ( $line =~ /position (.+)\(([+-])\):(\d+)-(\d+)\|(.+)\(([+-])\):(\d+)-(\d+), support (\d+)\./ );

            my $id = 0;  my $pos1 = 0;  my $pos2 = 0;
            if ( ( $chr1 eq $chr2 ) and ( $strand1 eq $strand2 ) ) {
                my $start = ( $start1 < $start2 ) ? $start1 : $start2;
                my $end = ( $end1 > $end2 ) ? $end1 : $end2;
                $pos1 = int ( ( $start1 + $end1 ) /2 );
                $pos2 = int ( ( $start2 + $end2 ) /2 );

                if ( not defined $readGroup_interval{$chr1}{$strand1} ) { $readGroup_interval{$chr1}{$strand1} = Set::IntervalTree->new; }
                my $overlap = $readGroup_interval{$chr1}{$strand1}->fetch ( $start, $end );
                my %overlapID = ();
                foreach my $o ( @{$overlap} ) { $overlapID{$o->{id}} = 1; }
                while ( defined $overlapID{$id} ) { $id++;  }
                $readGroup_interval{$chr1}{$strand1}->insert ( { 
                        id => $id, 
                        start => $start,
                        end => $end }, $start, $end );
            }
            else {
                if ( not defined $readGroup_interval{$chr1}{$strand1} ) { $readGroup_interval{$chr1}{$strand1} = Set::IntervalTree->new; }
                my $overlap1 = $readGroup_interval{$chr1}{$strand1}->fetch ( $start1, $end1 );
                if ( not defined $readGroup_interval{$chr2}{$strand2} ) { $readGroup_interval{$chr2}{$strand2} = Set::IntervalTree->new; }
                my $overlap2 = $readGroup_interval{$chr2}{$strand2}->fetch ( $start2, $end2 );
                my %overlapID = ();
                foreach my $o ( @{$overlap1} ) { $overlapID{$o->{id}} = 1; }
                foreach my $o ( @{$overlap2} ) { $overlapID{$o->{id}} = 1; }
                while ( defined $overlapID{$id} ) { $id++;  }
                $readGroup_interval{$chr1}{$strand1}->insert ( { 
                        id => $id, 
                        start => $start1,
                        end => $end1 }, $start1, $end1 );
            }
            $line = <RG>;
            INNER: while ( $line =<RG> ) {
                last INNER if ( not ( $line =~ /\S/ ) );
                my ( $readID ) = split ( /\t/, $line );
                $ref_read_group->{$readID} = $id;
                if ( $pos1 ) {
                    $ref_read_start->{$readID} = $pos1;
                    $ref_read_end->{$readID} = $pos2;
                }
            }
        }
    }
    close RG;
    print "in total $readGroupCount read groups read from $readGroupFile.\n";

    return $readGroupCount;
}

sub outputSam
{
    my $outputFile = shift;
    my $inputSamFileList = shift;
    my $ref_read_group = shift;
    my $ref_read_start = shift;
    my $ref_read_end = shift;

    open ( OUT, ">$outputFile" ) or die ( "Error in openning file $outputFile to output supporting reads!\n" );
    print "Output supporting reads to $outputFile...\n\tTime: ", `date`;
    my $lineCount = 0;
    my @samFiles = split ( /:/, $inputSamFileList );
    foreach my $samFile ( @samFiles ) { 
        open ( SAM, $samFile ) or die ( "Error in reading sam file $samFile!\n" );
        my $pseudoSam = $samFile; $pseudoSam =~ s/.sam//; $pseudoSam .= ".pseudo.sam";
        open ( PS, ">$pseudoSam" ) or die ( "Error in openning file $pseudoSam to output pseudo sam file!\n" );

        print "check for supporting reads from sam file $samFile...\n\tTime: ", `date`;
        while ( my $line = <SAM> ) {
            next if ( $line =~ /^#/ );
            if ( $line =~ /^@/ ) { }
            else {
                $lineCount++;
                if ( $lineCount % 1000000 == 0 ) { print "line: $lineCount\n", `date`; }

                chomp $line;
                my @data = split ( /\t/, $line );
                next if ( not defined $ref_read_group->{$data[0]} );

                print OUT $line, "\tNO:i:$ref_read_group->{$data[0]}\n";
                if ( defined $ref_read_start->{$data[0]} ) {
                    my $start = $ref_read_start->{$data[0]};
                    my $end = $ref_read_end->{$data[0]};

                    my $lenN = $end - $start;
                    $data[3] = $start;
                    $data[5] = "1M" . $lenN . "N1M";
                    $data[9] = "NN";
                    $data[10] = "BB";
                    print PS join "\t", @data, "\tNO:i:$ref_read_group->{$data[0]}\n";
                }
            }
        }
        close SAM;
        close OUT;
        close PS;
        print "in total $lineCount lines read from sam files.\n";

    }

    return ( $lineCount );
}


sub _reverseComplement
{
    my $inseq = shift;
    my $outseq = reverse ( $inseq );

    $outseq =~ tr/AGTCagtc/TCAGtcag/;
    return $outseq;
}

sub _parseCigar
{
    my $cigar = shift;
    my %parameters = @_;

    my @match = split ( /[0-9]+/, $cigar );             # CIGAR: \*|([0-9]+[MIDNSHPX=])+
    shift @match;
    my @matchSize = split ( /[MIDNSHPX=]/, $cigar );

    if ( $_debug ) {
        print STDERR "CIGAR\t", $cigar, "\nmatchOper";
        _printArray ( \@match );
        print STDERR "\nmatchSize";
        _printArray ( \@matchSize );
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

