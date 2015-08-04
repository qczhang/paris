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
my $readGroup = shift;
my $output = shift;

&main ( $readGroup, $output );
sub main
{
    my $readGroupFile = shift;
    my $outputFile = shift;
    my %parameters = @_;

    my $memoryUsage = Memory::Usage->new(); $memoryUsage->record('starting work');
    my $totalGroup = loadReadGroup ( $readGroupFile, $output );

    $memoryUsage->record('final memory usage'); $memoryUsage->dump();

    1;
}

## ----------------------------------
sub loadReadGroup
{
    my $readGroupFile = shift;
    my $outputFile = shift;

    my $readGroupCount = 0;
    my %readGroup_interval = ();
    open ( RG, $readGroupFile ) or die "Cannot open $readGroupFile for reading!\n";
    open ( OUT, ">$outputFile" ) or die "Cannot open $outputFile for writing!\n";
    OUTER: while ( my $line = <RG> ) {
        next if ( $line =~ /^#/ );
        if ( $line =~ /^Group/ ) {
            $readGroupCount++;
            my ( $chr1, $strand1, $start1, $end1, $chr2, $strand2, $start2, $end2, $support ) = ( $line =~ /position (.+)\(([+-])\):(\d+)-(\d+)\|(.+)\(([+-])\):(\d+)-(\d+), support (\d+)\./ );
            if ( $support < 2 ) {
                $line = <RG>;
                INNER: while ( $line =<RG> ) {
                    last INNER if ( not ( $line =~ /\S/ ) );
                }
            }
            elsif ( ( $start1 == $end1 ) or ( $start2 == $end2 ) ) {
                $line = <RG>;
                INNER: while ( $line =<RG> ) {
                    last INNER if ( not ( $line =~ /\S/ ) );
                }
            }
            else {
                print OUT $line;

                $line = <RG>;
                print OUT $line;

                INNER: while ( $line =<RG> ) {
                    print OUT $line;
                    last INNER if ( not ( $line =~ /\S/ ) );
                }
            }
        }
    }
    close RG;
    print "in total $readGroupCount read groups read from $readGroupFile.\n";

    return $readGroupCount;
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

