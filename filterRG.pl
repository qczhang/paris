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
my $supportSAM = shift;

&main ( $readGroup, $output, $supportSAM );
sub main
{
    my $readGroupFile = shift;
    my $outputFile = shift;
    my $supportSamFile = shift;
    my %parameters = @_;

    my %supportReads = ();
    my $memoryUsage = Memory::Usage->new(); $memoryUsage->record('starting work');
    my $totalGroup = loadReadGroup ( $readGroupFile, $output, \%supportReads );
    if ( defined $supportSamFile ) { 
	my $validReads = loadSupportSam ( $supportSamFile, \%supportReads ); 
	print STDERR "$validReads valid reads in total!\n";
    } 

    $memoryUsage->record('final memory usage'); $memoryUsage->dump();

    1;
}

## ----------------------------------
sub loadReadGroup
{
    my $readGroupFile = shift;
    my $outputFile = shift;
    my $ref_supportReads = shift;

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

		    my @support = split ( /\s+/, $line );
		    if ( $support[0] ) { $ref_supportReads->{$support[0]} = 1; }
                    last INNER if ( not ( $line =~ /\S/ ) );
                }
            }
        }
    }
    close RG;
    print "in total $readGroupCount read groups read from $readGroupFile.\n";

    return $readGroupCount;
}

sub loadSupportSam
{
    my $supportSamFile = shift;
    my $ref_supportReads = shift;

    my $validReads = 0;
    my $filterSamFile = $supportSamFile;
    $filterSamFile =~ s/.sam$/.filtered.sam/;
    open ( SP, $supportSamFile );
    open ( FL, ">$filterSamFile" );
    while ( my $line = <SP> ) {
	my @data = split ( /\t/, $line );
	if ( defined $ref_supportReads->{$data[0]} ) {
	    $validReads++;
	    print FL $line;
	}
    }
    close SP;
    close FL;

    return $validReads;
}

