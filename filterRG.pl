#! /usr/bin/perl
#
use strict;
use warnings;

use Data::Dumper;
use Time::HiRes;
use Memory::Usage;
#
##--------------------------------------------------
## filter read group with minimum support
#  			  remove PCR duplicates
#  			  and also overhang length
#
my $_debug = 0;

my %enviroment = (
);

my %global = (
    support => 2,
    overhang => 20,
    checkPCRduplicate => 1
);

#
##--------------------------------------------------
#
my $readGroup = shift;
my $supportSAM = shift;
my $output = shift;

&main ( $readGroup, $output, $supportSAM );
sub main
{
    my $readGroupFile = shift;
    my $outputFile = shift;
    my $supportSamFile = shift;
    my %parameters = @_;

    my %readGroup = ();
    my %supportReads = ();
    my $memoryUsage = Memory::Usage->new(); $memoryUsage->record('starting work');
    my $validGroup = loadReadGroup ( $readGroupFile, \%readGroup, \%supportReads );
    if ( defined $supportSamFile ) { 
	my $validReads = loadSupportSam ( $supportSamFile, \%supportReads, overhang => 20, checkPCRduplicates => 1 ); 

	if ( defined $global{checkPCRduplicate} ) {
	    $validGroup = removePCRduplicate ( \%readGroup, \%supportReads, minSupport => $global{support} );
	}
    } 

    printRG ( \%readGroup, $outputFile, minSupport => $global{support} );

    my $filteredSamFile = $supportSamFile; $filteredSamFile =~ s/.sam$/.filtered.sam/;
    printSupportSam ( \%supportReads, $filteredSamFile );

    $memoryUsage->record('final memory usage'); $memoryUsage->dump();

    1;
}

## ----------------------------------
sub loadReadGroup
{
    my $readGroupFile = shift;
    my $ref_readGroup = shift;
    my $ref_supportReads = shift;
    my %parameters = @_;

    my $totalReadGroupCount = 0; my $validReadGroupCount = 0;
    open ( RG, $readGroupFile ) or die "Cannot open $readGroupFile for reading!\n";
    OUTER: while ( my $line = <RG> ) {
        next if ( $line =~ /^#/ );
        if ( $line =~ /^Group/ ) {
            $totalReadGroupCount++;
            my ( $group, $chr1, $strand1, $start1, $end1, $chr2, $strand2, $start2, $end2, $support ) = ( $line =~ /^Group (\d+) == position (.+)\(([+-])\):(\d+)-(\d+)\|(.+)\(([+-])\):(\d+)-(\d+), support (\d+)\./ );
            if ( $support < $global{support} ) {
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
		$validReadGroupCount++;
		$ref_readGroup->{$group}{support} = $support;
		$ref_readGroup->{$group}{info} = $line;

                $line = <RG>;
                INNER: while ( $line =<RG> ) {
		    my ( $readID, @data ) = split ( /\s+/, $line );
		    if ( $readID ) {
			chomp $line;
			$ref_readGroup->{$group}{$readID} = $line;
			$ref_supportReads->{$readID}{group} = $group; 
		    }
                    last INNER if ( not ( $line =~ /\S/ ) );
                }
            }
        }
    }
    close RG;
    print "in total $totalReadGroupCount read groups read from $readGroupFile.\n";

    return $validReadGroupCount;
}

sub loadSupportSam
{
    my $supportSamFile = shift;
    my $ref_supportReads = shift;
    my %parameters = @_;

    my $totalReads = 0; my $validReads = 0;
    my %pos_reads = ();
    open ( SP, $supportSamFile );
    while ( my $line = <SP> ) {
	my @data = split ( /\t/, $line );
	if ( defined $ref_supportReads->{$data[0]} ) {
	    $validReads++;
	    $ref_supportReads->{$data[0]}{sam} = $line;
	    my $posMatch = $data[2] . ":" . $data[3] . ":" . $data[5];

	    if ( defined $parameters{checkPCRduplicates} ) {
		if ( defined $pos_reads{$posMatch} ) { 
		    $ref_supportReads->{$pos_reads{$posMatch}}{duplicate} = $posMatch; 
		    $ref_supportReads->{$data[0]}{duplicate} = $posMatch; 
		}
		else { $pos_reads{$posMatch} = $data[0]; }
	    }
	}
    }
    close SP;

    print STDERR "$validReads valid reads in total!\n";
    return $validReads;
}

sub removePCRduplicate
{
    my $ref_readGroup = shift;
    my $ref_supportReads = shift;
    my %parameters = @_;

    my $validGroupCount = 0;
    foreach my $group ( keys %{$ref_readGroup} ) {
	my %dup = ();
	foreach my $readID ( keys %{$ref_readGroup->{$group}} ) {
	    next if ( ( $readID eq "support" ) or ( $readID eq "info" ) );
	    if ( defined $ref_supportReads->{$readID}{duplicate} ) { 
		$ref_readGroup->{$group}{$readID} .= "\t*";
		$dup{$ref_supportReads->{$readID}{duplicate}}++; 
	    }
	}

	foreach my $dupPos ( keys %dup ) { $ref_readGroup->{$group}{support} = $ref_readGroup->{$group}{support} - $dup{$dupPos} + 1; }
	my $support = $ref_readGroup->{$group}{support};
	$ref_readGroup->{$group}{info} =~ s/support /support $support\//;
	if ( $ref_readGroup->{$group}{support} < $parameters{minSupport} ) {
	    foreach my $readID ( keys %{$ref_readGroup->{$group}} ) {
	    	next if ( ( $readID eq "support" ) or ( $readID eq "info" ) );
		$ref_supportReads->{$readID}{group} = -1;
	    }
	}
	else { $validGroupCount++; }
    }

    print "$validGroupCount valid groups remained after removing potential PCR duplicates!\n";
    return $validGroupCount;
}

sub printRG 
{
    my $ref_readGroup = shift;
    my $outputFile = shift;
    my %parameters = @_;

    open ( OUT, ">$outputFile" ) or die "Cannot open $outputFile for writing!\n";
    foreach my $group ( keys %{$ref_readGroup} ) {
	next if ( ( defined $parameters{minSupport} ) and ( $ref_readGroup->{$group}{support} < $parameters{minSupport} ) );
	print OUT $ref_readGroup->{$group}{info};
	print OUT "----\n";
	foreach my $readID ( keys %{$ref_readGroup->{$group}} ) {
	    next if ( ( $readID eq "info" ) or ( $readID eq "support" ) );
	    print OUT $ref_readGroup->{$group}{$readID}, "\n";
	}
	print OUT "\n";
    }
    close OUT;

    1;
}

sub printSupportSam
{
    my $ref_supportReads = shift;
    my $outputFile = shift;

    open ( OUT, ">$outputFile" ) or die "Cannot open $outputFile for writing!\n";
    foreach my $readID ( keys %{$ref_supportReads} ) {
	if ( $ref_supportReads->{$readID}{group} >=0 ) {
	    if ( not defined $ref_supportReads->{$readID}{sam} ) { print STDERR "Warning! Skipping unrecognized read $readID\n"; }
	    print OUT $ref_supportReads->{$readID}{sam};
	}
    }
    close OUT;

    1;
}
