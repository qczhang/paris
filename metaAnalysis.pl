#! /usr/bin/perl
#
use strict;
use warnings;

use Data::Dumper;
use Getopt::Std;

#
use lib "module";
use PARISutil qw( &readGTF_ensembl_new &getExonID &getBioType &get5primeLen &get3primeLen );

##--------------------------------------------------
#
my $_debug = 0;

my %environment = (
    annotationFile => "/Share/home/zhangqf/cliff/paris/data/ref/gencode.v21.chr_patch_hapl_scaff.annotation.gtf",
    exonHeatMapSize => 6,
    sizeRscript    => "scripts/sizeDist.violinplot.R",
    countRscript   => "scripts/countDist.pieChart.R",
    preDefinedRNAs => [ "5UTR", "start_codon", "CDS", "stop_codon", "3UTR", "protein coding", "lincRNA", "retained intron", "nonsense mediated decay", "processed pseudogene", "processed transcript", "snoRNA", "snRNA", "rRNA", "Mt rRNA", "misc RNA", "antisense" ]
);

use vars qw ($opt_h $opt_V $opt_D $opt_i $opt_o $opt_a $opt_f );
&getopts('hVDi:o:a:f:');

my $usage = <<_EOH_;
## --------------------------------------
Command:
count duplex group genomic features
$0 -i annotated_read_group_file -o tabulated_duplex_group_statistics_file
# what it is:
 -i     input read group file
 -o     output annotated file

# more options:
 -a     annotation GTF file
 -f     filter ("")
_EOH_
;

&main();

sub main {
    my %parameters = &init();
    print Dumper \%parameters if ( $opt_D );

#   my $ref_annotation = undef;
    #my $ref_annotation = readGTF_ensembl_new ( $parameters{annotationFile}, verbose => 1 );

    my %readGroup = ();
    my $readGroupFile = $parameters{readGroupFile};
    my $outputFile = $parameters{output};
    #my $totalGroup = loadReadGroup ( $readGroupFile, $outputFile, $ref_annotation, filter => $parameters{filter} );

    plotDistribution ( $outputFile, filter => $parameters{filter}, countCutoff => 100, suppressSingleExon => 1, exonHeatMapSize => 4 );

    1;
}

## ----------------------------------
sub init 
{
    my %parameters = ();

    die $usage if ( $opt_h || ( not $opt_i ) || ( not $opt_o ) );
    $opt_V = 0 if ( not defined $opt_V );
    $opt_D = 0 if ( not defined $opt_D );

    if ( defined $opt_i ) { $parameters{readGroupFile} = $opt_i; }
    if ( defined $opt_o ) { $parameters{output} = $opt_o; }

    if ( defined $opt_a ) { $parameters{annotationFile} = $opt_a; }
    else  {  $parameters{annotationFile} = $environment{annotationFile};  }
    if ( defined $opt_f ) { $parameters{filter} = $opt_f; }
    else  {  $parameters{filter} = "noFiltering";  }

    return ( %parameters );
}

sub plotDistribution
{ 
    my $dataFile = shift;
    my %parameters = @_;
    my $exonHeatMapSize = sprintf("%.0f", $parameters{exonHeatMapSize} );
    my $halfExonHeatMapSize = sprintf("%.0f", $parameters{exonHeatMapSize}/2 );

    my %class_gDist = ();
    my %class_tDist = ();
    my %class_id = ();

    my @heatMapExon = ();
    for ( my $idx1 = 0; $idx1 < $exonHeatMapSize; $idx1++ ) { for ( my $idx2 = 0; $idx2 < $exonHeatMapSize; $idx2++ ) { $heatMapExon[$idx1][$idx2] = 0; } }
    my @heatMapProtein = ();
    for ( my $idx1 = 0; $idx1 < 5; $idx1++ ) { for ( my $idx2 = 0; $idx2 < 5; $idx2++ ) { $heatMapProtein[$idx1][$idx2] = 0; } }
    my %proteinCode = ( "5UTR" => 0, "start_codon" => 1, "CDS" => 2, "stop_codon" => 3, "3UTR" => 4 );

    open ( DATA, $dataFile );
    my $lineCount = 0;
    while ( my $line = <DATA> ) {
        next if ( $line =~ /^#/ );
        chomp $line;
        my @data = split ( /\t/, $line );
	$lineCount++;
        next if ( ( not $data[20] ) or ( not $data[10] ) or ( not $data[23] ) );
        next if ( ( $data[11] eq "null" ) or ( $data[11] ne $data[17] ) );
        push ( @{$class_gDist{$data[20]}}, abs($data[10]) );
        push ( @{$class_tDist{$data[20]}}, abs($data[23]) );
        push ( @{$class_id{$data[20]}}, $data[11] );

	if ( ( $data[16] >= $exonHeatMapSize ) and ( $data[22] >= $exonHeatMapSize ) ) {
	    if ( $data[15] <= $halfExonHeatMapSize )  {
	    	if ( $data[21] <= $halfExonHeatMapSize ) {
	            $heatMapExon[$data[15]-1][$data[21]-1]++;
	            $heatMapExon[$data[21]-1][$data[15]-1]++ if ( $data[21] != $data[15] );
	    	}
	    my $reverse2 = $data[22] - $data[21];
	    if ( $reverse2 < $halfExonHeatMapSize ) {
		$heatMapExon[$data[15]-1][$exonHeatMapSize-1-$reverse2]++;
		$heatMapExon[$exonHeatMapSize-1-$reverse2][$data[15]-1]++;
	    }
	}
	my $reverse1 = $data[16] - $data[15];
	if ( $reverse1 < $halfExonHeatMapSize ) {
	    if ( $data[21] <= $halfExonHeatMapSize ) {
	        $heatMapExon[$exonHeatMapSize-1-$reverse1][$data[21]-1]++;
	        $heatMapExon[$data[21]-1][$exonHeatMapSize-1-$reverse1]++;
	    }
	    my $reverse2 = $data[22] - $data[21];
	    if ( $reverse2 < $halfExonHeatMapSize ) {
		$heatMapExon[$exonHeatMapSize-1-$reverse1][$exonHeatMapSize-1-$reverse2]++;
		$heatMapExon[$exonHeatMapSize-1-$reverse2][$exonHeatMapSize-1-$reverse1]++ if ( $reverse1 != $reverse2 );
	    }
	}
	}

	if ( ( defined $proteinCode{$data[14]} ) and ( defined $proteinCode{$data[20]} ) ) {
	    $heatMapProtein[$proteinCode{$data[14]}][$proteinCode{$data[20]}]++;
	    $heatMapProtein[$proteinCode{$data[20]}][$proteinCode{$data[14]}]++ if ( $data[14] ne $data[20] );
	}
    }
    close DATA;

    my $heatMapExonData = $dataFile;  $heatMapExonData =~ s/.txt//; $heatMapExonData .= ".exon.heatmap";
    open ( EH, ">$heatMapExonData" );
    for ( my $idx1 = 0; $idx1 < $exonHeatMapSize; $idx1++ ) { print EH join ( "\t", @{$heatMapExon[$idx1]} ), "\n"; }
    close EH;
    my $heatMapProteinData = $dataFile;  $heatMapProteinData =~ s/.txt//; $heatMapProteinData .= ".pc.heatmap";
    open ( PC, ">$heatMapProteinData" );
    for ( my $idx1 = 0; $idx1 < 5; $idx1++ ) { print PC join ( "\t", @{$heatMapProtein[$idx1]} ), "\n"; }
    close PC;

    my $distData = $dataFile;  $distData =~ s/.txt//; $distData .= ".dist.dat";
    open ( RD, ">$distData" );
    print RD "type\tgenomeSize\ttransSize\ttransID\n";
    foreach my $class ( keys %class_gDist ) {
        if ( scalar (@{$class_gDist{$class}}) < $parameters{countCutoff} ) {
            for ( my $idx = 0; $idx < scalar ( @{$class_gDist{$class}} ); $idx++ )  {
                print RD "others\t", $class_gDist{$class}[$idx], "\t", $class_tDist{$class}[$idx], "\t", $class_id{$class}[$idx], "\n";
            }
        }
        else {
            for ( my $idx = 0; $idx < scalar ( @{$class_gDist{$class}} ); $idx++ )  {
                print RD $class, "\t", $class_gDist{$class}[$idx], "\t", $class_tDist{$class}[$idx], "\t", $class_id{$class}[$idx], "\n";
            }
        }
    }
    close RD;
    my $sizeDist = $dataFile;  $sizeDist =~ s/.txt//; $sizeDist .= ".size.pdf";
    print STDERR `Rscript $environment{sizeRscript} $distData $sizeDist`;

    my $countData = $dataFile;  $countData =~ s/.txt//; $countData .= ".count.dat";
    open ( CT, ">$countData" );
    print CT "type\tcount\n";
    my %class_count = ();
    foreach my $class ( @{$environment{preDefinedRNAs}} ) {
        if ( $class eq "protein coding" ) { print CT "protein coding (not resolved)\t", scalar (@{$class_gDist{$class}} ), "\n"; }
        else { print CT $class, "\t", scalar (@{$class_gDist{$class}} ), "\n"; }
        $class_count{$class} = scalar (@{$class_gDist{$class}} );
    }
    my $otherCount = 0;
    foreach my $class ( keys %class_gDist ) {
        if ( not defined $class_count{$class} ) {
            $otherCount += scalar (@{$class_gDist{$class}} );
        }
    }
    print CT "others\t$otherCount\n";
    close CT;
    my $countDist = $dataFile;  $countDist =~ s/.txt//; $countDist .= ".count.pdf";
    print STDERR `Rscript $environment{countRscript} $countData $countDist`;

    1;
}

sub loadReadGroup
{
    my $readGroupFile = shift;
    my $outputFile = shift;
    my $ref_annotation = shift;
    my %parameters = @_;

    my $readGroupCount = 0;
    open ( RG, $readGroupFile ) or die "Cannot open $readGroupFile for reading!\n";
    open ( OUT, ">$outputFile" ) or die "Cannot open $outputFile for writing!\n";
    my $lastLine = ""; my $tmpLine = "";
    print "Read duplex group from $readGroupFile\n";
    OUTER: while ( my $line = <RG> ) {
        next if ( $line =~ /^#/ );
        if ( $lastLine ) {
            $tmpLine = $line;
            $line = $lastLine;
        }

        if ( $line =~ /^Group/ ) {
            $readGroupCount++;
            last OUTER if ( ( $readGroupCount > 5 ) and $opt_D );
            my ( $dgID, $chr1, $strand1, $start1, $end1, $chr2, $strand2, $start2, $end2, $support ) = ( $line =~ /^Group (\d+) == position (.+)\(([+-])\):(\d+)-(\d+)\|(.+)\(([+-])\):(\d+)-(\d+), support (\d+)/ );
            my $genomeSpan = ".";  
            if ( ( $chr1 eq $chr2 ) and ( $strand1 eq $strand2 ) )  {  $genomeSpan = ( $start2 > $start1 ) ? ( $end2 - $start1 ) : ( $end1 - $start2 );  }

            if ( not $lastLine )  { $line = <RG>; }
            else { $line = $tmpLine; }

            if ( $line =~ /ENS/ ) { 
		chomp $line;
                my @data = split ( /\t/, $line );  my $transSpan = ".";
                if ( $data[0] eq $data[6] ) { $transSpan = ( $data[7] > $data[1] ) ? ( $data[8] - $data[1] ) : ( $data[2] - $data[7] ); }
                my $bioType1 = "null";  my $bioType2 = "null";
                my $ref_exonID1 = [];  my $ref_exonID2 = [];
                $data[1] = "." if ( $data[1] eq "-" ); $data[2] = "." if ( $data[2] eq "-" ); $data[7] = "." if ( $data[7] eq "-" ); $data[8] = "." if ( $data[8] eq "-" );
                if ( $data[0] =~ /ENS/ ) { 
                    $bioType1 = parseFeature ( $ref_annotation, $data[0], $data[1], $data[2] ); 
                    $ref_exonID1 = getExonID ( $ref_annotation, $data[0], substr ( $data[3], 1 ), substr ( $data[4], 0, -1) ); 
                }
                if ( $data[6] =~ /ENS/ ) { 
                    $bioType2 = parseFeature ( $ref_annotation, $data[6], $data[7], $data[8] ); 
                    $ref_exonID2 = getExonID ( $ref_annotation, $data[6], substr ( $data[9], 1 ), substr ( $data[10], 0, -1) ); 
                }
                print OUT join ( "\t", $dgID, $support, $chr1, $strand1, $start1, $end1, $chr2, $strand2, $start2, $end2, $genomeSpan ), "\t"; 
                print OUT join ( "\t", $data[0], $data[1], $data[2], $bioType1, join ( ",", @{$ref_exonID1} ), $ref_annotation->{transcript_info}{$data[0]}{exonNum} ), "\t"; 
                print OUT join ( "\t", $data[6], $data[7], $data[8], $bioType2, join ( ",", @{$ref_exonID2} ), $ref_annotation->{transcript_info}{$data[6]}{exonNum}, $transSpan ), "\n"; 
            }

            INNER: while ( $line =<RG> ) {
                if ( $line =~ /ENS/ ) {
		    chomp $line;
                    my @data = split ( /\t/, $line );  my $transSpan = ".";
                    if ( $data[0] eq $data[6] ) { $transSpan = ( $data[7] > $data[1] ) ? ( $data[8] - $data[1] ) : ( $data[2] - $data[7] ); }
                    $data[1] = "." if ( $data[1] eq "-" ); $data[2] = "." if ( $data[2] eq "-" ); $data[7] = "." if ( $data[7] eq "-" ); $data[8] = "." if ( $data[8] eq "-" );
                    my $bioType1 = "null";  my $bioType2 = "null";
                    my $ref_exonID1 = [];  my $ref_exonID2 = [];
		    if ( $data[0] =~ /ENS/ ) { 
			$bioType1 = parseFeature ( $ref_annotation, $data[0], $data[1], $data[2] ); 
			$ref_exonID1 = getExonID ( $ref_annotation, $data[0], substr ( $data[3], 1 ), substr ( $data[4], 0, -1) ); 
		    }
		    if ( $data[6] =~ /ENS/ ) { 
			$bioType2 = parseFeature ( $ref_annotation, $data[6], $data[7], $data[8] ); 
			$ref_exonID2 = getExonID ( $ref_annotation, $data[6], substr ( $data[9], 1 ), substr ( $data[10], 0, -1) ); 
		    }
		    print OUT join ( "\t", $dgID, $support, $chr1, $strand1, $start1, $end1, $chr2, $strand2, $start2, $end2, $genomeSpan ), "\t"; 
		    print OUT join ( "\t", $data[0], $data[1], $data[2], $bioType1, join ( ",", @{$ref_exonID1} ), $ref_annotation->{transcript_info}{$data[0]}{exonNum} ), "\t"; 
		    print OUT join ( "\t", $data[6], $data[7], $data[8], $bioType2, join ( ",", @{$ref_exonID2} ), $ref_annotation->{transcript_info}{$data[6]}{exonNum}, $transSpan ), "\n"; 
                }
                elsif ( $line =~ /^Group/ )  {
                    $lastLine = $line;
                    last INNER;
                }
            }
        }
    }
    close RG;
    print "in total $readGroupCount read groups read from $readGroupFile.\n";

    return $readGroupCount;
}


## ------------------------------------

sub parseFeature
{
    my $ref_annotation = shift;
    my $transID = shift;
    my $start = shift;
    my $end = shift;

    my $bioType = getBioType ( $ref_annotation, $transID );
    if ( $bioType eq "protein coding" ) {
        my $fivePrimeLength = get5primeLen ( $ref_annotation, $transID );
        my $threePrimeLength = get3primeLen ( $ref_annotation, $transID );
        my $CDSlength = $ref_annotation->{transcript_info}{$transID}{length} - $fivePrimeLength - $threePrimeLength;

        if ( $fivePrimeLength and $threePrimeLength and $CDSlength ) {
            if ( $end <= $fivePrimeLength )  { $bioType = "5UTR"; }
            elsif ( ( $start <= ( $fivePrimeLength + 3 ) ) and ( $end > $fivePrimeLength ) )  { $bioType = "start_codon"; }
            elsif ( ( $start > ( $fivePrimeLength + 3 ) ) and ( $end <= ( $fivePrimeLength + $CDSlength ) ) )  { $bioType = "CDS"; }
            elsif ( ( $start <= ( $fivePrimeLength + $CDSlength + 3 ) ) and ( $end > ( $fivePrimeLength + $CDSlength ) ) )  { $bioType = "stop_codon"; }
            elsif ( $start > ( $fivePrimeLength + $CDSlength + 3 ) )  { $bioType = "3UTR"; }
        }
    }

    return $bioType;
}
