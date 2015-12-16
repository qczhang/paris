#! /usr/bin/perl
#
use strict;
use warnings;

use Data::Dumper;
use Getopt::Std;

#
use lib "module";
use PARISutil qw( &readGTF_ensembl_new &getExonID &getBioType &get5primeLen &get3primeLen &getCDSLen );

##--------------------------------------------------
#
my $_debug = 0;

my %environment = (
    annotationFile => "/Share/home/zhangqf/cliff/paris/data/ref/gencode.v21.chr_patch_hapl_scaff.annotation.gtf",
    exonHeatMapSize => 6,
    heatmapRscript    => "scripts/heatmap.2.R",
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
    my $ref_annotation = readGTF_ensembl_new ( $parameters{annotationFile}, verbose => 1 );

    my %readGroup = ();
    my $readGroupFile = $parameters{readGroupFile};
    my $outputFile = $parameters{output};
    my $totalGroup = loadReadGroup ( $readGroupFile, $outputFile, $ref_annotation, filter => $parameters{filter} );

    plotDistribution ( $outputFile, filter => $parameters{filter}, countCutoff => 100, suppressSingleExon => 1, exonHeatMapSize => $environment{exonHeatMapSize} );

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
    my $ref_annotation = shift;
    my %parameters = @_;

    my $exonHeatMapSize = sprintf("%.0f", $parameters{exonHeatMapSize} );
    my $halfExonHeatMapSize = sprintf("%.0f", $parameters{exonHeatMapSize}/2 );

    my %class_gDist = (); my %class_tDist = (); my %class_id = ();

    my %exon_label = (); my %pc_label = ();
    my @heatMapExon = (); for ( my $idx1 = 0; $idx1 < $exonHeatMapSize; $idx1++ ) { for ( my $idx2 = 0; $idx2 < $exonHeatMapSize; $idx2++ ) { $heatMapExon[$idx1][$idx2] = 0; } }
    my @avgHeatMapExon = (); for ( my $idx1 = 0; $idx1 < $exonHeatMapSize; $idx1++ ) { for ( my $idx2 = 0; $idx2 < $exonHeatMapSize; $idx2++ ) { $avgHeatMapExon[$idx1][$idx2] = 1; } }
    my @heatMapProtein = (); for ( my $idx1 = 0; $idx1 < 5; $idx1++ ) { for ( my $idx2 = 0; $idx2 < 5; $idx2++ ) { $heatMapProtein[$idx1][$idx2] = 0; } }
    my @avgHeatMapProtein = (); for ( my $idx1 = 0; $idx1 < 5; $idx1++ ) { for ( my $idx2 = 0; $idx2 < 5; $idx2++ ) { $avgHeatMapProtein[$idx1][$idx2] = 1; } }
    my %proteinCode = ( "5UTR" => 0, "start_codon" => 1, "CDS" => 2, "stop_codon" => 3, "3UTR" => 4 );

    open ( DATA, $dataFile );
    my $lineCount = 0;
    while ( my $line = <DATA> ) {
        next if ( $line =~ /^#/ );
        chomp $line;
        my @data = split ( /\t/, $line );
        $lineCount++;
        next if ( ( not defined $data[27] ) or ( not $data[27] ) );
        next if ( ( $data[11] eq "null" ) or ( $data[11] ne $data[19] ) );
        push ( @{$class_gDist{$data[22]}}, abs($data[10]) );
        push ( @{$class_tDist{$data[22]}}, abs($data[27]) );
        push ( @{$class_id{$data[22]}}, $data[11] );

        if ( $data[18] >= $exonHeatMapSize )  {
            if ( not defined $exon_label{$data[11]} ) {
                $exon_label{$data[11]} = 1;

                for ( my $idx1 = 0; $idx1 < $exonHeatMapSize; $idx1++ ) {
                    if ( $idx1 >= $halfExonHeatMapSize ) { $idx1 = $data[18] + $idx1 - $exonHeatMapSize; }
                    for ( my $idx2 = 0; $idx2 < $exonHeatMapSize; $idx2++ ) {
                        if ( $idx2 > $halfExonHeatMapSize ) { $idx2 = $data[18] + $idx2 - $exonHeatMapSize; }
                        if ( $idx1 == $idx2 ) { $avgHeatMapExon[$idx1][$idx1] += getExonLen ( $data[11], $idx1 + 1 ); }
                        else {
                            $avgHeatMapExon[$idx1][$idx2] += getExonLen ( $data[11], $idx1 ) + getExonLen ( $data[11], $idx2 + 1 );
                            $avgHeatMapExon[$idx2][$idx1] += getExonLen ( $data[11], $idx1 ) + getExonLen ( $data[11], $idx2 + 1 );
                        }
                    }
                }
            }
            else {
                $exon_label{$data[11]}++;
            }

            if ( $data[16] <= $halfExonHeatMapSize )  {
                if ( $data[24] <= $halfExonHeatMapSize ) {
                    $heatMapExon[$data[16]-1][$data[24]-1]++;
                    if ( $data[24] != $data[16] ) {
                        $heatMapExon[$data[24]-1][$data[16]-1]++;
                    }
                }
                my $reverse2 = $data[26] - $data[24];
                if ( $reverse2 < $halfExonHeatMapSize ) {
                    $heatMapExon[$data[16]-1][$exonHeatMapSize-1-$reverse2]++;
                    $heatMapExon[$exonHeatMapSize-1-$reverse2][$data[16]-1]++;
                }
            }
            my $reverse1 = $data[18] - $data[16];
            if ( $reverse1 < $halfExonHeatMapSize ) {
                if ( $data[24] <= $halfExonHeatMapSize ) {
                    $heatMapExon[$exonHeatMapSize-1-$reverse1][$data[24]-1]++;
                    $heatMapExon[$data[24]-1][$exonHeatMapSize-1-$reverse1]++;
                }
                my $reverse2 = $data[26] - $data[24];
                if ( $reverse2 < $halfExonHeatMapSize ) {
                    $heatMapExon[$exonHeatMapSize-1-$reverse1][$exonHeatMapSize-1-$reverse2]++;
                    if ( $reverse1 != $reverse2 ) {
                        $heatMapExon[$exonHeatMapSize-1-$reverse2][$exonHeatMapSize-1-$reverse1]++;
                    }
                }
            }
        }

        if ( ( defined $proteinCode{$data[14]} ) and ( defined $proteinCode{$data[22]} ) ) {
            if ( not defined $pc_label{$data[11]} ) {
                $pc_label{$data[11]} = 1;

                my @pcLens = ( get5primeLen ( $ref_annotation, $data[11] ), 3, getCDSLen ( $ref_annotation, $data[11] ), 3, get3primeLen ( $ref_annotation, $data[11] ) );
                for ( my $idx1 = 0; $idx1 < 5; $idx1++ ) {
                    for ( my $idx2 = 0; $idx2 < 5; $idx2++ ) {
                        if ( $idx1 == $idx2 )   {  $avgHeatMapProtein[$idx1][$idx2] += $pcLens[$idx1];  }
                        else {
                            $avgHeatMapProtein[$idx1][$idx2] += $pcLens[$idx1] + $pcLens[$idx2];
                            $avgHeatMapProtein[$idx2][$idx1] += $pcLens[$idx1] + $pcLens[$idx2];
                        }
                    }
                }
            }
            else {
                $pc_label{$data[11]}++;
            }

            $heatMapProtein[$proteinCode{$data[14]}][$proteinCode{$data[22]}]++;
            if ( $data[14] ne $data[22] ) { $heatMapProtein[$proteinCode{$data[22]}][$proteinCode{$data[14]}]++; }
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

    for ( my $idx1 = 0; $idx1 < $exonHeatMapSize; $idx1++ ) { for ( my $idx2 = 0; $idx2 < $exonHeatMapSize; $idx2++ ) { $avgHeatMapExon[$idx1][$idx2] = 1000 * $heatMapExon[$idx1][$idx2] / $avgHeatMapExon[$idx1][$idx2]; } }
    my $avgHeatMapExonData = $dataFile;  $avgHeatMapExonData =~ s/.txt//; $avgHeatMapExonData .= ".exon.norm.heatmap";
    open ( EH, ">$avgHeatMapExonData" );
    for ( my $idx1 = 0; $idx1 < $exonHeatMapSize; $idx1++ ) { print EH join ( "\t", @{$avgHeatMapExon[$idx1]} ), "\n"; }
    close EH;
    for ( my $idx1 = 0; $idx1 < 5; $idx1++ ) { for ( my $idx2 = 0; $idx2 < 5; $idx2++ ) { $avgHeatMapProtein[$idx1][$idx2] = 1000 * $heatMapProtein[$idx1][$idx2] / $avgHeatMapProtein[$idx1][$idx2]; } }
    $avgHeatMapProtein[1][3] = 0;
    $avgHeatMapProtein[3][1] = 0;
    my $avgHeatMapProteinData = $dataFile;  $avgHeatMapProteinData =~ s/.txt//; $avgHeatMapProteinData .= ".pc.norm.heatmap";
    open ( PC, ">$avgHeatMapProteinData" );
    for ( my $idx1 = 0; $idx1 < 5; $idx1++ ) { print PC join ( "\t", @{$avgHeatMapProtein[$idx1]} ), "\n"; }
    close PC;
    my $heatMapExonPDF = $dataFile;  $heatMapExonPDF =~ s/.txt//; $heatMapExonPDF .= ".exon.heatmap.pdf";
    my $heatMapPCPDF = $dataFile;  $heatMapPCPDF =~ s/.txt//; $heatMapPCPDF .= ".pc.heatmap.pdf";
    my $avgHeatMapExonPDF = $dataFile;  $avgHeatMapExonPDF =~ s/.txt//; $avgHeatMapExonPDF .= ".exon.norm.heatmap.pdf";
    my $avgHeatMapPCPDF = $dataFile;  $avgHeatMapPCPDF =~ s/.txt//; $avgHeatMapPCPDF .= ".pc.norm.heatmap.pdf";
    print STDERR `Rscript $environment{heatmapRscript} $heatMapExonData $heatMapExonPDF $avgHeatMapExonData $avgHeatMapExonPDF $heatMapProteinData $heatMapPCPDF $avgHeatMapProteinData $avgHeatMapPCPDF`;

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
                my $featureLen1 = ".";  my $featureLen2 = ".";  
                my @data = split ( /\t/, $line );  my $transSpan = ".";
                if ( $data[0] eq $data[6] ) { $transSpan = ( $data[7] > $data[1] ) ? ( $data[8] - $data[1] ) : ( $data[2] - $data[7] ); }
                my $bioType1 = "null";  my $bioType2 = "null";
                my $ref_exonID1 = [];  my $ref_exonID2 = [];
                my $ref_exonLen1 = [];  my $ref_exonLen2 = [];
                $data[1] = "." if ( $data[1] eq "-" ); $data[2] = "." if ( $data[2] eq "-" ); $data[7] = "." if ( $data[7] eq "-" ); $data[8] = "." if ( $data[8] eq "-" );
                if ( $data[0] =~ /ENS/ ) { 
                    $bioType1 = parseFeature ( $ref_annotation, $data[0], $data[1], $data[2] ); 
                    if ( $bioType1 eq "5UTR" )  {  $featureLen1 = get5primeLen ( $ref_annotation, $data[0] );  }
                    elsif ( $bioType1 eq "3UTR" )  {  $featureLen1 = get3primeLen ( $ref_annotation, $data[0] );  }
                    elsif ( $bioType1 eq "CDS" )  {  $featureLen1 = getCDSLen ( $ref_annotation, $data[0] );  }
                    elsif ( $bioType1 eq "start_codon" )  {  $featureLen1 = 3;  }
                    elsif ( $bioType1 eq "stop_codon" )  {  $featureLen1 = 3;  }
                    ( $ref_exonID1, $ref_exonLen1 ) = getExonID ( $ref_annotation, $data[0], substr ( $data[3], 1 ), substr ( $data[4], 0, -1) ); 
                }
                if ( $data[6] =~ /ENS/ ) { 
                    $bioType2 = parseFeature ( $ref_annotation, $data[6], $data[7], $data[8] ); 
                    if ( $bioType2 eq "5UTR" )  {  $featureLen2 = get5primeLen ( $ref_annotation, $data[6] );  }
                    elsif ( $bioType2 eq "3UTR" )  {  $featureLen2 = get3primeLen ( $ref_annotation, $data[6] );  }
                    elsif ( $bioType2 eq "CDS" )  {  $featureLen2 = getCDSLen ( $ref_annotation, $data[6] );  }
                    elsif ( $bioType2 eq "start_codon" )  {  $featureLen2 = 3;  }
                    elsif ( $bioType2 eq "stop_codon" )  {  $featureLen2 = 3;  }
                    ( $ref_exonID2, $ref_exonLen2 ) = getExonID ( $ref_annotation, $data[6], substr ( $data[9], 1 ), substr ( $data[10], 0, -1) ); 
                }
                print OUT join ( "\t", $dgID, $support, $chr1, $strand1, $start1, $end1, $chr2, $strand2, $start2, $end2, $genomeSpan ), "\t"; 
                print OUT join ( "\t", $data[0], $data[1], $data[2], $bioType1, $featureLen1, join ( ",", @{$ref_exonID1} ), join ( ",", @{$ref_exonLen1} ), $ref_annotation->{transcript_info}{$data[0]}{exonNum} ), "\t"; 
                print OUT join ( "\t", $data[6], $data[7], $data[8], $bioType2, $featureLen2, join ( ",", @{$ref_exonID2} ), join ( ",", @{$ref_exonLen2} ), $ref_annotation->{transcript_info}{$data[6]}{exonNum}, $transSpan ), "\n"; 
            }

            INNER: while ( $line =<RG> ) {
                if ( $line =~ /ENS/ ) {
                    chomp $line;
                    my $featureLen1 = ".";  my $featureLen2 = ".";  
                    my @data = split ( /\t/, $line );  my $transSpan = ".";
                    if ( $data[0] eq $data[6] ) { $transSpan = ( $data[7] > $data[1] ) ? ( $data[8] - $data[1] ) : ( $data[2] - $data[7] ); }
                    $data[1] = "." if ( $data[1] eq "-" ); $data[2] = "." if ( $data[2] eq "-" ); $data[7] = "." if ( $data[7] eq "-" ); $data[8] = "." if ( $data[8] eq "-" );
                    my $bioType1 = "null";  my $bioType2 = "null";
                    my $ref_exonID1 = [];  my $ref_exonID2 = [];
                    my $ref_exonLen1 = [];  my $ref_exonLen2 = [];
                    if ( $data[0] =~ /ENS/ ) { 
                        $bioType1 = parseFeature ( $ref_annotation, $data[0], $data[1], $data[2] ); 
                        if ( $bioType1 eq "5UTR" )  {  $featureLen1 = get5primeLen ( $ref_annotation, $data[0] );  }
                        elsif ( $bioType1 eq "3UTR" )  {  $featureLen1 = get3primeLen ( $ref_annotation, $data[0] );  }
                        elsif ( $bioType1 eq "CDS" )  {  $featureLen1 = getCDSLen ( $ref_annotation, $data[0] );  }
                        ( $ref_exonID1, $ref_exonLen1 ) = getExonID ( $ref_annotation, $data[0], substr ( $data[3], 1 ), substr ( $data[4], 0, -1) ); 
                    }
                    if ( $data[6] =~ /ENS/ ) { 
                        $bioType2 = parseFeature ( $ref_annotation, $data[6], $data[7], $data[8] ); 
                        if ( $bioType1 eq "5UTR" )  {  $featureLen2 = get5primeLen ( $ref_annotation, $data[6] );  }
                        elsif ( $bioType1 eq "3UTR" )  {  $featureLen2 = get3primeLen ( $ref_annotation, $data[6] );  }
                        elsif ( $bioType1 eq "CDS" )  {  $featureLen2 = getCDSLen ( $ref_annotation, $data[6] );  }
                        ( $ref_exonID2, $ref_exonLen2 ) = getExonID ( $ref_annotation, $data[6], substr ( $data[9], 1 ), substr ( $data[10], 0, -1) ); 
                    }
                    print OUT join ( "\t", $dgID, $support, $chr1, $strand1, $start1, $end1, $chr2, $strand2, $start2, $end2, $genomeSpan ), "\t"; 
                    print OUT join ( "\t", $data[0], $data[1], $data[2], $bioType1, $featureLen1, join ( ",", @{$ref_exonID1} ), join ( ",", @{$ref_exonLen1} ), $ref_annotation->{transcript_info}{$data[0]}{exonNum} ), "\t"; 
                    print OUT join ( "\t", $data[6], $data[7], $data[8], $bioType2, $featureLen2, join ( ",", @{$ref_exonID2} ), join ( ",", @{$ref_exonLen2} ), $ref_annotation->{transcript_info}{$data[6]}{exonNum}, $transSpan ), "\n"; 
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
