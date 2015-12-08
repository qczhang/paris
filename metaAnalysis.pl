#! /usr/bin/perl
#
use strict;
use warnings;

use Data::Dumper;
use Getopt::Std;

#
use lib "module";
use PARISutil qw( &readGTF_ensembl_new &getBioType &get5primeLen &get3primeLen &getExonLen );

##--------------------------------------------------
#
my $_debug = 0;

my %environment = (
    annotationFile => "/Share/home/zhangqf/cliff/paris/data/ref/gencode.v21.chr_patch_hapl_scaff.annotation.gtf",
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

    plotDistribution ( $outputFile, filter => $parameters{filter}, countCutoff => 100 );

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

    my $distData = $dataFile;  $distData =~ s/.txt//; $distData .= ".dist.dat";
    my %class_gDist = ();
    my %class_tDist = ();
    my %class_id = ();
    open ( DATA, $dataFile );
    my $lineCount = 0;
    while ( my $line = <DATA> ) {
        next if ( $line =~ /^#/ );
        chomp $line;
        my @data = split ( /\t/, $line );
	$lineCount++;
        next if ( ( not $data[18] ) or ( not $data[10] ) or ( not $data[19] ) );
        next if ( ( $data[11] eq "null" ) or ( $data[11] ne $data[15] ) );
        push ( @{$class_gDist{$data[18]}}, abs($data[10]) );
        push ( @{$class_tDist{$data[18]}}, abs($data[19]) );
        push ( @{$class_id{$data[18]}}, $data[11] );
    }
    close DATA;

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
                my @data = split ( /\t/, $line );  my $transSpan = ".";
                if ( $data[0] eq $data[6] ) { $transSpan = ( $data[7] > $data[1] ) ? ( $data[8] - $data[1] ) : ( $data[2] - $data[7] ); }
                my $bioType1 = "null";  my $bioType2 = "null";
                $data[1] = "." if ( $data[1] eq "-" ); $data[2] = "." if ( $data[2] eq "-" ); $data[7] = "." if ( $data[7] eq "-" ); $data[8] = "." if ( $data[8] eq "-" );
                if ( $data[0] =~ /ENS/ ) { $bioType1 = parseFeature ( $ref_annotation, $data[0], $data[1], $data[2] ); }
                if ( $data[6] =~ /ENS/ ) { $bioType2 = parseFeature ( $ref_annotation, $data[6], $data[7], $data[8] ); }
                print OUT join ( "\t", $dgID, $support, $chr1, $strand1, $start1, $end1, $chr2, $strand2, $start2, $end2, $genomeSpan ), "\t"; 
                print OUT join ( "\t", $data[0], $data[1], $data[2], $bioType1, $data[6], $data[7], $data[8], $bioType2, $transSpan ), "\n"; 
            }

            INNER: while ( $line =<RG> ) {
                if ( $line =~ /ENS/ ) {
                    my @data = split ( /\t/, $line );  my $transSpan = ".";
                    if ( $data[0] eq $data[6] ) { $transSpan = ( $data[7] > $data[1] ) ? ( $data[8] - $data[1] ) : ( $data[2] - $data[7] ); }
                    $data[1] = "." if ( $data[1] eq "-" ); $data[2] = "." if ( $data[2] eq "-" ); $data[7] = "." if ( $data[7] eq "-" ); $data[8] = "." if ( $data[8] eq "-" );
                    my $bioType1 = "null";  my $bioType2 = "null";
                    if ( $data[0] =~ /ENS/ ) { $bioType1 = parseFeature ( $ref_annotation, $data[0], $data[1], $data[2] ); }
                    if ( $data[6] =~ /ENS/ ) { $bioType2 = parseFeature ( $ref_annotation, $data[6], $data[7], $data[8] ); }
                    print OUT join ( "\t", $dgID, $support, $chr1, $strand1, $start1, $end1, $chr2, $strand2, $start2, $end2, $genomeSpan ), "\t"; 
                    print OUT join ( "\t", $data[0], $data[1], $data[2], $bioType1, $data[6], $data[7], $data[8], $bioType2, $transSpan ), "\n"; 
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
