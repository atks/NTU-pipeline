#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

zsplit

=head1 SYNOPSIS

 zsplit [options]

  -i     input file
  -o     output directory
  -l     line number
  -c     line count file
 
=head1 DESCRIPTION

Splits a gzip file by line.

=cut

my $help;
my $sampleID;
my $inputFile1;
my $inputFile2;
my $outputDir;
my $lineNo;
my $lineCountFile;

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help,
                's:s'=>\$sampleID,
                'o:s'=>\$outputDir,
                'l:i'=>\$lineNo,
                'c:s'=>\$lineCountFile
               )
  || (scalar(@ARGV) != 2)
  || !defined($outputDir)
  || !defined($lineNo)
  || !defined($lineCountFile))
{
    if ($help)
    {
        pod2usage(-verbose => 2);
    }
    else
    {
        pod2usage(1);
    }
}

$inputFile1 = $ARGV[0];
$inputFile2 = $ARGV[1];

printf("zsplit\n");
printf("\n");
printf("options: input file 1         %s\n", $inputFile1);
printf("         input file 2         %s\n", $inputFile2);
printf("         output directory     %s\n", $outputDir);
printf("         line number          %s\n", $lineNo);
printf("         line count file      %s\n", $lineCountFile);
printf("\n");

my ($file, $dir, $suffix) = fileparse($inputFile1);
my $inputFile1LineCount = 0;
my $subFileNo = 1;
my $outputFile = "$outputDir/$subFileNo" . "_" . "$file";

#print "$file : $dir : $suffix\n";
#print "$outputFile\n\n";

my $lines_written = 0;
open(OUT, "| gzip -c > $outputFile");
print "creating $outputFile\n";
open(IN,"zcat $inputFile1 |") || die "Cannot open $inputFile1\n";
while (<IN>)
{
    if ($lines_written==$lineNo)
    {
        close(OUT);
        ++$subFileNo;
        $outputFile = "$outputDir/$subFileNo" . "_" . "$file";
        print "creating $outputFile\n";
        $lines_written = 0;
        open(OUT, "| gzip -c > $outputFile");
    }

    print OUT $_;
    ++$lines_written;       
    ++$inputFile1LineCount;
}
close(IN);
close(OUT);

($file, $dir, $suffix) = fileparse($inputFile2);
my $inputFile2LineCount = 0;
$subFileNo = 1;
$outputFile = "$outputDir/$subFileNo" . "_" . "$file";
print "creating $outputFile\n";

$lines_written = 0;
open(OUT, "| gzip -c > $outputFile");
open(IN,"zcat $inputFile2 |") || die "Cannot open $inputFile2\n";
while (<IN>)
{
    if ($lines_written==$lineNo)
    {
        close(OUT);
        ++$subFileNo;
        $outputFile = "$outputDir/$subFileNo" . "_" . "$file";
        print "creating $outputFile\n";
        $lines_written = 0;
        open(OUT, "| gzip -c > $outputFile");
    }

    print OUT $_;
    ++$lines_written;       
    ++$inputFile2LineCount;
}
close(IN);
close(OUT);

if ($inputFile1LineCount!=$inputFile2LineCount)
{
    die "paired FASTQ files do not have the same number of lines";
}

if ($inputFile1LineCount%4!=0)
{
    die "FASTQ file lines : $inputFile1LineCount not a multiple of 4";
}

open(OUT,">$lineCountFile") || die "Cannot open $lineCountFile\n";
print OUT "$sampleID\t$inputFile1\t$inputFile1LineCount\n";
print "writing $inputFile1LineCount to $lineCountFile\n";
close(OUT);

