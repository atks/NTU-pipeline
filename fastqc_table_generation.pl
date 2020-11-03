#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

rearrange_chromosomes.pl

=head1 SYNOPSIS

 rearrange_chromosomes.pl [options]

  -r     reference genome file
  -l     sequence length file
  -o     output file
  -m     make file name

=head1 DESCRIPTION

Rearranges chromosomes in a FASTA file.

=cut

my $help;

my $sequenceLengthFile;
my $refGenomeFASTAFile;
my $outputDir;

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help,
                'o=s'=>\$outputDir,
                'l=s'=>\$sequenceLengthFile,
                'r:s'=>\$refGenomeFASTAFile
               )
  || !defined($outputDir)
  || !defined($sequenceLengthFile)
  || !defined($refGenomeFASTAFile))
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

printf("rearrange_chromosomes.pl\n");
printf("\n");
printf("options: output dir           %s\n", $outputDir);
printf("         make file            %s\n", $refGenomeFASTAFile);
printf("         split line no        %s\n", $sequenceLengthFile);
printf("\n");

##########################
#Read sequence length file 
##########################
my %CHROM = ();
my @orderedChrom = ();
open(CHROM,"$sequenceLengthFile") || die "Cannot open $sequenceLengthFile\n";
while (<CHROM>)
{
    s/\r?\n?$//;
    if(!/^#/)
    {
        my ($chrom, $length) = split("\t");
        print "SEQ: $chrom\n";
        push(@orderedChrom, $chrom);
        $CHROM{$chrom} = $length;
    }
}
close(CHROM);

print  "Expecting " . scalar(@orderedChrom) . " sequences\n";

mkpath($outputDir);

open(CHROM,"$refGenomeFASTAFile") || die "Cannot open $refGenomeFASTAFile\n";
while (<CHROM>)
{
    s/\r?\n?$//;
    if(/^>/)
    {
        s/>//;
        
        my $chrom = $_;
        if ($.!=1)
        {   
            close(NEW);
        }
        
        open(NEW,">$outputDir/$chrom.fa") || die "Cannot open $outputDir/$chrom.fa\n";
        print "opened $chrom.fa\n";
        print NEW ">$_\n";
    }
    else
    {
        print NEW uc($_) . "\n";
    }
}
close(NEW);
close(CHROM);

my $noFiles = 0;
open(ORDERED_FA,">$outputDir/ordered.fa") || die "Cannot open $outputDir/ordered.fa\n";
for my $chrom (@orderedChrom)
{
    print "reading $chrom.fa .. ";
    
    open(CHROM_FA,"$outputDir/$chrom.fa") || die "Cannot open $outputDir/$chrom.fa\n";
    $/ = undef;
    my $fa = <CHROM_FA>;
    close(CHROM_FA);
    print ORDERED_FA $fa;  
    print "written\n";
    ++$noFiles;
}
close(ORDERED_FA);

print  "Written in order :  " . scalar($noFiles) . " to $outputDir/ordered.fa\n";
