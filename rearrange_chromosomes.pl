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

Reorders chromosomes in a FASTA file.

=cut

my $help;

my $chromosomeList;
my $orderedRefGenomeFASTAFile;
my $refGenomeFASTAFile;


#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help,
                'o=s'=>\$outputDir,
                's=s'=>\$sampleFile,
                'd=i'=>\$splitLineNo,
                'n:s'=>\$pipelineName,
                'm:s'=>\$makeFile,
                'r:s'=>\$refGenomeDir
               )
  || !defined($outputDir)
  || !defined($splitLineNo)
  || !defined($sampleFile)
  || !defined($makeFile)
  || !defined($pipelineName)
  || !defined($refGenomeDir))
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

$makeFile = "$outputDir/$makeFile";

#programs
my $bismark_genome_preparation = "/home/users/ntu/adrianta/programs/bismark-0.19.0/bismark_genome_preparation";
my $bowtie2 = "/app/bowtie2/2.29/bowtie2";
my $bowtie2Path = "/app/bowtie2/2.29";
my $samtools = "/app/samtools/1.3/bin/samtools";
my $fastqc = "/home/users/ntu/adrianta/programs/fastqc-0.11.5/fastqc";
my $trimGalore = "/home/users/ntu/adrianta/programs/trimGalore-0.4.5/trim_galore";
my $cutAdaptPath = "/home/users/ntu/adrianta/programs/cutadapt-1.15";
my $cutAdapt = "/home/users/ntu/adrianta/programs/cutadapt-1.15/cutadapt";
my $zsplit = "/home/users/ntu/adrianta/programs/NTU-pipeline/zsplit.pl";
my $bismark = "/home/users/ntu/adrianta/programs/bismark-0.19.0/bismark";
my $bismarkPath = "/home/users/ntu/adrianta/programs/bismark-0.19.0";

printf("generate_fastqc_pipeline_makefile.pl\n");
printf("\n");
printf("options: output dir           %s\n", $outputDir);
printf("         make file            %s\n", $makeFile);
printf("         split line no        %s\n", $splitLineNo);
printf("         sample file          %s\n", $sampleFile);
printf("         reference            %s\n", $refGenomeDir);
printf("\n");


#################
#Bisulfite Genome
#################
#not handled by this pipeline, partly because I wanted a specific path 
#for it to be in which cannot be specified on the command line.  
#/home/users/ntu/adrianta/ref/hg19/Bisulfite_Genome
#$BISMARK_PATH/bismark_genome_preparation --path_to_bowtie $BOWTIE_PATH --verbose $GENOME_PATH

#################
#Read sample file 
#################
my %SAMPLE = ();
my @SAMPLE = ();
open(SA,"$sampleFile") || die "Cannot open $sampleFile\n";
while (<SA>)
{
    s/\r?\n?$//;
    if(!/^#/)
    {
        my ($sampleID, $fastq1Path, $fastq2Path) = split(/\s+/, $_);
        
        if (exists($SAMPLE{$sampleID}))
        {
            exit("$sampleID already exists. Please fix.");
        }
        
        $SAMPLE{$sampleID}{FASTQ1} = $fastq1Path;
        $SAMPLE{$sampleID}{FASTQ2} = $fastq2Path;
        
        push(@SAMPLE, $sampleID);
    }
}
close(SA);
