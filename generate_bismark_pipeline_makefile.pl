rtgrth
#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

generate_bismark_pipeline_makefile

=head1 SYNOPSIS

 generate_bismark_pipeline_makefile [options]

  -s     sample file list giving the location of each sample
         column 1: sample name
         column 2: path of bam file
  -r     reference genome file
  -l     sequence length file
  -w     interval width
  -o     output directory
  -m     make file name

=head1 DESCRIPTION

This script generates the make file to discovery and genotype a set of individuals.

=cut

my $help;

my $sampleFile = "";
 
my $outputDir;
my $slurmScriptsSubDir = "";
my $intervalWidth = 1000000;
my $refGenomeFASTAFile;

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help,
                'o:s'=>\$outputDir,
                'm:s'=>\$makeFile,
                'd:s'=>\$slurmScriptsSubDir,
                's:s'=>\$sampleFile,
                'i:s'=>\$intervalWidth,
                'r:s'=>\$refGenomeFASTAFile
               )
  || !defined($makeFile)
  || !defined($sampleFile)
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

#programs
my $bismark_genome_preparation = "/home/users/ntu/adrianta/programs/bismark-0.19.0/bismark_genome_preparation";

printf("generate_bismarck_pipeline_makefile.pl\n");
printf("\n");
printf("options: output dir           %s\n", $outputDir);
printf("         make file            %s\n", $makeFile);
printf("         sample file          %s\n", $sampleFile);
printf("         interval width       %s\n", $intervalWidth);
printf("         reference            %s\n", $refGenomeFASTAFile);
printf("\n");

my $vcfOutDir = "$outputDir/vcf";
#mkpath($vcfOutDir);
my $finalVCFOutDir = "$outputDir/final";
mkpath($finalVCFOutDir);
my $statsDir = "$outputDir/stats";
mkpath($statsDir);
my $logDir = "$outputDir/log";
mkpath($logDir);
my $auxDir = "$outputDir/aux";
mkpath($auxDir);
my $slurmScriptsDir = "$outputDir/slurm_scripts/$slurmScriptsSubDir";
mkpath($slurmScriptsDir);
my $slurmScriptNo = 0;
my $logFile = "$outputDir/run.log";

########################################
#Read file locations and name of samples
########################################
my %SAMPLE = ();
open(SA,"$sampleFile") || die "Cannot open $sampleFile\n";
my $bamFiles = "";
while (<SA>)
{
    s/\r?\n?$//;
    if(!/^#/)
    {
        my ($sampleID, $bamPath) = split(/\s+/, $_);
        $SAMPLE{$sampleID} = $bamPath;
        $bamFiles .= "$bamPath\n";
    }
}
close(SA);


exit;

#$BISMARK_PATH/bismark_genome_preparation --path_to_bowtie $BOWTIE_PATH --verbose $GENOME_PATH
#/home/users/ntu/adrianta/programs/bismark-0.19.0/bismark_genome_preparation --path_to_bowtie  /app/bowtie2/2.29 --verbose ~/ref/

my $bamListFile = "$auxDir/bam.list";
open(OUT,">$bamListFile") || die "Cannot open $bamListFile\n";
print OUT $bamFiles;
close(OUT);

print "read in " . scalar(keys(%SAMPLE)) . " samples\n";

###################
#Generate intervals
###################
my %intervalsByChrom = ();
my @intervals = ();
my @intervalNames = ();
my @intervalFiles = ();
my @CHROM = ();

open(SQ,"$refGenomeFASTAFile.fai") || die "Cannot open $refGenomeFASTAFile.fai\n";
while (<SQ>)
{
    s/\r?\n?$//;
    if(!/^#/)
    {
        my ($chrom, $len) = split('\t', $_);

        last if ($chrom=~/^GL/);

        print "processing $chrom\t$len ";

        push(@CHROM, $chrom);

        $intervalsByChrom{$chrom} = ();
        my $count = 0;
        for my $i (0 .. floor($len/$intervalWidth))
        {
            my $interval = "";
            my $intervalName = "";
            my $file = "";
            if ($i<floor($len/$intervalWidth))
            {
                $interval = $chrom . ":" . ($intervalWidth*$i+1) . "-" . ($intervalWidth*($i+1));
                $intervalName = $chrom . "_" . ($intervalWidth*$i+1) . "_" . ($intervalWidth*($i+1));
            }
            elsif ($i*$intervalWidth!=$len)
            {
                $interval = $chrom . ":" . ($intervalWidth*$i+1) . "-" . $len;
                $intervalName = $chrom . "_" . ($intervalWidth*$i+1) . "_" . $len;
            }
            else
            {
                last;
            }
            
            push(@{$intervalsByChrom{$chrom}}, "$intervalName");
            push(@intervals, $interval);
            push(@intervalNames, $intervalName);
            push(@intervalFiles, $file);

            $count++;
        }

        print "added $count intervals\n";
    }
}
close(SQ);

my @tgts = ();
my @deps = ();
my @cmds = ();
my $tgt;
my $dep;
my @cmd;
my $inputVCFFile;
my $outputVCFFile;

########
#Calling
########

#**************
#log start time
#**************
$tgt = "$logDir/start.calling.OK";
$dep = "";
@cmd = ("date | awk '{print \"samtools variant calling pipeline\\n\\nstart calling: \"\$\$0}' > $logFile");
makeLocalStep($tgt, $dep, @cmd);

my $intervalVCFFilesOK = "";
for my $i (0 .. $#intervals)
{
    $outputVCFFile = "$vcfOutDir/$intervalNames[$i].genotypes.vcf.gz";
    $tgt = "$outputVCFFile.OK";
    $dep = "";
    @cmd = ("mpileup -ugf $refGenomeFASTAFile -b $bamListFile -r $intervals[$i] |call -vmO z -o $outputVCFFile"),
    makeJob($partition, $tgt, $dep, @cmd);

    $intervalVCFFilesOK .= " $outputVCFFile.OK";
}

#************
#log end time
#************
$tgt = "$logDir/end.calling.OK";
$dep = "$intervalVCFFilesOK";
@cmd = ("date | awk '{print \"end: \"\$\$0}' >> $logFile");
makeLocalStep($tgt, $dep, @cmd);

###########################################
#Concatenate, normalize and drop duplicates
###########################################

#**************
#log start time
#**************
$tgt = "$logDir/start.concatenating.normalizing.OK";
$dep = "$logDir/end.calling.OK";
@cmd = ("date | awk '{print \"start concatenating and normalizing: \"\$\$0}' >> $logFile");
makeLocalStep($tgt, $dep, @cmd);


my $inputVCFFiles = join(" ", map {"$finalVCFOutDir/$_.genotypes.vcf.gz"} @CHROM);
my $inputVCFFilesOK = join(" ", map {"$finalVCFOutDir/$_.genotypes.vcf.gz.OK"} @CHROM);


#************
#log end time
#************
$tgt = "$logDir/end.concatenating.normalizing.OK";
$dep = "$inputVCFFile.tbi.OK";
@cmd = ("date | awk '{print \"end concatenating and normalizing: \"\$\$0}' >> $logFile");
makeLocalStep($tgt, $dep, @cmd);

#*******************
#Write out make file
#*******************
open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
print MAK ".DELETE_ON_ERROR:\n\n";
print MAK "all: @tgts\n\n";

#clean
push(@tgts, "clean");
push(@deps, "");
push(@cmds, "\t-rm -rf $outputDir/*.* $vcfOutDir/*.* $vcfOutDir/*/*.* $finalVCFOutDir/*.* $statsDir/* $logDir/* $outputDir/intervals/*.*");

for(my $i=0; $i < @tgts; ++$i) {
    print MAK "$tgts[$i] : $deps[$i]\n";
    print MAK "$cmds[$i]\n";
}
close MAK;

##########
#Functions
##########

#run a job either locally or by slurm
sub makeJob
{
    my ($method, $tgt, $dep, @cmd) = @_;

    if ($method eq "local")
    {
        makeLocalStep($tgt, $dep, @cmd);
    }
    else
    {
        makeSlurm($partition, $tgt, $dep, @cmd);
    }
}

#run slurm jobs
sub makeSlurm
{
    my ($partition, $tgt, $dep, @cmd) = @_;

    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        #contains pipe
        if ($c=~/\|/)
        {
            ++$slurmScriptNo;
            my $slurmScriptFile = "$slurmScriptsDir/$slurmScriptNo.sh";
            open(IN, ">$slurmScriptFile");
            print IN "#!/bin/bash\n";
            print IN "set -o pipefail; $c";
            close(IN);
            chmod(0755, $slurmScriptFile);

            $cmd .= "\techo '" . $c . "'\n";
            $cmd .= "\tsrun -p $partition $slurmScriptFile\n";
        }
        else
        {
            $cmd .= "\tsrun -p $partition " . $c . "\n";
        }
    }
    $cmd .= "\ttouch $tgt\n";
    push(@cmds, $cmd);
}

#run a local job
sub makeLocalStep
{
    my ($tgt, $dep, @cmd) = @_;

    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        $cmd .= "\tset -o pipefail; " . $c . "\n";
    }
    $cmd .= "\ttouch $tgt\n";
    push(@cmds, $cmd);
}

#run a local phony job
sub makePhonyJob
{
    my ($tgt, $dep, @cmd) = @_;

    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        $cmd .= "\t" . $c . "\n";
    }
    push(@cmds, $cmd);
}
