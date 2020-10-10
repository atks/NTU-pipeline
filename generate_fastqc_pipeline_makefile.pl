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

This script implements the pipeline for Whole Genome Bisulfite Sequencing using Bismark.

=cut

my $help;

my $sampleFile;
my $outputDir = "/home/users/ntu/adrianta/12000713/20200915_wgbs";
#my $outputDir = "/home/users/ntu/adrianta/12000713/20200807_wgbs_pilot";
my $splitLineNo = 40000000;
#my $splitLineNo = 2000000;
#my $makeFile = "$outputDir/temp_run_pipeline.mk";
my $makeFile = "$outputDir/run_pipeline.mk";
my $refWholeSulfiteGenomeFASTAFile = "/home/users/ntu/adrianta/ref/hg19/Bisulfite_Genome";

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help,
                'o=s'=>\$outputDir,
                's=s'=>\$sampleFile,
                'd=i'=>\$splitLineNo,
                'm:s'=>\$makeFile,
                'r:s'=>\$refWholeSulfiteGenomeFASTAFile
               )
  || !defined($sampleFile)
  || !defined($refWholeSulfiteGenomeFASTAFile))
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
my $bowtie2 = "/app/bowtie2/2.29/bowtie2";
my $samtools = "/app/samtools/1.3/bin/samtools";
my $fastqc = "/home/users/ntu/adrianta/programs/fastqc-0.11.5/fastqc";
my $trimGalore = "/home/users/ntu/adrianta/programs/trimGalore-0.4.5/trim_galore";
my $cutAdaptPath = "/home/users/ntu/adrianta/programs/cutadapt-1.15";
my $cutAdapt = "/home/users/ntu/adrianta/programs/cutadapt-1.15/cutadapt";
my $zsplit = "/home/users/ntu/adrianta/programs/NTU-pipeline/zsplit.pl";

printf("generate_bismarck_pipeline_makefile.pl\n");
printf("\n");
printf("options: output dir           %s\n", $outputDir);
printf("         make file            %s\n", $makeFile);
printf("         split line no        %s\n", $splitLineNo);
printf("         sample file          %s\n", $sampleFile);
printf("         reference            %s\n", $refWholeSulfiteGenomeFASTAFile);
printf("\n");

#this pipeline generator generates 2 makefiles
my $preprocessMakeFile = 0;

################################################
#Helper data structures for generating make file
################################################
my @tgts = ();
my @deps = ();
my @cmds = ();
my $tgt;
my $dep;
my $log;
my $err;
my @cmd;
my $inputVCFFile;
my $outputVCFFile;

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

##########################
#create sample directories
##########################
for my $sampleID (@SAMPLE)
{
    mkpath("$outputDir/samples/$sampleID");
    mkpath("$outputDir/samples/$sampleID/split");
    mkpath("$outputDir/samples/$sampleID/fastqc_output");
    mkpath("$outputDir/samples/$sampleID/trim_galore_output");    
    mkpath("$outputDir/samples/$sampleID/trimmed_fastqc_output");
}

################################################
#Create/update/read augmented sample file
#
#a. lines_counted.OK exists and
#   a. augmented.sa exists   
#      a. age(lines_counted.OK) >  age(augmented.sa)    => overwrite
#      b. age(lines_counted.OK) <= age(augmented.sa)    => read                       
#   b. augmented.sa !exists                             => create     
#b. lines_counted.OK !exists                            => skip
#
################################################

###read augmented sample file
##if (-e $linesCountedFile)
##{
##    my $linesCountedFile = "$outputDir/lines_counted.OK";
##    my $lastTimeLinesCounted = `stat -t --printf="%Y" $linesCountedFile`;
##    my $createOrWriteAugmentedSampleFile = 0;
##    
##    if (-e $augmentedSampleFile)
##    {
##        my $lastTimeAugmentedSampleFileModified = `stat -t --printf="%Y" $augmentedSampleFile`;
##        
##        if ($lastTimeLinesCounted>$lastTimeAugmentedSampleFileModified)
##        {
##            $createOrWriteAugmentedSampleFile = 1;
##        }
##        else
##        {
##            #read augmented sample file
##            open(SA,"$sampleFile") || die "Cannot open $sampleFile\n";
##            while (<SA>)
##            {
##                s/\r?\n?$//;
##                if(!/^#/)
##                {
##                    my ($sampleID, $fastq1Path, $fastq2Path, $noLines) = split(/\s+/, $_);
##                    
##                    if (exists($SAMPLE{$sampleID}))
##                    {
##                        $SAMPLE{$sampleID}{NO_LINES} = $noLines;
##                    }
##                    else
##                    {
##                        warn("$sampleID not in the sample file submitted. This sample is ignored.");
##                    }
##                }
##            }
##            close(SA);
##        }
##    }
##    else
##    {
##        $createOrWriteAugmentedSampleFile = 1;
##    }
##    
##    #create/overwrite augmented sample file
##    if ($createOrWriteAugmentedSampleFile)
##    {
##        open(SA, ">$augmentedSampleFile") || die "Cannot open $augmentedSampleFile\n";
##        for my $sampleID (@SAMPLE)
##        {
##            my $lines = `cat "$outputDir/samples/$sampleID/no_lines.txt"`;
##            print SA "$sampleID\t$SAMPLE{$sampleID}{FASTQ1}\t$SAMPLE{$sampleID}{FASTQ2}\t$lines\n";
##        }
##        close(SA);  
##    }
##}
##else
##{
##    
##}

#####################
#look for file counts
#####################
for my $sampleID (@SAMPLE)
{
    my $augmentedSampleFile = "$outputDir/samples/$sampleID/split/augmented.sa";
    if (-e "$outputDir/samples/$sampleID/split/split.OK" && -e $augmentedSampleFile)
    {
        my ($sampleID, $file, $lineCount) = split("\t", `cat $augmentedSampleFile`);
        $SAMPLE{$sampleID}{NO_LINES} = $lineCount;
    }
}

##################################################################
#count the samples that require their fastq file to be precounted.
##################################################################
my @SAMPLES_TO_COUNTLINES = ();
for my $sampleID (@SAMPLE)
{
#    print "$sampleID\n";
    if (!exists($SAMPLE{$sampleID}{NO_LINES}))
    {
        push(@SAMPLES_TO_COUNTLINES, $sampleID);
    }    
}

my $noFilesToBeCounted = scalar(@SAMPLES_TO_COUNTLINES);
print "$noFilesToBeCounted" . "/" . scalar(@SAMPLE) . " samples to have their FASTQ files counted.\n";

############################
#Split files and count lines
############################
my $noLinesOKFiles = "";
if ($noFilesToBeCounted)
{
    for my $sampleID (@SAMPLES_TO_COUNTLINES)
    {
        my $sampleDir = "$outputDir/samples/$sampleID";
        my $outputDir = "$sampleDir/split";
        my $outputFile = "$outputDir/augmented.sa";
        $noLinesOKFiles .= "$sampleDir/split.OK ";
        $tgt = "$outputDir/split.OK";
        $dep = "$SAMPLE{$sampleID}{FASTQ1} $SAMPLE{$sampleID}{FASTQ2}";
        $log = "$outputDir/split.log";
        $err = "$outputDir/split.err";
        @cmd = ("$zsplit -s $sampleID -o $outputDir -l 40000000 -c $outputFile $SAMPLE{$sampleID}{FASTQ1} $SAMPLE{$sampleID}{FASTQ2}");
        #@cmd = ("$zsplit -s $sampleID -o $outputDir -l 2000000 -c $outputFile $SAMPLE{$sampleID}{FASTQ2} $SAMPLE{$sampleID}{FASTQ2}");
        makePBSCommandLine($tgt, $dep, $log, $err, "24:00:00", @cmd);
    }
       
#    $tgt = "$outputDir/split.OK";
#    $dep = "$noLinesOKFiles";
#    @cmd = ("echo $noFilesToBeCounted files counted.");
#    makeLocalStep($tgt, $dep, @cmd);
    
    $makeFile = "$outputDir/count_lines.mk";
    print "Please run \"make -f count_lines.mk -j 8 -k\n";


    goto GENERATE_MAKEFILE;
}

print "generate commands for trim galore\n";

for my $sampleID (@SAMPLE)
{
    ###########
    #trimgalore
    ###########
    my $noFile = ceil($SAMPLE{$sampleID}{NO_LINES}/40000000);
    
    print "$sampleID \t no files = $noFile\n";
    my $splitTrimmedOKFiles = "";
    my $splitTrimmedR1FASTQFiles = "";
    my $splitTrimmedR2FASTQFiles = "";
    for my $i (1 .. $noFile)
    {
        my $splitDir = "$outputDir/samples/$sampleID/split";
        my $trimGaloreOutputDir = "$outputDir/samples/$sampleID/trim_galore_output/$i";
        mkpath("$trimGaloreOutputDir");
        $tgt = "$trimGaloreOutputDir/trim_galore.OK";
        $splitTrimmedOKFiles .= $i==0 ? "$tgt" : " $tgt";
        my ($file1, $dir1, $suffix1) = fileparse($SAMPLE{$sampleID}{FASTQ1}, (".fastq.gz"));
        my ($file2, $dir2, $suffix2) = fileparse($SAMPLE{$sampleID}{FASTQ2}, (".fastq.gz"));
        my $R1File = "$splitDir/$i" . "_$file1$suffix1";
        my $R2File = "$splitDir/$i" . "_$file2$suffix2";
        
        my $trimmedR1File = "$trimGaloreOutputDir/$i" . "_$file1" . "_val_1.fq.gz";
        my $trimmedR2File = "$trimGaloreOutputDir/$i" . "_$file2" . "_val_2.fq.gz";
        
#        print "CHECK: $file1, $dir1, $suffix1\n";
#        print "     : $R1File\n";
#        print "     : $R2File\n";
#        print "     : $trimmedR1File\n";
#        print "     : $trimmedR2File\n";
        
        $splitTrimmedR1FASTQFiles .= $i==0 ? "$trimmedR1File" : " $trimmedR1File";
        $splitTrimmedR2FASTQFiles .= $i==0 ? "$trimmedR2File" : " $trimmedR2File";
        $dep = "$splitDir/split.OK";
        $log = "$trimGaloreOutputDir/trim_galore.log";
        $err = "$trimGaloreOutputDir/trim_galore.err";
        @cmd = ("$trimGalore -o $trimGaloreOutputDir " .
                            "--path_to_cutadapt $cutAdapt " .
                            "--keep --illumina --clip_R2 18 --three_prime_clip_R1 18 --phred33 " .
                            "--paired $R1File $R2File");
        makePBSCommandLine($tgt, $dep, $log, $err, "24:00:00", @cmd);
    }

    ##########################
    #combine all gzipped files
    ##########################
    my $trimGaloreOutputDir = "$outputDir/samples/$sampleID/trim_galore_output";
    
    my $trimmedR1File = "$trimGaloreOutputDir/trimmed_$sampleID" . "_1.fastq.gz";
    $tgt = "$trimmedR1File.OK";
    $dep = "$splitTrimmedOKFiles";
    $log = "$trimGaloreOutputDir/zcat.log";
    $err = "$trimGaloreOutputDir/zcat.err";
    @cmd = ("zcat $splitTrimmedR1FASTQFiles | gzip -c > $trimmedR1File");
    makePBSCommandLine($tgt, $dep, $log, $err, "24:00:00", @cmd);
    
    my $trimmedR2File = "$trimGaloreOutputDir/trimmed_$sampleID" . "_2.fastq.gz";
    $tgt = "$trimmedR2File.OK";
    $dep = "$splitTrimmedOKFiles";
    $log = "$trimGaloreOutputDir/zcat.log";
    $err = "$trimGaloreOutputDir/zcat.err";
    @cmd = ("zcat $splitTrimmedR2FASTQFiles | gzip -c > $trimmedR2File");
    makePBSCommandLine($tgt, $dep, $log, $err, "24:00:00", @cmd);
        
    #####################
    #fastQC trimmed files
    #####################
    my $fastqcOutputDir = "$outputDir/samples/$sampleID/trimmed_fastqc_output";
    $tgt = "$fastqcOutputDir/fastqc1.OK";
    $dep = "$trimmedR1File.OK";
    $log = "$fastqcOutputDir/fastqc1.log";
    $err = "$fastqcOutputDir/fastqc1.err";
    @cmd = ("$fastqc $trimmedR1File -o $fastqcOutputDir");
    makePBSCommandLine($tgt, $dep, $log, $err, "24:00:00", @cmd);
    $tgt = "$fastqcOutputDir/fastqc2.OK";
    $dep = "$trimmedR2File.OK";
    $log = "$fastqcOutputDir/fastqc2.log";
    $err = "$fastqcOutputDir/fastqc2.err";
    @cmd = ("$fastqc $trimmedR2File -o $fastqcOutputDir");
    makePBSCommandLine($tgt, $dep, $log, $err, "24:00:00", @cmd);
}



goto GENERATE_MAKEFILE;

#######
#fastQC
#######
for my $sampleID (@SAMPLE)
{
#    print "$sampleID\n";
    if (!exists($SAMPLE{$sampleID}{NO_LINES}))
    {
        my $fastqcOutputDir = "$outputDir/samples/$sampleID/fastqc_output";
        $tgt = "$fastqcOutputDir/fastqc1.OK";
        $dep = "$SAMPLE{$sampleID}{FASTQ1}";
        $log = "$fastqcOutputDir/fastqc1.log";
        $err = "$fastqcOutputDir/fastqc1.err";
        @cmd = ("$fastqc $SAMPLE{$sampleID}{FASTQ1} -o $fastqcOutputDir");
        makePBSCommandLine($tgt, $dep, $log, $err, "03:00:00", @cmd);
        $tgt = "$fastqcOutputDir/fastqc2.OK";
        $dep = "$SAMPLE{$sampleID}{FASTQ2}";
        $log = "$fastqcOutputDir/fastqc2.log";
        $err = "$fastqcOutputDir/fastqc2.err";
        @cmd = ("$fastqc $SAMPLE{$sampleID}{FASTQ2} -o $fastqcOutputDir");
        makePBSCommandLine($tgt, $dep, $log, $err, "03:00:00", @cmd);
        
    }
}

##################
#trimgalore
##################
for my $sampleID (@SAMPLE)
{
    if (!exists($SAMPLE{$sampleID}{NO_LINES}))
    {
        my $trimGaloreOutputDir = "$outputDir/samples/$sampleID/trim_galore_output";
        $tgt = "$trimGaloreOutputDir/trim_galore.OK";
        $dep = "$SAMPLE{$sampleID}{FASTQ1} $SAMPLE{$sampleID}{FASTQ2}";
        $log = "$trimGaloreOutputDir/trimgalore.log";
        $err = "$trimGaloreOutputDir/trimgalore.err";
        @cmd = ("$trimGalore -o $trimGaloreOutputDir " .
                            "--path_to_cutadapt $cutAdapt " .
                            "--illumina --clip_R2 18 --three_prime_clip_R1 18 --phred33 " .
                            "--paired $SAMPLE{$sampleID}{FASTQ1} $SAMPLE{$sampleID}{FASTQ2}");
        makePBSCommandLine($tgt, $dep, $log, $err, "48:00:00", @cmd);
    }
}


#####################################################
#Read files and count lines for augmented sample list
#####################################################





#$BISMARK_PATH/bismark_genome_preparation --path_to_bowtie $BOWTIE_PATH --verbose $GENOME_PATH
#/home/users/ntu/adrianta/programs/bismark-0.19.0/bismark_genome_preparation --path_to_bowtie  /app/bowtie2/2.29 --verbose ~/ref/

#my $bamListFile = "$auxDir/bam.list";
#open(OUT,">$bamListFile") || die "Cannot open $bamListFile\n";
#print OUT $bamFiles;
#close(OUT);
#
#print "read in " . scalar(keys(%SAMPLE)) . " samples\n";

#qsub -v LIB=$LIB $PBS_O_WORKDIR/processfastqc-raw.sh
##########
#Alignment
##########

#print "adding aligning steps\n";

#**************
#log start time
#**************
#$tgt = "$logDir/start.calling.OK";
#$dep = "";
#@cmd = ("date | awk '{print \"samtools variant calling pipeline\\n\\nstart calling: \"\$\$0}' > $logFile");
#makeLocalStep($tgt, $dep, @cmd);
#


###################
#Generate intervals
###################


#************
#log end time
#************
#$tgt = "$logDir/end.calling.OK";
##$dep = "$intervalVCFFilesOK";
#@cmd = ("date | awk '{print \"end: \"\$\$0}' >> $logFile");
#makeLocalStep($tgt, $dep, @cmd);

###########################################
#Concatenate, normalize and drop duplicates
###########################################

#**************
#log start time
#**************
#$tgt = "$logDir/start.concatenating.normalizing.OK";
#$dep = "$logDir/end.calling.OK";
#@cmd = ("date | awk '{print \"start concatenating and normalizing: \"\$\$0}' >> $logFile");
#makeLocalStep($tgt, $dep, @cmd);


#my $inputVCFFiles = join(" ", map {"$finalVCFOutDir/$_.genotypes.vcf.gz"} @CHROM);
#my $inputVCFFilesOK = join(" ", map {"$finalVCFOutDir/$_.genotypes.vcf.gz.OK"} @CHROM);


#************
#log end time
#************
#$tgt = "$logDir/end.concatenating.normalizing.OK";
#$dep = "$inputVCFFile.tbi.OK";
#@cmd = ("date | awk '{print \"end concatenating and normalizing: \"\$\$0}' >> $logFile");
#makeLocalStep($tgt, $dep, @cmd);

#*******************
#Write out make file
#*******************
GENERATE_MAKEFILE:
print "\nwriting makefile\n";

open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
print MAK ".DELETE_ON_ERROR:\n\n";
print MAK "all: @tgts\n\n";

#clean
push(@tgts, "clean");
push(@deps, "");
push(@cmds, "\t-rm -rf $outputDir/*.* $outputDir/intervals/*.*");

for(my $i=0; $i < @tgts; ++$i) {
    print MAK "$tgts[$i] : $deps[$i]\n";
    print MAK "$cmds[$i]\n";
}
close MAK;

##########
#Functions
##########

#run a job either locally or by pbs
sub makeJob
{
    my ($method, $tgt, $dep, @cmd) = @_;

    if ($method eq "local")
    {
        makeLocalStep($tgt, $dep, @cmd);
    }
    else
    {
        makePBS($tgt, $dep, @cmd);
    }
}

#run PBS jobs
sub makePBSCommandLine
{
    my ($tgt, $dep, $log, $err, $walltime, @cmd) = @_;
    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    $cmd = "set -o pipefail;" . join(";", @cmd);
    $cmd = "\techo \"$cmd\" | qsub -q normal -P 12000713 -W block=true -o $log -e $err -l select=1:ncpus=1,walltime=$walltime\n";
    $cmd .= "\ttouch $tgt\n";
    push(@cmds, $cmd);
}

#echo "zcat /home/projects/12000713/common/WGBS_pilot_2018/Novogene_Novaseq/EAL011_Novo_Nova_Swiftbio_indexed_R1.fastq.gz  | wc -l" |

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


