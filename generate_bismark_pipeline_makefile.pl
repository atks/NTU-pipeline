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
my $outputDir;
my $splitLineNo;
my $makeFile;
my $refGenomeDir;
my $pipelineName;

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

printf("generate_bismarck_pipeline_makefile.pl\n");
printf("\n");
printf("options: output dir           %s\n", $outputDir);
printf("         make file            %s\n", $makeFile);
printf("         split line no        %s\n", $splitLineNo);
printf("         sample file          %s\n", $sampleFile);
printf("         reference            %s\n", $refGenomeDir);
printf("\n");

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
#bismark_genome_preparation --path_to_bowtie $bowtie2Path --verbose $refGGenomeDir

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
    mkpath("$outputDir/samples/$sampleID/aligned");   
    mkpath("$outputDir/samples/$sampleID/dedup");   
    mkpath("$outputDir/samples/$sampleID/methylation");   
    mkpath("$outputDir/samples/$sampleID/sorted");   
}

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
        @cmd = ("$zsplit -s $sampleID -o $outputDir -l $splitLineNo -c $outputFile $SAMPLE{$sampleID}{FASTQ1} $SAMPLE{$sampleID}{FASTQ2}");
        makeJob("namedPBS", $pipelineName, $tgt, $dep, $log, $err, "24:00:00", 1, "1G", @cmd);
    }
    
    $makeFile = "$outputDir/count_lines.mk";
    print "Please run \"make -f count_lines.mk -j 8 -k\n";
    goto GENERATE_MAKEFILE;
}

print "generate commands for trim galore\n";

for my $sampleID (@SAMPLE)
{
    #######
    #fastQC
    #######
    my $fastqcOutputDir = "$outputDir/samples/$sampleID/fastqc_output";
    $tgt = "$fastqcOutputDir/fastqc1.OK";
    $dep = "$SAMPLE{$sampleID}{FASTQ1}";
    $log = "$fastqcOutputDir/fastqc1.log";
    $err = "$fastqcOutputDir/fastqc1.err";
    @cmd = ("$fastqc $SAMPLE{$sampleID}{FASTQ1} -o $fastqcOutputDir");
    makeJob("namedPBS", $pipelineName, $tgt, $dep, $log, $err, "24:00:00", 1, "1G", @cmd);
    $tgt = "$fastqcOutputDir/fastqc2.OK";
    $dep = "$SAMPLE{$sampleID}{FASTQ2}";
    $log = "$fastqcOutputDir/fastqc2.log";
    $err = "$fastqcOutputDir/fastqc2.err";
    @cmd = ("$fastqc $SAMPLE{$sampleID}{FASTQ2} -o $fastqcOutputDir");
    makeJob("namedPBS", $pipelineName, $tgt, $dep, $log, $err, "24:00:00", 1, "1G", @cmd);
    
    ###########
    #trimgalore
    ###########
    my $noFile = ceil($SAMPLE{$sampleID}{NO_LINES}/$splitLineNo);
    
    print "$sampleID \t no files = $noFile\n";
    my $splitTrimmedOKFiles = "";
    my $splitTrimmedR1FASTQFiles = "";
    my $splitTrimmedR2FASTQFiles = "";
   
    my $splitTrimmedAlignedOKFiles = "";
    my $splitTrimmedAlignedBAMFiles = "";
    for my $i (1 .. $noFile)
    {
        ##############
        #trimmed reads
        ##############
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
        
        $splitTrimmedR1FASTQFiles .= $i==0 ? "$trimmedR1File" : " $trimmedR1File";
        $splitTrimmedR2FASTQFiles .= $i==0 ? "$trimmedR2File" : " $trimmedR2File";
        $dep = "$splitDir/split.OK";
        $log = "$trimGaloreOutputDir/trim_galore.log";
        $err = "$trimGaloreOutputDir/trim_galore.err";
        @cmd = ("$trimGalore -o $trimGaloreOutputDir " .
                            "--path_to_cutadapt $cutAdapt " .
                            "--keep --illumina --clip_R2 18 --three_prime_clip_R1 18 --phred33 " .
                            "--paired $R1File $R2File");
        makeJob("namedPBS", $pipelineName, $tgt, $dep, $log, $err, "24:00:00", 1, "1G", @cmd);
        
        ####################
        #align trimmed reads
        ####################
        my $alignedOutputDir = "$outputDir/samples/$sampleID/aligned/$i";
        my $tempAlignedOutputDir = "$alignedOutputDir/temp";
        $tgt = "$alignedOutputDir/aligned.OK";
        $splitTrimmedAlignedOKFiles .= $i==0 ? "$alignedOutputDir/aligned.OK" : " $alignedOutputDir/aligned.OK";
        
        my ($file, $dir, $suffix) = fileparse($trimmedR1File, (".fq.gz"));
        my $bamFile = $alignedOutputDir . "/" . $file . "_bismark_bt2_pe.bam";
        
#        print "BAM: $bamFile  \n";
#        13_JC1_Macr_Nova_Swiftbio_indexed_R1_val_1_bismark_bt2_pe.bam
#        /home/users/ntu/adrianta/12000713/20200807_wgbs_pilot/samples/JC1/trim_galore_output/13/
#        13_JC1_Macr_Nova_Swiftbio_indexed_R1_val_1.fq.gz
        $splitTrimmedAlignedBAMFiles .= $i==0 ? "$bamFile" : " $bamFile";
   
        $dep = "$trimGaloreOutputDir/trim_galore.OK";
        $log = "$alignedOutputDir/aligned.log";
        $err = "$alignedOutputDir/aligned.err";
        @cmd = ("$bismark --bowtie2 -p 4 --bam " . 
                    "--temp_dir $tempAlignedOutputDir --un --ambiguous " . 
                    "--score_min L,0,-0.2 " .
                    "--path_to_bowtie $bowtie2Path " .
                    "-o $alignedOutputDir " . 
                    "--genome_folder $refGenomeDir " .
                    "-1 $trimmedR1File -2 $trimmedR2File");
        makeJob("namedPBS", $pipelineName, $tgt, $dep, $log, $err, "24:00:00", 1, "1G", @cmd);    
        
        
    }

    ######################
    #combined aligned bams
    ######################
    my $alignedOutputDir = "$outputDir/samples/$sampleID/aligned";
    my $outputBamFile = "$alignedOutputDir/$sampleID.name_sorted.bam";  
    $tgt = "$outputBamFile.OK";
    $dep = $splitTrimmedAlignedOKFiles;
    $log = "$alignedOutputDir/merge_name_sorted.log";
    $err = "$alignedOutputDir/merge_name_sorted.err";
    @cmd = ("$samtools merge -nf $outputBamFile $splitTrimmedAlignedBAMFiles");
    makeJob("namedPBS", $pipelineName, $tgt, $dep, $log, $err, "24:00:00", 1, "1G", @cmd);

    ############
    #deduplicate
    ############
    my $dedupOutputDir = "$outputDir/samples/$sampleID/dedup";
    my $inputBAMFile = "$outputDir/samples/$sampleID/aligned/$sampleID.name_sorted.bam";  
    $tgt = "$dedupOutputDir/dedup.OK";
    $dep = "$inputBAMFile.OK";
    $log = "$dedupOutputDir/dedup.log";
    $err = "$dedupOutputDir/dedup.err";
    @cmd = ("$bismarkPath/deduplicate_bismark -p --bam $inputBAMFile --output_dir $dedupOutputDir");
    makeJob("namedPBS", $pipelineName, $tgt, $dep, $log, $err, "24:00:00", 1, "1G", @cmd);
     
    ####################
    #extract methylation
    ####################    
    my $methylationOutputDir = "$outputDir/samples/$sampleID/methylation";
    $inputBAMFile = "$outputDir/samples/$sampleID/dedup/$sampleID.name_sorted.deduplicated.bam";  
    $tgt = "$methylationOutputDir/methylation.OK";
    $dep = "$outputDir/samples/$sampleID/dedup/dedup.OK";
    $log = "$methylationOutputDir/methylation.log";
    $err = "$methylationOutputDir/methylation.err";
    @cmd = ("$bismarkPath/bismark_methylation_extractor " .
             "-p --multicore 4 --no_overlap " . 
             "-o $methylationOutputDir " . 
             "--comprehensive --merge_non_CpG " . 
             "--cutoff 1 --buffer_size 40G --zero_based --cytosine_report " .
             "--genome_folder $refGenomeDir " .
             "$inputBAMFile");
    makeJob("namedPBS", $pipelineName, $tgt, $dep, $log, $err, "24:00:00", 4, "96G", @cmd);
 	#$BISMARK_PATH/bismark_methylation_extractor 
 	#-p --multicore 4 --no_overlap 
 	#-o $METH_OUTPUT_DIR --comprehensive --merge_non_CpG 
 	#--cutoff 1 --buffer_size 40G --zero_based --cytosine_report --genome_folder $GENOME_PATH $bam_file
 
    #########
    #sort bam
    #########
    my $sortedBAMOutputDir = "$outputDir/samples/$sampleID/sorted";
    $inputBAMFile = "$outputDir/samples/$sampleID/dedup/$sampleID.name_sorted.deduplicated.bam";  
    my $outputBAMFile = "$sortedBAMOutputDir/$sampleID.bam";  
    $tgt = "$sortedBAMOutputDir/sorted.OK";
    $dep = "$outputDir/samples/$sampleID/dedup/dedup.OK";
    $log = "$sortedBAMOutputDir/sorted.log";
    $err = "$sortedBAMOutputDir/sorted.err";
    @cmd = ("$samtools sort $inputBAMFile -o $outputBAMFile");
    makeJob("namedPBS", $pipelineName, $tgt, $dep, $log, $err, "24:00:00", 1, "1G", @cmd);

    ############################
    #combine trimmed fastq files
    ############################
    my $trimGaloreOutputDir = "$outputDir/samples/$sampleID/trim_galore_output";
    
    my $trimmedR1File = "$trimGaloreOutputDir/trimmed_$sampleID" . "_1.fastq.gz";
    $tgt = "$trimmedR1File.OK";
    $dep = "$splitTrimmedOKFiles";
    $log = "$trimGaloreOutputDir/zcat.log";
    $err = "$trimGaloreOutputDir/zcat.err";
    @cmd = ("zcat $splitTrimmedR1FASTQFiles | gzip -c > $trimmedR1File");
    makeJob("namedPBS", $pipelineName, $tgt, $dep, $log, $err, "24:00:00", 1, "1G", @cmd);
    
    my $trimmedR2File = "$trimGaloreOutputDir/trimmed_$sampleID" . "_2.fastq.gz";
    $tgt = "$trimmedR2File.OK";
    $dep = "$splitTrimmedOKFiles";
    $log = "$trimGaloreOutputDir/zcat.log";
    $err = "$trimGaloreOutputDir/zcat.err";
    @cmd = ("zcat $splitTrimmedR2FASTQFiles | gzip -c > $trimmedR2File");
    makeJob("namedPBS", $pipelineName, $tgt, $dep, $log, $err, "24:00:00", 1, "1G", @cmd);
        
    #####################
    #fastQC trimmed files
    #####################
    $fastqcOutputDir = "$outputDir/samples/$sampleID/trimmed_fastqc_output";
    $tgt = "$fastqcOutputDir/fastqc1.OK";
    $dep = "$trimmedR1File.OK";
    $log = "$fastqcOutputDir/fastqc1.log";
    $err = "$fastqcOutputDir/fastqc1.err";
    @cmd = ("$fastqc $trimmedR1File -o $fastqcOutputDir");
    makeJob("namedPBS", $pipelineName, $tgt, $dep, $log, $err, "24:00:00", 1, "1G", @cmd);
    $tgt = "$fastqcOutputDir/fastqc2.OK";
    $dep = "$trimmedR2File.OK";
    $log = "$fastqcOutputDir/fastqc2.log";
    $err = "$fastqcOutputDir/fastqc2.err";
    @cmd = ("$fastqc $trimmedR2File -o $fastqcOutputDir");
    makeJob("namedPBS", $pipelineName, $tgt, $dep, $log, $err, "24:00:00", 1, "1G", @cmd);    
}

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
    my ($method, @others) = @_;
#    print "method: $method \n";
#    print "\@other: @others \n";

    if ($method eq "local")
    {
        my ($tgt, $dep, @rest) = @others;
        makeLocalStep($tgt, $dep, @rest);
    }
    elsif ($method eq "pbs")
    {
        my ($tgt, $dep, @rest) = @others;
        makePBSStep($tgt, $dep, @rest);
    }
    elsif ($method eq "namedPBS")
    {
        
        my ($name, $tgt, $dep, @rest) = @others;  
        
#        print "\t name: $name \n";
#        print "\t tgt: $tgt \n";
#        print "\t dep: $dep \n";
#        print "\t \@rest: @rest \n";
              
        makeNamedPBSStep($name, $tgt, $dep, @rest);
    }
    else
    {
        die "unrecognized method of job creation : $method\n";
    }
}

sub makePBSStep
{
    my ($tgt, $dep, $log, $err, $walltime, $ncpu, $mem, @cmd) = @_;
    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = join(";", @cmd);
    $cmd = "set -o pipefail;" . join(";", @cmd);
    $cmd = "\techo \"$cmd\" | qsub -q normal -P 12000713 -W block=true -o $log -e $err -l select=1:ncpus=$ncpu:mem=$mem,walltime=$walltime\n";
    $cmd .= "\ttouch $tgt\n";
    push(@cmds, $cmd);
}

sub makeNamedPBSStep
{
    my ($name, $tgt, $dep, $log, $err, $walltime, $ncpu, $mem, @cmd) = @_;

#    print "\t\t name: $name \n";
#    print "\t\t tgt: $tgt \n";
#    print "\t\t dep: $dep \n";
#    print "\t\t log: $log \n";
#    print "\t\t err: $err \n";
#    print "\t\t wt : $walltime \n";
#    print "\t\t \@cmd: @cmd \n";
        
    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = join(";", @cmd);
    $cmd = "set -o pipefail;" . join(";", @cmd);
    $cmd = "\techo \"$cmd\" | qsub -q normal -P 12000713 -W block=true -N $name -o $log -e $err -l select=1:ncpus=$ncpu:mem=$mem,walltime=$walltime\n";
    $cmd .= "\ttouch $tgt\n";
    
    push(@cmds, $cmd);
}

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
