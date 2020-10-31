#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

generate_coverage_analysis_makefile

=head1 SYNOPSIS

 generate_coverage_analysis_makefile [options]

  -s     sample file list giving the location of each sample
         column 1: sample name
         column 2: path of bam file
  -r     reference genome file
  -l     sequence length file
  -w     interval width
  -o     output directory
  -m     make file name

=head1 DESCRIPTION

Coverage analysis

=cut

my $help;

my $sampleFile;
my $outputDir;
my $makeFile;
my $pipelineName;

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help,
                'o=s'=>\$outputDir,
                's=s'=>\$sampleFile,
                'n:s'=>\$pipelineName,
                'm:s'=>\$makeFile
               )
  || !defined($outputDir)
  || !defined($sampleFile)
  || !defined($makeFile)
  || !defined($pipelineName))
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
my $samtools = "/app/samtools/1.3/bin/samtools";
my $picard = "java -Xmx2g -jar /home/users/ntu/adrianta/programs/picard-1.141/picard.jar";

#reference files
my $refGenomeFASTAFile = "/home/users/ntu/adrianta/12000713/ref/hg19/hg19.fa"; 

printf("generate_coverage_analysis_makefile.pl\n");
printf("\n");
printf("options: output dir           %s\n", $outputDir);
printf("         make file            %s\n", $makeFile);
printf("         sample file          %s\n", $sampleFile);
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
        my ($sampleID, $bamPath) = split(/\s+/, $_);

        if (exists($SAMPLE{$sampleID}))
        {
            exit("$sampleID already exists. Please fix.");
        }

        $SAMPLE{$sampleID}{BAM} = $bamPath;

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
}

for my $sampleID (@SAMPLE)
{
    ###############
    #samtools stats
    ###############
    my $bamFileOutputDir = "$outputDir/samples/$sampleID/";
    $tgt = "$bamFileOutputDir/st_stats.txt.OK";
    $dep = "$SAMPLE{$sampleID}{BAM}";
    $log = "$bamFileOutputDir/st_stats.log";
    $err = "$bamFileOutputDir/st_stats.err";
    @cmd = ("$samtools stats $SAMPLE{$sampleID}{BAM} > $bamFileOutputDir/st_stats.txt");
    makeJob("namedPBS", $pipelineName, $tgt, $dep, $log, $err, "24:00:00", 1, "1G", @cmd);
    
    ######
    #depth
    ######
    $tgt = "$bamFileOutputDir/st_depth.txt.OK";
    $dep = "$SAMPLE{$sampleID}{BAM}";
    $log = "$bamFileOutputDir/st_depth.log";
    $err = "$bamFileOutputDir/st_depth.err";
    @cmd = ("$samtools depth -q 20 -Q 20 $SAMPLE{$sampleID}{BAM} > $bamFileOutputDir/st_depth.txt");
    makeJob("namedPBS", $pipelineName, $tgt, $dep, $log, $err, "24:00:00", 1, "1G", @cmd);

#    $tgt = "$bamFileOutputDir/st_coverage.txt.OK";
#    $dep = "$bamFileOutputDir/st_depth.txt.OK";
#    $log = "$bamFileOutputDir/st_coverage.log";
#    $err = "$bamFileOutputDir/st_coverage.err";
#    @cmd = ("cat $bamFileOutputDir/st_depth.txt | perl -lane 'BEGIN{$$d=0;$$n=0}{$$d+=$$F[2]; $$n++}END{print \"$$n\t\" . $d/$n . \"\n\";}' > $bamFileOutputDir/st_coverage.txt");
#    makeJob("namedPBS", $pipelineName, $tgt, $dep, $log, $err, "24:00:00", 1, "1G", @cmd);

    $tgt = "$bamFileOutputDir/picard_wgs_metrics.txt.OK";
    $dep = "$SAMPLE{$sampleID}{BAM}";
    $log = "$bamFileOutputDir/picard_wgs_metrics.log";
    $err = "$bamFileOutputDir/picard_wgs_metrics.err";
    @cmd = ("$picard CollectWgsMetrics " . 
              "REFERENCE_SEQUENCE=$refGenomeFASTAFile " .
               "MINIMUM_MAPPING_QUALITY=0 " .
               "INPUT=$SAMPLE{$sampleID}{BAM} " .
               "OUTPUT=$bamFileOutputDir/picard_wgs_metrics.txt");
    makeJob("namedPBS", $pipelineName, $tgt, $dep, $log, $err, "24:00:00", 1, "4G", @cmd);

    $tgt = "$bamFileOutputDir/picard_wgs_metrics.map20.txt.OK";
    $dep = "$SAMPLE{$sampleID}{BAM}";
    $log = "$bamFileOutputDir/picard_wgs_metrics.map20.log";
    $err = "$bamFileOutputDir/picard_wgs_metrics.map20.err";
    @cmd = ("$picard CollectWgsMetrics " . 
              "REFERENCE_SEQUENCE=$refGenomeFASTAFile " .
               "MINIMUM_MAPPING_QUALITY=20 " .
               "INPUT=$SAMPLE{$sampleID}{BAM} " .
               "OUTPUT=$bamFileOutputDir/picard_wgs_metrics.map20.txt");
    makeJob("namedPBS", $pipelineName, $tgt, $dep, $log, $err, "24:00:00", 1, "4G", @cmd);





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
