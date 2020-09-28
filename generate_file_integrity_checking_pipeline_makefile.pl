#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

generate_file_integrity_checking_pipeline_makefile

=head1 SYNOPSIS

 generate_file_integrity_checking_pipeline_makefile [options]

  -d     directory containing files
  -e     type of file to compute md5 hash
  -c     file name of expected MD5 hashes with file names
  -o     output directory

=head1 DESCRIPTION

This script implements the pipeline for calculating md5sum.
1. A directory is input.
2. The directory is searched for files that have to be checked.
3. The directory is searched for files that contain the expected
4. The results are stored in the output directory

=cut

my $help;
my $outputDir; 
my $dataDir;
my $fileExtension;
my $md5File;
my $makeFile = "run_check_file_integrity_pipeline.mk";
                
#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help,
                'd:s'=>\$dataDir,
                'e:s'=>\$fileExtension,
                'c:s'=>\$md5File,
                'o:s'=>\$outputDir,
                'm:s'=>\$makeFile,
               )
  || !defined($dataDir)
  || !defined($fileExtension)
  || !defined($md5File)
  || !defined($outputDir))
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

printf("generate_file_integrity_checking_pipeline_makefile\n");
printf("\n");
printf("options: output dir         %s\n", $outputDir);
printf("         data dir           %s\n", $dataDir);
printf("         file extension     %s\n", $fileExtension);
printf("         md5File            %s\n", $md5File);
printf("         makeFile           %s\n", $makeFile);
printf("\n");

mkpath($outputDir);

################################################
#Helper data structures for generating make file
################################################
my @tgts = ();
my @deps = ();
my @cmds = ();
my $tgt;
my $dep;
my @cmd;

#################
#Search for files 
#################
mkpath($outputDir);
my @md5sum_files = split ("\n" ,`find $dataDir -name $md5File | xargs cat`);
#print @md5sum_files;
#create hash
my %FILE = ();
my $no = 0;
for my $line (@md5sum_files)
{
    my ($hash, $fileName) = split(" ", $line);
#    print "$no : $hash : $fileName\n";
    $FILE{$fileName} = $hash;
    $no++;
}

my $files = `find $dataDir -name "*.$fileExtension"`;
my @files = split("\n", $files);
$no = 0;
my $md5sumTxtFile = "";
for my $file (@files)
{
    $no++;
    my ($name, $path, $suffix) = fileparse($file);
    if (exists($FILE{$name}))
    {
        $md5sumTxtFile .=  ($no==1 ? "" : "\n") . "$FILE{$name}  $file";
        delete $FILE{$name};
    }    
    else
    {
        print "$file has no expected md5 hash\n";
    }    
}

if (scalar(keys(%FILE))!=0)
{
    print "Some expected files are missing.\n";
    print join("\n", keys(%FILE)) . "\n";
}
else
{
    print "$no files with expected md5 hash to check\n";
}

my $result = `echo \"$md5sumTxtFile\" | sort -k2 > $outputDir/expected_md5sum.txt`;

if ($?) 
{
    print "error with making expected md5 file\n";
}

#################
#Compute md5 hash 
#################
$no = 1;
for my $file (@files)
{
    my $outputFile = "$outputDir/$no.md5";
    $tgt = "$outputFile.OK";
    $dep = $file;
    @cmd = ("md5sum $file > $outputFile");
    makeLocalStep($tgt, $dep, @cmd);
    $no++;
}

###################
#Aggregate md5 hash
###################
my $outputFile = "$outputDir/computed_md5sum.txt";
$tgt = "$outputFile.OK";
$dep = join(" ", map {$_="$outputDir/$_.md5.OK"} (1 .. scalar(@files) ));
@cmd = ("cat " . join(" ", map {$_="$outputDir/$_.md5"} (1 .. scalar(@files) )) . "| sort -k2 > $outputFile");
makeLocalStep($tgt, $dep, @cmd);

#########################
#Compare old and new hash
#########################
my $outputFile = "$outputDir/computed_md5sum.txt";
$tgt = "$outputDir/diff.OK";
$dep = "$outputDir/computed_md5sum.txt.OK";
@cmd = ("diff $outputDir/cßomputed_md5sum.txt $outputDir/expected_md5sum.txt > diff.log 2> diff.err");
makeLocalStep($tgt, $dep, @cmd);

#*******************
#Write out make file
#*****************
GENERATE_MAKEFILE:
print "\nwriting makefile\n";

open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
print MAK ".DELETE_ON_ERROR:\n\n";
print MAK "all: @tgts\n\n";

#clean
push(@tgts, "clean");
push(@deps, "");
push(@cmds, "\t-rm -rf $outputDir/*.*");

for(my $i=0; $i < @tgts; ++$i) {
    print MAK "$tgts[$i] : $deps[$i]\n";
    print MAK "$cmds[$i]\n";
}
close MAK;

##########
#Functions
##########

#run PBS jobs
sub makePBSStep
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