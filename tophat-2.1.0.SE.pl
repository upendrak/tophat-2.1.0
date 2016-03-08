#!/usr/bin/perl -w
use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

use constant TOPHATP    => ('2.1.0'  => '/home/upendra_35/tophat-2.1.0/tophat-2.1.0.Linux_x86_64/');
use constant BOWTIEP   => ('2.2.5'  => '/home/upendra_35/tophat-2.1.0/bowtie2-2.2.5/');

use constant SAMTOOLSP  => '/usr/bin/samtools/';

# remember to remove this
#report_input_stack();

my (@file_query, $database_path, $user_database_path, $annotation_path, 
$user_annotation_path, $file_names, $root_names, $null);

my $format    = 'SE';
my $version   = '2.1.0';
my $btversion = '2.2.5';

GetOptions( "file_query=s"      => \@file_query,
	    "database=s"        => \$database_path,
	    "user_database=s"   => \$user_database_path,
            "annotation=s"      => \$annotation_path,
            "user_annotation=s" => \$user_annotation_path,
	    "file_names=s"      => \$file_names,
	    "root_names=s"      => \$root_names,
            "tophat_version=s"  => \$version,
	    "bowtie_version=s"  => \$btversion
	    );

if (!($user_database_path || $database_path)) {
    die "No reference genome was supplied\n";
}
if (@file_query < 1) {
    die "No FASTQ files were supplied\n";
}

my %tophat = TOPHATP;
my %bowtie = BOWTIEP;

my $tophatp   = $tophat{$version};
my $bowtiep   = $bowtie{$btversion};

chomp($ENV{PATH} = `echo \$PATH`);
$ENV{PATH} = join(':',$ENV{PATH},$tophatp,$bowtiep,SAMTOOLSP);


# Sanity check for input ref. genome
unless ($database_path || $user_database_path) {
  die "No reference genome was selected" 
}

# Allow over-ride of system-level database path with user
if ($user_database_path) {
  $database_path = $user_database_path;
  unless (`grep \\> $database_path`) {
      die "Error: $database_path  the user supplied file is not a FASTA file";
  }
  my $name = basename($database_path, qw/.fa .fas .fasta .fna/);
  print STDERR "bowtie-indexing $name\n";
  $bowtiep .= $btversion eq '2.2.5' ? 'bowtie2-build' : 'bowtie_build';
  system $bowtiep . " $database_path $name";
  if ($database_path !~ /$name\.fa$/) {
      my $new_path = $database_path;
      $new_path =~ s/$name\.\S+$/$name\.fa/;
      system "cp $database_path $new_path";
  }
  $database_path = $name;
}
if ($user_annotation_path) {
    $annotation_path = $user_annotation_path;
}

# Should be a temporary hack
#if ($btversion eq '2.0.0' && !$user_database_path) {
#    my $name = basename($database_path, qw/.fa .fas .fasta .fna/);
#    print STDERR "bowtie-indexing $name\n";
#    system $bowtiep . "bowtie2-build $database_path $name";
#    $database_path = $name;
#}

# If we are using new tophat with old bowtie
#if ($version eq '2.0.5' && $btversion eq '0.12.7') {
#    push @ARGV, '--bowtie1';
#}

my $success = undef;

my (@basenames,%sample);

if ($root_names) {
    my @names = split(',',$root_names);
    for my $name (@names) {
	next if grep {/^$name|\/$name/} @file_query;
        die "Root name $name does not match any query file in the list\n";
    }
}
elsif ($file_names) {
    my $idx;
    my @names = split (/\s+/,$file_names); 
    for my $sample (@names) {
	$idx++;
	my @files = split(',',$sample);
	for (@files) {
	    $sample{$_} = "sample$idx";
	}
    }
    for my $file (@file_query) {
	my $f = $file;
	$f =~ s!\S+/!!;
	next if $sample{$f};
	die "$file does not match any name in the list\n";
    }
}


my $nocount;
my $samples;
#my @to_move = ('bam');
for my $query_file (@file_query) {
    # Grab any flags or options we don't recognize and pass them as plain text
    # Need to filter out options that are handled by the GetOptions call
    my @args_to_reject = qw(-xxxx);
    my $TOPHAT_ARGS = join(" ", @ARGV);
    foreach my $a (@args_to_reject) {
	if ($TOPHAT_ARGS =~ /$a/) {
	    report("Most TopHat arguments are legal for use with this script, but $a is not. Please omit it and submit again");
	    exit 1;
	}
    }

    my $app  = $tophatp.'tophat';
    if ($annotation_path) {
	$TOPHAT_ARGS .= " -G $annotation_path";
    }
    my $align_command = "$app $TOPHAT_ARGS $database_path $query_file";
    
    chomp(my $basename = `basename $query_file`);
    $basename =~ s/\.\S+$//;
    #push @to_move, $basename;

    report("Executing: $align_command\n");
    system $align_command;
    system "mv tophat_out $basename";
    $success++ if -e "$basename/accepted_hits.bam";
    #my $bam = 'bam';
    #mkdir($bam) unless -d $bam;

    #my %matched;
    #if ($root_names) {
    #	my @names = split(',',$root_names);
    #	for my $root (@names) {
    #	    if ($query_file =~ /^$root|\/$root/) {
		#system "mkdir $bam/$root" unless -d "$bam/$root";
		#system "cp $basename/accepted_hits.bam $bam/$root/$basename.bam";
    #		$samples .= join("\t", $root, "$bam/$root/$basename.bam")."\n";
    #		last;
    #	    }
    #	}
    #}
    #elsif ($file_names) {
    #	my @samples = map {[split(',',$_)]} split(/\s+/,$file_names);
    #	my $idx;
    #	for (@samples) {
    #	    $idx++ unless $nocount;
    #	    next unless grep {/$query_file/} @$_;
    #	    system "mkdir $bam/sample$idx" unless -e "$bam/sample$idx";
    #	    system "cp $basename/accepted_hits.bam $bam/sample$idx/$basename.bam";
    # 	    $samples.= join("\t", "sample$idx", "$bam/sample$idx/$basename.bam")."\n";
    #	    $nocount = undef;
    #	    last;
    #	}
    #}
    #else {
    #	system "cp $basename/accepted_hits.bam $bam/$basename.bam";
    #}
}



#mkdir 'tophat_out' unless -d 'tophat_out';

#mkdir 'bam';
#for (@to_move) {
#    system "mv $_ tophat_out" or warn $!;
#    system "ln -s $_/accepted_hits.bam bam/$_.bam";
#}

system "rm -f *.ebwt *.bt2 2>/dev/null";

$success ? exit 0 : exit 1;

sub report {
    print STDERR "$_[0]\n";
}

sub report_input_stack {
    my @stack = @ARGV;
    report(Dumper \@stack);
}
