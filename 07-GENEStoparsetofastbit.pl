#!/usr/bin/perl
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL for extracting genes stuff from the database
#MOA 2017

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - U S E R  V A R I A B L E S- - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
use strict;
use Getopt::Long;
use Pod::Usage;
# DATABASE ATTRIBUTES

my $basepath = "/home/modupe/public_html/GenesAtlas";
my $finalpath = "/home/modupe/public_html/OtherTranscriptAtlas";
`mkdir -p $finalpath`;
my $statusfile = "$finalpath/myFBstatus";
open(STDOUT, '>', "$statusfile\.log") or die "Log file doesn't exist";
open(STDERR, '>', "$statusfile\.err") or die "Error file doesn't exist";
 
#open(STATUS,">$statusfile"); print STATUS "libraryid\ttotal\n"; close (STATUS);
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
my %HashDirectory;
my $library;
GetOptions ("lib|library=s" => \$library);
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
unless ($library) {
	opendir (DIR, $basepath) or die "Folder \"$basepath\" doesn't exist\n"; 
	my @Directory = readdir(DIR); close(DIR);

	#pushing each subfolder
	foreach (@Directory){
		if ($_ =~ /^\d*$/){$HashDirectory{$_}= $_;}
	}
	foreach my $key (sort {$a <=> $b} keys %HashDirectory) {
		#import to FastBit file
		my $execute = "ardea -d $finalpath/chickengenes -m \"
			chrom:char,
			geneid:char,
			genename:char,
			species:key,
			fpkmstatus:char,
			tissue:char,
			line:char,
			coverage:double,
			tpm:double,
			fpkm:double,
			fpkmlow:double,
			fpkmhigh:double,
			library:int,
			chromstart:int,
			chromstop:int\" -t $basepath/$key/$key\.txt";
		system($execute);
	} 
} else {
	my $execute = "ardea -d $finalpath/chickengenes -m \"
		chrom:char,
		geneid:char,
		genename:char,
		species:key,
		fpkmstatus:char,
		tissue:char,
		line:char,
		coverage:double,
		tpm:double,
		fpkm:double,
		fpkmlow:double,
		fpkmhigh:double,
		library:int,
		chromstart:int,
		chromstop:int\" -t $basepath/$library/$library\.txt";
	system($execute);
}
#GETTING ALL THE LIBRARIES FROM THE DATABASE.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
close STDOUT; close STDERR;
exit;

