#!/usr/bin/perl
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MOA 2017
# Extract Librarys and the genes expression details from fastbit to the atlas LIBFPKM dump.

# - - - - - - - - - - - - - - U S E R  V A R I A B L E S- - - - - - - - - - - - - - -

use strict;
use Getopt::Long;
use Pod::Usage;
use DBI;
use lib '/home/modupe/SCRIPTS/SUB';
use routine;
use passw;

# ATTRIBUTES
my ($statusfile, $libras, $libraries, $finalfile, $overlapfile, $finaloutput);
my $basepath = "/home/modupe/public_html/atlas/LIBFPKMDUMP";
our ($chickengenes, $mousegenes, $alligatorgenes) = FBGENES();
`mkdir -p $basepath $chickengenes $mousegenes $alligatorgenes`;

my %Libraries =  ("Gallus" => $chickengenes, "Alligator" => $alligatorgenes, "Mouse" => $mousegenes, );

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
foreach my $index (keys %Libraries) {
	`/home/modupe/.bin/bin/ibis -d $Libraries{$index} -q 'select library, countdistinct(library)' -o temp.txt`;
	my $total = `cat temp.txt | wc -l`;
	my $old = `head -n 1 $basepath/$index.no`;
	chomp ($total, $old); 
	if ($total != $old) { #if the old file doesn't match the new total
		open (IN, "<temp.txt");
		while (<IN>) {
			chomp;
			$libras = (split(/, /))[0];
			$libraries = $libraries.",".$libras;
		} close (IN);
		`rm -rf temp.txt`;
		$libraries =~ s/^,|,$//g;
		$statusfile = "$basepath/mystatus.txt";
		$finalfile = "$basepath/$index"."All";
		$overlapfile = "$basepath/$index"."Overlap";
		
		$finaloutput = "Processing $index,". scalar(localtime);
		`echo $finaloutput >> $statusfile`;
		`perl /home/modupe/public_html/atlas/SQLscripts/outputgenequery.pl -1 $libraries -2 $finalfile`; #recompute the new files
                $finaloutput = "Finished ALL Genes for  $index,". scalar(localtime);
                `echo $finaloutput >> $statusfile`;
                $finaloutput = "Moving on to overlapping genes,". scalar(localtime);
                `echo $finaloutput >> $statusfile`;
		`perl /home/modupe/public_html/atlas/SQLscripts/outputcommagenes.pl -1 $libraries -2 $overlapfile`; #recompute the overlap files.
		`gzip -f $finalfile.txt $overlapfile.txt`;
		`rm -rf $basepath/*v1 $basepath/*v2 $basepath/*v3`;
		`echo $total > $basepath/$index.no`;
		$finaloutput = "Completed $index,". scalar(localtime);
		`echo $finaloutput >> $statusfile`;
	}
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -

exit;

