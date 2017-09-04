#!/usr/bin/perl
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MOA 2017
# Redo and save them into fastbit format just in case the increments was buggy

# - - - - - - - - - - - - - - U S E R  V A R I A B L E S- - - - - - - - - - - - - - -

use strict;
use Getopt::Long;
use Pod::Usage;
use DBI;
use lib '/home/modupe/SCRIPTS/SUB';
use routine;
use passw;

# DATABASE ATTRIBUTES
my ($statusfile, $finalpath);
my $gbasepath = "/home/modupe/public_html/TAFiles/GenesAtlas";
my $vbasepath = "/home/modupe/public_html/TAFiles/VariantAtlas";
our ($chickengenes, $mousegenes, $alligatorgenes) = FBGENES();
our ($chickenpath, $mousepath, $alligatorpath) = FBPATHS();
`mkdir -p $vbasepath $gbasepath $chickengenes $mousegenes $alligatorgenes`;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

my ($dbh, $sth);
# CONNECT TO THE DATABASE
print "\n\n\tCONNECTING TO THE DATABASE\n\n";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#GETTING ALL THE LIBRARIES FROM THE DATABASE.
my ($glibras, %glibraries, @glibraries);
$dbh = mysql();
$glibras = "select a.library_id, a.species from bird_libraries a join genes_summary b on a.library_id = b.library_id where b.nosql = 'done'";
$sth = $dbh->prepare($glibras); $sth->execute or die "SQL Error: $DBI::errstr\n";

while ( my ($row, $species) = $sth->fetchrow_array() ) {
	$glibraries{$species} = $glibraries{$species}.",".$row;
}
$sth->finish();
$dbh->disconnect(); 
foreach my $species (keys %glibraries){ #getting the libraries for the different species
	$glibraries{$species} =~ s/^,|,$//g;
	$statusfile = "$gbasepath/mystatus.txt";
	@glibraries = split (",", $glibraries{$species});

	print "Working on Organism : $species \n\n";

# - - - - - - - - - - - - - - - -FASTBIT IMPORT - - - - - - - - - - - - - - - - - - -
	
	if ($species =~ /alligator/){ #making sure the final path is accurately specified
		$finalpath = $alligatorgenes;
	} elsif ($species =~ /mus/){
		$finalpath = $mousegenes;
	} elsif ($species =~ /gallus/){
		$finalpath = $chickengenes;
	} else {
		exit "$species is incorrect species name\n";
	}

	if ($#glibraries >= 0) {
		$statusfile = "$finalpath/myFBstatus";
		
		foreach my $file (@glibraries) {
			#import to FastBit file
			`rm -rf $finalpath; mkdir -p $finalpath;`;
			my $execute = "/home/modupe/.bin/bin/ardea -d $finalpath -m \"
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
					chromstop:int\" -t $gbasepath/$file/$file\.txt";
				print $execute,"\n";
				
				`$execute 1>>$statusfile\.log 2>>$statusfile\.err`;
		}
	}

} #finish working with the libraries of the given species.

# - - - - - - - - - - - - - - - - - FINISHED WITH GENES - - - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -

# - - - - - - - - - - - - WORKING WITH VARIANTS- - - - - - - - - - - - - -
#GETTING ALL THE LIBRARIES FROM THE DATABASE.
my ($vlibras, %vlibraries, @vlibraries);
$dbh = mysql();
$vlibras = "select a.library_id, a.species from vw_libraryinfo a join variants_summary b on a.library_id = b.library_id where b.nosql = 'done'";
$sth = $dbh->prepare($vlibras); $sth->execute or die "SQL Error: $DBI::errstr\n";

while ( my ($row, $species) = $sth->fetchrow_array() ) {
	$vlibraries{$species} = $vlibraries{$species}.",".$row;
}
$sth->finish();
$dbh->disconnect(); 
foreach my $species (keys %vlibraries){ #getting the libraries for the different species
	$vlibraries{$species} =~ s/^,|,$//g;
	$statusfile = "$vbasepath/mystatus.txt";
	@vlibraries = split (",", $vlibraries{$species});

#print @libraries; die;
	print "Working on Organism : $species \n\n";

# - - - - - - - - - - - - - - - -FASTBIT IMPORT - - - - - - - - - - - - - - - - - - -
	
	if ($species =~ /alligator/){ #making sure the final path is accurately specified
		$finalpath = $alligatorpath;
	} elsif ($species =~ /mus/){
		$finalpath = $mousepath;
	} elsif ($species =~ /gallus/){
		$finalpath = $chickenpath;
	} else {
		exit "$species is incorrect species name\n";
	}
	
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	
	$vlibraries{$species} =~ s/^,|,$//g;
	@vlibraries = split (",", $vlibraries{$species});
	
	if ($#vlibraries >= 0) {
		$statusfile = "$finalpath/myFBstatus";

		foreach my $file (@vlibraries) {
			#import to FastBit file
			`rm -rf $finalpath; mkdir -p $finalpath;`;
			my $execute = "/home/modupe/.bin/bin/ardea -d $finalpath -m \"
					class:char,
					zygosity:char,
					dbsnp:char,
					consequence:char,
					geneid:char,
					genename:char,
					transcript:char,
					feature:char,
					genetype:char,
					ref:char,
					alt:char,
					line:char,
					tissue:char,
					chrom:char,
					aachange:char,
					codon:char,
					species:key,
					notes:text, 
					quality:double,
					library:int,
					position:int,
					proteinposition:int\" -t $vbasepath/$file/$file\.txt";
				print $execute,"\n";
				
				`$execute 1>>$statusfile\.log 2>>$statusfile\.err`;
		}
	}

} #finish working with the libraries of the given species.

# - - - - - - - - - - - - - - - - - FINISHED WITH VARIANTS- - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -

exit;

