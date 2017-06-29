#!/usr/bin/perl
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -

# MANUAL for extracting variant stuff from the database
# and save to fastbit format
#MOA 2017

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
my $basepath = "/home/modupe/public_html/TAFiles/VariantAtlas";
our ($chickenpath, $mousepath, $alligatorpath) = FBPATHS();
`mkdir -p $basepath $chickenpath $mousepath $alligatorpath`;
# - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - - - - -
my ($dbh, $sth); my %HashDirectory;
# CONNECT TO THE DATABASE
print "\n\n\tCONNECTING TO THE DATABASE\n\n";

# - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - - - - - - -
#GETTING ALL THE LIBRARIES FROM THE DATABASE.
my ($libras, %libraries, @libraries);
$dbh = mysql();
$libras = "select a.library_id, a.species from vw_libraryinfo a join variants_summary b on a.library_id = b.library_id where b.status = 'done' and b.nosql is null";
$sth = $dbh->prepare($libras); $sth->execute or die "SQL Error: $DBI::errstr\n";

while ( my ($row, $species) = $sth->fetchrow_array() ) {
	$libraries{$species} = $libraries{$species}.",".$row;
}
$sth->finish();
$dbh->disconnect(); 
foreach my $species (keys %libraries){ #getting the libraries for the different species
	$libraries{$species} =~ s/^,|,$//g;
	$statusfile = "$basepath/mystatus.txt";
	@libraries = split (",", $libraries{$species});

#print @libraries; die;
	print "Working on Organism : $species \n\n";
#PARSING EACH FILE OUT OF THE DATABASE. 
	foreach my $file (@libraries){
		opendir (DIR, $basepath) or die "Folder \"$basepath\" doesn't exist\n"; 
		my @Directory = readdir(DIR); close(DIR);
	
		#pushing each subfolder
		foreach (@Directory){ if ($_ =~ /^\d*$/){ $HashDirectory{$_}= $_; } }
	
		unless (exists $HashDirectory{$file}){
			print "Working on Library $file...\n";
			system("mkdir -p $basepath/$file");
			$dbh = mysql();
			my $syntax= "select
				b.variant_class,b.zygosity,b.existing_variant,c.consequence,c.gene_id,
				c.gene_name,c.transcript,c.feature,c.gene_type,b.ref_allele,b.alt_allele,a.line,
				a.tissue,b.chrom,c.aminoacid_change,c.codon_change,a.species,
				a.notes,b.quality,a.library_id,b.position,
				c.protein_position
				from vw_libraryinfo a join variants_result b 
					on a.library_id = b.library_id 
				join variants_annotation c 
					on b.library_id = c.library_id and b.chrom = c.chrom and b.position = c.position
				where a.library_id = $file";
	 
			$sth = $dbh->prepare($syntax);			
			$sth->execute() or die "SQL Error: $DBI::errstr\n";

			open(OUT,">$basepath/$file/$file\.txt") or die "wtf";
		
			#TABLE FORMAT
			my $i = 0;
			while ( my @row2 = $sth->fetchrow_array() ) {
			$i++;
				foreach my $list (0..$#row2-1){
					if ((length($row2[$list]) < 1) || ($row2[$list] =~ /^\-$/) ){
						$row2[$list] = "NULL";
					}
					$row2[$list] =~ s/^'|'$//g;
					$row2[11] = uc($row2[11]); $row2[12] = uc($row2[12]); #uppercase line and tissue
				
					if ($list < 18) {
						print OUT "\'$row2[$list]\',";
					}
					else {
						print OUT "$row2[$list],";
					}
				}
				if ((length($row2[$#row2]) < 1) || ($row2[$#row2] =~ /^\-$/) ){
					$row2[$#row2] = "NULL";
				}
				print OUT "$row2[$#row2]\n";
			}
			$sth->finish();
			$dbh->disconnect(); 
			open(STATUS, ">>$statusfile");
			print STATUS "$file\t$i\t",scalar(localtime),"\n";
			close (OUT);
		}	 # end unless
		else {
			print "Already processed library $file\n";
		}
	} # end foreach
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

print "\n\nDONE with export from the database. \n\n";

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
	
	$libraries{$species} =~ s/^,|,$//g;
	@libraries = split (",", $libraries{$species});
	
	if ($#libraries >= 0) {
		$statusfile = "$finalpath/myFBstatus";

		foreach my $file (@libraries) {
			#import to FastBit file
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
					proteinposition:int\" -t $basepath/$file/$file\.txt";
				print $execute,"\n";
				
				`$execute 1>>$statusfile\.log 2>>$statusfile\.err`;
		}
		#changing status to done
		$dbh = mysql();
		my $syntax = "update variants_summary set nosql = 'done' where library_id in \($libraries{$species}\)";
		$sth = $dbh->prepare($syntax); $sth->execute() or die "$DBI::errstr Failed to update variant summary\n";
		print "\n\tFinished with nosql output \n";
		$sth->finish();
		$dbh->disconnect();
		
	}

} #finish working with the libraries of the given species.

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -

exit;

