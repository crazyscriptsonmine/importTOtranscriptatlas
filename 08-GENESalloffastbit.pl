#!/usr/bin/perl
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MOA 2017
# Extrct Genes from the Databas and save them into fastbit format
#This is only for chicken

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
my $basepath = "/home/modupe/public_html/TAFiles/GenesAtlas";
our ($chickengenes, $mousegenes, $alligatorgenes) = FBGENES();
`mkdir -p $basepath $chickengenes $mousegenes $alligatorgenes`;
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

my ($dbh, $sth); my %HashDirectory;
# CONNECT TO THE DATABASE
print "\n\n\tCONNECTING TO THE DATABASE\n\n";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#GETTING ALL THE LIBRARIES FROM THE DATABASE.
my ($libras, %libraries, @libraries);
$dbh = mysql();
$libras = "select a.library_id, a.species from bird_libraries a join genes_summary b on a.library_id = b.library_id where b.status = 'done' and b.nosql is null";
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
				c.chrom_no, c.gene_id, c.gene_short_name, a.species, c.fpkm_status, a.tissue, a.line,
				c.coverage, c.tpm, c.fpkm, c.fpkm_conf_low, c.fpkm_conf_high,
				a.library_id, c.chrom_start, c.chrom_stop
				from bird_libraries a join genes_fpkm c 
					on a.library_id = c.library_id
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
					$row2[5] = uc($row2[5]); $row2[6] = uc($row2[6]); #uppercase line and tissue
					if ($row2[7] =~ /NULL/){ $row2[7] = 0; } #change tpm from NULL to zero
					if ($list < 7) {
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
		} # end unless
		else {
			print "Already processed library $file\n";
		}
	} # end foreach
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	print "\n\nDONE with export from the database. \n\n";
	
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
	
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	$libraries{$species} =~ s/^,|,$//g;
	@libraries = split (",", $libraries{$species});
	
	if ($#libraries >= 0) {
		$statusfile = "$finalpath/myFBstatus";
		
		foreach my $file (@libraries) {
			#import to FastBit file
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
					chromstop:int\" -t $basepath/$file/$file\.txt";
				print $execute,"\n";
				
				`$execute 1>>$statusfile\.log 2>>$statusfile\.err`;
		}
		#changing status to done
		$dbh = mysql();
		my $syntax = "update genes_summary set nosql = 'done' where library_id in \($libraries{$species}\)";
		$sth = $dbh->prepare($syntax); $sth->execute() or die "$DBI::errstr Failed to update genes summary\n";
		print "\n\tFinished with nosql output \n";
		$sth->finish();
		$dbh->disconnect(); 
	}

} #finish working with the libraries of the given species.

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -

exit;

