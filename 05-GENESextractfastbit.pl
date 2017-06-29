#!/usr/bin/perl
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - U S E R  V A R I A B L E S- - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
use strict;
use DBI;
use lib '/home/modupe/SCRIPTS/SUB';
use routine;
use passw;

# DATABASE ATTRIBUTES
my $basepath = "/home/modupe/public_html/GenesAtlas";
`mkdir -p $basepath`;
my $statusfile = "$basepath/mystatus.txt";
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
my ($dbh, $sth); my %HashDirectory;
# CONNECT TO THE DATABASE
print "\n\n\tCONNECTING TO THE DATABASE\n\n";
$dbh = mysql();
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CONNECT
#GETTING ALL THE LIBRARIES FROM THE DATABASE.
my $libraries;
#for Chicken
my $libras = "select a.library_id from bird_libraries a join genes_summary b on a.library_id = b.library_id where a.species = \"gallus\" and b.status = 'done'";
$sth = $dbh->prepare($libras); $sth->execute or die "SQL Error: $DBI::errstr\n";

while ( my $row = $sth->fetchrow_array() ) {
	$libraries .= $row.",";
}
$libraries = substr($libraries,0,-1);
my @libraries = split (",", $libraries);
#print @libraries; die;
#PARSING EACH FILE OUT OF THE DATABASE. 
foreach my $file (@libraries){
	opendir (DIR, $basepath) or die "Folder \"$basepath\" doesn't exist\n"; 
	my @Directory = readdir(DIR); close(DIR);
	#pushing each subfolder
	foreach (@Directory){
		if ($_ =~ /^\d*$/){
		$HashDirectory{$_}= $_;}
	}
	unless (exists $HashDirectory{$file}){
		print "Working on $file...\n";
		system("mkdir -p $basepath/$file");
		my $syntax= "select
			c.chrom_no, c.gene_id, c.gene_short_name, a.species, c.fpkm_status, a.tissue, a.line,
			c.coverage, c.tpm, c.fpkm, c.fpkm_conf_low, c.fpkm_conf_high,
			a.library_id, c.chrom_start, c.chrom_stop
			from bird_libraries a join genes_fpkm c 
				on a.library_id = c.library_id
			where a.species = \"gallus\" and a.library_id = $file";

 
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
				if ($row2[7] eq "NULL"){ $row2[7] = 0; } #change tpm from NULL to zero
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
		open(STATUS, ">>$statusfile");
		print STATUS "$file\t$i\n";
		close (OUT);
	} # end unless
	else {
		print "Already processed library $file\n";
	}
} # end foreach
#changing status to done
if ($#libraries >= 1) {
	my $syntax = "update genes_summary set nosql = 'done' where library_id in \($libraries\)";
	$sth = $dbh->prepare($syntax); $sth->execute() or die "$DBI::errstr Failed to update genes summary\n";
	print "\n\tFinished with nosql output \n";
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
exit;

