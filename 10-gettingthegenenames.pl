#!/usr/bin/perl
use strict;
use DBI;
use File::Basename;
use Sort::Key::Natural qw(natsort);
use lib '/home/modupe/SCRIPTS/SUB';
use passw;
use routine;

my ($ibis, $syntax, $out, $outphp, %LineHash);
my $folder = "/home/modupe/public_html/FBTranscriptAtlas/nameofGENES";
our ($chickenpath, $mousepath, $alligatorpath) = FBPATHS();
our ($chickengene, $mousegene, $alligatorgene) = FBGENES();

my %allspecies = ("chickenvariant" => $chickenpath, "alligatorvariant" => $alligatorpath,"mousevariant" => $mousepath, "chickengene" => $chickengene, "alligatorgene" => $alligatorgene,"mousegene" => $mousegene,); #genomes

foreach my $towork (keys %allspecies) {
	$ibis = "/home/modupe/.bin/bin/ibis -d $allspecies{$towork} -q '";
	$syntax = "select genename,countdistinct(genename) where genename != \"NULL\"' -v -o temp.txt";
	$out = $folder."/".$towork.".txt";
	$outphp = $folder."/".$towork.".php";
	`$ibis$syntax`; #execute ibis

	open(IN, "<temp.txt");
	undef %LineHash;
	while (<IN>){
		chomp;
		my $linename = (split /,/)[0];
		$linename =~ s/^"|"$//g;
		$LineHash{$linename} = $linename;
	}
	close (IN);
	`rm -rf temp.txt`; 
	open (OUT, ">", $out);
	open (OUT2, ">", $outphp);
	print OUT2 '<?php
$genes = array (
//';
print OUT2 scalar keys %LineHash," genes \n";
	foreach (natsort keys %LineHash) {
		print OUT2 "\"",$_,"\",";
		print OUT $_,"\n";
	}
	print OUT2 '
);
echo json_encode($genes);
?>

';
`cp -f $outphp /home/modupe/public_html/atlas/names/`;
	close (OUT2); close (OUT);
}
