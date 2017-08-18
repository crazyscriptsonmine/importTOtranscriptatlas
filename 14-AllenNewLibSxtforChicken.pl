# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - U S E R  V A R I A B L E S- - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# CODE FOR
# Code for modified to fit Allen's new library structure for chicken 
use strict;
use Data::Dumper;
use File::Basename;
use DBI;
use Getopt::Long;
use Time::localtime;
use Pod::Usage;
use Time::Piece;
use File::stat;
use threads;
use Thread::Queue;
use DateTime;
use POSIX qw( ceil );
use lib '/home/modupe/SCRIPTS/SUB';
use routine;
use passw;

# #CREATING LOG FILES
my $std_out = '/home/modupe/.LOG/RavTAD-'.`date +%m-%d-%y_%T`; chomp $std_out; $std_out = $std_out.'.log';
my $std_err = '/home/modupe/.LOG/RavTAD-'.`date +%m-%d-%y_%T`; chomp $std_err; $std_err = $std_err.'.err';
my $jobid = "RavenTAD-".`date +%m-%d-%y_%T`;
my $progressnote = "/home/modupe/.LOG/progressnote".`date +%m-%d-%y_%T`; chomp $progressnote; $progressnote = $progressnote.'.txt'; 

open(STDOUT, '>', "$std_out") or die "Log file doesn't exist";
open(STDERR, '>', "$std_err") or die "Error file doesn't exist";
 
#ARGUMENTS
my($help,$manual,$deletenotdone,$in1);
GetOptions (	
          "delete" 			=> 	\$deletenotdone,
					"h|help"  		=>  \$help,
          "man|manual"	=>  \$manual );

# VALIDATE ARGS
pod2usage( -verbose => 2 )  if ($manual);
pod2usage( -verbose => 1 )  if ($help);
@ARGV<=1 or pod2usage("Syntax error");
#file path for input THIS SHOULD BE CONSTANT
$in1 = $ARGV[0]; #files to transfer

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

my $mystderror = "Contact Modupe Adetunji amodupe\@udel.edu\n";

# RESULTS_HASH
my (%Hashresults, %Birdresults, %Nullresults);

# DATABASE VARIABLES
my ($dbh, $sth, $syntax, $row, @row);

#DIRECTORY
my (@parse, @NewDirectory, @numberz);

# TABLE VARIABLES
my ($accepted, $samfile, $alignfile, $isoformsfile, $genesfile, $deletionsfile, $insertionsfile, $transcriptsgtf, $junctionsfile, $run_log, $htseqcount,$variantfile, $vepfile, $annovarfile);
my ($parsedinput, $len, $alignrate, $varianttool, $fblocation);
my ($lib_id, $total, $mapped, $unmapped, $deletions, $insertions, $junctions, $genes, $isoforms, $prep, $date); #transcripts_summary TABLE
my ($track, $class, $ref_id, $gene, $gene_name, $tss, $locus, $chrom_no, $chrom_start, $chrom_stop, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat); # GENES_FPKM &

#FILE VERSIONS
my ($gatk_version,$picard_version, $vep_version);
my ($found, $diffexpress);
#VARIANTS FOLDER & HASHES
my $Mfolder;
my (%VCFhash, %extra, %DBSNP, %VEPhash, %ExtraWork, %AMISG, %AMIST);
my (@allgeninfo, $mappingtool, $refgenome, $refgenomename, %ALL);
my ($stranded, $sequences, $annotation, $annotationfile, $annfileversion);
my (@foldercontent, @VAR, @threads, $queue);
my (%ARFPKM,%CHFPKM, %BEFPKM, %CFPKM, %DFPKM, %DHFPKM, %DLFPKM, %cfpkm, %dfpkm, %dhfpkm, %dlfpkm, %TPM,%tpm)= ();
my (%HASHDBVARIANT, %HASHDBVEP, %HASHRESULT, %HASHNEW) ;

#PARSABLE GENOMES FOR ANALYSIS
my $GENOMES="/home/modupe/.GENOMES/";
my $STORAGEPATH = "/home/modupe/CHICKENSNPS"; #variant files are stored in this directory
my %parsablegenomes = ("chicken" => 1, "alligator" => 2,"mouse" => 3, ); #genomes that work.
my %VEPparse = ("chicken" => 1,"mouse" => 2, ); #for VEP

#INDEPENDENT PROGRAMS TO RUN
my $PICARDDIR="/home/modupe/.software/picard-tools-1.136/picard.jar";
my $GATKDIR="/home/modupe/.software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar";
my $VEP="/home/modupe/.software/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl";
my $SNPdat="/home/modupe/.software/SNPdat_package_v1.0.5/SNPdat_v1.0.5.pl";
my $email = 'amodupe@udel.edu';

#OPENING FOLDER
opendir(DIR,$in1) or die "Folder \"$in1\" main doesn't exist\n"; 
my @Directory = readdir(DIR);
close(DIR);
#pushing each subfolder

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# -----------------------------------
#CREATING EMAIL NOTIFICATION
NOTIFICATION("Starting Job");
# -----------------------------------

# CONNECT TO THE DATABASE
print "\n\n\tCONNECTING TO THE DATABASE :".`date`."\n\n";
$dbh = mysql();
if ($deletenotdone) {DELETENOTDONE();}

if ($in1 =~ /\w.*_(\d+)_/){
	 my $libraryidnumber = $1;
  CHECKING();
  unless (exists $Hashresults{$libraryidnumber}){
    if (exists $Birdresults{$libraryidnumber}){
      $parsedinput = $in1;
			@foldercontent = split("\n", `find $parsedinput`); #get details of the folder
			foreach (grep /\.gtf/, @foldercontent) { unless (`head -n 3 $_ | wc -l` <= 0 && $_ =~ /skipped/) { $transcriptsgtf = $_; } }
			$accepted = (grep /sorted_bam.bam/, @foldercontent)[0];
			$alignfile = (grep /summary.txt/, @foldercontent)[0];
			$genesfile = (grep /genes.fpkm/, @foldercontent)[0];
			$isoformsfile = (grep /isoforms.fpkm/, @foldercontent)[0];
			$deletionsfile = (grep /deletions.bed/, @foldercontent)[0];
			$insertionsfile = (grep /insertions.bed/, @foldercontent)[0];
			$junctionsfile = (grep /junctions.bed/, @foldercontent)[0];
			$run_log = (grep /logs\/run.log/, @foldercontent)[0];
			$samfile = (grep /.sam$/, @foldercontent)[0];
			$variantfile = (grep /.vcf$/, @foldercontent)[0]; 
			$vepfile = (grep /vep.txt$/i, @foldercontent)[0];
			$annovarfile = (grep /anno.txt$/, @foldercontent)[0];
			$htseqcount = (grep /.counts$/, @foldercontent)[0];
			LOGFILE($libraryidnumber);
			my $verdict = PARSING($libraryidnumber,$parsedinput);
				
			#progress report
			if ($verdict == 1) {
				open (NOTE, ">>$progressnote");
				print NOTE "Subject: Update notes : $jobid\n\nCompleted library\t$1\n";
				system "sendmail $email < $progressnote"; close NOTE;
			} #end if
		} else {
			print "\nSkipping \"library_$libraryidnumber\" in \"$in1\" folder because it isn't in birdbase\n$mystderror\n";
		} #end if
	}else {print "\nLibrary => $libraryidnumber exists in the database\n";} #end unless
} #end if	 
#SUMMARYstmts(); 
system "rm -rf $progressnote";
# DISCONNECT FROM THE DATABASE
print "\n\tDISCONNECTING FROM THE DATABASE\n\n";
$dbh->disconnect();

# -----------------------------------
#send finish notification
NOTIFICATION("Job completed");
# -----------------------------------
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - -S U B R O U T I N E S- - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub DELETENOTDONE { #deleting the incomplete libraries (only when requested using option delete [-delete])
  print "\n\tDELETING NOT DONE\n";
  #CHECKING TO MAKE SURE NOT "done" FILES ARE REMOVED
  $syntax = "select library_id from transcripts_summary where status is NULL";
  $dbh->disconnect();
	$dbh = mysql();
  $sth = $dbh->prepare($syntax);
  $sth->execute or die "SQL Error: $DBI::errstr\n";
  my $incompletes = undef; my $count=0; my @columntoremove;
  while ($row = $sth->fetchrow_array() ) {
    $count++;
    $incompletes .= $row.",";
  }
  if ($count >= 1){
    $incompletes = substr($incompletes,0,-1);
    print "\tDeleted rows $incompletes\n";
    #DELETE FROM variants_annotation
    $sth = $dbh->prepare("delete from variants_annotation where library_id in ( $incompletes )"); $sth->execute();
    #DELETE FROM variants_result
    $syntax = "delete from variants_result where library_id in \( $incompletes \)";
    $sth = $dbh->prepare($syntax); $sth->execute();
    #DELETE FROM variants_summary
    $sth = $dbh->prepare("delete from variants_summary where library_id in ( $incompletes )"); $sth->execute();
    #DELETE FROM genes_fpkm
    $sth = $dbh->prepare("delete from genes_fpkm where library_id in ( $incompletes )"); $sth->execute();
    #DELETE FROM isoforms_fpkm
    $sth = $dbh->prepare("delete from isoforms_fpkm where library_id in ( $incompletes )"); $sth->execute();
    #DELETE FROM htseq
    $sth = $dbh->prepare("delete from htseq where library_id in ( $incompletes )"); $sth->execute();
    #DELETE FROM frnak_metadata
    $sth = $dbh->prepare("delete from frnak_metadata where library_id in ( $incompletes )"); $sth->execute();
    #DELETE FROM genes_summary
    $sth = $dbh->prepare("delete from genes_summary where library_id in ( $incompletes )"); $sth->execute();
		#DELETE FROM transcripts_summary
    $sth = $dbh->prepare("delete from transcripts_summary where library_id in ( $incompletes )"); $sth->execute();
  }
}
sub CHECKING { #subroutine for checking the libraries in the database and those processed
  #CHECKING THE LIBRARIES ALREADY IN THE DATABASE
  $syntax = "select library_id from transcripts_summary where status is not null";
  $sth = $dbh->prepare($syntax);
  $sth->execute or die "SQL Error: $DBI::errstr\n";
  my $number = 0;
  while ($row = $sth->fetchrow_array() ) {
    $Hashresults{$row} = $number; $number++;
  }
	$syntax = "select library_id from transcripts_summary where status is null";
  $sth = $dbh->prepare($syntax);
  $sth->execute or die "SQL Error: $DBI::errstr\n";
  $number = 0;
  while ($row = $sth->fetchrow_array() ) {
    $Nullresults{$row} = $number; $number++;
  }
  $syntax = "select library_id,date from bird_libraries";
  $sth = $dbh->prepare($syntax);
  $sth->execute or die "SQL Error: $DBI::errstr\n";
  $number = 0;
  while (my ($row1, $row2) = $sth->fetchrow_array() ) {
    $Birdresults{$row1} = $row2;
  }
}
sub LOGFILE { #subroutine for getting metadata
	#print "working on abc\tcc$samfile\tdd$run_log\n"; die;
	if ($samfile) {
		@allgeninfo = split('\s',`grep -m 1 "\@PG" $samfile | head -1`);
		foreach my $no (0..$#allgeninfo){ if ($allgeninfo[$no] =~ /^CL/) { ++$no; @allgeninfo = split('\s',`grep -m 1 "\@PG" $samfile | head -1`,$no);} }
		#getting metadata info
		if ($#allgeninfo > 1) {
			my $tool = (grep /^ID\:/, @allgeninfo)[0];
			my $ver_no = (grep /^VN\:/, @allgeninfo)[0];
			my $command = (grep /^CL\:/, @allgeninfo)[0];
			$mappingtool = ((split(':',((grep /^ID\:/, @allgeninfo)[0])))[-1])." v".((split(':',((grep /^VN\:/, @allgeninfo)[0])))[-1]);
			if ($allgeninfo[1] =~ /ID\:(\S*)$/){ $mappingtool = $1." v".(split(':',$allgeninfo[3]))[-1]; } #mapping tool name and version
			if ($mappingtool =~ /hisat/i) { 
				$command =~ /\-x\s(\w+).*\s/;
				$refgenome = $1;
				$refgenomename = (split('\/', $refgenome))[-1]; #reference genome name
				if ($command =~ /-1/){
					$command =~ /\-1\s(\S+)\s-2\s(\S+)"$/;
					
					unless ($1 =~ /pipe/) { #making sure .pipe isn't the name of the sequence
						my @nseq = split(",",$1); my @pseq = split(",",$2);
						foreach (@nseq){ $sequences .= ( (split('\/', $_))[-1] ).",";}
						foreach (@pseq){ $sequences .= ( (split('\/', $_))[-1] ).",";}
						chop $sequences;
					} else {
						$in1 =~ /\w.*_(\d+_.*R)\d.*/;
						my $path = `readlink -f $in1`; chop $path; my @path = split('\/', $path); undef $path; foreach (0..$#path-1) { $path .= $path[$_]."/"; }
						my $locale = `locate $path$1 | head -n 2`; chop $locale;
						my @seqall = split("\n", $locale, 2);
						my @nseq = split(",",$seqall[0]); my @pseq = split(",", $seqall[1]);
						foreach (@nseq){ $sequences .= ( (split('\/', $_))[-1] ).",";}
						foreach (@pseq){ $sequences .= ( (split('\/', $_))[-1] ).",";}
						
						unless ($alignfile){
							$Mfolder = "$STORAGEPATH/library_$_[0]"; `mkdir -p $Mfolder`; print "\tmade $Mfolder\n"; #decided to keep all the hisat details.
							$alignfile = "$Mfolder/align_summary.txt";
							unless (-e $alignfile) {
                        print "Job: Reference Mapping using HISAT2";
								`hisat2 -x $path/chicken/chicken -1 $seqall[0] -2 $seqall[1] -S $Mfolder/library_$_[0].sam 2>$alignfile`;
                        print ". . . Done\n";
							}
						}
					}
				}
				elsif ($command =~ /-U/){
					$command =~ /\-U\s(\S+)"$/;
					my @nseq = split(",",$1);
					foreach (@nseq){ $sequences .= ( (split('\/', $_))[-1] ).",";}
					chop $sequences;
				} #end if toggle for sequences
				$stranded = undef;
				$annotation = undef;
			} elsif ($mappingtool =~ /tophat/i) { 
				undef %ALL;
				my ($no, $number) = (0,1);
				@allgeninfo = split('\s',$command);
				while ($number <= $#allgeninfo){
				  unless ($allgeninfo[$number] =~ /-no-coverage-search/){
				    if ($allgeninfo[$number] =~ /^\-/){
				      my $old = $number++;
				      $ALL{$allgeninfo[$old]} = $allgeninfo[$number];
				    } else {
				      unless (exists $ALL{$no}){
				        $ALL{$no} = $allgeninfo[$number];
				        $no++;
				      }
				    }
					}
					$number++;
				}
				unless ((exists $ALL{"-G"}) || (exists $ALL{"--GTF"})) {
				  $annotation = undef;
				} else {
				  if (exists $ALL{"-G"}){ $annotationfile = $ALL{"-G"} ; } else { $annotationfile = $ALL{"--GTF"};}
				  $annotation = uc ( (split('\.',((split("\/", $annotationfile))[-1])))[-1] ); #(annotation file)
				}
				unless (exists $ALL{"--library-type"}) { $stranded = undef; } else { $stranded = $ALL{"--library-type"}; }
			
				$refgenome = $ALL{0}; my $seq = $ALL{1}; my $otherseq = $ALL{2};
				$refgenomename = (split('\/', $ALL{0}))[-1];
				unless(length($otherseq)<1){ #sequences
				  $sequences = ( ( split('\/', $seq) ) [-1]).",". ( ( split('\/', $otherseq) ) [-1]);
				} else {
				  $sequences = ( ( split('\/', $seq) ) [-1]);
				} #end if seq
			}
		} else {
			$annotation = undef;
			$stranded = undef; $sequences = undef;
		}
	} # end if samfile
	elsif ($run_log){
		@allgeninfo = split('\s',`head -n 1 $run_log`);
		#getting metadata info
		if ($#allgeninfo > 1){
			if ($allgeninfo[0] =~ /tophat/){ $mappingtool = "TopHat";}
			undef %ALL;
			my ($no, $number) = (0,1);
      while ($number <= $#allgeninfo){
        unless ($allgeninfo[$number] =~ /-no-coverage-search/){
          if ($allgeninfo[$number] =~ /^\-/){
            my $old = $number++;
            $ALL{$allgeninfo[$old]} = $allgeninfo[$number];
          } else {
            unless (exists $ALL{$no}){
              $ALL{$no} = $allgeninfo[$number];
              $no++;
            }
          }
        }
        $number++;
      }
      unless ((exists $ALL{"-G"}) || (exists $ALL{"--GTF"})) {
        $annotation = undef;
      } else {
        if (exists $ALL{"-G"}){ $annotationfile = $ALL{"-G"} ; } else { $annotationfile = $ALL{"--GTF"};}
        $annotation = uc ( (split('\.',((split("\/", $annotationfile))[-1])))[-1] ); #(annotation file)
      }
      unless (exists $ALL{"--library-type"}) { $stranded = undef; } else { $stranded = $ALL{"--library-type"}; }
			
      $refgenome = $ALL{0}; my $seq = $ALL{1}; my $otherseq = $ALL{2};
			$refgenomename = (split('\/', $ALL{0}))[-1];
      unless(length($otherseq)<1){ #sequences
        $sequences = ( ( split('\/', $seq) ) [-1]).",". ( ( split('\/', $otherseq) ) [-1]);
      } else {
        $sequences = ( ( split('\/', $seq) ) [-1]);
      } #end if seq
		}
	} else {
		print "ERROR: SAM file or TopHat LOG file is requiredt\n";
	}
}

sub GENES_FPKM { #subroutine for getting gene information
	#INSERT INTO DATABASE: #genes_summary table
	$sth = $dbh->prepare("select library_id from genes_summary where library_id = '$_[0]'"); $sth->execute(); $found = $sth->fetch();
	unless ($found) {
		print "NOTICE:\t $_[0] inserting $_[0] genes summary details into the database ...";
		$sth = $dbh->prepare("insert into genes_summary (library_id,date) values (?,?)");
		$sth ->execute($_[0], $date) or die "\nERROR:\t Complication in genes_summary table,\n";
		print " Done\n";
	} else {
		print "NOTICE:\t $_[0] Metadata already in genes_summary table... Moving on \n";
	}
	my ($genecount, $isoformcount) = (0,0);
	$sth = $dbh->prepare("select status from genes_summary where library_id = '$_[0]' and status ='done'"); $sth->execute(); $found = $sth->fetch();
	unless ($found) {
		$genecount = $dbh->selectrow_array("select count(*) from genes_fpkm where library_id = '$_[0]'");
		if ($genesfile){ #working with genes.fpkm_tracking file
			#cufflinks expression tool name
			$diffexpress = "Cufflinks";
			$genes = `cat $genesfile | wc -l`; if ($genes >=2){ $genes--;} else {$genes = 0;} #count the number of genes
			$sth = $dbh->prepare("update genes_summary set genes = $genes, diffexpresstool = '$diffexpress' where library_id= '$_[0]'"); $sth ->execute(); #updating genes_summary table.
			unless ($genes == $genecount) {
				unless ($genecount == 0 ) {
					print "NOTICE:\t Removed incomplete records for $_[0] in genes_fpkm table\n";
		      $sth = $dbh->prepare("delete from genes_fpkm where library_id = '$_[0]'"); $sth->execute();
				}
				print "NOTICE:\t Importing $diffexpress expression information for $_[0] to genes_fpkm table ...";
				#import into FPKM table;
				open(FPKM, "<", $genesfile) or die "\nERROR:\t Can not open file $genesfile\n";

				##change to thread
				my @fpkmdetails = <FPKM>; close(FPKM);
				shift @fpkmdetails;
				undef @VAR; undef @threads;
				push @VAR, [ splice @fpkmdetails, 0, 200 ] while @fpkmdetails; #sub the files to multiple subs

				$queue = new Thread::Queue();
				my $builder=threads->create(\&main); #create thread for each subarray into a thread
				push @threads, threads->create(\&cuffprocessor) for 1..5; #execute 10 threads
				$builder->join; #join threads
				foreach (@threads){$_->join;}

				print " Done\n";
			} else {
				print "NOTICE:\t $_[0] already in genes_fpkm table... Moving on \n";
			}
			if ($isoformsfile) {
				$isoforms = `cat $isoformsfile | wc -l`; if ($isoforms >=2){ $isoforms--;} else {$isoforms = 0;} #count the number of isoforms in file
				$sth = $dbh->prepare("update genes_summary set isoforms = $isoforms where library_id= '$_[0]'"); $sth ->execute(); #updating genes_summary table.
				$isoformcount = $dbh->selectrow_array("select count(*) from isoforms_fpkm where library_id = '$_[0]'");
				unless ($isoforms == $isoformcount) {
					unless ($isoformcount == 0 ) {
						print "NOTICE:\t Removed incomplete records for $_[0] in isoforms_fpkm table\n";
						$sth = $dbh->prepare("delete from isoforms_fpkm where library_id = $_[0]"); $sth->execute();
					}
					print "NOTICE:\t Importing $diffexpress expression information for $_[0] to isoforms_fpkm table ...";
					#import into ISOFORMSFPKM table;
					open(FPKM, "<", $isoformsfile) or die "\nERROR:\t Can not open file $isoformsfile\n";
					
					##change to thread
					my @fpkmdetails = <FPKM>; close(FPKM);
					shift @fpkmdetails;
					undef @VAR; undef @threads;
					push @VAR, [ splice @fpkmdetails, 0, 200 ] while @fpkmdetails; #sub the files to multiple subs
	
					$queue = new Thread::Queue();
					my $builder=threads->create(\&main); #create thread for each subarray into a thread
					push @threads, threads->create(\&isoprocessor) for 1..5; #execute 5 threads
					$builder->join; #join threads
					foreach (@threads){$_->join;}
					###
				
					print " Done\n";
				} else {
					print "NOTICE:\t $_[0] already in isoforms_fpkm table... Moving on \n";
				}
			}
			#set genes_summary to Done
			$sth = $dbh->prepare("update genes_summary set status = 'done' where library_id = '$_[0]'");
			$sth ->execute() or die "\nERROR:\t Complication in genes_summary table, consult documentation\n";
			
		} elsif ($transcriptsgtf){ #working with gtf transcripts file
			#differential expression tool names
			if (`head -n 1 $transcriptsgtf` =~ /cufflinks/i) { #working with cufflinks transcripts.gtf file
				$diffexpress = "Cufflinks";
				open(FPKM, "<", $transcriptsgtf) or die "\nERROR:\t Can not open file $transcriptsgtf\n";
				(%ARFPKM,%CHFPKM, %BEFPKM, %CFPKM, %DFPKM, %DHFPKM, %DLFPKM, %cfpkm, %dfpkm, %dhfpkm, %dlfpkm)= ();
				my $i=1;
				while (<FPKM>){
					chomp;
					my ($chrom_no, $tool, $typeid, $chrom_start, $chrom_stop, $qual, $orn, $idk, $therest ) = split /\t/;
					unless ($chrom_no =~ /^chr/) { $chrom_no = "chr".$chrom_no; }
               if ($typeid =~ /^transcript/){ #check to make sure only transcripts are inputed
						my %Drest = ();
						foreach (split("\";", $therest)) { $_ =~ s/\s+|\s+//g;my($a, $b) = split /\"/; $Drest{$a} = $b;}
						my $dstax;
						if (length $Drest{'gene_id'} > 1) {
							$dstax = "$Drest{'gene_id'}-$chrom_no";
                  } else {$dstax = "xxx".$i++."-$chrom_no";}
						if (exists $CHFPKM{$dstax}){ #chromsome stop
							if ($chrom_stop > $CHFPKM{$dstax}) {
								$CHFPKM{$dstax} = $chrom_stop;
							}
						}else {
							$CHFPKM{$dstax} = $chrom_stop;
						}
						if (exists $BEFPKM{$dstax}){ #chromsome start
							if ($chrom_start < $BEFPKM{$dstax}) {
								$BEFPKM{$dstax} = $chrom_start;
							}
						}else {
							$BEFPKM{$dstax} = $chrom_start;
						}
						unless (exists $CFPKM{$dstax}{$Drest{'cov'}}){ #coverage
							$CFPKM{$dstax}{$Drest{'cov'}}= $Drest{'cov'};
						}unless (exists $DFPKM{$dstax}{$Drest{'FPKM'}}){ #FPKM
							$DFPKM{$dstax}{$Drest{'FPKM'}}= $Drest{'FPKM'};
						}
						unless (exists $DHFPKM{$dstax}{$Drest{'conf_hi'}}){ #FPKM_hi
							$DHFPKM{$dstax}{$Drest{'conf_hi'}}= $Drest{'conf_hi'};
						}
						unless (exists $DLFPKM{$dstax}{$Drest{'conf_lo'}}){ #FPKM_lo
							$DLFPKM{$dstax}{$Drest{'conf_lo'}}= $Drest{'conf_lo'};
						}
						$ARFPKM{$dstax}= "$_[0],$Drest{'gene_id'},$chrom_no";
					}
				} close FPKM;
				#sorting the fpkm values and coverage results.
				foreach my $a (keys %DFPKM){
					my $total = 0;
					foreach my $b (keys %{$DFPKM{$a}}) { $total = $b+$total; }
					$dfpkm{$a} = $total;
				}
				foreach my $a (keys %CFPKM){
					my $total = 0;
					foreach my $b (keys %{$CFPKM{$a}}) { $total = $b+$total; }
					$cfpkm{$a} = $total;
				}
				foreach my $a (keys %DHFPKM){
					my $total = 0;
					foreach my $b (keys %{$DHFPKM{$a}}) { $total = $b+$total; }
					$dhfpkm{$a} = $total;
				}
				foreach my $a (keys %DLFPKM){
					my $total = 0;
					foreach my $b (keys %{$DLFPKM{$a}}) { $total = $b+$total; }
					$dlfpkm{$a} = $total;
				}
				#end of sort.
				#insert into database.
				$genes = scalar (keys %ARFPKM);
				$sth = $dbh->prepare("update genes_summary set genes = $genes, diffexpresstool = '$diffexpress' where library_id= '$_[0]'"); $sth ->execute(); #updating genes_summary table.
				unless ($genes == $genecount) {
					unless ($genecount == 0 ) {
						print "NOTICE:\t Removed incomplete records for $_[0] in genes_fpkm table\n";
						$sth = $dbh->prepare("delete from genes_fpkm where library_id = $_[0]"); $sth->execute();
					}
					print "NOTICE:\t Importing $diffexpress expression information for $_[0] to genes_fpkm table ...";
					
					##change to thread
					my @fpkmdetails = map {$_; } sort keys %ARFPKM;
					undef @VAR; undef @threads;
					push @VAR, [ splice @fpkmdetails, 0, 200 ] while @fpkmdetails; #sub the files to multiple subs
					$queue = new Thread::Queue();
					my $builder=threads->create(\&main); #create thread for each subarray into a thread
					push @threads, threads->create(\&gtfcuffprocessor) for 1..5; #execute 5 threads
					$builder->join; #join threads
					foreach (@threads){$_->join;}
					
					print " Done\n";
					#set genes_summary to Done
					$sth = $dbh->prepare("update genes_summary set status = 'done' where library_id = '$_[0]'");
					$sth ->execute() or die "\nERROR:\t Complication in genes_summary table, consult documentation\n";
				}	else {
						print "NOTICE:\t $_[0] already in genes_fpkm table... Moving on \n";	
				}	
			} # end if cufflinks
			elsif (`head -n 1 $transcriptsgtf` =~ /stringtie/i) { #working with stringtie output
				$diffexpress = substr( `head -n 2 $transcriptsgtf | tail -1`,2,-1);
				open(FPKM, "<", $transcriptsgtf) or die "\nERROR:\t Can not open file $transcriptsgtf\n";
				(%ARFPKM,%CHFPKM, %BEFPKM, %CFPKM, %DFPKM, %TPM, %cfpkm, %dfpkm, %tpm)= ();
				my $i=1;
				while (<FPKM>){
					chomp;
					my ($chrom_no, $tool, $typeid, $chrom_start, $chrom_stop, $qual, $orn, $idk, $therest ) = split /\t/;
					unless ($chrom_no =~ /^chr/) { $chrom_no = "chr".$chrom_no; }
               if ($typeid && $typeid =~ /^transcript/){ #check to make sure only transcripts are inputed
						my %Drest = ();
						foreach (split("\";", $therest)) { $_ =~ s/\s+|\s+//g; my($a, $b) = split /\"/; $Drest{$a} = $b;}
						my $dstax;
						if (length $Drest{'gene_id'} > 1) {
							$dstax = "$Drest{'gene_id'}-$chrom_no";
                  } else {$dstax = "xxx".$i++."-$chrom_no";}
						if (exists $CHFPKM{$dstax}){ #chromsome stop
							if ($chrom_stop > $CHFPKM{$dstax}) {
								$CHFPKM{$dstax} = $chrom_stop;
							}
						}else {
							$CHFPKM{$dstax} = $chrom_stop;
						}
						if (exists $BEFPKM{$dstax}){ #chromsome start
							if ($chrom_start < $BEFPKM{$dstax}) {
								$BEFPKM{$dstax} = $chrom_start;
							}
						}else {
							$BEFPKM{$dstax} = $chrom_start;
						}
						unless (exists $CFPKM{$dstax}{$Drest{'cov'}}){ #coverage
							$CFPKM{$dstax}{$Drest{'cov'}}= $Drest{'cov'};
						}unless (exists $DFPKM{$dstax}{$Drest{'FPKM'}}){ #FPKM
							$DFPKM{$dstax}{$Drest{'FPKM'}}= $Drest{'FPKM'};
						}
						unless (exists $TPM{$dstax}{$Drest{'TPM'}}){ #FPKM_hi
							$TPM{$dstax}{$Drest{'TPM'}}= $Drest{'TPM'};
						}
						unless ($Drest{'ref_gene_name'}){
							$ARFPKM{$dstax}= "$_[0],$Drest{'gene_id'}, ,$chrom_no";
						} else {
							$ARFPKM{$dstax}= "$_[0],$Drest{'gene_id'},$Drest{'ref_gene_name'},$chrom_no";
						}
					}
				} close FPKM;
				#sorting the fpkm values and coverage results.
				foreach my $a (keys %DFPKM){
					my $total = 0;
					foreach my $b (keys %{$DFPKM{$a}}) { $total = $b+$total; }
					$dfpkm{$a} = $total;
				}
				foreach my $a (keys %CFPKM){
					my $total = 0;
					foreach my $b (keys %{$CFPKM{$a}}) { $total = $b+$total; }
					$cfpkm{$a} = $total;
				}
				foreach my $a (keys %TPM){
					my $total = 0;
					foreach my $b (keys %{$TPM{$a}}) { $total = $b+$total; }
					$tpm{$a} = $total;
				}
				#end of sort.
				#insert into database.
				$genes = scalar (keys %ARFPKM);
				$sth = $dbh->prepare("update genes_summary set genes = $genes, diffexpresstool = '$diffexpress' where library_id= '$_[0]'"); $sth ->execute(); #updating genes_summary table.
			
				unless ($genes == $genecount) {
					unless ($genecount == 0 ) {
						print "NOTICE:\t Removed incomplete records for $_[0] in genes_fpkm table\n";
						$sth = $dbh->prepare("delete from genes_fpkm where library_id = '$_[0]'"); $sth->execute();
					}
					print "NOTICE:\t Importing StringTie expression information for $_[0] to genes_fpkm table ...";
					
					##change to thread
					my @fpkmdetails = map {$_; } sort keys %ARFPKM;
					undef @VAR; undef @threads;
					push @VAR, [ splice @fpkmdetails, 0, 200 ] while @fpkmdetails; #sub the files to multiple subs
					$queue = new Thread::Queue();
					my $builder=threads->create(\&main); #create thread for each subarray into a thread
					push @threads, threads->create(\&strprocessor) for 1..5; #execute 5 threads
					$builder->join; #join threads
					foreach (@threads){$_->join;}

					print " Done\n";
					#set genes_summary to Done
					$sth = $dbh->prepare("update genes_summary set status = 'done' where library_id = '$_[0]'");
					$sth ->execute() or die "\nERROR:\t Complication in genes_summary table, consult documentation\n";
				}	else {
						print "NOTICE:\t $_[0] already in genes_fpkm table... Moving on \n";
				}	
			} else {
				die "\nFAILED:\tCan not identify source of Genes Expression File '$transcriptsgtf', consult documentation.\n";
			}
		} else {
			die "\nERROR:\t Can not find gene expression file, making sure transcript files are present or StringTie file ends with .gtf, 'e.g. <xxx>.gtf'.\n";
		}
	} else {
		print "NOTICE:\t $_[0] already completed in genes_summary tables ... Moving on \n";
	}
}

sub PARSING {
  $dbh->disconnect();
	$dbh = mysql();
  print "\n\tINSERTING TRANSCRIPTS INTO THE DATABASE : \t library_$_[0]\n\n";
  $lib_id = $_[0]; my $librarydir = $_[1]; my $verdict = 0;
  #created a log file check because
  #making sure the input genome is chicken or gallus
	unless (($refgenomename =~ /chicken/i) || ($refgenomename =~ /galgal/i)){ die "$refgenomename is not chicken or gallus\n";} 
  #making sure I'm working on only the chicken files for now, need to find annotation of alligator
  my @checkerno = split('\/',$allgeninfo[$#allgeninfo]);
	unless ($checkerno[$#checkerno] =~ /pipe/) {
		@numberz = split('_', $checkerno[$#checkerno]);
	} else {
		$numberz[0] = $lib_id;
	}
	if ($numberz[0] == $lib_id){
    #making sure the arguments are accurately parsed
    if (exists $parsablegenomes{$refgenomename} || $refgenome =~ /Galgal4/){
      open(ALIGN, "<$alignfile") or die "Can't open file $alignfile\n";
			
      # PARSER FOR transcripts_summary TABLE
      if ($alignfile) {
				`head -n 1 $alignfile` =~ /^(\d+)\sreads/; $total = $1;
				open(ALIGN,"<", $alignfile) or die "\nFAILED:\t Can not open Alignment summary file '$alignfile'\n";
        while (<ALIGN>){
          chomp;
          if (/Input/){my $line = $_; $line =~ /Input.*:\s+(\d+)$/;$total = $1;}
					if (/overall/) {  my $line = $_; $line =~ /(\d+.\d+)%\s/; $alignrate = $1;}
					if (/overall read mapping rate/) {
						if ($mappingtool){
							unless ($mappingtool =~ /TopHat/i){
								die "\nERROR:\t Inconsistent Directory Structure, $mappingtool SAM file with TopHat align_summary.txt file found\n";
							}
						} else { $mappingtool = "TopHat"; }
					}
					if (/overall alignment rate/) {
						if ($mappingtool){
							unless ($mappingtool =~ /hisat/i){
								die "\nERROR:\t Inconsistent Directory Structure, $mappingtool LOG file with HISAT align_summary.txt file found\n";
							}
						} else { $mappingtool = "HISAT";}
					}
				} close ALIGN; 
				$mapped = ceil($total * $alignrate/100);
      } else {die "\nFAILED:\t Can not find Alignment summary file as 'align_summary.txt'\n";}
     	$deletions = undef; $insertions = undef; $junctions = undef;
			if ($deletionsfile){ $deletions = `cat $deletionsfile | wc -l`; $deletions--; } 
			if ($insertionsfile){ $insertions = `cat $insertionsfile | wc -l`; $insertions--; }
			if ($junctionsfile){ $junctions = `cat $junctionsfile | wc -l`; $junctions--; }
      $unmapped = $total-$mapped;
      $date = `date +%Y-%m-%d`;
      #PARSING FOR SNPanalysis
      @parse = split('\/\/',$accepted); $accepted = undef; $len = $#parse+1; foreach(@parse){$accepted .= $_; if($len>1){$accepted .="\/"; $len--;}};
      @parse = split('\/\/',$run_log); $run_log = undef; $len = $#parse+1; foreach(@parse){$run_log .= $_; if($len>1){$run_log .="\/"; $len--;}};
    
      #INSERT INTO DATABASE : transcriptatlas
			$sth = $dbh->prepare("select library_id from transcripts_summary where library_id = '$lib_id'"); $sth->execute(); $found = $sth->fetch();
			unless ($found) {
				print "NOTICE:\t $lib_id inserting $lib_id metadata details into the database ...";
				#transcripts_summary table
				$sth = $dbh->prepare("insert into transcripts_summary (library_id, date ) values (?,?)");
				$sth ->execute($lib_id, $date);
				
				#frnak_metadata table
				$annfileversion = substr(`head -n 1 $annotationfile`,2,-1); #annotation file version
				$sth = $dbh->prepare("insert into frnak_metadata (library_id,ref_genome, ann_file, ann_file_ver, stranded, sequences,user ) values (?,?,?,?,?,?,?)");
				$sth ->execute($lib_id, "Galgal5", $annotation, $annfileversion, $stranded,$sequences,"from raven" );
				
				print "Done \n";
			}else {
				print "NOTICE:\t $lib_id already in transcripts_summary tables... Moving on \n";
			}
			
			#update
			if ($insertions > 0) {
				$sth = $dbh->prepare("update transcripts_summary set total_reads = $total, mapped_reads = $mapped, unmapped_reads = $unmapped, deletions = $deletions, insertions = $insertions, junctions = $junctions where library_id = $lib_id");
				$sth ->execute();
			} else {
				$sth = $dbh->prepare("update transcripts_summary set total_reads = $total, mapped_reads = $mapped, unmapped_reads = $unmapped where library_id = $lib_id");
				$sth ->execute();
			}
			GENES_FPKM($lib_id);

      if ($htseqcount){ HTSEQ($htseqcount); } #htseqcount details to the database.

      #variant analysis
      if ($refgenome =~ /Galgal/){$refgenomename = "chicken";}
      VARIANTS($lib_id, $accepted, $refgenomename, $annotationfile);

      #Finally : the last update. transcripts_summary table updating status column with 'done'
      $dbh = mysql();
			$sth = $dbh->prepare("update transcripts_summary set status='done' where library_id = $lib_id");
      $sth ->execute();

      #TRY to implement nosql ###fix
      $verdict = 1;
    }
    else { 
      my $parsabletemp = 1;
      $verdict = 0;
      print "The reference genome isn't available, available genomes are : ";
      foreach my $pargenomes (keys %parsablegenomes){
        print "\"$pargenomes\"";
        if($parsabletemp < (keys %parsablegenomes)){ print ", "; $parsabletemp++; }
      }
      print " rather what you have is $refgenome\n$mystderror\n";
    }
  }
  else {
    print "library_id dont match $numberz[0] == $lib_id\n";
    $verdict = 0;
  }
  return $verdict;
}
sub HTSEQ { #importing Htseqcount details to the database
   print "\n\tSTORING HTSEQ IN THE DATABASE\n\n";
	my $htseqnumber = $dbh->selectrow_array("select count(*) from htseq where library_id = '$lib_id'");
	my $htseq = `cat $_[0] | wc -l`; if ($htseq >1){ $htseq -= 1;} else {$htseq = 0;} 
   unless ($htseq == $htseqnumber) {
		if ($htseqnumber > 1 ){ $sth = $dbh->prepare("delete from htseq where library_id = '$lib_id'"); $sth->execute(); }
      open(HTSEQ, "<$_[0]") or die "Can't open file $_[0]\n"; 
		##change to thread
		my @htseqdetails = <HTSEQ>; close(HTSEQ);
		undef @VAR; undef @threads;
		push @VAR, [ splice @htseqdetails, 0, 200 ] while @htseqdetails; #sub the files to multiple subs
		
		$queue = new Thread::Queue();
		my $builder=threads->create(\&main); #create thread for each subarray into a thread
		push @threads, threads->create(\&htseqprocessor) for 1..5; #execute 5 threads
		$builder->join; #join threads
		foreach (@threads){$_->join;}
	}
}

sub VARIANTS { #process variants
  $sth = $dbh->prepare("select library_id from variants_summary where library_id = '$_[0]' and status = 'done'"); $sth->execute(); $found = $sth->fetch();	
	unless ($found) {
		print "\n\tWORKING ON VARIANT ANALYSIS\n\n";
		my $libraryNO = "library_".$_[0];
      $Mfolder = "$STORAGEPATH/$libraryNO"; `mkdir -p $Mfolder`; print "\tmade $Mfolder\n"; #decided to keep all the variant folders.
		my $specie = $_[2];
		my $REF= "$GENOMES/$_[2]/$_[2]".".fa";
		my $ANN = $_[3];
		unless ($variantfile){ #check if variantfile exists if not, it will be generated using VAP
			print "NOTICE: Variant file isn't present, creating variantfile using VAP details & is stored in $STORAGEPATH\n";
			my $bamfile = $_[1];
         my $path = `readlink -f $in1`; chop $path; my @path = split('\/', $path); undef $path; foreach (0..$#path-1) { $path .= $path[$_]."/"; }

         #create .fna dictionary
         my ($DATA,$OUT1,$OUT2, %SEQ, %ORDER, %SEQnum, %SEQheader);
         open($OUT1, "> $Mfolder/chicken.fna") or die $!;
         open($OUT2, "> $Mfolder/chicken.fna.fai") or die $!;
         $/ = "\>";
         open ($DATA,"<$path/chicken/chicken.fa") or die $!;
         my @fastqfile = <$DATA>;
         shift(@fastqfile);
         my $ii = 0;
         foreach my $entry (@fastqfile){
            my @pieces = split(/\n/, $entry);
            $pieces[0] = (split(' ',$pieces[0]))[0];
            $ii++;
            $ORDER{$ii} = $pieces[0];
            my $seq = '';
            foreach my $num (1.. $#pieces-1){
               $seq .= $pieces[$num];
            }
            if($pieces[$#pieces] =~ />$/) { $seq .= substr($pieces[$#pieces],0,-1); }
            else { $seq .= $pieces[$#pieces]; }
            $SEQ{$pieces[0]} = $seq;
            $SEQnum{$pieces[0]} = length($seq);
            $SEQheader{$pieces[0]} = length($pieces[0]);
         }
         my ($check, $start, $newstart, $last);
         foreach my $number (sort {$a <=> $b} keys %ORDER){
            my $header = $ORDER{$number};
            if (length($header) >= 1) {
               print $OUT1 ">$header\n$SEQ{$header}\n";
               unless ($check){
                  $start = $SEQheader{$header}+2;
                  $last = $SEQnum{$header}+1;
                  print $OUT2 "$header\t$SEQnum{$header}\t$start\t$SEQnum{$header}\t$last\n";
                  $check = "yes";
               }
               else {
                  $newstart = $SEQheader{$header}+2+$last+$start;
                  $start = $newstart;
                  $last = $SEQnum{$header}+1;
                  print $OUT2 "$header\t$SEQnum{$header}\t$start\t$SEQnum{$header}\t$last\n";
                  $check = "yes";
               }
            }
         }
         close $DATA; close $OUT1; close $OUT2;
         $/ = "\n";
         
         my $REF = "$Mfolder/chicken.fna";
			my $ANN = $_[3];
  
			$Mfolder = "$STORAGEPATH/$libraryNO"; `mkdir -p $Mfolder`; print "\tmade $Mfolder\n"; #decided to keep all the variant folders.
			$gatk_version = ( ( split('\/',$GATKDIR)) [-2] );
			$picard_version = ( ( split('\/',$PICARDDIR)) [-2] );
	
			#VARIANT ANALYSIS
			#PICARD
			my $filename = "$Mfolder/$libraryNO.vcf";
			unless (-e $filename) {
            
            #create GATK sequence dictionary
            $filename = "$Mfolder/chicken.dict";
            unless (-e "$filename"){
               `java -jar $PICARDDIR CreateSequenceDictionary R=$REF O=$filename CREATE_INDEX=true`;
            }
         
				#SORT BAM
            $filename = "$Mfolder/$libraryNO.bam";
				unless (-e $filename){
					`java -jar $PICARDDIR SortSam INPUT=$bamfile OUTPUT=$filename SO=coordinate`;
				} else { print "NOTICE: $filename exists\n"; }
		
				#ADDREADGROUPS
				$filename = "$Mfolder/$libraryNO"."_add.bam";
				unless (-e "$filename"){
					my $addreadgroup = "java -jar $PICARDDIR AddOrReplaceReadGroups INPUT=$Mfolder/$libraryNO".".bam OUTPUT=$filename SO=coordinate RGID=LAbel RGLB=Label RGPL=illumina RGPU=Label 	RGSM=Label";
					`$addreadgroup`;
					print "NOTICE: Add read groups complete\n";
				} else { print "$Mfolder/$libraryNO"."_add.bam already exists\n"; }
				
				#MARKDUPLICATES
				unless (-e "$Mfolder/$libraryNO"."_mdup.bam" ) {
					my $markduplicates = "java -jar $PICARDDIR MarkDuplicates INPUT=$Mfolder/".$libraryNO."_add.bam OUTPUT=$Mfolder/".$libraryNO."_mdup.bam M=$Mfolder/".$libraryNO."_mdup.metrics 			CREATE_INDEX=true";
					`$markduplicates`;
					print "NOTICE: Mark duplicates complete\n";
				} else { print "$Mfolder/$libraryNO"."_mdup.bam already exists\n"; }
  
				#SPLIT&TRIM
				unless (-e "$Mfolder/$libraryNO"."_split.bam" ) {  
					my $splittrim = "java -jar $GATKDIR -T SplitNCigarReads -R $REF -I $Mfolder/".$libraryNO."_mdup.bam -o $Mfolder/".$libraryNO."_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 --filter_reads_with_N_cigar";
					`$splittrim`;
					print "NOTICE: Split N Cigar reads complete\n";
					} else { print "$Mfolder/$libraryNO"."_split.bam already exists \n"; }
			
				#GATK
				unless (-e "$Mfolder/$libraryNO.vcf" ) { 
					my $gatk = "java -jar $GATKDIR -T HaplotypeCaller -R $REF -I $Mfolder/".$libraryNO."_split.bam -o $Mfolder/$libraryNO.vcf";
					`$gatk`;
					print "NOTICE: Haplotype caller complete\n";
				} else { print "$Mfolder/$libraryNO".".vcf already exists \n"; }	
			} else {
				print "NOTICE: Variants VCF in $STORAGEPATH is already created\n";
			}
			#perl to select DP > 5 & get header information
			FILTERING($Mfolder, "$Mfolder/$libraryNO.vcf");
			DBVARIANTS("$Mfolder/$libraryNO"."_DP5.vcf", $libraryNO);
			
		} else {
			print "NOTICE:\t Variants VCF already created for $_[0]\n";
			#perl to select DP > 5 & get header information
			FILTERING($Mfolder, $variantfile);
			unless ($vepfile) { #filter if variant-annotation file doesn't already exists
				@foldercontent = split("\n", `find $parsedinput`); #get details of the main folder
				my $DP5file = (grep /_DP5.vcf$/, @foldercontent)[0];
				DBVARIANTS($DP5file, $libraryNO);
			} else {
				DBVARIANTS($variantfile, $libraryNO);
			}
		}
		#ANNOTATIONS : running VEP
		my ($DP5file, $VEPtxt);
		unless ($vepfile){ #check if vepfile exists if not, it will be generated using VEP
			unless ($variantfile){
				$vep_version = ( ( split('\/',$VEP)) [-4] );
				$DP5file = $Mfolder."/".$libraryNO."_DP5.vcf";
				$VEPtxt = $Mfolder."/".$libraryNO."_VEP.txt";
#				print "this is the species $specie\n"; #Remove this Modupe
			} else {
				$DP5file = (grep /_DP5.vcf$/, @foldercontent)[0]; #parsing DP5 file to VEPVARIANT
				$VEPtxt = $DP5file."_VEP.txt";
			}
			if (exists $VEPparse{$specie}){
				print "NOTICE:\t VEP Gene Variant Annotation being performed for $_[0] . . . ";
				my $veptxt = "perl $VEP -i $DP5file --fork 24 --species $specie  --dir /home/modupe/.vep/ --cache --merged --everything on --terms ensembl --output_file $VEPtxt"; `$veptxt`;
				#my $vepvcf = "perl $VEP -i $Mfolder/".$libraryNO."_DP5.vcf --fork 24 --species $specie  --dir /home/modupe/.vep/ --cache --vcf --merged --everything on --terms ensembl 		--output_file $Mfolder/".$libraryNO."_VEP.vcf"; `$vepvcf`;
				print "Done\n";
				VEPVARIANT($VEPtxt, $libraryNO);  #import variants to database
			} else {
				next "Unidentified genome\t$mystderror\n";
			}
		} else {
			print "NOTICE: VEP file already exists\n";
			VEPVARIANT($vepfile, $libraryNO);  #import variants to database
		}
	} #end unless the variants are already in the database.
	else {
		print "NOTICE: $_[0] Already exists in variants summary table\n";
	}
} 
sub FILTERING {
  my $input = $_[1];
  my $wkdir = $_[0];
  unless(open(FILE,$input)){
    print "File \'$input\' doesn't exist\n";
    exit;
  }
  my $out = fileparse($input, qr/(\.vcf)?$/);
  my $output = "$out"."_DP5.vcf";
  open(OUT,">$wkdir/$output");
  my $output2 = "$out"."_header.vcf";
  open(OUT2,">$wkdir/$output2");

  my @file = <FILE>; chomp @file; close (FILE);
  foreach my $chr (@file){
    unless ($chr =~ /^\#/){
      my @chrdetails = split('\t', $chr);
      my $chrIwant = $chrdetails[7];
      my @morechrsplit = split(';', $chrIwant);
      foreach my $Imptchr (@morechrsplit){
        if ($Imptchr =~ m/^DP/) {
          my @addchrsplit = split('=', $Imptchr);
          if ($addchrsplit[1] > 4){print OUT "$chr\n";}
        }
      }
    }
    else {
      print OUT "$chr\n"; print OUT2 "$chr\n";
    }
  }
  close (OUT); close (OUT2);
}

sub DBVARIANTS {
	#INSERT INTO DATABASE: #variants_summary table
	print "\n\tINSERTING VARIANTS INTO THE DATABASE\n\n";
  #disconnecting and connecting again to database just incase
  my ($toolvariant, $verd, $variantclass);
	$dbh->disconnect(); 
  $dbh = mysql();
	
	$_[1] =~ /^library_(\d*)$/;
  my $libnumber = $1;
  my $folder = undef;
	my ($itsnp,$itindel,$itvariants) = (0,0,0);
	#VEP file
  my @splitinput = split('\/', $_[0]);
  foreach my $i (0..$#splitinput-1){$folder.="$splitinput[$i]/";$i++;}
  my $information = fileparse($_[0], qr/(\.vcf)?$/);
	
	undef %VCFhash;
	if($_[0]){ open(VARVCF,$_[0]) or die ("\nERROR:\t Can not open variant file $_[0]\n"); } else { die ("\nERROR:\t Can not find variant file. make sure variant file with suffix '.vcf' is present\n"); }
	while (<VARVCF>) {
		chomp;
		if (/^\#/) {
			if (/^\#\#GATK/) {
				$_ =~ /ID\=(.*)\,.*Version\=(.*)\,Date/;
				$toolvariant = "GATK v.$2,$1";
				$varianttool = "GATK";
			} elsif (/^\#\#samtoolsVersion/){
				$_ =~ /Version\=(.*)\+./;
				$toolvariant = "samtools v.$1";
				$varianttool = "samtools";
			}
		} else {
			my @chrdetails = split "\t";
			my @morechrsplit = split(';', $chrdetails[7]);
         unless ($chrdetails[0] =~ /^chr/) { $chrdetails[0] = "chr".$chrdetails[0]; }
			if (((split(':', $chrdetails[9]))[0]) eq '0/1'){$verd = "heterozygous";}
			elsif (((split(':', $chrdetails[9]))[0]) eq '1/1'){$verd = "homozygous";}
			elsif (((split(':', $chrdetails[9]))[0]) eq '1/2'){$verd = "heterozygous alternate";}
			$VCFhash{$chrdetails[0]}{$chrdetails[1]} = "$chrdetails[3]|$chrdetails[4]|$chrdetails[5]|$verd";
		}
	} close VARVCF;
	$sth = $dbh->prepare("select library_id from variants_summary where library_id = '$libnumber' and status = 'done'"); $sth->execute(); $found = $sth->fetch();	
	unless ($found) {
		$sth = $dbh->prepare("select library_id from variants_summary where library_id = '$libnumber'"); $sth->execute(); $found = $sth->fetch();	
		if ($found) { #deleting previous records if there are
			$sth = $dbh->prepare("delete from variants_annotation where library_id = '$libnumber'"); $sth->execute();
			$sth = $dbh->prepare("delete from variants_result where library_id = '$libnumber'"); $sth->execute();
			$sth = $dbh->prepare("delete from variants_summary where library_id = '$libnumber'"); $sth->execute();
		} # end unless
		$vep_version = ( ( split('\/',$VEP)) [-4] ); 
		$sth = $dbh->prepare("insert into variants_summary ( library_id, ANN_version, Picard_version, GATK_version, variant_tool, date ) values (?,?,?,?,?,?)");
		$sth ->execute($libnumber, $vep_version, $picard_version, $gatk_version, $toolvariant, $date);
    
		#VARIANT_RESULTS
		print "NOTICE:\t Importing $varianttool variant information for $libnumber to variants_result table ...";
		my $ii = 0;
		foreach my $abc (sort keys %VCFhash) {
			foreach my $def (sort {$a <=> $b} keys %{ $VCFhash{$abc} }) {
				my @vcf = split('\|', $VCFhash{$abc}{$def});
				if ($vcf[3] =~ /,/){
					my $first = split(",",$vcf[1]);
					if (length $vcf[0] == length $first){ $itvariants++; $itsnp++; $variantclass = "SNV"; }
					elsif (length $vcf[0] < length $first) { $itvariants++; $itindel++; $variantclass = "insertion"; }
					else { $itvariants++; $itindel++; $variantclass = "deletion"; }
				}
				elsif (length $vcf[0] == length $vcf[1]){ $itvariants++; $itsnp++; $variantclass = "SNV"; }
				elsif (length $vcf[0] < length $vcf[1]) { $itvariants++; $itindel++; $variantclass = "insertion"; }
				else { $itvariants++; $itindel++; $variantclass = "deletion"; }
		
				#putting variants info into a hash table
				my @hashdbvariant = ($libnumber, $abc, $def, $vcf[0], $vcf[1], $vcf[2], $variantclass, $vcf[3]); 
				$HASHDBVARIANT{$ii++} = [@hashdbvariant];
				
				#to variant_result
				#$sth = $dbh->prepare("insert into variants_result ( library_id, chrom, position, ref_allele, alt_allele, quality, variant_class, zygosity ) values (?,?,?,?,?,?,?,?)");
				#$sth ->execute($libnumber, $abc, $def, $vcf[0], $vcf[1], $vcf[2], $variantclass, $vcf[3]) or die "\nERROR:\t Complication in variants_result table, consult documentation\n";
			}
		}
		#update variantsummary with counts
		$sth = $dbh->prepare("update variants_summary set total_VARIANTS = $itvariants, total_SNPS = $itsnp, total_INDELS = $itindel where library_id= '$libnumber'"); 
		$sth ->execute() or die "$DBI::errstr Error in updating the Variants Summary table\n";
		$sth->finish(); 
		
		`date`; 
		#threads to import variants
		#print Data::Dumper->Dump( [ \%HASHDBVARIANT ], [ qw(*thehash) ] );
		my @hashdetails = keys %HASHDBVARIANT; #print "First $#hashdetails\n"; die;
		undef @VAR; undef @threads;
		push @VAR, [ splice @hashdetails, 0, 200 ] while @hashdetails; #sub the files to multiple subs
		$queue = new Thread::Queue();
		my $builder=threads->create(\&main); #create thread for each subarray into a thread 
		push @threads, threads->create(\&dbvarprocessor) for 1..5; #execute 5 threads
		$builder->join; #join threads print "I'm here 2\n";
		foreach (@threads){$_->join;}
		`date`;
		print "Done\n";
	} else { print "NOTICE: $libnumber Already exists in variants summary table\n"; }
}

sub VEPVARIANT {
  #INSERT INTO DATABASE: #variants_annotation table
	print "\n\tINSERTING VARIANTS - ANNOTATION INTO THE DATABASE\n\n";
  #disconnecting and connecting again to database just incase
	$dbh->disconnect(); 
  $dbh = mysql();
	undef %extra;
	undef %DBSNP;
	$_[1] =~ /^library_(\d*)$/;
  print "NOTICE:\t Importing VEP variant annotation for $1 to variants_result table ...";
	my $libnumber = $1; my $ii = 0;
	my ($chrom, $position);
	if($_[0]){ open(VEP,$_[0]) or die ("\nERROR:\t Can not open vep file $_[0]\n"); } else { die ("\nERROR:\t Can not find VEP file. make sure vep file with suffix '.vep.txt' is present\n"); }
	while (<VEP>) {
		chomp;
		unless (/^\#/) {
			unless (/within_non_coding_gene/i || /coding_unknown/i) {
				my @veparray = split "\t"; #14 columns
				my @extraarray = split(";", $veparray[13]);
				foreach (@extraarray) { my @earray = split "\="; $extra{$earray[0]}=$earray[1]; }
				my @indentation = split("_", $veparray[0]);
				if ($#indentation > 2) { $chrom = $indentation[0]."_".$indentation[1]; $position = $indentation[2]; }
				else { $chrom = $indentation[0]; $position = $indentation[1]; }
				unless ($chrom =~ /^chr/) { $chrom = "chr".$chrom; }
				unless ( $extra{'VARIANT_CLASS'} =~ "SNV" or $extra{'VARIANT_CLASS'} =~ "substitution" ){ $position--; }
				else {
					my @poly = split("/",$indentation[$#indentation]);
					unless ($#poly > 1){ unless (length ($poly[0]) == length($poly[1])){ $position--; } }
				}
				my $geneid = $veparray[3];
				my $transcriptid = $veparray[4];
				my $featuretype = $veparray[5];
				my $consequence = $veparray[6]; 
				if ($consequence =~ /NON_(.*)$/){ $consequence = "NON".$1; } elsif ($consequence =~ /STOP_(.*)$/) {$consequence = "STOP".$1; }
				my $pposition = $veparray[9];
				my $aminoacid = $veparray[10];
				my $codons = $veparray[11];
				my $dbsnp = $veparray[12];
				my $locate = "$_[1],$chrom,$position,$consequence,$geneid,$pposition";
				if ( exists $VEPhash{$locate} ) {
					unless ( $VEPhash{$locate} eq $locate ){ die "\nERROR:\t Duplicate annotation in VEP file, consult documentation\n"; }
				} else {
					$VEPhash{$locate} = $locate;
					#$sth = $dbh->prepare("insert into variants_annotation ( library_id, chrom, position, consequence, gene_id, gene_name, transcript, feature, gene_type,protein_position, aminoacid_change, codon_change ) values (?,?,?,?,?,?,?,?,?,?,?,?)");
					if (exists $extra{'SYMBOL'}) { $extra{'SYMBOL'} = uc($extra{'SYMBOL'}); } else { $extra{'SYMBOL'} = "NULL"; }
					
					my @hashdbvariant = ($libnumber, $chrom, $position, $consequence, $geneid, $extra{'SYMBOL'}, $transcriptid, $featuretype, $extra{'BIOTYPE'} , $pposition, $aminoacid, $codons);
					$HASHDBVEP{$ii++} = [@hashdbvariant];
					#$sth ->execute($libnumber, $chrom, $position, $consequence, $geneid, $extra{'SYMBOL'}, $transcriptid, $featuretype, $extra{'BIOTYPE'} , $pposition, $aminoacid, $codons) or die "\nERROR:\t Complication in variants_annotation table, consult documentation\n";
					$HASHRESULT{$chrom}{$position} = $extra{'VARIANT_CLASS'};
					#$HASHRESULT{$kk++} = ($libnumber, $chrom, $position, $extra{'VARIANT_CLASS'});
					
					#$sth = $dbh->prepare("update variants_result set variant_class = '$extra{'VARIANT_CLASS'}' where library_id = '$libnumber' and chrom = '$chrom' and position = $position"); $sth ->execute() or die "\nERROR:\t Complication in updating VarResult table, consult documentation\n";
					
					$DBSNP{$chrom}{$position} = $dbsnp; #updating dbsnp
					#$DBSNP{$jj++} = ($libnumber,$chrom,$position,$dbsnp); #updating dbsnp	
				}
			}
		} else {
			if (/API (version \d+)/){
				unless($vep_version) {
					$vep_version = $1;
					$sth = $dbh->prepare("update variants_summary set ANN_version = 'VEP $vep_version' where library_id = '$libnumber'"); $sth ->execute();
				}
			} #getting VEP version
		}
	} close VEP; #end of processing vep file
	print "\n\nTotal number of annotations is $ii\n\n";
	
	#adding the dbsnp annotation and variant_class
	$ii = 0; my @hashdb;
	foreach my $chrom (sort keys %DBSNP) { #updating existing_variant
		foreach my $position (sort keys %{ $DBSNP{$chrom} }) {
			if (exists $HASHRESULT{$chrom}{$position}){
				@hashdb = ($libnumber,$chrom,$position,$DBSNP{$chrom}{$position},$HASHRESULT{$chrom}{$position});
				delete $HASHRESULT{$chrom}{$position};
			} else {
				@hashdb = ($libnumber,$chrom,$position,$DBSNP{$chrom}{$position},"NULL");
			}
			$HASHNEW{$ii++} = [@hashdb];
			delete $DBSNP{$chrom}{$position};
		}
	}
	foreach my $chrom (sort keys %HASHRESULT) { #updating variant class
		foreach my $position (sort keys %{ $HASHRESULT{$chrom} }) {
			if (exists $DBSNP{$chrom}{$position}){
				@hashdb = ($libnumber,$chrom,$position,$DBSNP{$chrom}{$position},$HASHRESULT{$chrom}{$position});
				delete $DBSNP{$chrom}{$position};
			} else {
				@hashdb = ($libnumber,$chrom,$position,"NULL", $HASHRESULT{$chrom}{$position});
			}
			$HASHNEW{$ii++} = [@hashdb];
			delete $HASHRESULT{$chrom}{$position};
		}
	}
	
	`date`;
		#threads to import variantsvep
		my @hashdetails = keys %HASHDBVEP;
		undef @VAR; undef @threads;
		push @VAR, [ splice @hashdetails, 0, 200 ] while @hashdetails; #sub the files to multiple subs
		$queue = new Thread::Queue();
		my $builder=threads->create(\&main); #create thread for each subarray into a thread
		push @threads, threads->create(\&dbvepprocessor) for 1..5; #execute 5 threads
		$builder->join; #join threads
		foreach (@threads){$_->join;}
	`date`;
		#threads to update dbsnp & EXISTING VARIANT
		print "done with the vep\n";
		@hashdetails = keys %HASHNEW;
		undef @VAR; undef @threads;
		push @VAR, [ splice @hashdetails, 0, 200 ] while @hashdetails; #sub the files to multiple subs
		$queue = new Thread::Queue();
		$builder=threads->create(\&main); #create thread for each subarray into a thread
		push @threads, threads->create(\&dbdetprocessor) for 1..5; #execute 5 threads
		$builder->join; #join threads
		foreach (@threads){$_->join;}
	`date`;
	
	$dbh = mysql();
	$sth = $dbh->prepare("update variants_summary set status = 'done' where library_id= '$libnumber'"); #set variants_summary status as done
	$sth ->execute();
}
sub NOTIFICATION {
  my $notification = '/home/modupe/.LOG/note.txt';
  open (NOTE, ">$notification");
  print NOTE "Subject: ". $_[0] .": $jobid\n\nName of log files\n\t$std_out\n\t$std_err\n";
  system "sendmail $email < $notification";
  close NOTE;
  system "rm -rf $notification";
}

sub main {
  foreach my $count (0..$#VAR) {
		while(1) {
			if ($queue->pending() < 100) {
				$queue->enqueue($VAR[$count]);
				last;
			}
		}
	}
	foreach(1..5) { $queue-> enqueue(undef); }
}

sub cuffprocessor { my $query; while ($query = $queue->dequeue()){ geneparseinput(@$query); } }
sub gtfcuffprocessor { my $query; while ($query = $queue->dequeue()){ gtfcuffparseinput(@$query); } }
sub isoprocessor { my $query; while ($query = $queue->dequeue()){ isoparseinput(@$query); } }
sub strprocessor { my $query; while ($query = $queue->dequeue()){ strparseinput(@$query); } }
sub htseqprocessor { my $query; while ($query = $queue->dequeue()){	htseqparseinput(@$query); } }
sub dbvarprocessor { my $query; while ($query = $queue->dequeue()){ dbvarinput(@$query); } }
sub dbvepprocessor { my $query; while ($query = $queue->dequeue()){ dbvepinput(@$query); } }
sub dbdetprocessor { my $query; while ($query = $queue->dequeue()){ dbdetinput(@$query); } }

sub gtfcuffparseinput {
	$syntax = "insert into genes_fpkm (library_id, gene_id, chrom_no, chrom_start, chrom_stop, coverage, fpkm, fpkm_conf_low, fpkm_conf_high ) values (?,?,?,?,?,?,?,?,?)";
	foreach my $a (@_) {
		$dbh = mysql();
		$sth = $dbh->prepare($syntax);
		my @array = split(",",$ARFPKM{$a});
		$sth -> execute(@array, $BEFPKM{$a}, $CHFPKM{$a}, $cfpkm{$a}, $dfpkm{$a}, $dlfpkm{$a}, $dhfpkm{$a}) or print "\nERROR:\t Complication in genes_fpkm table, consult documentation\n";
		$sth->finish;
	}
}
sub geneparseinput{
	$syntax = "insert into genes_fpkm (library_id, gene_id, gene_short_name, chrom_no, chrom_start, chrom_stop, coverage, fpkm, fpkm_conf_low, fpkm_conf_high, fpkm_status ) values (?,?,?,?,?,?,?,?,?,?,?)";
  foreach (@_){
		chomp;
		$dbh = mysql();
		$sth = $dbh->prepare($syntax);
		my ($track, $class, $ref_id, $gene, $gene_name, $tss, $locus, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat ) = split /\t/;
		unless ($track eq "tracking_id"){ #check & specifying undefined variables to null
			if($coverage =~ /-/){$coverage = undef;}
			my ($chrom_no, $chrom_start, $chrom_stop) = $locus =~ /^(.+)\:(.+)\-(.+)$/; $chrom_start++;
			$sth ->execute($lib_id, $gene, $gene_name, $chrom_no, $chrom_start, $chrom_stop, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat ) or print "\nERROR:\t Complication in genes_fpkm table, consult documentation\n";
			$sth->finish;
		}
	}
}
sub strparseinput{
	$syntax = "insert into genes_fpkm (library_id, gene_id, gene_short_name, chrom_no, chrom_start, chrom_stop, coverage, fpkm, tpm ) values (?,?,?,?,?,?,?,?,?)";
	foreach my $a (@_) {
		$dbh = mysql();
		$sth = $dbh->prepare($syntax);
		my @array = split(",",$ARFPKM{$a});
		$sth -> execute(@array, $BEFPKM{$a}, $CHFPKM{$a}, $cfpkm{$a}, $dfpkm{$a}, $tpm{$a}) or print "\nERROR:\t Complication in genes_fpkm table, consult documentation\n";
		$sth->finish;
	}
}
sub isoparseinput{
	$syntax = "insert into isoforms_fpkm (library_id, tracking_id, gene_id, gene_short_name, chrom_no, chrom_start, chrom_stop, coverage, fpkm, fpkm_conf_low, fpkm_conf_high, fpkm_status ) values (?,?,?,?,?,?,?,?,?,?,?,?)";
  foreach (@_){
		chomp;
		$dbh = mysql();
		$sth = $dbh->prepare($syntax);
		my ($track, $class, $ref_id, $gene, $gene_name, $tss, $locus, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat ) = split /\t/;
		unless ($track eq "tracking_id"){ #check & specifying undefined variables to null
			if($coverage =~ /-/){$coverage = undef;}
			my ($chrom_no, $chrom_start, $chrom_stop) = $locus =~ /^(.+)\:(.+)\-(.+)$/; $chrom_start++;
			$sth ->execute($lib_id, $track, $gene, $gene_name, $chrom_no, $chrom_start, $chrom_stop, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat ) or print "\nERROR:\t Complication in isoforms_fpkm table, consult documentation\n";
			$sth->finish;
		}
	}
}				
sub htseqparseinput {
	$syntax = "insert into htseq ( library_id, gene_name, count) values (?,?,?)";
	foreach my $a (@_) {
		chomp;
		$dbh = mysql();
		$sth = $dbh->prepare($syntax);
		my ($NAME, $VALUE) = split (/\t/, $a);
      if ($NAME =~ /^[a-z0-9A-Z]/i) {
         if ($VALUE >= 0) {
#           if ($NAME =~ /ENSGALG/) {
#			   	my $path = `readlink -f $in1`; chop $path; my @path = split('\/', $path); undef $path; foreach (0..$#path-1) { $path .= $path[$_]."/"; }
#			   	my $tobesyntax = 'j=$(grep "gene" '.$path.'chicken/*gff3 | grep "'.$NAME.'" | grep "Name"| head -n1); k=$(echo $j | awk -F\';\' \'{print $2}\'); echo $k | awk -F\'=\' \'{print $2}\'';
#			   	my $newgene = `$tobesyntax`; chop $newgene;
#			   	if ($newgene) { $NAME = $newgene;}
#			     	$sth->execute($lib_id, $NAME, $VALUE);
#           } else {
              $sth->execute($lib_id, $NAME, $VALUE);
            #}
         }
      }
		$sth->finish;
	}
}
sub dbvarinput {
	$syntax = "insert into variants_result ( library_id, chrom, position, ref_allele, alt_allele, quality, variant_class, zygosity ) values (?,?,?,?,?,?,?,?)";
	foreach my $a (@_) { 
		$dbh = mysql();
		$sth = $dbh->prepare($syntax);
		$sth -> execute(@{$HASHDBVARIANT{$a}}) or print "\nERROR:\t Complication in variants_result table, consult documentation\n";
	}
}
sub dbvepinput {
	$syntax = "insert into variants_annotation ( library_id, chrom, position, consequence, gene_id, gene_name, transcript, feature, gene_type,protein_position, aminoacid_change, codon_change ) values (?,?,?,?,?,?,?,?,?,?,?,?)";
	foreach my $a (@_) {
		$dbh = mysql();
		$sth = $dbh->prepare($syntax);
		#print "==> $a\t", @{$HASHDBVEP{$a}},"\n\n";
		$sth -> execute(@{$HASHDBVEP{$a}}) or print "\nERROR:\t Complication in variants_annotation table, consult documentation\n";
	}
}
sub dbdetinput {
	foreach my $a (@_) {
		$dbh = mysql();
		my ($liz, $chz, $poz, $dbz, $exz) = @{$HASHNEW{$a}};
		$sth = $dbh->prepare("update variants_result set existing_variant = '$dbz', variant_class = '$exz' where library_id = '$liz' and chrom = '$chz' and position = $poz");
		$sth ->execute() or print "\nERROR:\t Complication in updating variants_result table, consult documentation\n";
	}
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
close STDOUT; close STDERR;
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
exit;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR RAVENinserttranscriptome.pl

=pod

=head1 NAME

$0 -- Comprehensive pipeline : Inputs frnakenstein results from tophat and cufflinks and generates a metadata which are all stored in the database : transcriptatlas
: Performs variant analysis using a suite of tools from the output of frnakenstein and input them into the database.

=head1 SYNOPSIS

RAVENinserttranscriptome.pl [--help] [--manual] <directory of files>

=head1 DESCRIPTION

Accepts all folders from frnakenstein output.
 
=head1 OPTIONS

=over 3

=item B<--delete>

Delete incomplete libraries.  (Optional) 

=item B<-h, --help>

Displays the usage message.  (Optional) 

=item B<-man, --manual>

Displays full manual.  (Optional) 

=back

=head1 DEPENDENCIES

Requires the following Perl libraries (all standard in most Perl installs).
   DBI
   DBD::mysql
   Getopt::Long
   Pod::Usage

=head1 AUTHOR

Written by Modupe Adetunji, 
Center for Bioinformatics and Computational Biology Core Facility, University of Delaware.

=head1 REPORTING BUGS

Report bugs to amodupe@udel.edu

=head1 COPYRIGHT

Copyright 2017 MOA.  
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.  
This is free software: you are free to change and redistribute it.  
There is NO WARRANTY, to the extent permitted by law.  

Please acknowledge author and affiliation in published work arising from this script's usage
=cut


   
