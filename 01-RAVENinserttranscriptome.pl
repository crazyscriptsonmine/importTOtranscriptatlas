#!/usr/bin/perl

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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - U S E R  V A R I A B L E S- - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# CODE FOR 
use strict;
use File::Basename;
use DBI;
use Getopt::Long;
use Time::localtime;
use Pod::Usage;
use Time::Piece;
use File::stat;
use DateTime;
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
                                "delete" 	=> 	\$deletenotdone,
				"h|help"        =>      \$help,
                                "man|manual"	=>      \$manual );

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

# FOLDER VARIABLES
#my ($top_folder, $cuff_folder, $htseq_folder);

# RESULTS_HASH
my %Hashresults; my %Birdresults; my %HTSEQ;

# DATABASE VARIABLES
my ($dbh, $sth, $syntax, $row, @row);

#DIRECTORY
my (@parse, @NewDirectory);

# TABLE VARIABLES
my ($accepted, $samfile, $alignfile, $isoformsfile, $genesfile, $deletionsfile, $insertionsfile, $transcriptsgtf, $junctionsfile, $prepfile, $run_log, $htseqcount,$variantfile, $vepfile, $annofile);
my ($parsedinput, $len);
my ($lib_id, $total, $mapped, $unmapped, $deletions, $insertions, $junctions, $genes, $isoforms, $prep, $date); #transcripts_summary TABLE
my ($track, $class, $ref_id, $gene, $gene_name, $tss, $locus, $chrom_no, $chrom_start, $chrom_stop, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat); # GENES_FPKM & ISOFORMS_FPKM TABLE

#VARIANTS FOLDER & HASHES
my $Mfolder; my %VCFhash; my %VEPhash; my %ExtraWork; my %AMISG; my %AMIST;

#PARSABLE GENOMES FOR ANALYSIS
my $GENOMES="/home/modupe/.GENOMES/";
my $STORAGEPATH = "/home/modupe/CHICKENSNPS"; #variant files are stored in this directory
my %parsablegenomes = ("chicken" => 1, "alligator" => 2,"mouse" => 3, ); #genomes that work.
my %VEPparse = ("chicken" => 1,"mouse" => 2, ); #for VEP
my $executemouse; #to reworking the mouse assembly.

#INDEPENDENT PROGRAMS TO RUN
my $PICARDDIR="/home/modupe/.software/picard-tools-1.136/picard.jar";
my $GATKDIR="/home/modupe/.software/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar";
my $VEP="/home/modupe/.software/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl";
my $SNPdat="/home/modupe/.software/SNPdat_package_v1.0.5/SNPdat_v1.0.5.pl";
my $email = 'amodupe@udel.edu';

#OPENING FOLDER
opendir(DIR,$in1) or die "Folder \"$in1\" main doesn't exist\n"; 
my @Directory = readdir(DIR);
close(DIR);
#pushing each subfolder
foreach (@Directory){
  if ($_ =~ /^\w*_\d*$/){
    push (@NewDirectory, $_);
  }
}

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

foreach my $SubNewFolder (@NewDirectory) {
  if ($SubNewFolder =~ /^\w.*_(\d.*)$/){
    CHECKING();
    unless (exists $Hashresults{$1}){
      if (exists $Birdresults{$1}){
        $parsedinput = "$in1/$SubNewFolder";
        $Mfolder = "$parsedinput/variant_output";
        my $verdict = PARSING($1,$parsedinput);
        #progress report
        if ($verdict == 1) {
          open (NOTE, ">>$progressnote");
          print NOTE "Subject: Update notes : $jobid\n\nCompleted library\t$1\n";
          system "sendmail $email < $progressnote"; close NOTE;
        } #end if
      }
      else {
        print "\nSkipping \"library_$1\" in \"$SubNewFolder\" folder because it isn't in birdbase\n$mystderror\n";
      } #end if
    }else {print "\nLibrary => $1 exists in the database\n";} #end unless
  } #end if	 
} #end foreach
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
sub DELETENOTDONE {
  print "\n\tDELETING NOT DONE\n";
  #CHECKING TO MAKE SURE NOT "done" FILES ARE REMOVED
  $syntax = "select library_id from transcripts_summary where status is NULL";
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
    #DELETE FROM transcripts_summary
    $sth = $dbh->prepare("delete from transcripts_summary where library_id in ( $incompletes )"); $sth->execute();
  }
}
sub CHECKING {
  #CHECKING THE LIBRARIES ALREADY IN THE DATABASE
  $syntax = "select library_id from transcripts_summary";
  $sth = $dbh->prepare($syntax);
  $sth->execute or die "SQL Error: $DBI::errstr\n";
  my $number = 0;
  while ($row = $sth->fetchrow_array() ) {
    $Hashresults{$row} = $number; $number++;
  }
  $syntax = "select library_id,date from bird_libraries";
  $sth = $dbh->prepare($syntax);
  $sth->execute or die "SQL Error: $DBI::errstr\n";
  $number = 0;
  while (my ($row1, $row2) = $sth->fetchrow_array() ) {
    $Birdresults{$row1} = $row2;
  }
}
sub PARSING {
  $dbh = mysql();
  print "\n\tINSERTING TRANSCRIPTS INTO THE DATABASE : \t library_$_[0]\n\n";
  $lib_id = $_[0]; my $librarydir = $_[1]; my $verdict = 0;
  my @foldercontent = split("\n", `find $librarydir`); #get details of the folder
  foreach (grep /\.gtf/, @foldercontent) { unless (`head -n 3 $_ | wc -l` <= 0 && $_ =~ /skipped/) { $transcriptsgtf = $_; } }
	$accepted = (grep /accepted_hits.bam/, @foldercontent)[0];
  $alignfile = (grep /summary.txt/, @foldercontent)[0];
  $genesfile = (grep /genes.fpkm/, @foldercontent)[0];
  $isoformsfile = (grep /isoforms.fpkm/, @foldercontent)[0];
  $deletionsfile = (grep /deletions.bed/, @foldercontent)[0];
  $insertionsfile = (grep /insertions.bed/, @foldercontent)[0];
  $junctionsfile = (grep /junctions.bed/, @foldercontent)[0];
  $run_log = (grep /logs\/run.log/, @foldercontent)[0];
	$samfile = (grep /.sam$/, @foldercontent)[0];
  $variantfile = (grep /.vcf$/, @foldercontent)[0]; 
  $vepfile = (grep /.vep.txt$/, @foldercontent)[0];
  $annofile = (grep /anno.txt$/, @foldercontent)[0];
  $htseqcount = (grep /.counts$/, @foldercontent)[0];
  $prepfile = (grep /prep_reads.info/, @foldercontent)[0];

  #$top_folder = "$second/tophat_out"; $cuff_folder = "$second/cufflinks_out"; $htseq_folder = "$second/htseq_output";
  
  #created a log file check because
  #I'm having issues with the authenticity of the results caz some of the log files is different
  my $logfile = `head -n 1 $run_log`;
  my @allgeninfo = split('\s',$logfile);

  #also getting metadata info
  my ($str, $ann, $ref, $seq,$allstart, $allend) = 0; #defining calling variables
  #making sure I'm working on only the chicken files for now, need to find annotation of alligator
  my @checkerno = split('\/',$allgeninfo[$#allgeninfo]);
  my @numberz =split('_', $checkerno[$#checkerno]);
  if ($numberz[0] == $lib_id){
    #making sure the arguments are accurately parsed
    if ($allgeninfo[1] =~ /.*library-type$/ && $allgeninfo[3] =~ /.*no-coverage-search$/){$str = 2; $ann = 5; $ref = 10; $seq = 11; $allstart = 4; $allend = 7;}
    elsif ($allgeninfo[1] =~ /.*library-type$/ && $allgeninfo[3] =~ /.*G$/ ){$str = 2; $ann = 4; $ref = 9; $seq = 10; $allstart = 3; $allend = 6;}
    elsif($allgeninfo[3] =~ /\-o$/){$str=99; $ann=99; $ref = 5; $seq = 6; $allstart = 3; $allend = 6;}
    else {print "File format doesn't match what was encoded for frnak_metadata\nSyntax error:\n\t@allgeninfo\n$mystderror\n";next;}
    #assuring we are working on available genomes
    my $refgenome = (split('\/', $allgeninfo[$ref]))[-1]; #reference genome name
    if (exists $parsablegenomes{$refgenome} || $refgenome =~ /Galgal4/){
      open(ALIGN, "<$alignfile") or die "Can't open file $alignfile\n";
      open(GENES, "<$genesfile") or die "Can't open file $genesfile\n";
      open(ISOFORMS, "<$isoformsfile") or die "Can't open file $isoformsfile\n";
    
      # PARSER FOR transcripts_summary TABLE
      while (<ALIGN>){
        chomp;
        if (/Input/){my $line = $_; $line =~ /Input.*:\s+(\d+)$/;$total = $1;}
        if (/Mapped/){my $line = $_; $line =~ /Mapped.*:\s+(\d+).*$/;$mapped = $1;}
      } close ALIGN;
      $unmapped = $total-$mapped;
      $deletions = `cat $deletionsfile | wc -l`; $deletions--;
      $insertions = `cat $insertionsfile | wc -l`; $insertions--;
      $junctions = `cat $junctionsfile| wc -l`; $junctions--;
      $genes = `cat $genesfile| wc -l`; $genes--;
      $isoforms = `cat $isoformsfile | wc -l`; $isoforms--;
      $prep = `cat $prepfile`;
      $date = `date +%Y-%m-%d`;

      #PARSING FOR SNPanalysis
      @parse = split('\/\/',$accepted); $accepted = undef; $len = $#parse+1; foreach(@parse){$accepted .= $_; if($len>1){$accepted .="\/"; $len--;}};
      @parse = split('\/\/',$run_log); $run_log = undef; $len = $#parse+1; foreach(@parse){$run_log .= $_; if($len>1){$run_log .="\/"; $len--;}};
    
      #INSERT INTO DATABASE : transcriptatlas
      #transcripts_summary table
      $sth = $dbh->prepare("insert into transcripts_summary (library_id, total_reads, mapped_reads, unmapped_reads, deletions, insertions, junctions, isoforms, genes, info_prep_reads, date ) values (?,?,?,?,?,?,?,?,?,?,?)");
      $sth ->execute($lib_id, $total, $mapped, $unmapped, $deletions, $insertions, $junctions, $isoforms, $genes, $prep, $date);
      
      #GENES_FPKM table
      $sth = $dbh->prepare("insert into genes_fpkm (library_id, tracking_id, class_code, nearest_ref_id, gene_id, gene_short_name, tss_id, chrom_no, chrom_start, chrom_stop, length, coverage, fpkm, fpkm_conf_low, fpkm_conf_high, fpkm_status ) values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
      while (<GENES>){
        chomp;
        my ($track, $class, $ref_id, $gene, $gene_name, $tss, $locus, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat ) = split /\t/;
        unless ($track eq "tracking_id"){ #check & specifying undefined variables to null
          if($class =~ /-/){$class = undef;} if ($ref_id =~ /-/){$ref_id = undef;}
          if ($length =~ /-/){$length = undef;} if($coverage =~ /-/){$coverage = undef;}
            $locus =~ /^(.+)\:(.+)\-(.+)$/;
            $chrom_no = $1; $chrom_start = $2; $chrom_stop = $3;
            $sth ->execute($lib_id, $track, $class, $ref_id, $gene, $gene_name, $tss, $chrom_no, $chrom_start, $chrom_stop, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat );
        }
      } close GENES;

      #ISOFORMS_FPKM table
      $sth = $dbh->prepare("insert into isoforms_fpkm (library_id, tracking_id, class_code, nearest_ref_id, gene_id, gene_short_name, tss_id, chrom_no, chrom_start, chrom_stop, length, coverage, fpkm, fpkm_conf_low, fpkm_conf_high, fpkm_status ) values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
      while (<ISOFORMS>){
        chomp;
        my ($track, $class, $ref_id, $gene, $gene_name, $tss, $locus, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat ) = split /\t/;
        unless ($track eq "tracking_id"){
          if ($class =~ /-/){$class = undef;} if ($ref_id =~ /-/){$ref_id = undef;}
          if ($length =~ /-/){$length = undef;} if($coverage =~ /-/){$coverage = undef;}
          $locus =~ /^(.+)\:(.+)\-(.+)$/;
          $chrom_no = $1; $chrom_start = $2; $chrom_stop = $3;
          $sth ->execute($lib_id, $track, $class, $ref_id, $gene, $gene_name, $tss, $chrom_no, $chrom_start, $chrom_stop, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat );
        }
      } close ISOFORMS;
  
      #extracting the syntax used for METADATA 
      #temporary replace the file path
      my ($replace, $temptest,$stranded, $sequences, $annotation, $annotationfile, $annfileversion);
      unless ($ann ==99){
        $temptest = `head -n 1 $allgeninfo[$ann]`;
        if (length($temptest) < 1 ){
          $replace = "/home/modupe/.big_ten"; foreach my $change (0..$#allgeninfo){if ($allgeninfo[$change] =~ /^(.*big_ten).*$/){$allgeninfo[$change] =~ s/$1/$replace/g;}}
        }
        if ($refgenome =~ /mouse/) {
          foreach (0..2){ #adding the first commands for allignment
            $executemouse .= $allgeninfo[$_]." ";
          }
          $executemouse .= "--no-coverage-search ";
          foreach ($allstart..$allend){
            $executemouse .= "$allgeninfo[$_] ";
          }
          $executemouse .= "-o $STORAGEPATH/library_".$lib_id." $GENOMES/$refgenome/$refgenome ";
        }
        $annotation = $allgeninfo[$ann];
        $annotationfile = uc ( (split('\.',((split("\/", $allgeninfo[$ann]))[-1])))[-1] ); #(annotation file)
        $annfileversion = substr(`head -n 1 $allgeninfo[$ann]`,2,-1); #annotation file version
      }
      else { $annotation = undef; $annotationfile = undef; $annfileversion = undef; }
      if ($str == 99){ $stranded = undef; }else { $stranded = $allgeninfo[$str]; } # (stranded or not)	
      my $otherseq = $seq++;
      unless(length($allgeninfo[$otherseq])<1){ #sequences 
        if ($refgenome =~ /mouse/) {
          my $seqtest = `ls $allgeninfo[$seq]`; # to check if the file is zipped or not
          if (length($seqtest) < 1 ){$executemouse .= $allgeninfo[$seq].".gz ".$allgeninfo[$otherseq].".gz";}
          else {$executemouse .= "$allgeninfo[$seq] $allgeninfo[$otherseq]";}
        }
        $sequences = ( ( split('\/', $allgeninfo[$seq]) ) [-1]).",". ( ( split('\/', $allgeninfo[$otherseq]) ) [-1]);
      } else {
        if ($refgenome =~ /mouse/) {
          my $seqtest = `ls $allgeninfo[$seq]`; # to check if the file is zipped or not
          if (length($seqtest) < 1 ){$executemouse .= $allgeninfo[$seq];}
          else {$executemouse .= $allgeninfo[$seq].".gz";}
        }
        $sequences = ( ( split('\/', $allgeninfo[$seq]) ) [-1]);}
   
      #frnak_metadata table
      $sth = $dbh->prepare("insert into frnak_metadata (library_id,ref_genome, ann_file, ann_file_ver, stranded, sequences,user ) values (?,?,?,?,?,?,?)");
      $sth ->execute($lib_id, $refgenome, $annotationfile, $annfileversion, $stranded,$sequences,"from raven" );
      
      # -----------------------
      #not used anymore
      #htseq table: Has to be created on the fly for different genomes; NOT WORKING (08/20/2015) too many columns
      #new idea, like genes_fpkm table.
      HTSEQ($htseqcount);
      # ----------------------
      
      #variant analysis
      if ($refgenome =~ /Galgal4/){$refgenome = "chicken";}
      VARIANTS($lib_id, $accepted, $refgenome, $annotation);

      #Finally : the last update. transcripts_summary table updating status column with 'done'
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
      print " rather what you have is $allgeninfo[$ref]\n$mystderror\n";
    }
  }
  else {
    print "library_id dont match $numberz[0] == $lib_id\n";
    $verdict = 0;
  }
  return $verdict;
}
sub HTSEQ {
  print "\n\tSTORING HTSEQ IN THE DATABASE\n\n";
# too many columns, I have to figure out another way to import this table: maybe NoSQL.
  #my %HTSEQ = undef; my $syntaxforhtseq = undef; my %dbtables = undef;
  open(HTSEQ, "<$_[0]") or die "Can't open file $_[0]\n";
  $sth = $dbh->prepare("insert into htseq ( library_id, gene_name, count) values (?,?,?)");
  while (<HTSEQ>){
    chomp;
    my ($NAME, $VALUE) = split /\t/;
    if ($NAME =~ /^[a-z0-9A-Z]/i) {
      $sth->execute($lib_id, $NAME, $VALUE);
      #$HTSEQ{uc($NAME)} = $VALUE;
    }
  } close HTSEQ;
}
sub VARIANTS{
  print "\n\tWORKING ON VARIANT ANALYSIS\n\n";
  my $libraryNO = "library_".$_[0]; my $bamfile = $_[1]; 
  my $REF= "$GENOMES/$_[2]/$_[2]".".fa";
  my $specie = $_[2];
  my $geneLIST = "$GENOMES/$_[2]/$_[2]".".txt";
  my $DICT = "$GENOMES/$_[2]/$_[2]".".dict";
  my $ANN = $_[3];
  
  $Mfolder = "$STORAGEPATH/$libraryNO"; `mkdir -p $Mfolder`; print "made $Mfolder\n"; #decided to keep all the variant folders.
  
  # MAYBE RUN ALIGNMENT PROCESS BECAUSE OF DOWNSTREAM CONFLICT -- but for only mouse genome
  if ($specie =~ /mouse/) {
    `$executemouse`;
    $bamfile = "$Mfolder/accepted_hits.bam"
  }
  
  #VARIANT ANALYSIS
  #PICARD
  my $filename = "$Mfolder/$libraryNO"."_DP5.vcf";
  unless (-e $filename) {
	$filename = "$Mfolder/$libraryNO".".bam";
	unless (-e $filename){
		`java -jar $PICARDDIR SortSam INPUT=$bamfile OUTPUT=$filename SO=coordinate`;
  	}else { print "$filename exists\n"; }
  
	#ADDREADGROUPS
	$filename = "$Mfolder/$libraryNO"."_add.bam";
  	unless (-e "$filename"){
		my $addreadgroup = "java -jar $PICARDDIR AddOrReplaceReadGroups INPUT=$Mfolder/$libraryNO".".bam OUTPUT=$filename SO=coordinate RGID=LAbel RGLB=Label RGPL=illumina RGPU=Label RGSM=Label";
  		`$addreadgroup`;
  	}else { print "$Mfolder/$libraryNO"."_add.bam exists\n "; }
  	#MARKDUPLICATES
  	unless (-e "$Mfolder/".$libraryNO."_mdup.bam" ) {
		my $markduplicates = "java -jar $PICARDDIR MarkDuplicates INPUT=$Mfolder/".$libraryNO."_add.bam OUTPUT=$Mfolder/".$libraryNO."_mdup.bam M=$Mfolder/".$libraryNO."_mdup.metrics CREATE_INDEX=true";
  		`$markduplicates`;
  	} else { print "$Mfolder/$libraryNO"."_mdup.bam exists\n"; }
  
	#SPLIT&TRIM
  	unless (-e "$Mfolder/".$libraryNO."_split.bam" ) {  
		my $splittrim = "java -jar $GATKDIR -T SplitNCigarReads -R $REF -I $Mfolder/".$libraryNO."_mdup.bam -o $Mfolder/".$libraryNO."_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 --filter_reads_with_N_cigar";
  		`$splittrim`;
  	} else { print "$Mfolder/$libraryNO"."_split.bam exists \n"; }
  
	#GATK
  	my $gatk = "java -jar $GATKDIR -T HaplotypeCaller -R $REF -I $Mfolder/".$libraryNO."_split.bam -o $Mfolder/$libraryNO.vcf";
  	`$gatk`;
  
  	#perl to select DP > 5 & get header information
  	FILTERING($Mfolder, "$Mfolder/$libraryNO.vcf");
  } else {print "Variants VCF already created\n";}

  #ANNOTATIONS : running VEP
  print "this is the species $specie\n"; #Remove this Modupe
  if ($specie =~ /alligator/) {
    my $snpdat = "perl $SNPdat -i $Mfolder/".$libraryNO."_DP5.vcf -f $REF -g $ANN"; 
    `$snpdat`;
            
    #DATABASE INSERT
    ALLIGATOR($geneLIST,$ANN);
    DBSNPDATVARIANTS($Mfolder."/".$libraryNO."_DP5.vcf", $Mfolder."/".$libraryNO."_DP5.vcf.output", $libraryNO);
    
  } elsif (exists $VEPparse{$specie}){
    #my $veptxt = "perl $VEP -i $Mfolder/".$libraryNO."_DP5.vcf --fork 24 --species $specie  --dir /home/modupe/.vep/ --cache --merged --everything on --terms ensembl --output_file $Mfolder/".$libraryNO."_VEP.txt"; `$veptxt`;
    my $vepvcf = "perl $VEP -i $Mfolder/".$libraryNO."_DP5.vcf --fork 24 --species $specie  --dir /home/modupe/.vep/ --cache --vcf --merged --everything on --terms ensembl --output_file $Mfolder/".$libraryNO."_VEP.vcf"; `$vepvcf`;
        
    #DATABASE INSERT
    DBVARIANTS($Mfolder."/".$libraryNO."_VEP.vcf", $libraryNO);
    
  } else {
    next "Unidentified genome\t$mystderror\n";
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
sub DBVARIANTS{	
  print "\n\tINSERTING VARIANTS INTO THE DATABASE\n\n";
  #disconnecting and connecting again to database just incase
  $dbh->disconnect(); 
  $dbh = mysql();

  $_[1] =~ /^library_(\d*)$/;
  my $libnumber = "$1";
  my $folder = undef;

  #VEP file
  my @splitinput = split('\/', $_[0]);
  foreach my $i (0..$#splitinput-1){$folder.="$splitinput[$i]/";$i++;}
  my $information = fileparse($_[0], qr/(\.vcf)?$/);
  
  #variant_metadata
  my $gatk_version = ( ( split('\/',$GATKDIR)) [-2] ); my $vep_version = ( ( split('\/',$VEP)) [-4] ); my $picard_version = ( ( split('\/',$PICARDDIR)) [-2] );

  #OUTPUT file
  my $output = "$folder$information".".v1table";
  my $output2 = "$folder$information".".v2table";
  open(OUTDBVAR,">$output");
  open(OUTDBVAR2,">$output2");
  
  print OUTDBVAR "#CHR\tPOS\tREF\tALT\tQUAL\tCLASS\t";
  print OUTDBVAR "ZYGOSITY\tdbSNP\n";
  
  print OUTDBVAR2 "#CHR\tPOS\tCONQ\t";
  print OUTDBVAR2 "GENE\tSYMBOL\tTRANSCRIPT\t";
  print OUTDBVAR2 "FEATURE\tTYPE OF GENE\tPROTEIN POSITION\t";
  print OUTDBVAR2 "AA CHANGE\tCODON CHANGE\n";
  my $date = `date +%Y-%m-%d`;

  #initializing the hash tables . . .
  %VCFhash = ();
  %VEPhash = ();
  %ExtraWork = ();
  #running through subroutines . . . 
  VEPVARIANT($_[0]);
  my ($itsnp,$itindel,$itvariants) = 0;
  #printing to output table & variant_results table
  #VARIANT_SUMMARY
  $sth = $dbh->prepare("insert into variants_summary ( library_id, ANN_version, Picard_version, GATK_version, date ) values (?,?,?,?,?)");
  $sth ->execute($libnumber, $vep_version, $picard_version, $gatk_version, $date);
  #VARIANT_RESULTS
  foreach my $abc (sort keys %VCFhash) {
    foreach my $def (sort {$a <=> $b} keys %{ $VCFhash{$abc} }) {
      my @vcf = split('\|', $VCFhash{$abc}{$def});
      my @ext = split('\|', $ExtraWork{$abc}{$def});
      if ($vcf[3] =~ /,/){
        my $first = split(",",$vcf[1]);
        if (length $vcf[0] == length $first){
          $itvariants++; $itsnp++;
        }
        else {
          $itvariants++; $itindel++;
        }
      }
      elsif (length $vcf[0] == length $vcf[1]){
        $itvariants++; $itsnp++;
      }
      else {
        $itvariants++; $itindel++;
      }
              
      print OUTDBVAR "$abc\t$def\t$vcf[0]\t$vcf[1]\t$vcf[2]\t$ext[0]\t";
      print OUTDBVAR "$vcf[3]\t$ext[1]\n";
              
      #to variant_result
      $sth = $dbh->prepare("insert into variants_result ( library_id, chrom, position, ref_allele, alt_allele, quality, variant_class,
                           zygosity, existing_variant ) values (?,?,?,?,?,?,?,?,?)");
      $sth ->execute($libnumber, $abc, $def, $vcf[0], $vcf[1], $vcf[2], $ext[0], $vcf[3], $ext[1]);
    }
  }	
  close (OUTDBVAR);	
  foreach my $abc (sort keys %VEPhash) {
    foreach my $def (sort {$a <=> $b} keys %{ $VEPhash{$abc} }) {
      foreach my $ghi (sort keys %{ $VEPhash{$abc}{$def} }) {
        foreach my $jkl (sort keys %{ $VEPhash{$abc}{$def}{$ghi} }) {
          foreach my $mno (sort keys %{ $VEPhash{$abc}{$def}{$ghi}{$jkl} }){
            my @vep = split('\|', $VEPhash{$abc}{$def}{$ghi}{$jkl}{$mno});
            if(length($vep[4]) < 1){$vep[4] = "-";}
            if (length($vep[1]) < 1) {$vep[1] = "-";}
            print OUTDBVAR2 "$abc\t$def\t$ghi\t";
            print OUTDBVAR2 "$vep[1]\t$vep[8]\t$vep[2]\t";
            print OUTDBVAR2 "$vep[0]\t$vep[7]\t$vep[4]\t";
            print OUTDBVAR2 "$vep[5]\t$vep[6]\n";
            
            #to variants_annotation
            $sth = $dbh->prepare("insert into variants_annotation ( library_id, chrom, position, consequence, gene_id, gene_name,
                                transcript, feature, gene_type,protein_position, aminoacid_change, codon_change ) values
                                (?,?,?,?,?,?,?,?,?,?,?,?)");
            $sth ->execute($libnumber, $abc, $def, $ghi, $vep[1], $vep[8], $vep[2], $vep[0], $vep[7], $vep[4], $vep[5], $vep[6]);
          }
        }
      }
    }
  }
  close (OUTDBVAR2);
  
  #VARIANT_SUMMARY
  $syntax = "update variants_summary set total_VARIANTS = $itvariants, total_SNPS = $itsnp, 
  total_INDELS = $itindel, status = \'done\' where library_id like \"$libnumber\"";
  $sth = $dbh->prepare($syntax);
  $sth ->execute();
}
sub VEPVARIANT {
  #working on VEP variants
  my %Geneinfo = ''; my %Transcriptinfo = ''; my %Conqinfo = ''; my %Varclass = '';
  my %Featureinfo = ''; my %Proinfo = ''; my %Prochangeinfo = ''; my %Codoninfo = '';
  my %dbSNPinfo = ''; my %locinfo = ''; my %Entrezinfo = ''; my %GENEtype = '';
  my %GENEname = ''; my $position; my %location = '';
  unless(open(FILE,$_[0])){print "File \'$_[0]\' doesn't exist\n";exit;}
  my $verd;
  my @file = <FILE>;
  chomp @file;
  close (FILE);
  foreach my $chr (@file){
    unless ($chr =~ /^#/){
      my @chrdetails = split('\t', $chr);
      
      #removing the random chromosomes (scaffolds) - because no vital information can be found for them.
      my @morechrsplit = split(';', $chrdetails[7]);
      if (((split(':', $chrdetails[9]))[0]) eq '0/1'){$verd = "heterozygous";}
      elsif (((split(':', $chrdetails[9]))[0]) eq '1/1'){$verd = "homozygous";}
      elsif (((split(':', $chrdetails[9]))[0]) eq '1/2'){$verd = "heterozygous alternate";}
      
      #Processing the VEP section
      my @finalchrsplit = split("\,",(((split('=',((split(';',$chrdetails[7]))[-1])))[1]))); 
      foreach my $FCR (0..$#finalchrsplit){
        my @vepdetails = split('\|', $finalchrsplit[$FCR]);
        if ($vepdetails[1] !~ /WITHIN_NON_CODING_GENE/){
      
          #VCFhash information
          $VCFhash{$chrdetails[0]}{$chrdetails[1]} = "$chrdetails[3]|$chrdetails[4]|$chrdetails[5]|$verd";
      
          if ($vepdetails[14] < 1) {
            $vepdetails[14] = "-";
          }
          $location{"$chrdetails[0]|$chrdetails[1]|$vepdetails[1]|$vepdetails[14]|$vepdetails[4]"} = "$chrdetails[0]|$chrdetails[1]|$vepdetails[1]|$vepdetails[14]|$vepdetails[4]"; #specifying location with consequence.
          
          #GENE - 1
          if (exists $Geneinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}){
            unless ($Geneinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} eq $vepdetails[4]){
              print "something is wrong with the code >> GENE";
              print "\na $chrdetails[0]\tb $chrdetails[1]\tc $vepdetails[14]\td $vepdetails[4]\te $vepdetails[4]\tf $Geneinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}\n";		
            }
          }
          else {$Geneinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $vepdetails[4];}
          
          #TRANSCRIPT - 2
          if (exists $Transcriptinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}){
            unless ($Transcriptinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} eq $vepdetails[6]){
              my $temp1details = "$Transcriptinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]},$vepdetails[6]";
              $Transcriptinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $temp1details;
            }
          }
          else {$Transcriptinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $vepdetails[6];}
          
          #CONSEQUENCE - 3
          if (exists $Conqinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} ){
            unless ($Conqinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} eq $vepdetails[1]){
              print "something is wrong with the code >> consequence";
              print "\na $chrdetails[0]\tb $chrdetails[1]\tc $vepdetails[14]\td $vepdetails[4]\te $vepdetails[1]\tf $Conqinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}\n";		
              #exit;
            }
          }
          else {$Conqinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $vepdetails[1];}
          
          #VARIANT CLASS - 4
          if (exists $Varclass{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}){
            unless ($Varclass{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} eq $vepdetails[20]){
              print "something is wrong with the code >> variantclass";
              print "\na $chrdetails[0]\tb $chrdetails[1]\tc $vepdetails[14]\td $vepdetails[4]\te $vepdetails[20]\tf $Varclass{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}\n";		
              #exit;
            }
          }
          else {$Varclass{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $vepdetails[20];}
         
          #FEATURE TYPE - 5
          if (exists $Featureinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}){
            unless ($Featureinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} eq $vepdetails[5]){
              print "something is wrong with the code >> feature type";
              print "\na $chrdetails[0]\tb $chrdetails[1]\tc $vepdetails[14]\td $vepdetails[4]\te $vepdetails[5]\tf $Featureinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}\n";		
            }
          }
          else {$Featureinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $vepdetails[5];}
          
          #PROTEIN POSITION - 6
          if (exists $Proinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}){
            unless ($Proinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} eq $vepdetails[14]){
              print "something is wrong with the code >> protein position";
              print "\na $chrdetails[0]\tb $chrdetails[1]\tc $vepdetails[14]\td $vepdetails[4]\te $vepdetails[14]\tf $Proinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}\n";		
            }		
          }
          else {$Proinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $vepdetails[14];}
         
          #PROTEIN CHANGE - 7 (or Amino Acid)
          if (exists $Prochangeinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}){
            unless ($Prochangeinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} eq $vepdetails[15]){
              my $temp2details = "$Prochangeinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]},$vepdetails[15]";
              $Prochangeinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $temp2details;
            }
          }
          else {$Prochangeinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $vepdetails[15];}
          
          #CODON CHANGE - 8
          if (exists $Codoninfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}){
            unless ($Codoninfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} eq $vepdetails[16]){
              my $temp3details = "$Codoninfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]},$vepdetails[16]";
              $Codoninfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $temp3details;
            }		
          }
          else {$Codoninfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $vepdetails[16];}
          
          #dbSNP - 9
          if (exists $dbSNPinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}){
            unless ($dbSNPinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} eq $vepdetails[17]){
              print "something is wrong with the code >> dbSNP";
              print "\na $chrdetails[0]\tb $chrdetails[1]\tc $vepdetails[14]\td $vepdetails[4]\te $vepdetails[17]\tf $dbSNPinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}\n";		
            }		
          }
          else {$dbSNPinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $vepdetails[17];}
          
          #GENE name - 10
          if (exists $GENEname{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}){
            unless ($GENEname{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} eq $vepdetails[3]){
              print "something is wrong with the code >> genename";
              print "\na $chrdetails[0]\tb $chrdetails[1]\tc $vepdetails[14]\td $vepdetails[4]\te $vepdetails[3]\tf $GENEname{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}\n";		
            }	
          }
          else {$GENEname{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $vepdetails[3];}
          
          #GENE type - 11
          if (exists $GENEtype{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}){
            unless ($GENEtype{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} eq $vepdetails[7]){
              my $temp4details = "$GENEtype{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]},$vepdetails[7]";
              $GENEtype{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $temp4details;
            }	
          }
          else {$GENEtype{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $vepdetails[7];}
        }
      }
    }
  }

  foreach my $alldetails (keys %location){
        my ($chrdetails1, $chrdetails2, $consequences, $pposition, $genename) = split('\|', $alldetails);
        #cleaning up the text
        my $clean1 = CLEANUP($Varclass{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename});
        my $clean2 = CLEANUP($Featureinfo{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename});
        my $clean3 = CLEANUP($Geneinfo{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename});
        my $clean4 = CLEANUP($Transcriptinfo{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename});
        my $clean5 = CLEANUP($Conqinfo{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename});
        my $clean6 = CLEANUP($Proinfo{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename});
        my $clean7 = CLEANUP($Prochangeinfo{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename});
        my $clean8 = CLEANUP($Codoninfo{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename});
        my $clean9 = CLEANUP($dbSNPinfo{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename});
        my $clean10 = CLEANUP($GENEtype{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename});
        my $clean11 = CLEANUP($GENEname{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename});
        $VEPhash{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename} = "$clean2|$clean3|$clean4|$clean5|$clean6|$clean7|$clean8|$clean10|$clean11";
        $ExtraWork{$chrdetails1}{$chrdetails2} = "$clean1|$clean9";
  }   
}
sub DBSNPDATVARIANTS{	
  print "\n\tINSERTING SNPDAT VARIANTS INTO THE DATABASE\n\n";
  #disconnecting and connecting again to database just incase
  $dbh->disconnect(); 
  $dbh = mysql();

  $_[2] =~ /^library_(\d*)$/;
  my $libnumber = $1;
  my $folder = undef;

  #VCF file
  my @splitinput = split('\/', $_[1]);
  foreach my $i (0..$#splitinput-1){$folder.="$splitinput[$i]/";$i++;}
  my $information = fileparse($_[0], qr/(\.vcf)?$/);
  
  #variant_metadata
   my $gatk_version = ( ( split('\/',$GATKDIR)) [-2] ); my $ann_version = ( ( split('\/',$SNPdat)) [-2] ); my $picard_version = ( ( split('\/',$PICARDDIR)) [-2] );

  #OUTPUT file
  my $output = "$folder$information".".v1table";
  my $output2 = "$folder$information".".v2table";
  
  open(OUTDBVAR,">$output");
  open(OUTDBVAR2,">$output2");
  
  print OUTDBVAR "#CHR\tPOS\tREF\tALT\tQUAL\tCLASS\t";
  print OUTDBVAR "ZYGOSITY\tdbSNP\n";
  
  print OUTDBVAR2 "#CHR\tPOS\tCONQ\t";
  print OUTDBVAR2 "GENE\tSYMBOL\tTRANSCRIPT\t";
  print OUTDBVAR2 "FEATURE\tTYPE OF GENE\tPROTEIN POSITION\t";
  print OUTDBVAR2 "AA CHANGE\tCODON CHANGE\n";
  my $date = `date +%Y-%m-%d`;

  #initializing the hash tables . . .
  %VCFhash = ();
  %VEPhash = ();
  %ExtraWork = ();
  #running through subroutines . . . 
  SNPdatVARIANT($_[0], $_[1]);
  my ($itsnp,$itindel,$itvariants) = 0;
  #printing to output table & variant_results table
  #VARIANT_SUMMARY
  $sth = $dbh->prepare("insert into variants_summary ( library_id, ANN_version, Picard_version, GATK_version, date ) values (?,?,?,?,?)");
  $sth ->execute($libnumber, $ann_version, $picard_version, $gatk_version, $date);

  #VARIANT_RESULTS
  foreach my $abc (sort keys %VCFhash) {
    foreach my $def (sort {$a <=> $b} keys %{ $VCFhash{$abc} }) {
      my @vcf = split('\|', $VCFhash{$abc}{$def});
      my $variantclass;
      my $first = split(",",$vcf[1]);
      if ($vcf[3] =~ /,/){                   
        if (length($vcf[0]) == length($first)){
          $itvariants++; $itsnp++; $variantclass = "SNV";
        }
        elsif (length($vcf[0]) > length($first)) {
          $itvariants++; $itindel++; $variantclass = "deletion";
        }
        elsif (length($vcf[0]) < length($first)) {
          $itvariants++; $itindel++; $variantclass = "insertion";
        }
        else {
          print "first False variant, cross-check asap"; exit; 
        }
      }
      elsif (length $vcf[0] == length $vcf[1]){
        $itvariants++; $itsnp++; $variantclass = "SNV";
      }
      elsif (length $vcf[0] > length ($vcf[1])) {
        $itvariants++; $itindel++; $variantclass = "deletion";
      }
      elsif (length $vcf[0] < length ($vcf[1])) {
        $itvariants++; $itindel++; $variantclass = "insertion";
      }
      else {
        print "second False variant, cross-check asap\n"; exit;
      }
      print OUTDBVAR "$abc\t$def\t$vcf[0]\t$vcf[1]\t$vcf[2]\t$variantclass\t";
      print OUTDBVAR "$vcf[3]\t$ExtraWork{$abc}{$def}\n";
              
      $sth = $dbh->prepare("insert into variants_result ( library_id, chrom, position, ref_allele, alt_allele, quality, variant_class,
                          zygosity, existing_variant ) values (?,?,?,?,?,?,?,?,?)");
      $sth ->execute($libnumber, $abc, $def, $vcf[0], $vcf[1], $vcf[2], $variantclass, $vcf[3], $ExtraWork{$abc}{$def});                  
    }
  }
  close (OUTDBVAR);
  
  foreach my $abc (sort keys %VEPhash) {
    foreach my $def (sort {$a <=> $b} keys %{ $VEPhash{$abc} }) {
      foreach my $ghi (sort keys %{ $VEPhash{$abc}{$def} }) {
          foreach my $jkl (sort keys %{ $VEPhash{$abc}{$def}{$ghi} }) {
            my @vep = split('\|', $VEPhash{$abc}{$def}{$ghi}{$jkl});
            
            if(length($vep[4]) < 1){
                $vep[4] = "-";
            }
            if (length($vep[1]) < 1) {
                $vep[1] = "-";
            }
              print OUTDBVAR2 "$abc\t$def\t$ghi\t";
              print OUTDBVAR2 "$vep[1]\t$vep[8]\t$vep[2]\t";
              print OUTDBVAR2 "$vep[0]\t$vep[7]\t$vep[4]\t";
              print OUTDBVAR2 "$vep[5]\t$vep[6]\n";
              
              #to variants_annotation
              $sth = $dbh->prepare("insert into variants_annotation ( library_id, chrom, position, consequence, gene_id, gene_name,
                                 transcript, feature, gene_type,protein_position, aminoacid_change, codon_change ) values
                                (?,?,?,?,?,?,?,?,?,?,?,?)");
              $sth ->execute($libnumber, $abc, $def, $ghi, $vep[1], $vep[8], $vep[2], $vep[0], $vep[7], $vep[4], $vep[5], $vep[6]);
          }
      }
    }
  }
  close (OUTDBVAR2);
  
  #VARIANT_SUMMARY
  $syntax = "update variants_summary set total_VARIANTS = $itvariants, total_SNPS = $itsnp, 
  total_INDELS = $itindel, status = \'done\' where library_id like \"$libnumber\"";
  $sth = $dbh->prepare($syntax);
  $sth ->execute();
}
sub SNPdatVARIANT {
  #working on VEP variants
  my %Geneinfo = ''; my %Transcriptinfo = ''; my %Conqinfo = '';
  my %Featureinfo = ''; my %Proinfo = ''; my %Prochangeinfo = ''; my %Codoninfo = '';
  my %dbSNPinfo = ''; my %locinfo = ''; my %Entrezinfo = ''; my %GENEtype = '';
  my %GENEname = ''; my $position; my %location = ''; my ($consqofgene,$idofgene, $codonchange, $aminoacid);
	
  unless(open(VAR,$_[0])){print "File \'$_[0]\' doesn't exist\n";exit;}
  my $verd; my @variantfile = <VAR>; chomp @variantfile; close (VAR);
  foreach my $chr (@variantfile){
    unless ($chr =~ /^#/){
      my @chrdetail = split('\t', $chr);
      #removing the random chromosomes (scaffolds) - because no vital information can be found for them.
      my @morechrsplit = split(';', $chrdetail[7]);
      if (((split(':', $chrdetail[9]))[0]) eq '0/1'){$verd = "heterozygous";}
      elsif (((split(':', $chrdetail[9]))[0]) eq '1/1'){$verd = "homozygous";}
      elsif (((split(':', $chrdetail[9]))[0]) eq '1/2'){$verd = "heterozygous almyternate";}
      #VCFhash information
      $VCFhash{$chrdetail[0]}{$chrdetail[1]} = "$chrdetail[3]|$chrdetail[4]|$chrdetail[5]|$verd";
    }
  }
  unless(open(FILE,$_[1])){print "File \'$_[0]\' doesn't exist\n";exit;}
  my @file = <FILE>; chomp @file; close (FILE); shift(@file);
  foreach my $chr2 (@file){
    #Processing the ANNotation section
    my @chrdetails = split('\t', $chr2);
    $chrdetails[0] = "\L$chrdetails[0]";
    if ($chrdetails[10] =~ /ID=AMISG/){
      #getting gene id and name of gene
      my @geneidlist = split("\;", $chrdetails[10]);
      foreach (@geneidlist){
        if($_ =~ /^ID=AMISG/){$idofgene = substr($_,3);}
        last;
      }
              
      #change codon format to be consistent
      unless ($chrdetails[19] eq "NA"){
        $chrdetails[19] = "\L$chrdetails[19]";
    
        my $aq = ((split('\[', $chrdetails[19],2))[0]); my $bq = ((split('\]', $chrdetails[19],2))[1]); 
        my $dq = uc ((split('\[', ((split('\/', $chrdetails[19],2))[0]),2)) [1]); 
        my $eq = uc ((split('\]',((split('\/', $chrdetails[19],2))[1]),2)) [0]);
        $codonchange = "$aq$dq$bq/$aq$eq$bq";
        $aminoacid = substr(substr($chrdetails[20],1),0,-1);
      }else {
        $codonchange = undef; $aminoacid = undef;
      }
      #determine the consequence of the SNP change
      if($chrdetails[21] eq "N") {
        if ($aminoacid =~ /^M/){$consqofgene = "START_LOST";}
        elsif($aminoacid =~ /M$/){$consqofgene = "START_GAINED";}
        elsif($aminoacid =~/^\-/){$consqofgene = "STOP_LOST";}
        elsif($aminoacid =~ /\-$/){$consqofgene = "STOP_GAINED";}	
        else {$consqofgene = "NON_SYNONYMOUS_CODING";}
      }
      elsif ($chrdetails[21] eq "Y"){
        $consqofgene = "SYNONYMOUS_CODING";
        $aminoacid = ((split('\/', $aminoacid,2))[0]);
      }
      elsif($chrdetails[21] eq "NA") {
        if ($chrdetails[3] =~ /Intergenic/){$consqofgene = "INTERGENIC";}
        elsif($chrdetails[3] =~ /Intronic/){$consqofgene = "INTRONIC";}
        elsif($chrdetails[3] =~ /Exonic/){$consqofgene = "EXON";}
        else {$consqofgene = "NA";}
      }
      else {
        print "There's a problem with the consequence of the gene\n\n"; exit;
      }
      foreach (0..$#chrdetails){if ($chrdetails[$_] eq "NA"){ $chrdetails[$_] = undef;}} #converting NA to undef
              
      $location{"$chrdetails[0]|$chrdetails[1]|$consqofgene|$idofgene"} = "$chrdetails[0]|$chrdetails[1]|$consqofgene|$idofgene"; #specifying location with consequence.
              
      #GENE - 1
      if (exists $Geneinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene}){
        unless ($Geneinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} eq $idofgene){
          print "something is wrong with the code >> Geneinfo\n";
          exit;
        }
      }
      else {	
        $Geneinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = $idofgene;
      }
      
      #TRANSCRIPT - 2
      if (exists $Transcriptinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene}){
        unless ($Transcriptinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} eq $AMIST{$idofgene}){
          print "something is wrong with the code >> transcript\n";
          exit;
        }
      }
      else {	
        $Transcriptinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = $AMIST{$idofgene};
      } 	
      
      #CONSEQUENCE - 3
      if (exists $Conqinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} ){
        unless ($Conqinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} eq $consqofgene){
          print "something is wrong with the code >> consequence\n";
          exit;
        }			
      }
      else {	
        $Conqinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = $consqofgene;
      }	
        
      #FEATURE TYPE - 5 #SNPdat
      if (exists $Featureinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene}){
        unless ($Featureinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} eq $chrdetails[5]){
          print "something is wrong with the code >> feature\n";
          exit;
        }			
      }
      else {	
        $Featureinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = $chrdetails[5];
      }

      #PROTEIN POSITION - 6 #SNPdat this isn't provided
      $Proinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = undef;

      #PROTEIN CHANGE - 7 (or Amino Acid) #SNPdat
      if (exists $Prochangeinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene}){
        unless ($Prochangeinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} eq $aminoacid){
          my $temp1details = "$Prochangeinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene},$aminoacid";
          $Prochangeinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = $temp1details;
        }			
      }
      else {	
        $Prochangeinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = $aminoacid;
      }
	
      #CODON CHANGE - 8 #SNPdat
      if (exists $Codoninfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene}){
        unless ($Codoninfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} eq $codonchange){
          my $temp2details = "$Codoninfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene},$codonchange";
          $Codoninfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = $temp2details;
        }			
      }
      else {	
        $Codoninfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = $codonchange;
      }
      
	#dbSNP - 9 #SNPdat # no dbSNP info for alligator
	$dbSNPinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = undef;			
	
	#GENE name - 10 #name of gene from the alligator gene list. #SNPdat
      if (exists $GENEname{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene}){
        unless ($GENEname{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} eq $AMISG{$idofgene}){
          print "something is wrong with the code >> gene name\n";
          exit;
        }
	}
	else {	
	  $GENEname{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = $AMISG{$idofgene};
	}
	
	#GENE type - 11 #SNPdat
      if (exists $GENEtype{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene}){
        unless ($GENEtype{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} eq $chrdetails[3]){
          print "something is wrong with the code >> gene type\n";
          exit;
        }
	}
	else {	
	  $GENEtype{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = $chrdetails[3];
	}
    }
  }

  foreach my $alldetails (keys %location){
    my ($chrdetails1, $chrdetails2, $consequences, $genename) = split('\|', $alldetails);
    #cleaning up the text
    my $clean2 = CLEANUP($Featureinfo{$chrdetails1}{$chrdetails2}{$consequences}{$genename});
    my $clean3 = CLEANUP($Geneinfo{$chrdetails1}{$chrdetails2}{$consequences}{$genename});
    my $clean4 = CLEANUP($Transcriptinfo{$chrdetails1}{$chrdetails2}{$consequences}{$genename});
    my $clean5 = CLEANUP($Conqinfo{$chrdetails1}{$chrdetails2}{$consequences}{$genename});
    my $clean6 = CLEANUP($Proinfo{$chrdetails1}{$chrdetails2}{$consequences}{$genename});
    my $clean7 = CLEANUP($Prochangeinfo{$chrdetails1}{$chrdetails2}{$consequences}{$genename});
    my $clean8 = CLEANUP($Codoninfo{$chrdetails1}{$chrdetails2}{$consequences}{$genename});
    my $clean9 = CLEANUP($dbSNPinfo{$chrdetails1}{$chrdetails2}{$consequences}{$genename});
    my $clean10 = CLEANUP($GENEtype{$chrdetails1}{$chrdetails2}{$consequences}{$genename});
    my $clean11 = CLEANUP($GENEname{$chrdetails1}{$chrdetails2}{$consequences}{$genename});
    $VEPhash{$chrdetails1}{$chrdetails2}{$consequences}{$genename} = "$clean2|$clean3|$clean4|$clean5|$clean6|$clean7|$clean8|$clean10|$clean11";
    $ExtraWork{$chrdetails1}{$chrdetails2} = $clean9;
  } 
}


sub CLEANUP {
  #cleaning up the VEP variants so that it doesn't have repetitions in the output
  my @unclean = split(',', $_[0]);
  my ($cleansyntax, %Hashdetails, %Hash2) = undef;
  foreach (0..$#unclean){
    if ($unclean[$_] =~ /^[a-zA-Z0-9]/){
      $Hashdetails{$unclean[$_]} = $_;
    }
  }
  foreach my $unique (keys %Hashdetails){
    if ($unique =~ /^[a-zA-Z0-9]/){
      $Hash2{$Hashdetails{$unique}} = $unique;
    }
  }
  foreach my $final (sort keys %Hash2){
    if ($Hash2{$final} =~ /^[a-zA-Z0-9]/){
      $cleansyntax .= ",$Hash2{$final}";
    }
  }
  my $returnclean = substr($cleansyntax, 1);
  return $returnclean;
}
sub ALLIGATOR {
  my ($theparent, $theid); 
  unless(open(FILE,$_[0])){exit "File \'$_[0]\' doesn't exist\n";}
  my @file = <FILE>;
  chomp @file; shift (@file);
  close (FILE);
  foreach my $line (@file){
    my @details = split(',', $line);
    $AMISG{$details[0]} = $details[1];
  } 
  unless(open(FILER,$_[1])){exit "File \'$_[1]\' doesn't exist\n";}
  my @filer = <FILER>;
  chomp @filer; shift (@filer);
  close (FILER);
  foreach (@filer){
    my @newdetails = split('\t', $_);
    if ($newdetails[8] =~ /Parent=AMISG/){
	my @whatIwant = split("\;", $newdetails[8]);
	foreach (@whatIwant){
	  if ($_ =~ /ID=/){
	    $theid = substr($_, 3);
        }
	  elsif ($_ =~ /Parent=/){
	    $theparent = substr($_, 7);
        }
	}
	$AMIST{$theparent} = $theid;
    }
  }    
}
sub SUMMARYstmts{
  #transcripts_summary from the database
  $dbh = mysql();
  print "\n\tEXECUTING SELECT STATEMENT ON THE DATABASE TABLES \n";
  $syntax = "select count(*) from transcripts_summary";
  $sth = $dbh->prepare($syntax);
  $sth->execute or die "SQL Error: $DBI::errstr\n";
  while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"transcripts_summary\" table \t:\t @row\n";}
  #genes_fpkm
  $syntax = "select count(*) from genes_fpkm";
  $sth = $dbh->prepare($syntax);
  $sth->execute or die "SQL Error: $DBI::errstr\n";
  while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"genes_fpkm\" table \t\t:\t @row\n";}
  #isoforms_fpkm
  $syntax = "select count(*) from isoforms_fpkm";
  $sth = $dbh->prepare($syntax);
  $sth->execute or die "SQL Error: $DBI::errstr\n";
  while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"isoforms_fpkm\" table \t:\t @row\n";}
  #variant_summary
  $syntax = "select count(*) from variants_summary";
  $sth = $dbh->prepare($syntax);
  $sth->execute or die "SQL Error: $DBI::errstr\n";
  while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"variants_summary\" table \t:\t @row\n";}
  #variant_list
  $syntax = "select count(*) from variants_result";
  $sth = $dbh->prepare($syntax);
  $sth->execute or die "SQL Error: $DBI::errstr\n";
  while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"variants_results\" table \t:\t @row\n";}
  #variant_annotion
  $syntax = "select count(*) from variants_annotation";
  $sth = $dbh->prepare($syntax);
  $sth->execute or die "SQL Error: $DBI::errstr\n";
  while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"variants_annotation\" table \t:\t @row\n";}
  #frnak_metadata
  $syntax = "select count(*) from frnak_metadata";
  $sth = $dbh->prepare($syntax);
  $sth->execute or die "SQL Error: $DBI::errstr\n";
  while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"frnak_metadata\" table \t:\t @row\n";}
}
sub NOTIFICATION {
  my $notification = '/home/modupe/.LOG/note.txt';
  open (NOTE, ">$notification");
  print NOTE "Subject: ". $_[0] .": $jobid\n\nName of log files\n\t$std_out\n\t$std_err\n";
  system "sendmail $email < $notification";
  close NOTE;
  system "rm -rf $notification";
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
close STDOUT; close STDERR;
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
exit;
