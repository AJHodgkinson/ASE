use strict; 
use Parallel::Loops;
use Getopt::Long;
use Pod::Usage;

#name: ase_normalisation_sim.pl
#description: A program to obtain reads and pairs for heterozygote positions in RNA data and simulate both haplotypes
#date: 22nd April 2014
#Written By: Alan Hodgkinson

###Defaults, Input and file checks

###Required:
my $exome_id = undef;
my $rna_id = undef;
my $exome_vcf = undef;
my $rna_file = undef;

my $maxProcs = 1; #set max processes to 1 by default
my $sample_depth = 2000; #set sample depth 2000 by default

my $add_vcf = undef;
my $add_id = undef;

my $cmd_help = 0;

GetOptions (
	    'DNAid=s' => \$exome_id,
	    'OUTid=s' => \$rna_id,
	    'DNAvcf=s' => \$exome_vcf,
	    'RNAbam=s' => \$rna_file,
	    'MaxProcs=i' => \$maxProcs,
	    'SampleDepth:i' => \$sample_depth,
	    'ADDid:s' => \$add_id,
	    'ADDvcf:s' => \$add_vcf,
	    'help' => \$cmd_help
) or exit(1);

pod2usage(1) if $cmd_help;

print "--DNAid $exome_id
--OUTid $rna_id
--DNAvcf $exome_vcf
--RNAbam $rna_file
--MaxProcs $maxProcs
--SampleDepth $sample_depth
--ADDid $add_id
--ADDvcf $add_vcf\n";

print "Option --DNAid required: Please supply name of Individual as defined in VCF file\n" if not defined ($exome_id);
print "Option --OUTid required: Please supply name Output Directory to write to\n" if not defined ($rna_id);
print "Option --DNAvcf required: Please supply name of VCF file\n" if not defined ($exome_vcf);
print "Option --RNAbam required: Please supply name of BAM file containing mapped RNA reads\n" if not defined ($rna_file);

if ((defined $add_vcf)&&(not defined $add_id)) {
  die "You have supplied an additional vcf file to add to the simulation, you must also supply the name of the ID in that file using --ADDid\n";
}

my $rna_dir = $rna_id."_simdir";
if (-d "$rna_dir") {
  die "Folder $rna_id already exists - please rename the --OUTid variable, or delete/rename the existing folder\n";
}
if (!(-d "$rna_dir")) {
  mkdir $rna_dir, 0755;
}

if ($maxProcs>8) {
  warn "Warning: Using more than 8 cpus on the same node will probably be very slow as it is likely that system calls to samtools within the simulation will overrun the node.  During testing we found that speed increased up to 8 cpus, but not beyond\n";
}

#Set output file names (no file checks since muts be new folder)

my $final = $rna_dir."/ref_allele_count_".$rna_id.".txt"; #Allele count output file
my $final_non = $rna_dir."/ref_allele_count_".$rna_id."_non.txt"; #Allele count output file

my $fastq1 = $rna_dir."/fastq_sim_".$rna_id."_R1.fastq";
my $fastq1gz = $rna_dir."/fastq_sim_".$rna_id."_R1.fastq.gz";
my $fastq2 = $rna_dir."/fastq_sim_".$rna_id."_R2.fastq";
my $fastq2gz = $rna_dir."/fastq_sim_".$rna_id."_R2.fastq.gz";

for (1..24) {
  my $f1_chr = $rna_dir."/fastq_sim_".$rna_id."_".$_."_R1.fastq";
  my $f2_chr = $rna_dir."/fastq_sim_".$rna_id."_".$_."_R2.fastq";
}

###MAIN PROGRAM

$SIG{INT} = sub { die "Caught a signal interuption $!" };

my %hets = (); #Het SNPs to use in normalisation
my %insertion = (); #insertions to insert
my %deletion = (); #deletions to remove
my %rna_hets = (); #optional RNA variation 
my @het_array = (); #array containing SNPs for normalisation
my %snp_allele_coverage = (); #hash to store the coverage of each heterozygous position in simulated data - should be at least 2000X 
my %snp_allele_coverage_non = (); #hash to store the coverage of each heterozygous position in simulated data (non overkapping pairs) - should be at least 2000X 


&CollectDNASNPs($exome_id,$exome_vcf); #Collect DNA Heterozygote Positions
if (defined $add_vcf) { #if extra vcf supplied
  &CollectRNASNPs($add_id,$add_vcf); #Collect Additional Heterozygote Positions
}

my $pl = Parallel::Loops->new($maxProcs);
$pl->share(\@het_array);
$pl->share(\%snp_allele_coverage);
$pl->share(\%snp_allele_coverage_non);
my @chrs = (1..24);
$pl->foreach( \@chrs, sub {
		my $chrs = $_; 
		&InsertSNPs($chrs);
	      });

for (my $chrs=1;$chrs<25;$chrs++) {
  my $fastq1_chr = $rna_dir."/fastq_sim_".$rna_id."_".$chrs."_R1.fastq";
  my $fastq2_chr = $rna_dir."/fastq_sim_".$rna_id."_".$chrs."_R2.fastq";
  system ("cat $fastq1_chr >> $fastq1");
  system ("cat $fastq2_chr >> $fastq2");
  system ("rm $fastq1_chr");
  system ("rm $fastq2_chr");
}
system ("pigz $fastq1");
system ("pigz $fastq2");

open (FINAL, ">$final") || die "Unable to open counts file to write to: $!\n";
print FINAL "#Position\tRef_alleles\n";
while (my ($a,$b) = each %snp_allele_coverage) {
  print FINAL "$a\t$b\n";
}

open (FINALN, ">$final_non") || die "Unable to open non-overlapping counts file to write to: $!\n";
print FINALN "#Position\tRef_alleles\n";
while (my ($a,$b) = each %snp_allele_coverage_non) {
  print FINALN "$a\t$b\n";
}

################Subroutines

sub CollectDNASNPs ($$) {
  my ($exome_id_local,$exome_vcf_local) = @_;
  my @head = (); my $size = 0;
  open (EXOME, "$exome_vcf_local") || die "Unable to open file $exome_vcf: $!\n";
  while (<EXOME>) {  #get edit positions and edit types
    if ($_ =~ /^#CHROM/) {
      @head = split;
    } 
    if (!($_ =~ /^#CHROM/)) {
      my @array = split;
      for (my $i=0;$i<@array;$i++) {
	if ($head[$i] eq $exome_id_local) {
	  if (($array[$i] =~ /0[\/\|]1/)||($array[$i] =~ /1[\/\|]0/)) {
	    if (!($array[4] =~ /,/)) {  #just look at single allele sites (may need to improve this later)
	      my $chr = $array[0]; $chr =~ s/chr//; $chr = "chr".$chr; #Make sure that vcf contains "chr";
	      my $tag = $chr."_".$array[1];
	      my $len1 = length($array[3]); my $len2 = length($array[4]);
	      if (($len1==1)&&($len2==1)) { #just for snps
		$hets{$tag} = $array[3]."_".$array[4]; #store het positions with base change
		push @het_array, $tag;
		$size++;
	      }
	      if (($len1!=1)||($len2!=1)) { #store insertions and deletions
		if ($len1>$len2) {
		  $deletion{$tag} = $array[3]."_".$array[4];
		}
		if ($len2>$len1) {
		  $insertion{$tag} = $array[3]."_".$array[4];
		}
	      }
	    }
	  }
	}
      }
    }
  }
  close (EXOME);
  if ($size==0) {
    die "ERROR: No Heterozygous SNPs (DNA) found for ID selected ($exome_id_local) - Please check VCF and/or input parameters\n"; 
  }
}

sub CollectRNASNPs ($$) {
  my ($rna_id_local,$rna_vcf_local) = @_;
  my @head = (); my $size = 0;
  open (RNA, "$rna_vcf_local") || die "Unable to open file $rna_vcf_local: $!\n";
  while (<RNA>) {  #get edit positions and edit types
    chomp;
    if ($_ =~ /^#CHROM/) {
      @head = split;
    } 
    if (!($_ =~ /^#CHROM/)) {
      my @array = split;
      for (my $i=0;$i<@array;$i++) {
	if ($head[$i] eq $rna_id_local) {
	  if (($array[$i] =~ /0[\/\|]1/)||($array[$i] =~ /1[\/\|]0/)) {
	    if (!($array[4] =~ /,/)) {
	      my $chr = $array[0]; $chr =~ s/chr//; $chr = "chr".$chr; #Make sure that vcf contains "chr";
	      my $tag = $chr."_".$array[1];
	      my $len1 = length($array[3]); my $len2 = length($array[4]);
	      if (($len1==1)&&($len2==1)) { #just for snps
		if (!($hets{$tag})) { #if not already stored in DNA hash
		  $rna_hets{$tag} = $array[3]."_".$array[4]; #store het positions with base change
		  $size++;
		}
	      }
	      if (($len1!=1)||($len2!=1)) { #store insertions and deletions
		if ($len1>$len2) {
		  if (!(exists $deletion{$tag})) {
		    $deletion{$tag} = $array[3]."_".$array[4];
		  }
		}
		if ($len2>$len1) {
		  if (!(exists $insertion{$tag})) {
		    $insertion{$tag} = $array[3]."_".$array[4];
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  close (RNA);
  if ($size==0) {
    print "WARNING: No Additional Heterozygous SNPs detected in 2nd VCF file for ID selected ($rna_id_local): Please check VCF and/or input parameters\n"; 
  }
}

sub InsertSNPs ($) {
  my ($chrs_local) = @_[0];
  my $chromosome = "chr".$chrs_local; $chromosome =~ s/23/X/; $chromosome =~ s/24/Y/;
  my $chromosome_number = $chromosome; $chromosome_number =~ s/chr//;
  
  my $fastq1_chr = $rna_dir."/fastq_sim_".$rna_id."_".$chrs_local."_R1.fastq";
  my $fastq2_chr = $rna_dir."/fastq_sim_".$rna_id."_".$chrs_local."_R2.fastq";
  open (FAST1, ">$fastq1_chr") || die "Unable to open fastq1 file to write to: $!\n";
  open (FAST2, ">$fastq2_chr") || die "Unable to open fastq2 file to write to: $!\n";
  
  my @het_array_chr = ();
  for (my $i=0;$i<@het_array;$i++) { #collect hets just in each chromosome
    if ($het_array[$i] =~ /^$chromosome\_/) {
      push @het_array_chr, $het_array[$i];
    }
  }
 
  my @capture = ();
  my %tags = ();
  for (my $i=0;$i<@het_array_chr;$i++) { #for each het pos - speed up using BED file with all positions?
    my $pos = $het_array_chr[$i];
    my @a = split (/_/,$pos);
    my $tag = $a[0].":".$a[1]."-".$a[1]; #create region identifier
    my @temp_out = `samtools view $rna_file $tag`;
    push @capture, $_ foreach @temp_out;
    @temp_out = sort (map {my @a = split/\s/; $a[0]} @temp_out); #get first element then sort
    my @temp_unique = do { my %seen; grep { !$seen{$_}++ } @temp_out }; #get unique elements
    my $temp_string = join (',',@temp_unique); $tags{$pos} = $temp_string;
  }
  @capture = sort (map {my @a = split/\s/; $a[0]} @capture); #get first element then sort
  my @unique = do { my %seen; grep { !$seen{$_}++ } @capture }; #get unique elements
  my $tempout_tags = $rna_dir."/hetbam_".$rna_id."_".$chromosome_number."_tags.txt"; #file for read identifiers
  open (TRIAL, ">$tempout_tags") || die "Unable\n";
  print TRIAL "$_\n" foreach @unique;
  close (TRIAL);
  
  my @collected_sam = `samtools view $rna_file $chromosome | grep -Fwf $tempout_tags`; #for all read identifiers overlapping hets, get read and read pair
  @collected_sam = sort(@collected_sam);
  my %read_data = ();
  for (my $i=0;$i<@collected_sam;$i=$i+2) {
    my @array = split (/\s/,$collected_sam[$i]);
    my @array_pair = split (/\s/,$collected_sam[$i+1]);
    if ($array[0] eq $array_pair[0]) {
      $read_data{$array[0]} = $collected_sam[$i]."second".$collected_sam[$i+1]; #create tag and store each pair of reads in hash
    }
    if ($array[0] ne $array_pair[0]) {
      warn "Either $array[0] or $array_pair[0] are not paired - this may cause biases in the final output: please ensure reads are paired\n";
      $i--;
    }
  }

  ###########################Now circle through hets on this chr and use master SAM file:
  
  for (my $i=0;$i<@het_array_chr;$i++) { #foreach het, create site specific SAM
    my $master_pos = $het_array_chr[$i];
    my @a = split (/_/,$master_pos); 
    my $covered_tag = $a[1];  #store position for covered check later on (read only included if it actually covers SNP with nucleotide (not N in cigar)
    my $hits = $tags{$master_pos};    
    if ($hits) { #if site has reads covering it
      my @hits = split (/,/,$hits);
      my @reads = ();
      foreach (@hits) { #foreach read
	my @pairs = split (/second/,$read_data{$_}); #get the pairs from read hash and push into site specific array
	push @reads, $_ foreach @pairs;
      }
      @reads = reverse(sort(@reads));
      
      my @new_reads = (); my @snps_contained = (); my @snps_contained_non = (); # store how many alleles covering each position exists in each read
      my @snp_weights = (); #array for number of snps contained in read pair - important to balance coverage across sites
      for (my $j=0;$j<@reads;$j=$j+2) { #for each pair of reads, look for snps and replace with alt hap
	my $ref_reads; my $covered_check=0; #this is to see if SNP position has actually been covered by a base (fine for indels - not counted)
	my $snps_contained_string = ""; #add if SNP position included
	my $snp_weight = 0; #variable to store snp numbers per pair
	my %snp_catch = (); #store snps found in read pairs - but don't count same snp twice
	for (my $k=0;$k<2;$k++) { #for each pair of reads:
	  my $read = $reads[$j+$k];
	  my @array = split(/\t/,$read);
	  my $flag = $array[1];
	  my $binary = sprintf ("%08b",$flag);
	  (($k==0)&&($binary=~/^0/)) || (($k==1)&&($binary=~/^1/)) || die "Error: Either reads are not paired, or they are out of order\n"; #check for correct order and pairing of reads
	  my $cigar =  $array[5];
	  my $sequence = $array[9]; my @seq = split (//,$sequence); #get sequence
	  my $quality = $array[10]; my @qual = split (//,$quality); #get quality
	  my @alt = @seq; #create alternative sequence
	  my @alt_qual = @qual; #create alternative quality sequence
	  my @c = split (//,$cigar);
	  my $nums; my $pos = $array[3]; my $match = 0; my $in_got = 0;
	  for (my $l=0;$l<@c;$l++) { #read through cigar
	    if (!($c[$l] =~ /[MDIHNS]/)) {
	      $nums .= $c[$l];
	    }
	    if ($c[$l] eq "M") {  #for matching sequences check if there is a het SNP present from exome data
	      for (my $m=$pos;$m<($pos+$nums);$m++) {
		my $tag = $chromosome."_".$m; #generate position
		if ($m==$covered_tag) { 
		  $covered_check++;
		}
		if (exists $hets{$tag}) { #if snp is in DNA
		  my @a = split (/_/,$hets{$tag});
		  my $change;
		  foreach my $base (@a) {
		    if (!($seq[$match] =~ /$base/i)) { #for each base in SNP, if the sequence contains one of these, change the alt seq to the other base
		      if ($alt[$match] =~ /[AGTC]/i) { #if base exists (not removed by indel)
			$alt[$match] = $base;
			$snps_contained_string .= ",".$tag; #count alleles covered by reads
			$snp_weight++ if (!(exists $snp_catch{$tag})); #add weight when snp has not yet been counted
			$snp_catch{$tag} = 1; #store that snp has been counted
		      }
		    }
		  }
		} 
		if (exists $rna_hets{$tag}) { #look for snps in RNA and do same as above
		  my @a = split (/_/,$rna_hets{$tag});
		  my $change;
		  foreach my $base (@a) {
		    if (!($seq[$match] =~ /$base/i)) { #for each base in SNP, if the sequence contains one of these, change the alt seq to the other base
		      if ($alt[$match] =~ /[AGTC]/i) { #if base exists (not removed by indel)
			$alt[$match] = $base;
		      }
		    }
		  }
		}
		if (exists $insertion{$tag}) { #look for insertion
		  if ($m<($pos+$nums-1)) { #if not at the end of the matching sequence, add in the indel (if at the end it means that the read has finished or and indel will likely follow
		    my @a = split (/_/,$insertion{$tag});
		    my $base = $a[0];
		    if ($seq[$match] =~ /$base/i) {
		      $in_got = 1;
		      $alt[$match] = $a[1]; #change the sequence for the insertion
		      my $ins_length = length($a[1]);
		      my $aq;
		      for (my $qq=0;$qq<$ins_length;$qq++) {
			$aq .= $qual[$match];
		      }
		      $alt_qual[$match] = $aq;
		    }
		  }
		}
		if (exists $deletion{$tag}) { #look for insertion
		  if ($m<($pos+$nums-1)) { #if not at the end of the matching sequence, add in the indel (if at the end it means that the read has finished or and indel will likely follow
		    my @a = split (/_/,$deletion{$tag});
		    my $change = $a[0]; my $change_length = length($change); #get deletion from vcf
		    my $check_read;
		    for (my $loc=($match);$loc<($match+$change_length);$loc++) { #check read to see if same insertion is present
		      $check_read .= $seq[$loc];
		    }
		    if ($check_read eq $change) { #if deletion matches vcf deletion
		      for (my $loc=($match+1);$loc<($match+$change_length);$loc++) { #remove it from alternative read
			$alt[$loc] = "";
			$alt_qual[$loc] = "";
		      }
		    }
		  }
		}
		$match++; #move along sequence array
	      }
	      $pos+=$nums;
	      $nums="";
	    }
	    if ($c[$l] =~ /[N]/) {  #if skip reference (gap), reposition location to check by increasing position along reference
	      $pos+=$nums;
	      $nums="";
	    }
	    if ($c[$l] =~ /[D]/) {  #if delete, reposition location to check by increasing position along reference
	      my $back_one = $pos-1; #go back on reference by one
	      my $tag = $chromosome."_".$back_one; #generate position
	      if (exists $deletion{$tag}) { #if deletion has been called
		my @a = split (/_/,$deletion{$tag});
		my $base = $a[1]; my $new_match=$match-1;
		if ($seq[$new_match] =~ /$base/i) {
		  $alt[$new_match] = $a[0]; #change the sequence for the insertion
		  my $del_length = length($a[0]);
		  my $aq;
		  for (my $qq=0;$qq<$del_length;$qq++) {
		    $aq .= $qual[$new_match];
		  }
		  $alt_qual[$new_match] = $aq;
		}
	      }
	      $pos+=$nums;
	      $nums="";
	    }
	    if ($c[$l] eq "I") {  #if insertion, reposition location to check by increasing position along read
	      my $trial = $pos-1; my $tag = $chromosome."_".$trial; #generate position before base to look for indel
	      if (exists $insertion{$tag}) { #if indel is found as variation
		my @a = split (/_/,$insertion{$tag});
		my $change = $a[1]; my $change_length = length($change); #get insertion from vcf
		my $check_read;
		for (my $loc=($match-1);$loc<($match-1+$change_length);$loc++) { #check read to see if same insertion is present
		  $check_read .= $seq[$loc];
		}
		if ($check_read eq $change) { #if insertion matches vcf insertion
		  $in_got = 1;
		  for (my $loc=($match);$loc<($match-1+$change_length);$loc++) { #remove it from alternative read
		    $alt[$loc] = "";
		    $alt_qual[$loc] = "";
		  }
		}
	      }
	      $match+=$nums;
	      $nums="";
	    }
	    if ($c[$l] eq "S") {  #if softclipped, move along the sequence array by the number softclipped (start refers to seq after clipping)
	      $match+=$nums;
	      $nums="";
	    }
	  }
	  my $altseq = join('',@alt); #remake alt sequencing read and store both original and alternative reads in string
	  my $altqual = join ('',@alt_qual); #remake alternative quality scores
	  if($flag & hex("0x10")) { #if read is in reverse direction, flip over - required to get properly paired reads
	    $sequence =reverse($sequence);
	    $sequence =~ tr/ACGTRYSWKMBDHV/TGCAYRSWMKVHDB/;
	    $altseq =reverse($altseq);
	    $altseq =~ tr/ACGTRYSWKMBDHV/TGCAYRSWMKVHDB/;
	    $quality = reverse($quality);
	    $altqual = reverse($altqual);
	  }
	  $ref_reads .= "\t".$array[0]."_0\t".$sequence."\t".$quality."\t".$array[0]."_1\t".$altseq."\t".$altqual; ###EXTRA and remove below
	}
	$ref_reads =~ s/\t//;
	if ($covered_check>0) { #check if position has been covered by read or pair  ####EXTRA
	  push @new_reads, $ref_reads; #push original and alternative read pairs into array
	  $snps_contained_string =~ s/,//; #remove spare comma from snp string
	  push @snps_contained, $snps_contained_string; #store the snps contained in read and pair
	  push @snp_weights, $snp_weight; #store number of heterozygotes covered by the read
	  my $snps_contained_string_non;
	  while (my ($ptag,$irrel) = each %snp_catch) { #store weight of non-overlapping SNPs (only count position once if read pairs overlap)
	    $snps_contained_string_non .= ",".$ptag;
	  }
	  $snps_contained_string_non =~ s/,//; #remove spare comma from snp string
	  push @snps_contained_non, $snps_contained_string_non;
	}
      }
      
      if (@snp_weights != @new_reads) {
	print "$master_pos: Weights do not match reads\n"; #Check that a weight has been recorded for each read pair
      }
      my $read_number = @new_reads;
      if ($read_number>0) { #ADDED
	my @read_count = ();
	for (my $n=0;$n<$read_number;$n++) { #build index array to name each selected read pair
	  push @read_count, 0;
	}

	for (my $n=0;$n<$sample_depth;$n++) { #randomly select new reads to obtain sample depth coverage (2X - default 1000) at het position (one ref, one alt)
	  my $r = int(rand($read_number));
	  my $w = $snp_weights[$r]; #check the snp weight of the read (if greater than 1, will be sampled elsewhere too)
	  my $rw = int(rand($w));
	  if ($rw==0) { #sample according to weight: e.g reads/pairs with 2 snps will be sampled again for the other snp, so only add haf of the time here
	    my $read_data = $new_reads[$r];
	    my @a = split(/\t/,$read_data);
	    print FAST1 "\@"."$a[0]\_$read_count[$r]\n$a[1]\n+\n$a[2]\n"; #generate fastq1 file with this data - gives 50:50 split of alt and ref at snp pos
	    print FAST1 "\@"."$a[3]\_$read_count[$r]\n$a[4]\n+\n$a[5]\n";
	    print FAST2 "\@"."$a[6]\_$read_count[$r]\n$a[7]\n+\n$a[8]\n"; #generate fastq2 file with this data - gives 50:50 split of alt and ref at snp pos
	    print FAST2 "\@"."$a[9]\_$read_count[$r]\n$a[10]\n+\n$a[11]\n";
	    $read_count[$r]++;
	    my @snps_collected = split(/,/,$snps_contained[$r]);
	    foreach my $snp_item (@snps_collected) {
	      $snp_allele_coverage{$snp_item}++;
	    }
	    my @snps_collected_non = split(/,/,$snps_contained_non[$r]);
	    foreach my $snp_item_non (@snps_collected_non) {
	      $snp_allele_coverage_non{$snp_item_non}++;
	    }
	  }
	}
      }
    }
  }
  system ("rm $tempout_tags");
  close (FAST1);
  close (FAST2);
}

__END__

=head1 NAME
 
ASE Program
 
=head1 USAGE

perl ase_normalisation_sim.pl [options] --DNAvcf <VCF_file_name> --DNAid <ind_id_from_vcf> --OUTid <Outfolder_name> --RNAbam <Mapped_RNA_BAM>

=head1 OPTIONS

--DNAvcf <VCF_file_name> : Name and path to VCF file containing variation that will be used within the simulation. (Required)

--DNAid <ind_id_from_vcf> : Name of the individual as shown in the VCF file. (Required)

--OUTid <Outfolder_name> : Name of Output folder (also used in output file names).  Directory will be created with this name and simulation files will be produced within this directory. (Required)

--RNAbam <Mapped_RNA_BAM> : Name of the BAM file containing mapped RNA sequencing data for the individual of interest.  This file must contain properly paired, uniquely mapped reads. (Required)

--MaxProcs <number_of_cpus> : The simulation can be run in parallel on multiple cpus. We recommend using no more than 8 cpus on a node containing 12 cpus. [default: 1] (Optional)

--SampleDepth <depth_of_reads> : The depth of reads to be simulated across each site. Increasing this may improve accuracy during normalisation, decreasing will improve speed of simulation and subsequent mapping. [default: 2000] (Optional)

--ADDvcf <additional_vcf> : It is possible to supply an additional VCF file, and additional variation not present in the required VCF will be simulated in surrounding regions.  This variation will not be used for read collection at heterozygous sites, but will be used if it falls within selected reads.  One example would be to supply a VCF from exome sequencing in order to find heterozygous SNVs and simulate reads across them, and then provide a whole-genome VCF to add extra SNPs around these sites for more realistic simulations. (Optional)

--ADDid <id_additional_vcf> : If an additional VCF file is provided, use this option to define the name of the individual as specified in that VCF (can be different from ID in required VCF file). (Optional)

Alan Hodgkinson, E<lt>alan.j.hodgkinson@gmail.comE<gt>.
Sainte-Justine CHU Research Center, Montreal University.

Perl Module Requirements: 

Parallel::Loops
Getopt::Long
Pod::Usage

Program requirements (in path):

Samtools (1.0 or newer)
Pigz

=head1 AUTHOR

Alan Hodgkinson, E<lt>alan.j.hodgkinson@gmail.comE<gt>.
Sainte-Justine UHC Research Center, Montreal University.

=head1 DATE 

25-March-2015

=cut  


