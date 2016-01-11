use strict;
use Parallel::Loops; 
use Statistics::R;
use Getopt::Long;
use Pod::Usage;

#name: ase_normalisation_compare.pl
#description: A program to obtain reads and pairs for heterozygote positions in RNA data and simulate both haplotypes
#date: 22nd April 2014
#Written By: Alan Hodgkinson

###8th October: adjusted so alt and ref alleles can be flipped (if vcf has been called to a different reference)
###9th October: chrX added

###Required:
my $exome_id = undef; #exome id
my $rna_id = undef; #rna id
my $exome_vcf = undef; #exome vcf file
my $rna_file = undef; #original rna bam file
my $rna_file_sim = undef; #simulated rna bam file
my $reference = undef; #path to reference sequence for use by samtools
my $map_null = undef; #path to null mapped allele file from simulations
my $maxProcs = 1; #set max processes to 1 by default
my $baq = 0; #Include BAQ realignment in samtools mpileup

my $cmd_help = 0;

GetOptions (
	    'DNAid=s' => \$exome_id,
	    'OUTid=s' => \$rna_id,
	    'DNAvcf=s' => \$exome_vcf,
	    'RNAbam=s' => \$rna_file,
	    'RNASIMbam=s' => \$rna_file_sim,
	    'Ref=s' => \$reference,
	    'MapNull=s' => \$map_null,
	    'MaxProcs=i' => \$maxProcs,
	    'Baq' => \$baq,
	    'help' => \$cmd_help
) or exit(1);

pod2usage(1) if $cmd_help;

print "Option --DNAids required: Please supply name of Individual as defined in VCF file\n" if not defined ($exome_id);
print "Option --OUTid required: Please supply name of Output Directory to write to\n" if not defined ($rna_id);
print "Option --DNAvcf required: Please supply name of VCF file\n" if not defined ($exome_vcf);
print "Option --RNAbam required: Please supply name of BAM file containing original mapped RNA reads\n" if not defined ($rna_file);
print "Option --RNASIMbam required: Please supply name of BAM file containing the simulated mapped RNA reads\n" if not defined ($rna_file_sim);
print "Option --Ref required: Please supply name and path of reference fasta files to be used in mpileup\n" if not defined ($reference);
print "Option --MapNull required: Please supply name of Map Null file produced by simulation program\n" if not defined ($map_null);

if ($maxProcs>8) {
  warn "Warning: Using more than 8 cpus on the same node will probably be very slow as it is likely that system calls to samtools within the simulation will overrun the node.  During testing we found that speed increased up to 8 cpus, but not beyond\n";
}

my $rna_res = $rna_id."_results";
if (-d "$rna_res") {
  die "Folder $rna_id already exists - please rename the --OUTid variable, or delete/rename the existing folder\n";
}
if (!(-d "$rna_res")) {
  mkdir $rna_res, 0755;
}

####MAIN Program
$SIG{INT} = sub { die "Caught a signal interuption $!" };

my %hets = (); 
my %null = ();
my @het_array = (); 

my $R = Statistics::R->new();

&CollectDNASNPs($exome_id,$exome_vcf); #Collect DNA Heterozygote Positions
&CollectNull($map_null);

my $pl = Parallel::Loops->new($maxProcs);
$pl->share(\@het_array);
$pl->share(\%hets);
$pl->share(\%null);
my @chrs = (1..23);
$pl->foreach( \@chrs, sub {
		my $chrs = $_; 
		&PileupCompare($chrs,$rna_file,$rna_file_sim,$reference,$rna_id);
	      });

################Collect Heterozygote Positions:

sub CollectDNASNPs ($$) {
  my ($exome_id_local,$exome_vcf_local) = @_;
  my @head = ();
  open (EXOME, "$exome_vcf_local") || die "Unable to open snp file to read to: $!\n";
  while (<EXOME>) {  #get edit positions and edit types
    if ($_ =~ /^#CHROM/) {
      @head = split;
    } 
    if (!($_ =~ /^#CHROM/)) {
      my @array = split;
      for (my $i=0;$i<@array;$i++) {
	if ($head[$i] eq $exome_id) {
	  if (($array[$i] =~ /0[\/\|]1/)||($array[$i] =~ /1[\/\|]0/)) {
	    my $chr = $array[0]; $chr =~ s/chr//; $chr = "chr".$chr; #Make sure that vcf contains "chr";
	    my $tag = $chr."_".$array[1];
	    if (($array[3] =~ /^[AGTC]$/i)&&($array[4] =~ /^[AGTC]$/i)) { #Single bp SNPs only
	      if (!($array[4] =~ /,/)) {  #just look at single allele sites (may need to improve this later)
		$hets{$tag} = $array[3]."_".$array[4]; #store het positions with base change
		push @het_array, $tag;
	      }
	    }
	  }
	}
      }
    }
  }
  close (EXOME);
}

sub CollectNull ($$) {
  my ($mapped_null_local) = @_;
  open (ALS, "$mapped_null_local") || die "Unable to open Null Allele file to read: $!\n";
  while (<ALS>) {
    if (!($_ =~ /^#/)) {
      my @array = split;
      my $tag = $array[0];
      $tag =~ s/chr//; $tag = "chr".$tag;
      $null{$tag} = $array[1];
    }
  }
  close (ALS);
}

sub PileupCompare ($$$$$) {
  my ($chromosome_number,$rna_file,$rna_file_sim,$reference,$rna_id_local) = @_;
  $chromosome_number =~ s/23/X/;
  my $chromosome = "chr".$chromosome_number;
  my $bed_file = $rna_res."/snps_".$rna_id_local."_".$chromosome.".bed";
  open (BED, ">$bed_file") || die "Unable to open BED file to write to: $!\n";
  my @chr_het_array = ();
  foreach (@het_array) {
    if ($_ =~ /^$chromosome\_(\d{1,10})/) {
      my $site = $1;
      my $tag = $chromosome."_".$site;
      push @chr_het_array, $tag;
      my $b=$site-1; my $a=$site+1; #create SNP regions for BED file (mpileup)
      print BED "$chromosome\t$b\t$a\n";
    }
  }
  close (BED);

  my %real_ref = (); 
  my %real_cov = (); 
  my %sim_ref = (); 
  my %sim_cov = ();

  my $pileup1 = $rna_res."/pileup_real_".$rna_id."_temp_".$chromosome.".pu";
  if ($baq) {
    system ("samtools mpileup -x -d 10000000 -f $reference -r $chromosome -l $bed_file $rna_file > $pileup1");
  }
  if (!($baq)) {
    system ("samtools mpileup -x -B -d 10000000 -f $reference -r $chromosome -l $bed_file $rna_file > $pileup1");
  }
  open (PILEUP1, "$pileup1") || die "Unable to open real pileup file to read: $!\n";
  while (<PILEUP1>) {
    my @array = split;
    my $tag = $array[0]."_".$array[1];
    $tag =~ s/chr//; $tag = "chr".$tag;
    if (exists $hets{$tag}) {
      my @als = split (/_/,$hets{$tag});
      my $alt = $als[1];
      if ($alt =~ /$array[2]/i) {
	$alt=$als[0];
      }
      my ($ref,$prop,$cov) = &countpos($array[4],$array[5],$alt);
      if ($cov>=1) {
	$real_ref{$tag} = $ref;
	$real_cov{$tag} = $cov;
      }
    }
  }
  close (PILEUP1);


  my $pileup2 = $rna_res."/pileup_sim_".$rna_id."_temp_".$chromosome.".pu";
  if ($baq) {
    system ("samtools mpileup -x -d 10000000 -f $reference -r $chromosome -l $bed_file $rna_file_sim > $pileup2");
  }
  if (!($baq)) {
    system ("samtools mpileup -x -B -d 10000000 -f $reference -r $chromosome -l $bed_file $rna_file_sim > $pileup2");
  }
  open (PILEUP2, "$pileup2") || die "Unable to open sim pileup file to read: $!\n";
  while (<PILEUP2>) {
    my @array = split;
    my $tag = $array[0]."_".$array[1];
    $tag =~ s/chr//; $tag = "chr".$tag;
    if (exists $hets{$tag}) {
      my @als = split (/_/,$hets{$tag});
      my $alt = $als[1];
      if ($alt =~ /$array[2]/i) {
	$alt=$als[0];
      }
      my ($ref,$prop,$cov) = &countpos($array[4],$array[5],$alt);
      if ($cov>=1) {
	$sim_ref{$tag} = $ref;
	$sim_cov{$tag} = $cov;
      }
    }
  }
  close (PILEUP2);

  system ("rm $pileup1");
  system ("rm $pileup2");
  system ("rm $bed_file");

  my $outfile = $rna_res."/results_".$rna_id."_".$chromosome.".txt";
  open (RES, ">$outfile") || die "Unable to open results file to write to: $!\n";
  print RES "#Chr_Pos\tMapRef\tMapCov\tNull\tSimRef\tSimCov\tSimAdjRef\tSimAdjCov\tMapAseP\tSimAdjP\n";

  foreach my $tag (@chr_het_array) {
    if ($real_cov{$tag}>0) { #&&(exists $null{$chromosome."_".$tag})) {
      my $cov = $real_cov{$tag}; #get cov
      my $ref = $real_ref{$tag}; #get mapped reference
      my $alt = $cov-$ref;
      my $ref_prop = $ref/$cov;
      my $adj_ref = 0; my $adj_alt = 0; my $adj_cov = 0; my $null = "NA"; my $ref_sim = "NA"; my $cov_sim = "NA"; my $sim_ref_prop = "NA"; my $adj_ref_prop = "NA";
      if (exists $sim_cov{$tag}) {
	$null = $null{$tag} if (exists $null{$tag});
	$cov_sim = $sim_cov{$tag};
        $ref_sim = $sim_ref{$tag};
	my $alt_sim = $cov_sim-$ref_sim;
	$adj_ref = int((((($null-$ref_sim)/$cov_sim)*$cov)+$ref)+0.5);
	$adj_alt = int((((($null-$alt_sim)/$cov_sim)*$cov)+$alt)+0.5);
	$adj_ref = 0 if $adj_ref<0;
	$adj_alt = 0 if $adj_alt<0;
	$adj_cov = $adj_ref+$adj_alt;
	$sim_ref_prop = $ref_sim/$cov_sim;
      }

      my $adj_test = "NA";
      if ($adj_cov>0) {

	$adj_ref_prop = $adj_ref/$adj_cov;
	my $make ="binom.test($adj_ref,$adj_cov,0.5,alternative=\"two.sided\")\$p.value";
	$R->send($make);
	$adj_test = $R->read;
	$adj_test =~ s/\[1\]//;
	$adj_test =~ s/\s//g;
      }

      my $map_test = "NA";
      if ($cov>0) {
	my $make ="binom.test($ref,$cov,0.5,alternative=\"two.sided\")\$p.value";
	$R->send($make);
	$map_test = $R->read;
	$map_test =~ s/\[1\]//;
	$map_test =~ s/\s//g;
      }
 
      #print "N:$null\tC:$cov_sim\tR:$ref_sim\tAR:$adj_ref\tAA:$adj_alt\tAC:$adj_cov\n";
      #print "$chromosome\t$tag\t$ref\t$cov\t$null\t$ref_sim\t$cov_sim\t$adj_ref\t$adj_cov\n";
    
      print RES "$tag\t$ref_prop\t$cov\t$null\t$sim_ref_prop\t$cov_sim\t$adj_ref_prop\t$adj_cov\t$map_test\t$adj_test\n";
    }
  }
  close (RES);
}

sub countpos ($$$) {
  my ($seq,$score,$alt) = @_;
  $seq =~ s/\$//g; #end of read tag - not associated with any quality score
  $seq =~ s/\^.//g; #start of read tag - always followed by a mapping quality for the read, this is different from the base quality in the final column
  my @s = split (//,$seq); my $r = 0; my $c = 0;
  for (my $i=0;$i<@s;$i++) {
    if ($s[$i] =~ /[\.\,$alt]/i) {
      $c++;
      if ($s[$i] =~ /[\.\,]/) {
	$r++;
      }
    }
    if ($s[$i] =~ /[\d]/) {
      $i+=$s[$i];
    }
  }
  my $p = 0;
  if ($c>0) {
    $p = $r/$c;
  }
  return($r,$p,$c);
}


__END__

=head1 NAME
 
ASE Calc Program
 
=head1 USAGE

perl ase_normalisation_compare.pl [options] --DNAvcf <VCF_file_name> --DNAid <ind_id_from_vcf> --OUTid <Outfolder_name> --RNAbam <Mapped_RNA_BAM> --RNASIMbam <Mapped_RNASIM_BAM> --Ref <Reference_fasta> --MapNull <Null_sim_file>

=head1 OPTIONS

--DNAvcf <VCF_file_name> : Name and path to VCF file containing variation that will be used within the simulation. (Required)

--DNAid <ind_id_from_vcf> : Name of the individual as shown in the VCF file. (Required)

--OUTid <Outfolder_name> : Name of Output folder (also used in output file names).  Directory will be created with this name and simulation files will be produced within this directory. (Required)

--RNAbam <Mapped_RNA_BAM> : Name of the BAM file containing mapped RNA sequencing data for the individual of interest. (Required)

--RNASIMbam <Mapped_RNA_BAM> : Name of the BAM file containing mapped RNA sequencing data for the individual of interest, generated from the simulation program. (Required)

--Ref <Reference_fasta> : Full path to Reference fasta file for use in Samtools Mpileup.  This should be the same reference as used to call SNVs present in the supplied VCF file. (Required)

--MapNull <Null_sim_file> : The path to null allele file generated by the simulation program. (Required)

--MaxProcs <number_of_cpus> : The simulation can be run in parallel on multiple cpus. We recommend using no more than 8 cpus on a node containing 12 cpus. [default: 1] (Optional)

--Baq : Use to turn on probabilistic realignment for the computation of base alignment quality in Samtools Mpileup (BAQ). [default: off] (Optional)

Alan Hodgkinson, E<lt>alan.j.hodgkinson@gmail.comE<gt>.
Sainte-Justine CHU Research Center, Montreal University.

Perl Module Requirements: 

Parallel::Loops
Statistics::R
Getopt::Long
Pod::Usage

Program requirements (in path):

Samtools (1.0 or newer)

=head1 AUTHOR

Alan Hodgkinson, E<lt>alan.j.hodgkinson@gmail.comE<gt>.
Sainte-Justine UHC Research Center, Montreal University.

=head1 DATE 

25-March-2015

=cut  

