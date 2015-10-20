# ASE

Programs: ase_normalisation_sim.pl and ase_normalisation_compare.pl

These two programs can be used to normalise mapped Paired End RNA sequencing files (bam) and produce an output file giving a P-value for the liklihood that the site is deviating from an allele ratio of 50:50.

Details for running the programs can be obtained through the --help options.

## Dependencies

Perl Module Requirements: 

Parallel::Loops;
Getopt::Long;
Pod::Usage;

Software Requirements:

[Samtools (1.0 or newer)] (http://www.htslib.org);
[Pigz] (http://zlib.net/pigz/);

## Usage Pipeline

Step 1: 
-------

Map RNA Sequencing Data with software of choice.

### Examples:
	STAR --genomeDir STARgenome2pass --runThreadN threads --alignEndsType EndToEnd --readFilesIn fastq_1.fq.gz
	fastq_2.fq.gz --readFilesCommand zcat --outFileNamePrefix Outname --outFilterMismatchNmax 10
	
	or
	
	tophat2 -o Outname -N 5 --read-edit-dist 5 -p threads --transcriptome-index=index bowtie2_index
	fastq_1.fq.gz fastq_2.fq.gz

This can then be followed by desired filtering strategies such as PCR duplicate removal and/or keeping properly paired and uniquely mapped reads. Note that in order to perform the next step, only paired reads should be kept (non-paired reads will be discarded and may thus bias the simulation).


Step 2: 
-------

Use the mapped data as the basis to create the simulated dataset:

	perl ase_normalisation_sim.pl [options] --DNAvcf <VCF_file_name> --DNAid <ind_id_from_vcf> 
	--OUTid <Outfolder_name> --RNAbam <Mapped_RNA_BAM>

--DNAvcf <VCF_file_name> : Name and path to VCF file containing variation that will be used within the simulation. (Required)

--DNAid <ind_id_from_vcf> : Name of the individual as shown in the VCF file. (Required)

--OUTid <Outfolder_name> : Name of Output folder (also used in output file names).  Directory will be created with this name and simulation files will be produced within this directory. (Required)

--RNAbam <Mapped_RNA_BAM> : Name of the BAM file containing mapped RNA sequencing data for the individual of interest.  This file must contain properly paired, uniquely mapped reads. (Required)

--MaxProcs <number_of_cpus> : The simulation can be run in parallel on multiple cpus. We recommend using no more than 8 cpus on a node containing 12 cpus. \[default: 1\] (Optional)

--SampleDepth <depth_of_reads> : The depth of reads to be simulated across each site. Increasing this may improve accuracy during normalisation, decreasing will improve speed of simulation and subsequent mapping. \[default: 2000\] (Optional)

--ADDvcf <additional_vcf> : It is possible to supply an additional VCF file, and additional variation not present in the required VCF will be simulated in surrounding regions.  This variation will not be used for read collection at heterozygous sites, but will be used if it falls within selected reads.  One example would be to supply a VCF from exome sequencing in order to find heterozygous SNVs and simulate reads across them, and then provide a whole-genome VCF to add extra SNPs around these sites for more realistic simulations. (Optional)

--ADDid <id_additional_vcf> : If an additional VCF file is provided, use this option to define the name of the individual as specified in that VCF (can be different from ID in required VCF file). (Optional)

Step 3: 
-------

The simulation program will generate two fastq.gz files with paired reads.  Map these files using the exact same approach that you used in step 1. 

The simulation program will also produce two count files showing the number of reference (and thus alternative) reads covering each site for both trimmed (ignoring read pair overlaps) and non-trimmed data.  These are for use in the ase_normalisation_compare.pl program to generate results.

Step 4: 
-------

Use the original mapped data, the simulated mapped data and the count files (Null_sim_file below) to generate P-values for each heterozygous site:


	perl ase_normalisation_compare.pl [options] --DNAvcf <VCF_file_name> --DNAid <ind_id_from_vcf> --OUTid <Outfolder_name> 
	--RNAbam <Mapped_RNA_BAM> --RNASIMbam <Mapped_RNASIM_BAM> --Ref <Reference_fasta> --MapNull <Null_sim_file>
	
--DNAvcf <VCF_file_name> : Name and path to VCF file containing variation that will be used within the simulation. (Required)

--DNAid <ind_id_from_vcf> : Name of the individual as shown in the VCF file. (Required)

--OUTid <Outfolder_name> : Name of Output folder (also used in output file names).  Directory will be created with this name and simulation files will be produced within this directory. (Required)

--RNAbam <Mapped_RNA_BAM> : Name of the BAM file containing mapped RNA sequencing data for the individual of interest. (Required)

--RNASIMbam <Mapped_RNA_BAM> : Name of the BAM file containing mapped RNA sequencing data for the individual of interest, generated from the simulation program. (Required)

--Ref <Reference_fasta> : Full path to Reference fasta file for use in Samtools Mpileup.  This should be the same reference as used to call SNVs present in the supplied VCF file. (Required)

--MapNull <Null_sim_file> : The path to null allele file generated by the simulation program. (Required)

--MaxProcs <number_of_cpus> : The simulation can be run in parallel on multiple cpus. We recommend using no more than 8 cpus on a node containing 12 cpus. \[default: 1\] (Optional)

--Baq : Use to turn on probabilistic realignment for the computation of base alignment quality in Samtools Mpileup (BAQ). \[default: off\] (Optional)

## Results File

The results file will detail allele counts before and after normalisation, as well as a P-value from a binomial test:

	#Chr_Pos	  MapRef  MapCov	Null	SimRef	SimCov	SimAdjRef	SimAdjCov	MapAseP	SimAdjP
	chr1_881627	0.483695652173913	184	2000	0.5	      4000	0.483695652173913	184	0.7125308   0.7125308
	chr1_900505	0.772727272727273	22	2000	0.525129  3283	0.703703703703704	27	0.01690054  0.05223899
	

Column1 (Chr_Pos): Position of heterozygous site

Column2 (MapRef): Proportion of reference alleles in original mapped data

Column3 (MapCov): Coverage in original mapped data

Column4 (Null): Number of reference alleles in null simulated dataset

Column5 (SimRef): Proportion of reference alleles in mapped simulated data

Column6 (SimCov): Coverage in mapped simulated data

Column7 (SimAdjRef): Proportion of reference alleles heterozygous site after normalisation

Column8 (SimAdjCov): Coverage at heterozygous site after normalisation

Column9 (MapAseP): P value for allelic bias in original mapped data

Column10 (SimAdjP): P value for allelic bias after normalisation





	
