# SSRgenotyper

Many programs can identify SSR motifs in genomic data. SSRgenotyper extends SSR identification to genotype calling across multiple individuals in diversity panels and mapping populations. SSRgenotyper will find simple sequence repeats (SSRs) of lengths 2, 3 and 4 from SAM files and a modified reference fasta. Mononucleotide SSRs are excluded. Soft masked sequences are excluded during the analysis process. SSRgenotyper has only been tested on diploid organisms - use with caution on polyploids. Several outputs are possible including a simple table with the SSR marker name, position and SSR alleles (defined by the repeat number of the repeat motif). Specific output files for genetic diversity analysis include a Genepop formated file and a traditional A, H, B mapping files output can be selected - phased to the parents of the population for bi-parental linkage mapping populations. Questions regarding useage can be sent to ssrgenotyperhelp@gmail.com.

SSRgenotyper requires a modified reference which lists each targeted SSR with ~100 bp of flanking sequence. The modified reference can be easily created using a combination of easy to use bioinformatic tools, specifically [MISA](https://webblast.ipk-gatersleben.de/misa/misa_sourcecode_22092015.zip) and Bedtools. MISA is run against the reference genome of the species of interest to identify the location of the targeted SSRs. The reference genome can be a gold standard reference genome or a simple draft reference. Bedtools is then used to extract the targeted SSRs and their flanking sequences. Flanking sequence of ~100 bp upstream and downstream are needed for mapping/genotyping purposes. We refer to the MISA/Bedtools output as the modified reference (below it is referred to as "my_modified_Reference.fasta").

## Make the modified reference with MISA and Bedtools:

### MISA:
perl misa.pl my_Reference.fasta

### MISA requires a misa.ini file that should looks like this:

definition(unit_size,min_repeats):                   2-6 3-4\
interruptions(max_difference_between_2_SSRs):        100\
GFF:                                                     true

### Modify the gff files:
for i in *.gff; do grep -v "compound" $i | awk '{if ($5-$4 >10 && $5-$4 <50) print $1 "\t" $4-100 "\t" $5+100}' > $i.mod.gff; echo "processing $i"; done

*this process removes any compound SSRs and calculates how much flanking sequence is available.

### Concatenated the modified gff files:
for i in *.mod.gff; do cat $i >> cat.gff; echo "processing $i"; done

### Remove SSRs that do not have sufficient flanking seqeunce:

awk '($2 >= 0)' cat.gff > cat_filter1.gff 

### Use bedtools getfasta to make the modified reference:

bedtools getfasta -fi my_Reference.fasta -bed cat_filter1.gff -fo my_modified_Reference.fasta

## Map the Illumina reads to the modified reference

Map trimmed and quality-controlled Illumina reads (FASTQ) to the modified reference. We provide an example using BWA mem below, however any mapping software should work (minimap2, bowtie2, etc.). The illumina reads should be quality controlled (e.g., Trimmomatic) with PCR duplicates marked (e.g., Samtools -markdups). To genotype multiple individuals in a population, each individual should have its own FASTQ file (i.e., one FASTQ per individual, producing one SAM file per individual). To map Illumina reads with BWA:

### Index the modified reference file:

bwa index my_modified_Reference.fasta my_modified_Reference.fasta

### Map the Illumina reads to the modified reference.fasta (single end reads process shown - can be done for paired-end reads):

for i in \*.fq; do bwa mem myReferenceForSSRgenotyper.fasta $i > $i.sam; done 

*each individual is represented by a different .fq file

### Quality control SAM files for SSRgenotyper
While the whole SAM file can be passed to SSRgenotyper we encourage users to first filter the sam file with samtools to improve performance:

for i in *.sam; do samtools view $i -q 45 > $i.Q45; done

This will remove reads with mapping quality less than 45 as well as unnessary header information. SSRgenotyper provides further filtering (option -Q) that can be used for additional filter stringency. The SAM files do not need to be sorted or indexed. Lastly a file listing all SAM files to be processed is required by SSRgenotyper and can be produced with:

ls *.Q45 > samFiles.txt

## How SSRgenotyper Works

SSRgenotyper identifies reads mapping to each of the SSR in the modified reference file for each of the SAM files and makes a genotypic call based on the number of SSR units.  Processing of the SAM files is very fast with minimal memory utilization.

usage: python3 SSRgenotyper.py <my_modified_reference.fasta> <SamFiles.txt> <OutputFileName>
  
## Options

positional arguments:

**SsrReferenceFile** - This is the my_modified_reference.fasta referred to previously (FASTA)\
**SamFiles.txt** - Text document with the SAM file names separated by newlines\
**OutputFileName** - Output file name (".ssr" extension will be added)

Optional arguments:
  
**-M --MinorAlleleHet** The minimum percentage of the minor allele for a genotype to be considered heterozygous. For diploid organisms, a maximum of two alleles should be present at the SSR locus (a minor and major allele). We use this parameter to set the cut off for calling a heterozygote. If two alleles are found but the percentage of minor allele reads is below this cutoff the genotype will be called as a homozygote for the major allele. Between 0 and 1. [0.20]\
**-S --Support** The minimum number supporting reads needed for a genotype to be called. If the total number of reads is less than this threshold, the genotype will not be called and will be reported as "0,-4" (see below). For example, if "-S 6" and there are three reads supporting an allele of five repeat units and two reads supporting an allele with seven repeat units, no alleles will be called. [3]\
**-R --RefUnitsMin** The minimum SSR repeat number used to identify repeats using MISA. Set this to the "min_repeats" used to identify the SSRs in the misa.ini file.  In our example above, the smallest repeat we searched for was a trinucleotide motif repeated four times (3-4) - thus we would set this to 4. [4]\
**-P --PopUnitsMin** The minimum SSR repeat number allowed within the population. It is recommended that this be smaller than the --RefUnits. This parameter allows individuals within the population to have a smaller repeat number than the reference repeat number. [3]\
**-F --FlankSize** The number of flanking nucleotides on each side of the SSR that must match the reference for a read to be used to support an allelic call. If there is a concern with non-specific read mapping (allopolyploidy, paralogous/duplicated loci), increasing this parameter can add specficity. [20] \
**-W --WindowOffset** Window search offset. We do not recommend changing the default offset. [1]\
**-r --MissingDataFilter** The maximum missing data threshold for reporting an SSR locus. For example, if "-r 0.3" then any SSR loci with greater than 30% missing data will not be reported. Between 0 and 1. [0.30]\
**-Q --QualityFilter** Only Reads with equal to or greater than the specified mapping quality are used to support genotype calling. [45]\
**-A --AlignmentShow** Provides a text display of the read alignments for a specific SAM file and a specific SSR locus. The user specifies a marker name and a specific SAM file name separated by "," (e.g. SSR1,File1.sam). This is useful for spot checking alignments and specific genotypic calls.
**-L --LinkageMapFile** Outputs an "A,B,H" genotype file for downstream linkage map construction of biparental populations. The first two SAM files in the SAMFiles.txt should correspond to two parents of the biparental population and are necessary for correctly phasing the genotypes. If the parental genotypes are monomorphic (same genotype) or if one of the parents is missing a genotypic call the SSR locus is excluded. Similarly, if one or both of the parents have a heterozygous genotype, phase cannot be determined, and the marker is also be excluded.\
**-s --spuriousAlleleRemoval** SSR loci are known to produce spurious (stutter) alleles in both PCR and sequencing which can result in ambigous genotypes (multi-allelic). These spurious alleles are normally minimal in frequency but need to be removed for accurate genotyping. SSR loci that have additional alelles, beyond the expected two (major/minor) alleles, are marked as ambigous if the additional allele is supported by more than this threshold percentage of reads, otherwise the allele is considered to be a spurious and excluded.  This parameter only affects alleles beyond the major and minor allele. [0.10]\
**-m mismatch** The number of mismatches (substitutions and indels) allowed in the flanking region between the reads and the reference. For more distantly related samples the number can be modestly increased. [0]\
**-N --NameSize** Used to shorten output names. Sam files often have excessively long file names as they often include long sequencing library names. This option will limit the number of characters printed in the output table for each sam file. When using this option make sure that new names will still be unique as sam files with the same name will result in inaccurate results. [100]

## Example
python3 SSRgenotyper.py my_modified_reference.fasta SamFiles.txt OutputFileName -H 30 -S 3

## Output files

### The .ssr file 
The basic output file is a tab delimited table with SSR names (and position) in rows and the names of the sam files (individuals) as column names. The genotypes are called as alleles that reflect the number of repeat units found at each SSR locus.  For example, a "9,9" genotypic call reflects a homozygote where all the reads mapping to the SSR locus had a repeat unit of “9”.  Similarly, a “9,7” genotypic call reflects a heterozygote where reads with both repeat numbers of “9” and “7” were identified in the read mapping.  

If the first number is "0" followed by a negative number, then no alleles were called. The second number codes for why no alleles were called. For example "0,-2" indicates that no alleles were called because no reads mapped to this marker. All no call codes are defined here:

-1: No SSR motif was identified in the modified reference sequence.\
-2: No reads mapped to the SSR locus.\
-3: More than expected two alleles were found.  See “—SpuriousAlleleRemoval”. \
-4: Insufficient reads were identified to support calling an allele. The minimum number of supporting reads need to support a call is determined by “-S”.

### The .ssrstat file
This file includes basic genotyping statistics.
 
### The .pop file
This is a genepop file.  Need to add more info here!

