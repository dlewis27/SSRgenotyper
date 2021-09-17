# SSRgenotyper

Many programs can identify simple sequence repeats (SSRs) in genomic data. SSRgenotyper extends SSR identification to genotype calling across multiple individuals in diversity panels and mapping populations. SSRgenotyper will find SSR motifs of lengths 2, 3 and 4 from SAM files and an SSR reference fasta file (SsrReferenceFile.fasta). Mononucleotide SSRs and soft masked sequences (lower case sequences) are excluded during the analysis process. SSRgenotyper has only been tested on diploid and allopolyploid organisms - use with caution on autopolyploids. Several outputs are possible including a simple table with the SSR marker name, position and SSR alleles (defined by the repeat number of the repeated motif in the individual). Specific output files include a Genepop formatted file for genetic diversity analyses and a traditional A, H, B mapping file output - phased to the parents of the population for bi-parental linkage mapping populations. Problems and questions regarding usage can be sent to ssrgenotyperhelp@gmail.com.

SSRgenotyper requires a SSR reference fasta file which lists each targeted SSR with ~100 bp of flanking sequence. The SSR reference fasta file can be easily created using a combination of easy to use bioinformatic tools, specifically [MISA](https://webblast.ipk-gatersleben.de/misa/misa_sourcecode_22092015.zip), [Samtools](https://github.com/samtools/samtools) and [Bedtools](https://bedtools.readthedocs.io/en/latest/). MISA is run against the genome assembly of the species of interest to identify the location of the putative SSR loci. The reference genome, referred to below as my_reference.fasta, can be a gold standard, chromosomal-scale reference genome or a simple draft assembly consisting of several thousand contigs. Bedtools is then used to extract the targeted SSRs and their flanking sequences. Flanking sequences of ~100 bp upstream and downstream are needed for mapping/genotyping purposes. We refer to the MISA/Bedtools output as the SsrReferenceFile.fasta.

SSRgenotyper requires python libraries [regex](https://anaconda.org/conda-forge/regex) and [biopython](https://anaconda.org/anaconda/biopython).  An easy way to install these is through [conda](https://docs.conda.io/en/latest/). SSRgenotyper uses python 3.7. Using earlier versions may result in errors. 

## Make the SSR reference fasta (SsrReferenceFile.fasta) with MISA and Bedtools:

## 1. MISA
`perl misa.pl my_Reference.fasta`

#### MISA requires a misa.ini file in the directory where MISA is being executed that should look like this

definition(unit_size,min_repeats):                   2-6 3-4 4-4\
interruptions(max_difference_between_2_SSRs):        100\
GFF:                                                     true

## 2. Modify the MISA produced gff files as follows
### Remove any compound SSRs and calculate how much flanking sequence is available at each SSR locus
`for i in *.gff; do grep -v "compound" $i | awk '{if ($5-$4 >10 && $5-$4 <50) print $1 "\t" $4-100 "\t" $5+100}' > $i.mod.gff; echo "processing $i"; done`

### Concatenate the modified gff files
`for i in *.mod.gff; do cat $i >> cat.gff; echo "processing $i"; done`

### Remove SSRs that do not have sufficient flanking sequence

`awk '($2 >= 0)' cat.gff > cat_filter1.gff` 

### Use bedtools getfasta to make the SsrReferenceFile.fasta

`bedtools getfasta -fi my_Reference.fasta -bed cat_filter1.gff -fo SsrReferenceFile.fasta`

## 3. Map the Illumina reads to the SsrReferenceFile.fasta

We trim and quality control our reads with [Trimmomatic](https://github.com/timflutre/trimmomatic), which produces paired forward (_1P.fq.gz) and a reverse reads files (_2P.fq.gz) for each sample. In our example code below the trimmed reads are then mapped to the SsrReferenceFile.fasta using [BWA](https://github.com/lh3/bwa), however any short read mapping software should work (minimap2, bowtie2, etc.). After mapping the reads, PCR duplicates are removed using [Samtools](https://github.com/samtools/samtools) -markdup. To genotype multiple individuals in a population, each individual should have its own FASTQ file (i.e., one FASTQ per individual, producing one SAM file per individual). The basic steps are as follows:

### Index the SsrReferenceFile.fasta file

`bwa index SsrReferenceFile.fasta`

### Map the Illumina reads to the SsrReferenceFile.fasta (paired-end reads process shown)

`for forward_file in *_1P.fq.gz; do name=echo $forward_file | sed 's/_1P.fq.gz//\'; bwa mem -M ../reference/SsrReferenceFile.fasta ${name}_1P.fq.gz ${name}_2P.fq.gz -o $name.sam; done`

*after trimming with Trimmomatic each sample has a *_1P.fq.gz and *_2P.fq.gz file

### Remove PCR duplicate reads
Sort SAM files by read name:
`for i in *.sam; do samtools sort -n -o $i.sorted $i; done`

Identify mate coordinates:
`for i in *.sorted; do samtools fixmate -m $i $i.fixmate; done` 

Re-sort SAM files by genomic location:
`for i in *.fixmate; do samtools sort -o $i.position $i; done`

Mark and remove duplicates:
`for i in *.position; do samtools markdup -r $i $i.markdup; done`

*each sample will have a markdup file where PCR duplicates have been removed (the other files can be deleted)

### Quality control SAM files for SSRgenotyper
While the whole SAM file can be passed to SSRgenotyper we encourage users to first filter the markdup file with Samtools to improve performance and to remove errantly mapped reads:

`for i in *.markdup; do samtools view $i -q 45 > $i.Q45.sam; done`

This will remove reads with mapping quality less than 45 as well as unnessary header information. SSRgenotyper provides further filtering (option -Q) that can be used for additional filter stringency if required. 

A file listing all SAM files to be processed is required by SSRgenotyper and can be produced with:

`ls *.Q45.sam > samFiles.txt`

## How SSRgenotyper Works

SSRgenotyper identifies reads mapping to each of the SSRs identified in the SsrReferenceFile.fasta for each of the SAM files and makes a genotypic call based on the number of SSR units. Processing of the SAM files is very fast with minimal memory utilization.

*usage: python3 SSRgenotyper.py <SsrReferenceFile.fasta> <SamFiles.txt> <OutputFileName*>
  
## Options

Positional arguments:

**SsrReferenceFile.fasta** - This is a fasta file consisting of the SSRs and flanking sequence\
**SamFiles.txt** - Text document with the SAM file names separated by newlines\
**OutputFileName** - Output file name (".ssr" extension will be added)

Optional arguments:
  
**-M --MinorAlleleHet** The minimum percentage of the minor allele for a genotype to be considered heterozygous. For diploid organisms, a maximum of two alleles should be present at the SSR locus (a minor and major allele). We use this parameter to set the cutoff for calling a heterozygote. If two alleles are found, but the percentage of minor allele reads is below this cutoff, the genotype will be called as a homozygote for the major allele. If the species or population is expected to be homozygous (i.e., cleistogamous, recombinant inbred lines, etc.), setting this argument to 0.51 will produce only homozygous genotypic calls. Between 0 and 1. [0.20]\
**-S --Support** The minimum number of supporting reads needed for a genotype to be called at an SSR locus. If the total number of reads is less than this threshold, the genotype will not be called and will be reported as "0,-4" (see below). For example, if "-S 6" and there are three reads supporting an allele of five repeat units and two reads supporting an allele with seven repeat units, no alleles will be called as the minimum threshold of 6 total reads was not met. [3]\
**-R --RefUnitsMin** The minimum SSR repeat number used to identify repeats using MISA. Set this to the "min_repeats" used to identify the SSRs in the misa.ini file.  In our example above, the smallest repeat we searched for was a trinucleotide motif repeated four times (3-4) - thus we would set this to 4. [4]\
**-P --PopUnitsMin** The minimum SSR repeat number allowed within the population. It is recommended that this be smaller than the --RefUnits. This parameter allows individuals within the population to have a smaller repeat number than the reference repeat number. [3]\
**-B --BorderSeq** The number of nucleotides on each side of the reference SSR motif that will be evaluated to determine whether a read supports a genotypic call. If there is a concern with non-specific read mapping (allopolyploidy, paralogous/duplicated loci), increasing this parameter can add specificity to the allele calling.  See also –mismatch. [20] \
**-s --spuriousAlleleRemoval** SSR loci are known to produce spurious (stutter) alleles in both PCR and sequencing which can result in ambiguous genotypes (multi-allelic). These spurious alleles are normally minimal in frequency but need to be removed for accurate genotyping. Diploid species should only have a maximum of two alleles possible. SSR loci that have more than the expected two (major/minor) alleles are marked as ambiguous (missing) if the additional allele(s) is supported by more than this threshold percentage of reads. Otherwise the allele(s) is considered to be spurious and excluded. This parameter only affects alleles beyond the major and minor allele. [0.10]\
**-m --mismatch** The number of mismatches (substitutions and indels) allowed in the border sequence (--BorderSeq) region between the reads and the reference. For more distantly related individuals (e.g., difference species or subspecies) the number can be modestly increased to accommodate expected sequence differences. [0]\
**-A --AlignmentShow** Provides a text file of the read alignments for a specific individual (SAM file) and a specific SSR locus. The user specifies a specific SSR locus name and a specific SAM file name separated by "," (e.g. SSR1,samFile.markdup.Q45.sam). This is useful for spot checking alignments and specific genotypic calls.\
**-W --WindowOffset** Window search offset. We do not recommend changing the default offset. [1]

Data filtering parameters\
**-F --FilterDataLoci** The maximum missing data threshold for reporting an SSR locus. For example, if "--FilterDataLoci 0.3" then any SSR locus with more than 30% missing data is removed from the output files. Between 0 and 1. [0.30]\
**-f --filterDataSam** The maximum missing data threshold for reporting an individual (SAM file). For example, if "--filterDataSam 0.3" then any individual with more than 30% missing data is removed from the output files. Between 0 and 1. [0.30]\
**-Q --QualityFilter** Only reads with equal to or greater than the specified mapping quality (MAPQ) threshold are utilized in the genotype calling. MAPQ (SAM file format, column 5) defines the probability that a read is wrongly mapped: MAPQ = -10 * long10(p), where p is the probability that the read is incorrectly mapped. [45]\
**-N --NameSize** Used to shorten output names. SAM files often have excessively long file names as they often include long sequencing library names. This option will limit the number of characters printed in the output table for each SAM file. When using this option make sure that new names will still be unique as SAM files with the same name will result in inaccurate results. [100]

Output file types\
**-G --Genepop**. Ouputs a genepop formated file. [GENEPOP](https://kimura.univ-montp2.fr/~rousset/Genepop.htm) is a population genetics software package originally developed by M. Raymond and F. Rousset at the Laboratiore de Genetique et Environment, Montpellier, France.  See output file descriptions for more information.\
**-L --LinkageMapFile** Outputs an "A,B,H" genotype file for downstream linkage map construction of biparental populations. The first two SAM files in the SamFiles.txt should correspond to two parents of the biparental population and are necessary for correctly phasing the genotypes. If the parental genotypes are monomorphic (same genotype) or if one of the parents is missing a genotypic call, the SSR locus is excluded. Similarly, if one or both of the parents have a heterozygous genotype, phase cannot be determined, and the marker is also be excluded. SSRgenotyper can impute the genotype of a missing parent (the other parent must have a genotype) based on the progeny by providing the -L argument with a with a specific numeric threshold value (e.g., -L 0.30). If a threshold is provided and two alleles (major and minor) are found in the progeny with the minor allele frequency above the provided threshold and where one allele agrees with the genotyped parent, SSRgenotyper will internally assign the second allele as the genotype of the missing parent, however the imputed genotype will still be reported in the .map file as missing - thus allowing the user to identify which markers had an imputed parental genotype. [.30]


## Example
`python3 SSRgenotyper.py SsrReferenceFile.fasta SamFiles.txt OutputFileName -M 0.30 -S 6 -B 3 -F 0.20 -f 0.30 -Q 60 -s 0.10 -G `

## Output files

### The .ssr file 
The basic output file is a tab delimited table with SSR names (and position) in rows and the names of the SAM files (individuals) as column names. The genotypes are called as alleles that reflect the number of repeat units found at each SSR locus.  For example, a "9,9" genotypic call reflects a homozygote where all the reads mapping to the SSR locus had a repeat unit of “9”.  Similarly, a “7,9” genotypic call reflects a heterozygote where reads with both repeat numbers of “7” and “9” were identified in the read mapping. 

The missing data in this file provide addition information as to why the genotype was coded as missing.  Specifically missing data is coded as either "0,-1", "0,-2", "0,-3", "0,-4", where the second number codes for why no genotypes were called:

0,-1: The reference SSR was not evaluated by SSRgenotyper for one of several possible reasons, including the reference sequence was soft masked or the repeat was not a di-, tri- or tetranucleotide repeat , see also --RefUnitsMin. Note that these errors will not show up unless --FilterDataLoci is 0.\
0,-2: No reads mapped to the SSR locus.\
0,-3: More than the expected two alleles were identified in the mapped reads, see also “--spuriousAlleleRemoval”. \
0,-4: Insufficient reads were identified to support calling an allele. The minimum number of supporting reads needed to support a call is determined by “--Support”.

### The .ssrstat file
This file includes basic run and genotyping statistics, including a listing of the positional and optional arguments selected as well as various run statistics. The statistics calculated for each of the various output files types (.ssr, .pop, .map) take into account the missing data and spurious allele filters. The stats reported include total number of SAM files (individual) reported, total number of SSR loci reported, total number of genotypes called (separated as homozygous and heterozygous calls) and the percentage of missing calls.

### The .pop file
SSRgenotyper can output genotypes in a format accepted by [GENEPOP](https://kimura.univ-montp2.fr/~rousset/Genepop.htm). GENEPOP is popular population genetics program that can be used to test for Hardy-Weinberg Exact tests, linkage disequilibrium, population differentiation, Nm estimates as well as file conversion to other population genetics platforms (i.e., FSTAT, BIOSYS, and LINKDOS). GENEPOP files can easily be converted to other file formats used by other population genetics platforms (e.g., [ARLEQUIN](http://cmpg.unibe.ch/software/arlequin35/), [STRUCTURE](https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/html/structure.html), [MSA](http://i122server.vu-wien.ac.at/MSA/info.html/MSA_info.html), etc) using [PGDSpider](http://www.cmpg.unibe.ch/software/PGDSpider/). All missing data is coded as "000000". Specific missing codes are provided in the .ssr file. Monomorphic loci are excluded in the output.

### The .map file
SSRgenotyper can output an "A,B,H" genotype file to facilitate downstream linkage map construction of biparental populations using linkage map construction software (e.g., [JoinMap](https://www.kyazma.nl/index.php/JoinMap/)). The .map file is a tab delimited table with SSR names in rows and the names of the SAM files (individuals) as column names. The first two SAM files in the SAMFiles.txt should correspond to the two parents of the biparental population and are necessary for correctly phasing the genotypes. The .map file is solicited using the optional argument "--LinkageMapFile". By default if the parental genotypes are monomorphic or if one of the parents is missing a genotypic call the SSR locus is excluded from the output. SSRgenotyper can impute a missing genotype for a parental line based on the segregation patterns observed in the progeny (see --LinkageMapFile documentation above). If one or both of the parents have a heterozygous genotype, phase cannot be determined, and the marker is excluded. Similarly, any progeny genotypic calls that do not agree with the alleles identified in the parents are deemed problematic and are assigned a missing value. All missing data is coded as "-". Specific missing codes are provided in the .ssr file. 

## High Performance Computing (HPC) scripts

For users wishing to prepare their SAM files using an HPC cluster we provide an example job scripts for SLURM based systems that can be easily manipulated for other HPC scheduling systems. The following script will perform bwa mapping, samtools sort and markdup as well as remove read with poor mapping quality (q<45):

#!/bin/bash

#SBATCH --time=8:00:00   # walltime\
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)\
#SBATCH --nodes=1   # number of nodes\
#SBATCH --mem-per-cpu=4G   # memory per CPU core\
#SBATCH -J "bwa_mem"   # job name

module load bwa_0.7.17\
module load samtools/1.6

bwa mem -M -t $SLURM_NPROCS $BWA_INDEX ${file}_1P.fq.gz ${file}_2P.fq.gz -o ${file}.sam && samtools sort -n --threads $SLURM_NPROCS -o ${file}.sorted.sam ${file}.sam && samtools fixmate -m ${file}.sorted.sam ${file}.sorted.fixmate.sam && samtools sort --threads $SLURM_NPROCS -o ${file}.sorted.fixmate.position.sam ${file}.sorted.fixmate.sam && samtools markdup -r ${file}.sorted.fixmate.position.sam ${file}.sorted.fixmate.position.markdup.sam && samtools view ${file}.sorted.fixmate.position.markdup.sam -q 45 --threads $SLURM_NPROCS > ${file}.sorted.fixmate.position.markdup.Q45.sam

