# SSRgenotyper

SSRgenotyper will find simple sequence repeats (SSRs) of lengths 2, 3 and 4 from SAM files and a modified reference fasta. Mono-nucleotide SSRs are excluded. Soft masked sequences are excluded during the analysis process. SSRgenotyper has only been tested on diploid organisms - use with caution on polyploids. Several output are possible including a simple table with the SSR marker name, position and SSR alleles (defined by the repeat number of the repeat motif). Specific output files for genetic diversity analysis include a Genepop formated file and a traditional A, H, B mapping files output can be selected - phased to the parents of the population for bi-parental linkage mapping populations. Questions regarding useage can be sent to ssrgenotyperhelp@gmail.com.

SSRgenotyper requires a modified reference which lists each targeted SSR with ~100 bp of flanking sequence. The modified reference can be easily created using a combination of easy to use bioinformatic tools, specifically [MISA](https://webblast.ipk-gatersleben.de/misa/misa_sourcecode_22092015.zip) and Bedtools. MISA is run against the reference genome of the species of interest to identify the location of the targeted SSRs. The reference genome can be a gold standard reference genome or a simple draft reference. Bedtools is then used to extract the targeted SSRs and their flanking sequences. Flanking sequence of ~100 bp upstream and downstream are needed for mapping/genotyping purposes. We refer to the MISA/Bedtools output as the modified reference (below it is referred to as "my_modified_Reference.fasta").

## Make the modified reference with MISA and Bedtools:

### MISA:
perl misa.pl my_Reference.fasta

### MISA requires a misa.ini file that should looks like this:

definition(unit_size,min_repeats):                   2-6 3-4\
interruptions(max_difference_between_2_SSRs):        100\
GFF:                                                     true

### Modify the gff files:
for i in \*.gff; do grep -v "compound" $i | awk '{if ($5-$4 >10 && $5-$4 <50) print $1 "\t" $4-100 "\t" $5+100}' > $i.mod.gff; echo "processing $i"; done

*this process removes any compound SSRs and calculates how much flanking sequence is available.

### Concatenated the modified gff files:
for i in \*.mod.gff; do cat $i >> cat.gff; echo "processing $i"; done

### Remove SSRs that do not have sufficient flanking seqeunce:

awk '($2 >= 0)' cat.gff > cat_filter1.gff 

### Use bedtools getfasta to make the modified reference:

bedtools getfasta -fi my_Reference.fasta -bed cat_filter1.gff -fo my_modified_Reference.fasta

## Map the Illumina reads to the modified reference

Map trimmed and quality controlled Illumina reads (FASTQ) to the modified reference. We provide an example using BWA mem below, however any mapping software should work (minimap2, bowtie2, etc.). The illumina reads should be quality controlled (e.g., Trimmomatic) with PCR duplicates marked (e.g., Samtools -markdups). To map Illumina reads with BWA:

### Index the modified reference file:

bwa index my_modified_Reference.fasta my_modified_Reference.fasta

### Map the Illumina reads to the modified reference.fasta (single end reads process shown - can be done for paired-end reads):

for i in \*.fq; do bwa mem myReferenceForSSRgenotyper.fasta $i > $i.sam; done 

### Prepare the SAM files for SSRGenotyper
while the whole SAM file can be passed to SSRGenotyper, we encourage users to first filter the sam file with samtools to improve performance:

for i in *.sam; do samtools view $i -q 45 > $i.Q45; done

This will remove reads with mapping quality less than 45. SSRgenotyper provides further filtering arguments (-Q) that can be used for additional filter stringency. The SAM files do not need to be sorted or indexed. A file listing all SAM file to be processed is required by SSRGenotyper and can be produced with:

ls *.Q45 > samFiles.txt

## How it Works

SSRgenotyper will go through the reference file, finding the SSRs and the specified number of flanking nucleotides on both sides of the SSR. It then goes through each SAM file, finding reads that mapped to the reference sequence, and looks the SSR pattern and matching flanking nucleotides. A call is then made based on the number of SSR units.

## Output

A file with a tab separated table that includes the name of the reference sequence as the names of the rows, and the names of the SAM files for the column names. The elements of the table show the number of SSR units found. For example, "9,9" means that the SSR unit was found 9 times. This is likely a homogenous allele. If this were to show "9,8" then reads supporting an allele with 9 SSR units were found as well as reads supporting 8 SSR units. This is a heterozygous allele. If the first number is "0" followed by a negative number, then no alleles were called. The negative number is the code for why no alleles were called. For example "0,-2" shows no alleles were called because no reads in the SAM file mapped to this marker. The following codes are:

-1: No SSR was found in the reference sequence.\
-2: No reads from the accession were mapped to this marker.\
-3: More than 2 alleles were found so the marker in this SAM file. It is considered is ambiguous\
-4: Not enough reads with the target SSR were found to support calling an allele. The minimum number of supporting reads is determined by option -S

Three files are generated. The .ssr file has the tab-delimited table and the .ssrstat shows the genotyping statistics. The genepop file ends in .pop

## Options

usage: python3 SSRgenotyper.py <ReferenceFile> <SamFiles> <OutputFile>

positional arguments:

**ReferenceFile** - The reference file (FASTA)\
**SamFiles** - Text document with the SAM file names separated by newlines\                       
**OutputFile** - Output file name ( ".ssr" will be added onto it)

optional arguments:
  
**-A --AlleleRatio** The minimum ratio of reads for 2 alleles to be considered heterozygous. If 2 alleles are found but the ratio of reads for each does not meet this threshold, it will be reported as a homozygous allele. (default = .2)\
**-N --NameSize** The number of characters, starting from the first character, from the name of the SAM file to be listed in the output table. Make sure that this is still unique. If this causes SAM files to have have the same name, the results will be inaccurate. (default = 100)\
**-R --RefUnits** The minimum number of SSR units in a reference SSR. For example, if the parameter for this is 4, "ATGATGATG" will not be found but "ATGATGATGATG" will.(default = 4)\
**-P --PopUnits** The minimum number of SSR units in an accession SSR. This is the same as the REFUNITS parameter, but for the SAM files (default = 3)\
**-F --FlankSize** The number of flanking bases on each side of the SSR that must match the reference. If there is a high amount of ambiguous non-calls, then increasing this may help (default = 20).\
**-S --Support** Then minimum number of total supporting reads for alleles to be called. For example, if this parameter is 6 and there are 3 reads that show an allele with 5 SSR units and 2 reads that show an allele with 7 SSR units, no allele will not be called and then reported as "0,-4". (default = 3)\
**-W --WindowOffset** Offset on each side of the reference sequence, making a window for searching for the SSR. It is recommended that this not be changed. (default = 1)\
**-r --refFilter** If the porportion of the population that had no call for a marker meets this threshold, then the marker will not be reported (i.e. if this is .8 and 90% of the population had no call at this marker, then the marker will be omitted from the output table.). In other words, how much of data will you allow to be missing. The smaller this number, the more data in the rows. Should be between 0 and 1 (default = 1)\
**-Q --QualityFilter** filters the reads from the SAM file. Only reads above this threshold will be considered in SSRgenotyper. If there are a lot of ambiguous non-calls, increasing this can help. (default = 45)\
**-X --Xdebug** Provide marker name and SAM file name separated by ",". This will output the reads from the SAM file that mapped to the marker. If this option is not "" then the main program will not run. The output will be in debug.txt (default = "")\
**-M --Map** Output a map for a biparental population. The first 2 SAM files in the SAM file list should be the two parents. Non-informative (parent 1 and parent 2 have the same allele) markers and markers where there is no call for one or both the parents will be excluded. If one or both of the parents are heterozygous, the marker will also be excluded.\
**-a ambiguoussalvage** If the reads supporting the 3rd most supported allele divided by the total reads supporting the first 2 alleles is equal to or greater than this, the call will be ambiguous (default = .1).\
**-m mismatch** The number of mismatch allowance for each flanking region. Insertions, deletions, and substitutions considered (default = 0)
## Example
python3 SSRgenotyperV2.py myReferenceForSSRgenotyper.fasta samFiles.txt myOutput -F 20 -S 1
