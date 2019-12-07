# SSRgenotyper

SSRgenotyper will find simple sequence repeats (SSRs) of length 2 and 3 from given SAM files and a reference. Homogenous SSRs with only one letter (i.e. GGGGGGG) will be excluded. Lowercase characters are not recognized and will be skipped. Currently it only works for diploid organisms. The output is a table that with the marker names and SSR alleles. An optional output is available that shows alleles based on two parents.

## Making the Reference

The reference is expected to be a FASTA file with a target SSR on each sequence surrounded by flanking nucleotides. This can be created using [MISA](https://webblast.ipk-gatersleben.de/misa/misa_sourcecode_22092015.zip) on a reference genome to find the location of the SSRs, and then using Bedtools to extend the sequence on both sides for mapping purposes. Extending it by 75 bp upstream and downstream seems to work well. SSRgenotyper will find which SSR is on each sequence so this does not need to be provided.

To make the reference with MISA and Bedtools:

#### running misa:
perl misa.pl myReference.fasta

#### and misa.ini looks like this:

definition(unit_size,min_repeats):                   2-6 3-4\
interruptions(max_difference_between_2_SSRs):        100\
GFF:                                                     true

#### modify the resulting gff files and concatenated them
for i in \*.gff; do grep -v "compound" $i | awk '{if ($5-$4 >10 && $5-$4 <50) print $1 "\t" $4-100 "\t" $5+100}' > $i.mod.gff; echo "processing $i"; done

for i in \*.mod.gff; do cat $i >> cat.gff; echo "processing $i"; done

#### remove lines where the range is less than 0

awk '($2 >= 0)' cat.gff > cat_f1.gff 

#### use bedtools to make the reference sequence

bedtools getfasta -fi myReference.fasta -bed cat_f1.gff -fo myReferenceForSSRgenotyper.fasta

#### Mapping

We used BWA mem to map FASTQ files to the newly made reference though any mapping software can be used. To map FASTQ files with BWA mem, first index the reference file with:

bwa index myReferenceForSSRgenotyper.fasta

then map the files with:

for i in \*.fq; do bwa mem myReferenceForSSRgenotyper.fasta $i > $i.sam; done 

## Prepping the SAM Files
while the whole SAM file can be passed in, this is strongly discouraged as it will slow down the run time. Filtering the file with samtools will speed it up. For example, run:

samtools view -q 45 <samFile> > samFile.f1

This will remove reads with mapping quality less than 45. SSRgenotyper has an option, -Q, to further filter the reads. There is no need to convert the SAM files to BAM format or to sort or index. The sam files should be provided as a list in a text document with each file name on a new line. This can be done with:

ls *.f1 > samFiles.txt

## How it Works

SSRgenotyper will go through the reference file, finding the SSRs and the specified number of flanking nucleotides on both sides of the SSR. It then goes through each SAM file, finding reads that mapped to the reference sequence, and looks the SSR pattern and matching flanking nucleotides. A call is then made based on the number of SSR units.

## Output

A file with a tab separated table that includes the name of the reference sequence as the names of the rows, and the names of the SAM files for the column names. The elements of the table show the number of SSR units found. For example, "9,9" means that the SSR unit was found 9 times. This is likely a homogenous allele. If this were to show "9,8" then reads supporting an allele with 9 SSR units were found as well as reads supporting 8 SSR units. This is a heterozygous allele. If the first number is "0" followed by a negative number, then no alleles were called. The negative number is the code for why no alleles were called. For example "0,-2" shows no alleles were called because no reads in the SAM file mapped to this marker. The following codes are:

-1: No SSR was found in the reference sequence.\
-2: No reads from the accession were mapped to this marker.\
-3: More than 2 alleles were found so the marker in this SAM file. It is considered is ambiguous\
-4: Not enough reads were found to support calling an allele. The minimum number of supporting reads is determined by option -S

Two files are generated. The .ssr file has the tab-delimited table and the .ssrstat shows the genotyping statistics

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
**-r --refFilter** If the porportion of the population that had no call for a marker meets this threshold, then the marker will not be reported (i.e. if this is .8 and 90% of the population had no call at this marker, then the marker will be omitted from the output table.) Should be between 0 and 1 (default = 1)\
**-Q --QualityFilter** filters the reads from the SAM file. Only reads above this threshold will be considered in SSRgenotyper. If there are a lot of ambiguous non-calls, increasing this can help. (default = 45)\
**-X --Xdebug** Provide marker name and SAM file name separated by ",". This will output the reads from the SAM file that mapped to the marker. If this option is not "" then the main program will not run. The output will be in debug.txt (default = "")\
**-M --Map** Output a map for a biparental population. The first 2 SAM files in the SAM file list should be the two parents. Non-informative (parent 1 and parent 2 have the same allele) markers and markers where there is no call for one or both the parents will be excluded. If one or both of the parents are heterozygous, the marker will also be excluded.

## Example
python3 SSRgenotyper.py myReferenceForSSRgenotyper.fasta samFiles.txt myOutput -F 20 -S 1
