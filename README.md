# SSRgenotyper
find SSRs in a population

SSRgenotyper will find simple sequence repeats of length 2 and 3 from given sam files and a reference. SSR with the same letters (i.e. GGG) or are lowercase will be excluded. Currently it only works for diploid organisms. The output is a table that with the marker names and SSR alleles. 

## Making the Reference

The reference is expected to be a .fasta file with a target SSR on each sequence surrounded by flanking nucleotides. This can be created using MISA on a reference and then using Bedtools to extend the sequence on both sides for mapping purposes. Extending it by 75 bp upstream and downstream seems to work well. SSRgenotyper will find which SSR is on each sequence so this does not need to be provided.

For example:

#### run misa
perl misa.pl myReference.fasta

#### and misa.ini looks like this:

definition(unit_size,min_repeats):                   2-6 3-4

interruptions(max_difference_between_2_SSRs):        100

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

for i in \*.fq; do BWA mem myReferenceForSSRgenotyper.fasta $i > $i.sam; done 

## Prepping the SAM Files
while the whole SAM file can be passed in, this is strongly discouraged as it will slow down the run time. Filtering the file with samtools will speed it up. For example, run:

samtools view -q 45 <samFile> > samFile.f1

There is no need to convert the file to BAM format or to sort or index the file. The sam files should be provided in as a list in a text document with each file name on a new line. This can be done with:
  ls *.f1 > samFiles.txt

## How it Works

The script go through the reference file, finding the SSRs and the specified number of flanking nucleotides on both sides of the SSR. It then goes through each sam file, looks at the reads that mapped to the reference sequence, and looks for reads with the SSR pattern and matching flanking nucleotides.

## Output

A file with a tab separated table that includes the name of the reference sequence for the row names, and the first 7 characters (unless changed by option -N) in the names of the sam files for the column names. The elements of the table show the number of SSR units found. For example, "9,9" means that the SSR unit was found 9 times. This is likely a homogenous allele. If this were to show "9,8" then reads supporting an allele with 9 SSR units were found as well as reads supporting 8 SSR units. This is a hetrozygous allele. If the first number is "0" followed by a negative number, then no alleles were called and the negative number is the code for whay no alleles were called. For example "0,-2" shows no alleles were called becauses no reads in the SAM file mapped to this marker. The codes corrispond to the following reasons:

-1: No SSR was found in the reference marker. 

-2: No reads from the accession were mapped to this marker.

-3: More than 2 alleles were found so the marker for this asseccion is ambiguous

-4: Not enough reads were found to support calling an allele. The minimum number of supporting reads is determined by option -S.

Two files are generated. The .ssr file has the tab-delemited table and the .ssrstat has a log file that shows the genotyping statistics

## Options

usage: python3 SSRgenotyperV2.py <ReferenceFile> <SamFiles> <OutputFile>

positional arguments:

**ReferenceFile** - The reference file (FASTA)

**SamFiles** - Text document with the SAM file names seperated by
                        newline
                        
**OutputFile** - Output file name ( ".ssr" will be added onto it)

optional arguments:
  
  -A --AlleleRatio
                        The minimum ratio of reads for 2 alleles to be considered heterozygous. If 2 alleles are found but the ratio of reads for each does not meet this threshold, it will be reported as a homozygous allele. (default = .2)

-N --NameSize 
  The number of characters, starting from the first character, from the name of the SAM file to be listed in the output table
  
  -R REFUNITS, --RefUnits REFUNITS
                        The minimum number of SSR units in a reference SSR. For example, if the parameter for this is 4, "ATGATGATG" will not be found but "ATGATGATGATG" will.
                        (default = 4)
  
  -P POPUNITS, --PopUnits POPUNITS
                        The minimum number of SSR units in an accession SSR. This is the same as the REFUNITS parameter, but for the population
                        (default = 3)
  
  -F FLANKSIZE, --FlankSize FLANKSIZE
                        The number of flanking bases on each side of the SSR
                        that must match the reference (default = 10)
  
  -S SUPPORT, --Support SUPPORT
                        Then minimum number of total supporting reads for an allele
                        to be called. For example, if this parameter is 6 and there are 3 reads that show an allele with 5 SSR units and 2 reads that show an allele with 7 SSR units, then this allele will not be called and then reported as "0,-4" (default = 2)
  
  -W WINDOWOFFSET, --WindowOffset WINDOWOFFSET
                        Offset on each side of the reference sequence, making
                        a window for searching for the SSR. It is recomended that this not be changed. (default = 1)
                        
-r --refFilter If the porportion of the population that had no call for a marker meets this threshold, then the marker will not be reported (i.e. if this is .3 and 40% of the population had no call at this marker, then the marker will be ommitted from the output table.) Should be between 0 and 1 (default = 0)

-Q --QualityFilter, filters the reads from the SAM file. Only reads above this threshold will be considered in SSRgenotyper (default = 45)

-X, --Xdebug Provide marker name and SAM file name seperated by ?%?. This will output the reads from the SAM file that mapped to the marker. If this option is not "" then the main program will not run. The output will be in debug.txt (default = "").

## Example
python3 SSRgenotyperV2.py myReferenceForSSRgenotyper.fasta samFiles.txt myOutput -F 20 -S 1
