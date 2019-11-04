# SSRgenotyper
find SSRs in a population

SSRmulti will find simple sequence repeats of length 2,3, and 4 from given sam files and a reference. SSR with the same letters (i.e. GGG) or are lowercase will be excluded. The output is a table that with the marker names and SSR alleles. 

## Making the Reference

The reference is expected to be a .fasta file with a target SSR on each sequence. This can be created using MISA on a reference and then using Bedtools to extend the sequence on both sides for mapping purposes. Extending it by 75 bp upstream and downstream seems to work well. SSRmulti will find which SSR is on each sequence so this does not need to be provided.


## Prepping the SAM Files
The sam files should be filtered for quality, the headers removed and only the first 10 row should be kept. This can be done with samtools and cut:
  samtools view -q 45 <samFile> | cut -f1-10 > <samFile>.f1
The sam files should be provided in as a list in a text document with each file name on a new line. This can be done with:
  ls *.f1 > samFiles.txt

## How it Works

The script go through the reference file, finding the SSRs and the specified number of flanking nucleotides on both sides of the SSR. It then goes through each sam file, looks at the reads that mapped to the reference sequence, and looks for reads with the SSR pattern and matching flanking nucleotides.

## Output

A file with a tab separated table that includes the name of the reference sequence for the row names, and the first 7 characters in the names of the sam files for the column names. The elements of the table show the number of SSR units found. For example, "7-7" means that the SSR unit was found 7 times. This is likely a homogenous allele. If this were to show "7-6" then reads supporting an allele with 7 SSR units were found as well as reads supporting 6 SSR units. This is likely a hetrozygous allele. If "0-0" shows, then there was no call.

## Options

usage: python3 SSRFinder.py ReferenceFile SamFiles OutputFile

positional arguments:
  ReferenceFile         The refrence file (FASTA)
  SamFiles              Text document with the SAM file names seperated by
                        newline
  OutputFile            Output file name ( ".ssr" will be added onto it)

optional arguments:
  -A ALLELERATIO, --AlleleRatio ALLELERATIO
                        The minmum ration of major to minor alleles 
                        (default = .2)
  -R REFUNITS, --RefUnits REFUNITS
                        The minimum number of SSR units in a reference SSR
                        (default = 4)
  -P POPUNITS, --PopUnits POPUNITS
                        The minimum number of SSR units in an accession SSR
                        (default = 3)
  -F FLANKSIZE, --FlankSize FLANKSIZE
                        The number of flanking bases on each side of the SSR
                        that must match the reference (default = 10)
  -S SUPPORT, --Support SUPPORT
                        Then minimum number of supporting reads for an allele
                        to be called (default = 2)
  -W WINDOWOFFSET, --WindowOffset WINDOWOFFSET
                        Offset on each side of the reference sequence, making
                        a window for searching for the SSR (default = 1)

## Example
python SSRFinder <referenceFile.fa> <mySamFiles.txt> <myOutput>
