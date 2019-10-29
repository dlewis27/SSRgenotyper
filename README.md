# SSRFinder
find SSRs in a population

SSRmulti will find simple sequence repeats of length 2,3, and 4 from given sam files and a reference. SSR with the same letters (i.e. GGG) or are lowercase will be excluded.

MAKING THE REFERENCE:

The reference is expected to be a .fasta file with a target SSR on each sequence. This can be created using MISA on a reference and then using Bedtools to extend the sequence on both sides for mapping purposes. Extending it by 75 bp upstream and downstream seems to work well. SSRmulti will find which SSR is on each sequence so this does not need to be provided. 

PREPPING THE SAM FILES:
The sam files should be filtered for quality, the headers removed and only the first 10 row should be kept. This can be done with samtools and cut:
  samtools view -q 45 <samFile> | cut -f1-10 > <samFile>.f1
The sam files should be provided in as a list in a text document with each file name on a new line. This can be done with:
  ls *.f1 > samFiles.txt

HOW IT WORKS:

The script go through the reference file, finding the SSRs and the 2 flanking nucleotides on both sides of the SSR. It then goes through each sam file, looks at the reads that mapped to the refeerence sequence, and looks for reads with the SSR pattern and a matching 2 flanking nucleotides.

OUTPUT:

A file "SSRoutput.ssr" with a tab separated table that includes the name of the reference sequence for the row names, and the first 7 characters in the names of the sam files for the column names. The elements of the table show the number of SSR units found. For example, "7-7" means that the SSR unit was found 7 times. This is likely a homogenous allele. If this were to show "7-6" then reads supporting an allele with 7 SSR units were found as well as reads supporting 6 SSR units. This is likely a hetrozygous allele. If "0-0" shows, then there was no call.

EXAMPLE:
python SSRmulti <referenceFile.fa> <mySamFiles.txt>
