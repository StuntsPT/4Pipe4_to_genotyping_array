#4Pipe4_to_sequenom

###This program will transform a sorted and indexed bam file, plus a FASTA file generated by 4Pipe4 into SNPs and their flanking regions ready for a sequenom Mass Array (c).

###Usage:

    python2 4Pipe4_to_sequenom.py infile.bam infile.fasta variation_treshold(float: 0-1) max_variations(int)

The "infile.bam.bai" has to be in the same directory as the "infile.bam".

###Requirements:

* Python 2.x (tested with 2.7)
* pysam (code.google.com/p/pysam/‎)

###License:
GPLv3.