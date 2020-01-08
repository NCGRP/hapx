hapx extracts molecular haplotypes from bam files at user-specified sites

In working folder:
1) A fasta-formatted reference sequence
2) A bam file containing reads aligned to the reference in (1)

Requirements (in path):
1) samtools (hapx calls samtools sort/view/mpileup)
2) bwa (mem/index)
3) sambamba (index)
