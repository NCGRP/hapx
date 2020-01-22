hapx extracts molecular haplotypes from bam files at user-specified sites

In working folder:
1) A fasta-formatted reference sequence
2) A sorted, deduped bam file containing reads aligned to the reference in (1)

Requirements (in path):
1) samtools (hapx calls samtools sort/view/mpileup)
2) bwa (mem/index)
3) sambamba (index)
4) GNU parallel

Usage: hapx ref bam sites
where,
ref = path to reference genome sequence in multi-fasta format
bam = path to bam file of reads aligned to ref
sites = genomic positions to use. Provide a comma-delimited list of the form:
     jcf7180008454378:303-304,jcf7180008531951:103-495
     which specifies bp 303-304 of the contig named "jcf7180008454378" and bps 103-495 of contig "jcf7180008531951".
     Samtools view will not allow you to specify a single nucleotide position, so just use a 2bp range.

Examples: ./hapx.sh /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta /share/space/reevesp/patellifolia/xtr/AllP.merged.bam jcf7180008454378:303-304,jcf7180008531951:103-495;
