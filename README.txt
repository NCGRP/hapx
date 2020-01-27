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
Usage: hapx -r ref -b bam -a alnr [-f inc] [-F exc] [-q qual] -s sites
where,
ref = path to reference genome sequence in multi-fasta format [required]
bam = path to bam file of reads aligned to ref [required]
sites = genomic positions to use [required]
     Provide a comma-delimited list of the form:
         jcf7180008454378:303-304,jcf7180008531951:103-495
     which specifies bp 303-304 of the contig named "jcf7180008454378" and bps 103-495 of contig "jcf7180008531951".
     Samtools view will not allow you to specify a single nucleotide position, so just use a 2bp range.
alnr = aligner used to create bam file, options: gem bwamem
inc = integer flag value for Samtools view -f option (properties of reads to include), see https://broadinstitute.github.io/picard/explain-flags.html [default=3, include reads that are paired and are mapped in a proper pair]
exc = integer flag value for Samtools view -F option (properties of reads to exclude) [default=3852, exclude unmapped reads && reads whose mate or pair is unmapped && not primary alignment && read fails platform/vendor quality checks && read is PCR or optical duplicate && supplementary alignment]
qual = Samtools view -q option (minimum mapping quality of included reads) [default=60]

Examples: ./hapx.sh /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta /share/space/reevesp/patellifolia/xtr/AllP.merged.bam bwamem jcf7180008454378:303-304,jcf7180008531951:103-495,jcf7180008395354:294-295,jcf7180008827236:1031-1032,jcf7180008378511:277-278,jcf7180008637475:7-8,jcf7180008587925:106-107,jcf7180008527965:177-178,jcf7180008578969:84-85,jcf7180008484650:470-471,jcf7180008856767:98710-98886;

./hapx.sh /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam gem jcf7180008454378:303-304,jcf7180008531951:103-495,jcf7180008395354:294-295,jcf7180008827236:1031-1032,jcf7180008378511:277-278,jcf7180008637475:7-8,jcf7180008587925:106-107,jcf7180008527965:177-178,jcf7180008578969:84-85,jcf7180008484650:470-471,jcf7180008856767:98710-98886;



https://broadinstitute.github.io/picard/explain-flags.html