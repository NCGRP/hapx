hapx extracts molecular haplotypes from bam files at user-specified sites

In working folder:
1) A fasta-formatted reference sequence
2) A sorted, deduped bam file containing reads aligned to the reference in (1)
3) A line delimited list of target sites

Requirements (in path):
1) samtools (hapx calls samtools sort/view/mpileup)
2) bwa (mem/index)
3) sambamba (index)
4) GNU parallel
5) muscle

Usage: hapx ref bam sites
Usage: hapx -r ref -b bam -a alnr [-f inc] [-F exc] [-q qual] [-d -m] -s sites
where,
ref = path to reference genome sequence in multi-fasta format [required]
bam = path to bam file of reads aligned to ref [required]
sites = path to file containing genomic positions to use [required]
     Provide a line delimited list of the form:
         jcf7180008454378:303-304
         jcf7180008531951:103-495
     which specifies bp 303-304 of the contig named "jcf7180008454378" and bps 103-495 of contig "jcf7180008531951".
     Samtools view will not allow you to specify a single nucleotide position, so just use a 2bp range in that case.
alnr = aligner used to create bam file, options: gem bwamem
inc = integer flag value for Samtools view -f option (properties of reads to include), see https://broadinstitute.github.io/picard/explain-flags.html [default=1, include paired reads. Avoid including proper pairs, bwa mem requires many reads to calculate a distribution from which "proper pairs" are determined. In hapx, individual read pairs are aligned to their contig, thus no such distribution can be calculated and bwa mem will not return any proper pairs]
exc = integer flag value for Samtools view -F option (properties of reads to exclude) [default=3852, exclude unmapped reads && reads whose mate or pair is unmapped && not primary alignment && read fails platform/vendor quality checks && read is PCR or optical duplicate && supplementary alignment]
qual = Samtools view -q option (minimum mapping quality of included reads) [default=60]

-d = delete duplicate sequences and subsequences
-m = map (bwa mem) and align (muscle) extracted haploblocks

Examples: ./hapx.sh /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta /share/space/reevesp/patellifolia/xtr/AllP.merged.bam bwamem jcf7180008454378:303-304,jcf7180008531951:103-495,jcf7180008395354:294-295,jcf7180008827236:1031-1032,jcf7180008378511:277-278,jcf7180008637475:7-8,jcf7180008587925:106-107,jcf7180008527965:177-178,jcf7180008578969:84-85,jcf7180008484650:470-471,jcf7180008856767:98710-98886;

./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -a gem -s jcf7180008454378:303-304,jcf7180008531951:103-495,jcf7180008395354:294-295,jcf7180008827236:1031-1032,jcf7180008378511:277-278,jcf7180008637475:7-8,jcf7180008587925:106-107,jcf7180008527965:177-178,jcf7180008578969:84-85,jcf7180008484650:470-471,jcf7180008856767:98710-98886;
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -a gem -d -s jcf7180008454378:303-304,jcf7180008531951:103-495,jcf7180008856767:98710-98886;

./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -a gem -d -s jcf7180008454378:204-304,jcf7180008531951:395-495,jcf7180008395354:195-295,jcf7180008827236:932-1032,jcf7180008378511:178-278,jcf7180008637475:7-107,jcf7180008587925:7-107,jcf7180008527965:78-178,jcf7180008578969:84-184,jcf7180008484650:371-471,jcf7180008856767:98710-98810;
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -a -m gem --s sites.txt;

#use a function to generate target sites (effectively a sliding window)
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -a gem -m -s <(for i in $(seq 1 1 305); do echo jcf7180008587925:"$i"-$(( $i + 1 )); done;);
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -a gem -m -s <(for i in $(seq 1 25 544); do echo jcf7180008531951:"$i"-$(( $i + 25 )); done;)
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -a gem -m -s <(for i in $(seq 1 100 157500); do echo jcf7180008856767:"$i"-$(( $i + 100 )); done;)

<(for i in $(seq 1 1 305); do echo jcf7180008587925:"$i"-$(( $i + 1 )); done;)
<(for i in $(seq 1 25 544); do echo jcf7180008531951:"$i"-$(( $i + 25 )); done;)
<(for i in $(seq 1 100 157500); do echo jcf7180008856767:"$i"-$(( $i + 100 )); done;)




#some postprocessing
#postprocess haploblock counts in log.txt into something plottable
myq() {
      i=$1; #a position on the reference is incoming
      for j in $b;
      do a=$(grep _"$i"\\."$j" numblocks.txt | cut -d: -f4); #number of haploblock read pairs mapped to position i
        if [[ "$a" == "" ]]; then a=0; fi; #if no data at position i set number of haploblock read pairs to 0
        echo "$i $a" >> "$j".txt; #echo to output file specific for the readgroup
      done;
      
}
export -f myq;

grep ^# log.txt | tr '-' '_' > numblocks.txt;
a=$(cut -d_ -f3 numblocks.txt | cut -d. -f1 | sort -n -u); #determine set of possible sites
b=$(cut -d. -f2 numblocks.txt | cut -d$'\t' -f1 | sort -u); export b; #determine set of possible read groups

for i in $b; do echo "position $i"hblox > $i.txt; done; #initialize output file with a header
echo "$a" | parallel --bar --keep-order --env myq --env b myq;

#sort and tab delimit output
for i in $b;
  do sort -t' ' -k1,1n "$i".txt | tr ' ' '\t' > "$i".2.txt;
  done;
#swap back to original filename
for i in $b; 
  do mv "$i".2.txt "$i".txt;
  done;
  
  


