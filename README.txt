hapx extracts molecular haplotypes from bam files at user-specified sites

In working folder:
1) A fasta-formatted reference sequence
2) A sorted, deduped bam file containing reads aligned to the reference in (1)
3) A line delimited list of target sites

Requirements (in path):
1) samtools (hapx calls samtools sort/view/mpileup)
2) bwa (mem/index)
4) GNU parallel
5) muscle

Usage: hapx ref bam sites
Usage: hapx -r ref -b bam -o out -a alnr [-f inc] [-F exc] [-q qual] [-p maxp] [-d -m -x] -s sites
where,
ref = path to reference genome sequence in multi-fasta format [required]
bam = path to bam file of reads aligned to ref [required]
out = name of directory for output files
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
maxp = allow no more than maxp Ns between read pairs. This prevents read pairs from being assembled into an NNN-padded haploblock when they are too far apart.  You may want to set -p maxp according to average library insert size [default 1000];

-d = delete duplicate sequences and subsequences from output
-m = map (bwa mem) and align (muscle) extracted haploblocks
-x = do not write any output files except log (makes -d and -m irrelevant)

Examples: ./hapx.sh /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta /share/space/reevesp/patellifolia/xtr/AllP.merged.bam bwamem jcf7180008454378:303-304,jcf7180008531951:103-495,jcf7180008395354:294-295,jcf7180008827236:1031-1032,jcf7180008378511:277-278,jcf7180008637475:7-8,jcf7180008587925:106-107,jcf7180008527965:177-178,jcf7180008578969:84-85,jcf7180008484650:470-471,jcf7180008856767:98710-98886;

./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -a gem -s jcf7180008454378:303-304,jcf7180008531951:103-495,jcf7180008395354:294-295,jcf7180008827236:1031-1032,jcf7180008378511:277-278,jcf7180008637475:7-8,jcf7180008587925:106-107,jcf7180008527965:177-178,jcf7180008578969:84-85,jcf7180008484650:470-471,jcf7180008856767:98710-98886;
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -a gem -d -s jcf7180008454378:303-304,jcf7180008531951:103-495,jcf7180008856767:98710-98886;

./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -a gem -d -s jcf7180008454378:204-304,jcf7180008531951:395-495,jcf7180008395354:195-295,jcf7180008827236:932-1032,jcf7180008378511:178-278,jcf7180008637475:7-107,jcf7180008587925:7-107,jcf7180008527965:78-178,jcf7180008578969:84-184,jcf7180008484650:371-471,jcf7180008856767:98710-98810;
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o ss -a gem -d -m -s sites.txt;

#use a function to generate target sites (effectively a sliding window)
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o 1305 -a gem -m -s <(for i in $(seq 1 1 305); do echo jcf7180008587925:"$i"-$(( $i + 1 )); done;);
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o 25544 -a gem -m -s <(for i in $(seq 1 25 544); do echo jcf7180008531951:"$i"-$(( $i + 25 )); done;)
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o 157500 -a gem -m -s <(for i in $(seq 1 100 157500); do echo jcf7180008856767:"$i"-$(( $i + 100 )); done;)
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o 157500d -a gem -d -s <(for i in $(seq 1 100 157500); do echo jcf7180008856767:"$i"-$(( $i + 100 )); done;)
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o 157500dm -a gem -d -m -s <(for i in $(seq 1 100 157500); do echo jcf7180008856767:"$i"-$(( $i + 100 )); done;)
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o 157500step1 -a gem -d -s <(for i in $(seq 10000 1 15750); do echo jcf7180008856767:"$i"-$(( $i + 1 )); done;)

#takes about 4 hrs for 50kb (for 1GB then, 9132 years on 232 cores)
#now takes about 1.5 hrs for 50kb (for 1GB then, 3424 years on 232 cores)
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o BvFl1_50kb -a gem -d -s <(for i in $(seq 73000 1 123000); do echo jcf7180008856767:"$i"-$(( $i + 1 )); done;)
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o BvFl1_50kb2 -a gem -d -s <(for i in $(seq 73000 1 123000); do echo jcf7180008856767:"$i"-$(( $i + 1 )); done;)
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o BvFl1_MADS2 -a gem -d -m -s <(for i in $(seq 97315 1 97494); do echo jcf7180008856767:"$i"-$(( $i + 1 )); done;)
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o BvFl1_MADS -a gem -d -x -s <(for i in $(seq 97315 1 97494); do echo jcf7180008856767:"$i"-$(( $i + 1 )); done;)

./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o s1 -a gem -s sites1.txt;

<(for i in $(seq 1 1 305); do echo jcf7180008587925:"$i"-$(( $i + 1 )); done;)
<(for i in $(seq 1 25 544); do echo jcf7180008531951:"$i"-$(( $i + 25 )); done;)
<(for i in $(seq 1 100 157500); do echo jcf7180008856767:"$i"-$(( $i + 100 )); done;)




#some postprocessing
#postprocess haploblock counts in log.txt into something plottable
myq() {
      i=$1; #a position on the reference is incoming
      for j in $b;
      do a=$(grep _"$i"\\."$j" numblocks.txt | cut -d: -f4); #number of unique haploblock read pairs mapped to position i
        if [[ "$a" == "" ]]; then a=0; fi; #if no data at position i set number of haploblock read pairs to 0
        echo "$i $a" >> "$j".txt; #echo to output file specific for the readgroup
        b=$(grep _"$i"\\."$j" numblocks.txt | cut -d$'\t' -f2 | cut -d: -f1); #total number haploblock read pairs mapped to position i
        if [[ "$b" == "" ]]; then b=0; fi; #if no data at position i set to 0
        echo "$i $b" >> "$j".total.txt; #echo to output file specific for the readgroup
      done;
}
export -f myq;

grep ^# log.txt | tr '-' '_' > numblocks.txt;
a=$(cut -d_ -f3 numblocks.txt | cut -d. -f1 | sort -n -u); #determine set of possible sites
b=$(cut -d. -f2 numblocks.txt | cut -d$'\t' -f1 | sort -u); export b; #determine set of possible read groups

for i in $b;
  do echo "position $i"hblox > $i.txt;
  echo "position $i"totalhblox > $i.total.txt;
  done; #initialize output file with a header
echo "$a" | parallel --bar --keep-order --env myq --env b myq;

#sort and tab delimit output
for i in $b;
  do sort -t' ' -k1,1n "$i".txt | tr ' ' '\t' > "$i".2.txt;
    sort -t' ' -k1,1n "$i".total.txt | tr ' ' '\t' > "$i".2.total.txt;
  done;
#swap back to original filename
for i in $b; 
  do mv "$i".2.txt "$i".txt;
    mv "$i".2.total.txt "$i".total.txt;
  done;
  
#consolidate files into 1
paste -d$'\t' global.txt RG_Z_5[0-5].txt global.total.txt RG_Z_*.total.txt | cut -d$'\t' -f1,2,4,6,8,10,12,14,16,18,20,22,24,26,28 > summary.txt;

