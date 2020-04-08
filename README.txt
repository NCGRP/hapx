hapx extracts molecular haplotypes from bam files at user-specified sites

In working folder:
1) A fasta-formatted reference sequence
2) A sorted, deduped bam file containing reads aligned to the reference in (1)
3) A line delimited list of target sites

Requirements (in path):
1) samtools (hapx calls samtools sort/view/mpileup)
2) bwa (mem/index)
3) GNU parallel
4) muscle

Usage: hapx -r ref -b bam -o out -a alnr [-f inc] [-F exc] [-q qual] [-p maxp] [-d -mm -mb -x -db] -s sites
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
-mm = align (muscle) extracted haploblocks
-mb = map (bwa mem) extracted haploblocks
-x = do not write any output files except log (makes -d and -m irrelevant)
-db = debugging mode, save internal data structures as files

Examples: ./hapx.sh /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta /share/space/reevesp/patellifolia/xtr/AllP.merged.bam bwamem jcf7180008454378:303-304,jcf7180008531951:103-495,jcf7180008395354:294-295,jcf7180008827236:1031-1032,jcf7180008378511:277-278,jcf7180008637475:7-8,jcf7180008587925:106-107,jcf7180008527965:177-178,jcf7180008578969:84-85,jcf7180008484650:470-471,jcf7180008856767:98710-98886;

./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -a gem -s jcf7180008454378:303-304,jcf7180008531951:103-495,jcf7180008395354:294-295,jcf7180008827236:1031-1032,jcf7180008378511:277-278,jcf7180008637475:7-8,jcf7180008587925:106-107,jcf7180008527965:177-178,jcf7180008578969:84-85,jcf7180008484650:470-471,jcf7180008856767:98710-98886;
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -a gem -d -s jcf7180008454378:303-304,jcf7180008531951:103-495,jcf7180008856767:98710-98886;

./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -a gem -d -s jcf7180008454378:204-304,jcf7180008531951:395-495,jcf7180008395354:195-295,jcf7180008827236:932-1032,jcf7180008378511:178-278,jcf7180008637475:7-107,jcf7180008587925:7-107,jcf7180008527965:78-178,jcf7180008578969:84-184,jcf7180008484650:371-471,jcf7180008856767:98710-98810;
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o s1 -a gem -d -m -s sites.txt;

#use a function to generate target sites (effectively a sliding window)
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o 1305 -a gem -m -s <(for i in $(seq 1 1 305); do echo jcf7180008587925:"$i"-$(( $i + 1 )); done;);
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o 25544 -a gem -m -s <(for i in $(seq 1 25 544); do echo jcf7180008531951:"$i"-$(( $i + 25 )); done;)
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o 157500 -a gem -m -s <(for i in $(seq 1 100 157500); do echo jcf7180008856767:"$i"-$(( $i + 100 )); done;)
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o 157500d -a gem -d -s <(for i in $(seq 1 100 157500); do echo jcf7180008856767:"$i"-$(( $i + 100 )); done;)
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o 157500dx -a gem -d -x -s <(for i in $(seq 1 100 157500); do echo jcf7180008856767:"$i"-$(( $i + 100 )); done;)
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o 157500step1 -a gem -d -x -s <(for i in $(seq 1 1 157500); do echo jcf7180008856767:"$i"-$(( $i + 1 )); done;)
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o 157500step2 -a gem -d -x -s <(for i in $(seq 1 1 157500); do echo jcf7180008856767:"$i"-$(( $i + 1 )); done;)

#takes about 4 hrs for 50kb (for 1GB then, 9.132 years on 232 cores)
#after improvements, now takes about 1.5 hrs for 50kb (for 1GB then, 3.424 years on 232 cores), 5 hrs for 150kb (for 1GB then, 3.8 yrs on 232 cores)
#after adding freqs, 0.3hrs for 10kb (1GB 3.4 years), 5hrs for 150kb --no noticeable increase in time
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o BvFl1_50kb -a gem -x -s <(for i in $(seq 73000 1 123000); do echo jcf7180008856767:"$i"-$(( $i + 1 )); done;)
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o BvFl1_10kb -a gem -d -s <(for i in $(seq 73000 1 83000); do echo jcf7180008856767:"$i"-$(( $i + 1 )); done;)
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o BvFl1_50kb2 -a gem -d -s <(for i in $(seq 73000 1 123000); do echo jcf7180008856767:"$i"-$(( $i + 1 )); done;)
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o BvFl1_MADS2 -a gem -d -m -s <(for i in $(seq 97315 1 97494); do echo jcf7180008856767:"$i"-$(( $i + 1 )); done;)
./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o BvFl1_MADS -a gem -d -x -s <(for i in $(seq 97315 1 97494); do echo jcf7180008856767:"$i"-$(( $i + 1 )); done;)

./hapx.sh -r /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -b /share/space/reevesp/patellifolia/xtr/AllP.merged.gem.bam -o s1 -a gem -s sites1.txt;

<(for i in $(seq 1 1 305); do echo jcf7180008587925:"$i"-$(( $i + 1 )); done;)
<(for i in $(seq 1 25 544); do echo jcf7180008531951:"$i"-$(( $i + 25 )); done;)
<(for i in $(seq 1 100 157500); do echo jcf7180008856767:"$i"-$(( $i + 100 )); done;)




#some postprocessing
#postprocess haploblock counts in log.txt into something plottable
mytd1() {
        i=$1; #a position on the reference is incoming
        a=$(grep "$i"\\."$k" "$pd"/numblocks.txt | cut -d: -f4); #number of unique haploblock read pairs mapped to position i in the readgroup
        if [[ "$a" == "" ]]; then a="?"; fi; #if no data at position i set number of haploblock read pairs to 0
        echo "$i"."$k $a"; #report result to parallel statement
}
export -f mytd1;

mytd2() {
        i=$1; #a position on the reference is incoming
        b=$(grep "$i"\\."$k" "$pd"/numblocks.txt | cut -d$'\t' -f2 | cut -d: -f1); #total number haploblock read pairs mapped to position i in the readgroup
        if [[ "$b" == "" ]]; then b="?"; fi; #if no data at position i set to 0
        echo "$i"."$k $b"; #report result to parallel statement
}
export -f mytd2;

mypa() {
       i=$1; #contig:site-range
       j=$(grep "$i" "$pd"/counts.txt | grep -v global);
       names=$(cut -d$'\t' -f1 <<<"$j" | sed 's/$/ 0/');
       counts=$(cut -d$'\t' -f2 <<<"$j");
       nc=$(head -1 <<<"$counts" | awk -F: '{print NF}'); #number of alleles
       nr=$(wc -l <<<"$counts"); #number of read groups
       nz=$(( $nr - 1 )); #number of zeroes needed in a column of allele counts for a private allele to exist
       for m in $(seq 1 1 $nc);
         do col=$(cut -d: -f$m <<<"$counts"); #extract the column of allele counts
           fl=$(grep -n -v 0 <<<"$col"); #show line numbers where alleles are present
           
           #if only one line of the column of data has alleles, it is a private allele
           if [[ $(echo "$fl" | wc -l) == 1 ]];
           then r=$(cut -d: -f1 <<<"$fl"); #row number of private allele
             names=$(awk -F' ' -v r=$r '{if (NR==r) $2++; print}' <<<"$names"); #index up by one the row with the private allele
           fi;
         done;
         
       #report result to parallel statement
       echo "$names";
}
export -f mypa;

mypadpa() {
          aa=$1;
          cc=$(grep "$aa.$dd" "$pd"/names.txt);
          if [[ "$cc" == "" ]];
          then echo "$aa.$dd ?"; #if readgroup not found print a ?
          else echo "$cc";
          fi;
}
export -f mypadpa;


pd=$(pwd); export pd;

#count total and distinct haploblocks
grep ^# log.txt | tr '-' '_' > numblocks.txt;
a=$(cut -d. -f1 numblocks.txt | uniq); #list of contig:site-ranges
b=$(cut -d. -f2 numblocks.txt | cut -d$'\t' -f1 | sort -u); export b; #determine set of possible read groups

#iterate over read groups to count haploblocks
for j in $b;
  do k="$j"; export k; #do this so iterator can be sent to nodes using --sshloginfile
    echo "$j: counting distinct alleles";
    >"$pd"/"$j".distallel.txt; #initialize output file with a header for count of distinct alleles within populations
    echo "$a" | parallel --sshloginfile ~/machines --jobs 1 --pipe -N960 --env mytd1 --env pd --env k /home/reevesp/bin/parallel --jobs 96 --env mytd1 --env pd --env k mytd1 >> "$pd"/"$j".distallel.txt;
    #echo "$a" | parallel --jobs 1 --pipe -N960 --env mytd1 --env pd --env j /home/reevesp/bin/parallel --jobs 96 --env mytd1 --env pd --env j mytd1 >> "$pd"/"$j".distallel.txt;
    
    echo "$j: counting all alleles";
    >"$pd"/"$j".totallel.txt; #initialize output file with a header for count of all alleles within populations
    echo "$a" | parallel --sshloginfile ~/machines --jobs 1 --pipe -N960 --env mytd2 --env pd --env k /home/reevesp/bin/parallel --jobs 96 --env mytd2 --env pd --env k mytd2 >> "$pd"/"$j".totallel.txt;
  done;


#count private alleles
grep ^'@' log.txt | tr '-' '_' | cut -d$'\t' -f1,3 | sort -t_ -k2,2n > freqs.txt;
grep ^'@' log.txt | tr '-' '_' | cut -d$'\t' -f1,2 | sort -t_ -k2,2n > counts.txt;

c=$(cut -d. -f1 counts.txt | uniq); #list of contig:site-ranges
>names.txt;
echo "$c" | parallel --sshloginfile ~/machines --jobs 1 --pipe -N960 --env pd --env mypa /home/reevesp/bin/parallel --jobs 96 --env pd --env mypa mypa >> names.txt; #counting sub

#pad the file containing counts of private alleles with "?" for contig:site-range_readgroup combinations that weren't found
for bb in $b;
do >"$bb".privallel.txt; #create an output file for the readgroup
  dd="$bb"; export dd; #transfer iterator to variable that can be exported for gnu parallel
  echo "$a" | sed 's/#/@/g' | parallel --sshloginfile ~/machines --jobs 1 --pipe -N960 --env pd --env dd --env mypadpa \
                              /home/reevesp/bin/parallel --jobs 96 --env pd --env dd --env mypadpa mypadpa >> "$bb".privallel.txt;
done;

#sort and tab delimit output
for i in $b;
  do sort -t'_' -k1,1 "$i".distallel.txt | sort -t'_' -k2,2n | sed 's/^#//' | sed "1i position $i.distallel" | tr ' ' '\t' > "$i".2.distallel.txt;
    sort -t'_' -k1,1 "$i".totallel.txt | sort -t'_' -k2,2n | sed 's/^#//' | sed "1i position $i.totallel" | tr ' ' '\t'  > "$i".2.totallel.txt;
    sort -t'_' -k1,1 "$i".privallel.txt | sort -t'_' -k2,2n | sed 's/^@//' | sed "1i position $i.privallel" | tr ' ' '\t'  > "$i".2.privallel.txt;
  done;
#swap back to original filename
for i in $b; 
  do mv "$i".2.distallel.txt "$i".distallel.txt;
    mv "$i".2.totallel.txt "$i".totallel.txt;
    mv "$i".2.privallel.txt "$i".privallel.txt;
  done;


#consolidate files into 1
paste -d$'\t' global.distallel.txt RG_Z_*.distallel.txt global.totallel.txt RG_Z_*.totallel.txt global.privallel.txt RG_Z_*.privallel.txt\
  | cut -d$'\t' -f1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42 | sed 's/\.global//'> summary.txt;

