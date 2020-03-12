Notes on using the software bazam to extract read pairs or full length pacbio/10x reads
covering a user-defined site of interest.

#create paired end alignment
rsync -aP admin@10.177.9.14:"/share/Public/Data/PatReeves/PatellifoliaIlluminaData/*.fastq.gz" /scratch/reevesp/patellifolia/data; #get the raw data for webbiana


#trim reads of adapters and initial quality trimming. Use popoolation recommended Q20 filter, and minlen=50
#AdapterRead1 is found in *_R1_001.fastq reads, AdapterRead2 is found in *_R2_001.fastq reads.
cd /scratch/reevesp/patellifolia/data;
echo ">AdapterRead1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>AdapterRead2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" > adapters.fa;
	
mypp() {
        i=$1;
        echo "$i";
        java -jar ~/bin/trimmomatic.jar PE -threads 40 17134D-01-$i*_L001_R1_001.fastq.gz 17134D-01-$i*_L001_R2_001.fastq.gz \
          "$i"output_forward_paired.fq.gz "$i"output_forward_unpaired.fq.gz "$i"output_reverse_paired.fq.gz "$i"output_reverse_unpaired.fq.gz \
          ILLUMINACLIP:adapters.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50;
}
export -f mypp;
seq 50 55 | parallel --env mypp mypp; #use parallel because multithreading doesn't work.

#clean up unpaired reads
rm data/5[012345]output_*_unpaired.fq.gz;


#Perform paired end alignments for each sample using bwa-mem, and markdups using sambamba
#This saves all output, intermediate and otherwise, to the local node. Things can get messy if it fails.
mypp() { 
       i=$1;
       f="$wd"/"$i"output_forward_paired.fq.gz; #forward reads
       r="$wd"/"$i"output_reverse_paired.fq.gz; #reverse reads
       thr=$(lscpu | grep "^CPU(s):" | awk '{print $2}'); #max threads
       tmpd="/state/partition1/tmpbwa$i";
       if [ ! -d "$tmpd" ]; then mkdir "$tmpd"; fi; #make a local tmp directory on the main disk
       readgroup="@RG\tID:$i\tSM:$i\tPL:illumina\tLB:na\tPU:na";

       #do bwa-mem alignment, save result to local directory
       #stdout from bwa mem goes to samtools sort, a coordinate sorted bam file is the result
       # -k 19 is the default seed length, it can be varied below
       /home/reevesp/bin/bwa mem -R "$readgroup" -t "$thr" -k 19 "$v" "$f" "$r" | /home/reevesp/bin/samtools sort -O BAM --threads "$thr" -T "$tmpd"/tmp.$i.bam -o "$tmpd"/$i"_aligned_reads.bam"; 

       #markdups
       #x=$(ls /state/partition1/tmpbwa$i/"$i"_aligned_reads.bam | rev | cut -d_ -f3 | cut -c1-2 | rev); #get the sample number of the aligned reads file on this node (in case nodes have changed due to failure and repeat)
       /home/reevesp/bin/sambamba markdup -t "$thr" --tmpdir="$tmpd" --overflow-list-size 6000000 --remove-duplicates "$tmpd"/$i"_aligned_reads.bam" "$tmpd"/$i"_markdups.bam";
}
export -f mypp;
cd /share/space/reevesp/patellifolia/data;
v="/share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta"; export v;
#v="/share/space/reevesp/patellifolia/ref/PatellifoliaSampleRef.fasta"; export v;
wd=$(pwd); export wd;
seq 50 1 55 | parallel --sshloginfile ~/machinesbwa --env v --env wd --env mypp mypp;

#retrieve results from tmp directories of remote nodes
for i in {0..9};
  do echo -n "compute-0-$i ";
    ssh compute-0-$i 'ls -l /state/partition1/tmpbwa*/*_markdups.bam 2>/dev/null';
  done;
#results present on compute-0-[012347], rsync them to head node
for i in 0 1 2 3 4 7;
  do rsync -aP compute-0-$i:"/state/partition1/tmpbwa*/5[012345]_*.bam*" /share/space/reevesp/patellifolia/map;
  done;
##then go in and delete the tmp directory contents### 
seq 0 1 9 | parallel ssh compute-0-{} "'rm -r /state/partition1/tmpbwa*'";


#Combine sorted, deduped bam files, index, then calculate depth:
cd /share/space/reevesp/patellifolia/map;
time(sambamba merge -t40 -p AllP.merged.bam 5[012345]_markdups.bam;)
#time(sambamba index -t40 -p AllP.merged.bam AllP.merged.bai;)
time(sambamba depth base -t40 AllP.merged.bam -o AllP.merged.coverage;) #~ 4hrs



#use bazam to extract read pairs at site (bazam doesn't work right)
cd /share/space/reevesp/patellifolia/map;
time(java -jar /home/reevesp/bin/bazam/build/libs/bazam.jar -bam AllP.merged.bam --regions jcf7180008368139:50-51 > bazam1.fa;) #bazam is doing something other than advertised. it returns >12E6 read pairs for a span that doesn't exist in the alignment

#use pierre lindenbaum's samviewwithmate.jar to extract read pairs at site (doesn't work right either)
java -jar /home/reevesp/bin/jvarkit/dist/samviewwithmate.jar --region "jcf7180008359438:50-51" AllP.merged.bam > samvm.sam; #Lindenbaums samviewwithmate.jar produces some weird results inconsistent with samtools view and igv, namely sequences that don't show up with either, probably these are sequences that have a primary alignment somewhere else.  Anyway, don't use samviewwithmate.jar.
grep -v ^'@' samvm.sam | wc -l; #39 reads
grep -v ^'@' samvm.sam | sort | cut -c1-70

#do some test extractions with samtools to figure out the proper flags to use
samtools view  53filPproseDedup.bam "jcf7180008368139:50-51" > samout2.tmp; #include all reads (85 reads)
samtools view -f 0x2 53filPproseDedup.bam "jcf7180008368139:50-51" > samout2.tmp; #include only proper pairs (no result
samtools view -f 0x1 53filPproseDedup.bam "jcf7180008368139:50-51" > samout1.tmp; #include only paired reads (no result
samtools view  -f 0x8 53filPproseDedup.bam "jcf7180008368139:50-51" > samout8.tmp; #include reads whose mate is unmapped (0 reads)
samtools view  -f 0x9 53filPproseDedup.bam "jcf7180008368139:50-51" > samout9.tmp; #include paired reads and reads whose mate is unmapped (0 reads)
samtools view  -f 4 53filPproseDedup.bam "jcf7180008368139:50-51" > samout4.tmp; #include only unmapped reads (0 reads)

samtools view 53_markdups.bam "jcf7180008359438:50" > samout0.tmp; #include all reads (33 reads)
samtools view -f 0x2 53_markdups.bam "jcf7180008359438:50"  > samout2.tmp; #include only proper pairs (8 reads)
samtools view -f 0x1 53_markdups.bam "jcf7180008359438:50"  > samout1.tmp; #include only paired reads (33 reads)
samtools view  -f 0x8 53_markdups.bam "jcf7180008359438:50"  > samout8.tmp; #include reads whose mate is unmapped (0 reads)
samtools view  -f 0x9 53_markdups.bam "jcf7180008359438:50"  > samout9.tmp; #include paired reads and reads whose mate is unmapped (0 reads)
samtools view  -f 4 53_markdups.bam "jcf7180008359438:50"  > samout4.tmp; #include only unmapped reads (0 reads)
samtools view  -f 11 53_markdups.bam "jcf7180008359438:50"  > samout11.tmp; #include paired end reads, proper pairs and pairs where next segment is unmapped (0 reads)
samtools view  -f 15 53_markdups.bam "jcf7180008359438:50"  > samout15.tmp; #include paired end reads, proper pairs and pairs where next segment is unmapped (0 reads)

samtools view -F 4 AllP.merged.bam "jcf7180008359438:50" > asamout0.tmp; #include all reads (57 reads <- this includes "improper" pairs)
samtools view -f 2 -F 4 AllP.merged.bam "jcf7180008359438:50" > asamout2.tmp; #proper pairs (23 reads) In bwa a proper pair -probably- is one that both ends are mapped within a reasonable distance of one another (ie surely not on different chromosomes)
samtools view -f 8 -F4 AllP.merged.bam "jcf7180008359438:50" > asamout8.tmp; #reads whose mate is unmapped (1 reads, also included in asamout0.tmp

#settle on the following samtools command (-f2 include only proper pairs, -f4 exclude unmapped reads, Q20 cutoff (retains ~80% of reads alignment wide)
samtools view -f 2 -F 4 -q 20 AllP.merged.bam "jcf7180008359438:50" > asamout2q20.tmp; #proper pairs with mapping quality > q (pr >0.999 Q30 7 reads, 3 proper pairs) (Q20, 9 reads, 4 proper pairs)

#get a sample of mapping quality values from the bam file (important values are 0,27,40,60.  I like 20 as the cutoff, retains ~80% of aligned reads:
samtools view AllP.merged.bam | head -1000000 | awk -F$'\t' '{print $5}' | sort | uniq -c | sed 's/^ *//g' | sort -t' ' -k2,2n | awk -F' ' '{print $2,$1}' | tr ' ' '\t';  
0	142220
1	4011
2	3145
3	5213
4	4232
5	3512
6	4434
7	4183
8	2641
9	3699
10	1697
11	1300
12	2092
13	2383
14	1659
15	3509
16	1887
17	2117
18	2916
19	2885
20	1816
21	4564
22	3097
23	2562
24	3348
25	5666
26	1333
27	8499
28	1432
29	1506
30	1926
31	2241
32	1324
33	2554
34	1193
35	992
36	1713
37	1882
38	692
39	3223
40	35094
41	1819
42	1835
43	2436
44	1874
45	3137
46	2754
47	2979
48	5131
49	1795
50	1790
51	1880
52	1791
53	1313
54	1600
55	1834
56	1377
57	3188
58	2077
59	2170
60	670828

###BELOW CODE MIGRATED TO GITHUB PROJECT 'hapx' to build executable script

#randomly sample haplotypes at a few positions
#get list of contigs and length, shuffle it for random sampling (samtools idxstats output: ref,length,mapped reads,unmapped reads)
samtools idxstats AllP.merged.bam | shuf > AllP.idxstats.txt; #plot distribution in ReadsPerContig.png, chloroplast sequences have very high reads per contig (>1E6 reads <20kb contig)

nl=10; #number of loci desired, these will be on different contigs
a=$(head -$nl AllP.idxstats.txt);

#select a site at random from each contig
b=$(for i in $(echo "$a" | cut -d$'\t' -f2 | tr "\n" " ");
  do shuf -i 1-"$i" -n 1; #use shuf to emit one random number from the range
  done;
);

c=$(paste -d$'\t' <(echo "$a") <(echo "$b")); #place random site in column 5
echo "$c" > AllP.randomsites"$nl".txt; #save the randomly chosen sites in case you need them later

#get list like contig name:random site-random site +1 (e.g. jcf7180008527965:177-178), substitute file for "$c" AllP.randomsites"$nl".txt input if variable has been lost
e=$(while read j;
do d=$(echo "$j" | cut -d$'\t' -f1,5 | tr "\t" ":");
  dd=$(( $(echo "$d" | cut -d: -f2) + 1 ));
  echo "$d"-"$dd"; #format like contig1:167-168
done <<<"$c";)




#extract read pairs at site using samtools. Extract proper pairs with mapping quality > 60 excluding unmapped pairs
echo "$e" | parallel 'samtools view -f 2 -F 4 -q 60 AllP.merged.bam "{}" | sort > {}.tmp; \
  mv {}.tmp $(echo {} | tr ":" "_").tmp';
  
#describe number of read pairs and single ends of a read pair extracted in a tabular format
f=$(while read j;
  do
    fn=$(echo "$j" | tr ":" "_")".tmp"; #name of input file
    p=$(cut -d$'\t' -f1 "$fn" | sort | uniq -c | grep ^" \+2 "); #paired reads
    s=$(cut -d$'\t' -f1 "$fn" | sort | uniq -c | grep ^" \+1 "); #single ends of paired reads
    pp=$(echo "$p" | awk 'NF' | wc -l); #number of paired reads (awk 'NF' removes empty lines)
    ss=$(echo "$s" | awk 'NF' | wc -l); #number of single ends of paired reads
    echo "$j $pp $ss";
  done <<<"$e";)
echo "contig:site pairedreads singlereads"$'\n'"$f";

			jcf7180008454378:303-304 0 8
			jcf7180008531951:103-104 4 45
			jcf7180008395354:294-295 3 55
			jcf7180008827236:1031-1032 0 19
			jcf7180008378511:277-278 4 39
			jcf7180008637475:7-8 0 2
			jcf7180008587925:106-107 2 7
			jcf7180008527965:177-178 1 17
			jcf7180008578969:84-85 7 39
			jcf7180008484650:470-471 0 7


g=$(echo "$f" | awk -F' ' '$2+$3!=0 {print $0}'); #remove any contigs from consideration when there are 0 reads (no paired reads:$2, no unpaired reads:$3) that align to them (this could also be a variable that supports a cutoff)

#extract each contig from the reference genome to be used by itself as a reference for realignment with the extracted read pairs, then bwa index
echo "$g" | cut -d: -f1 | parallel 'sed -n -e "/^>{}$/,/>/p" Ppanfinal.genome.scf.fasta | sed "$ d" > {}_ref.txt; bwa index {}_ref.txt'; #sed extracts lines from name of contig of interest until next contig name line, second sed deletes the last line

#reconstruct fastq file for reads that aligned to each contig, adding readgroup and f/r designation, r is for reads marked with a negative alignment length in column 9
#then sort reads into individual read pairs, resulting files are like {contigname_readname.rp.fq}
#rm *rp.fq; #you should remove the old read pair output files before proceeding here
for i in $(echo "$g" | cut -d' ' -f1 | tr ':' '_' | tr '\n' ' ');
  do awk -F$'\t' '{print "@"$1,$9,$10,"+",$11,$17}' "$i".tmp | awk -F' ' '{if ($2 > 0) {print $1"_"$6":f",$3,$4,$5} else if ($2 < 0) {print $1"_"$6":r",$3,$4,$5}}' | tr " " "\n" > "$i"_reads.fq;

    m=$(grep ^'@' "$i"_reads.fq | cut -d: -f1-9 | sort -u); #capture list of unique read pair names (includes also unpaired reads with unique names) 
    for j in $m;
      do echo "$i / $j";
        fn=$(echo "$j" | sed 's/@//g' | tr ':' '_'); #filename for output
        grep -A3 "$j" "$i"_reads.fq > "$i"_"$fn".rp.fq; #save read pair into its own file
      done;
    
    rm "$i"_reads.fq; #clean up
    rm "$i".tmp;
  done;

			Here is a sample reconstructed *.rp.fq file, note readgroup in col9 of name, and "orientation" or "endedness of read pair" in col10:
			@E00558:144:HHGCMCCXY:5:2224:6654:40055_RG:Z:55:r
			CAAAAATATTTTGGTAATTATTCTCAACAAAATGATTTGAAAGGTGTTCATACAACACAAATCGCCTAAGAGACTATGACGGTTTTATCCTCTGATTTGAATTGAGTTTGATCCAAGGGCTTCATATGATTGAAATATAC
			+
			FFJJJJAFJ<JJFF7FJAJAA-A7AJJJJJJ<JJJJJJJJJJFJJJJJJJJJJJ<JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFFFAA
			@E00558:144:HHGCMCCXY:5:2224:6654:40055_RG:Z:55:f
			ATAAATATGCAAGGAGTCAAAAATATTTGGTAATTAAGCACAAAAACTGATTTGAAATGTGTTCGTACACCACAAATCACCTAAGAGACTATGACGGTTTTACCCTTTGATTTGAATTGAGTTTGATCCAAGGGCTTCATATGCTTTAAA
			+
			AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJ

#align *.rp.fq to *_ref.txt, output as bam file, then sambamba index (output files have names like jcf7180008378511_277.bam (contigname_site.bam)
#this step aligns, then verifies mapping quality since it changes from the original values, then realigns
#rm *_RG_Z_5[012345].[bs]am*; #clean up old files before repeating this step
#rm *2x*;
#This step should probably be parallelized
for i in $(find . -name "*.rp.fq" | cut -d'/' -f2 | sed 's/\.rp\.fq//g' | tr "\n" " ");
  do echo "$i";
    h=$(echo "$i" | cut -d_ -f1); #just the contig name
    bwa mem "$h"_ref.txt "$i".rp.fq > "$i".sam; #write out human readable sam file
    samtools sort "$i.sam" -O BAM -o "$i".bam; #convert sam to bam
    sambamba index "$i".bam "$i".bam.bai;

    #mapping quality changes upon realignment, extract only reads with new mapping quality > 60, excluding unmapped, convert those to proper fastq file
    samtools view -F 4 -q 60 "$i".bam | awk -F$'\t' '{print "@"$1,$10,"+",$11}' | tr " " "\n" > "$i".2xrp.fq;
    if [[ $(cat "$i".2xrp.fq) != "" ]]; #only re-realign if there are reads left after latest quality filter
    then bwa mem "$h"_ref.txt "$i".2xrp.fq > "$i".2x.sam; #perform second realignment
      #samtools sort "$i.2x.sam" -O BAM -o "$i".2x.bam; #convert sam to bam
      #sambamba index "$i".2x.bam "$i".2x.bam.bai;
    fi;
    
    #clean up
    rm "$i".sam "$i".bam "$i".bam.bai "$i".2xrp.fq "$i".rp.fq; #clean up
  done;


#print mapping quality after 1st realignment (it changes from original mapping quality)
awk -F$'\t' '{print $5}' *.sam | grep -v ^CL | grep -v ^$ | sort | uniq -c;
#print mapping quality after RE-realignment (all 60)
awk -F$'\t' '{print $5}' *.2x.sam | grep -v ^CL | grep -v ^$ | sort | uniq -c;

#figure out the possibilities for coding in the nucleotide column ($5)
> xxx.tmp;
for i in $(find . -name "jcf*_RG_Z_5[012345].2x.sam");
  do echo "$i" >> xxx.tmp;
    samtools mpileup "$i" >> xxx.tmp;
  done;
awk -F$'\t' '{print $4,$5}' xxx.tmp | sort | uniq -c; #print all unique possibilities for nucleotide coding, these have to be accounted for when assembling the consensus





#compute a consensus sequence for each aligned read pair
conseq() {
          echo "$1" | awk '$1!=p+1{print p+1"-"$1-1}{p=$1}';
}
export -f conseq;

#iupac() replaces character strings with their iupac equivalent or proper insertion and deletion coding.
#deletions relative to reference are given an 'x', insertions are given their bases
myiupac() {
        read i; #input is a piped row of data
        a=$(echo "$i" | tr " " "\n" | cut -d$'\t' -f1); #get the position numbering
        b=$(echo "$i" | tr " " "\n" | cut -d$'\t' -f2 \
          | sed 's/^AA$/A/g' | sed 's/^CC$/C/g' | sed 's/^GG$/G/g' | sed 's/^TT$/T/g' \
          | sed 's/^AC$/M/g' | sed 's/^CA$/M/g' \
          | sed 's/^AG$/R/g' | sed 's/^GA$/R/g' \
          | sed 's/^AT$/W/g' | sed 's/^TA$/W/g' \
          | sed 's/^CG$/S/g' | sed 's/^GC$/S/g' \
          | sed 's/^CT$/Y/g' | sed 's/^TC$/Y/g' \
          | sed 's/^GT$/K/g' | sed 's/^TG$/K/g' \
          | sed 's/^\*$/x/g' | sed 's/^\**$/x/g' \
          | sed 's/^A\*$/A/g' | sed 's/^\*A$/A/g' \
          | sed 's/^C\*$/C/g' | sed 's/^\*C$/C/g' \
          | sed 's/^G\*$/G/g' | sed 's/^\*G$/G/g' \
          | sed 's/^T\*$/T/g' | sed 's/^\*T$/T/g');
          paste -d' ' <(echo "$a") <(echo "$b"); #recombine position numbers with consensus calls
}
export -f myiupac;

myinsertion() {
              read i;
              i=$(echo "$i" | tr ' ' '\n'); #delimit with line breaks
              j=$(echo "$i" | grep "+"); #collect all lines with insertions
              ij=$(echo "$i" | grep -v "+"); #collect all lines without insertions

              #process lines with insertions ("$j"), save in a variable $lwi
              if [[ "$j" == "" ]];
              then lwi=""; #there are no insertions so create a dummy, empty variable
              else
                lwi=$(while read str;
                  do pos=$(echo "$str" | cut -d$'\t' -f1); #get the starting position of the insertion
                    m=$(echo "$str" | cut -d$'\t' -f2 | tr '+' '\n'); #break into lines on +
                    l1len=$(echo "$m" | head -1 | awk '{print length}'); #length of line one
                    nl=$(echo "$m" | wc -l); #number of lines
                    l2=$(echo "$m" | head -2 | tail -1); #get line 2 input
                    l2inslen=$(echo "$l2" | tr -d -c 0-9); #deletes all non-numeric characters from line 2 to read the length of the insertion
                    l2nchar=$(echo "$l2" | tr -d 0-9 | awk '{ print length }'); #number of chars in line 2, counts all nun numeric characters
                    
                    if [[ $l1len == 1 ]]; #case where there in an insertion in the first read there is one character on line 1. here, you need to get character from first line plus last character of the second line, if present
                    then if [[ $nl == 2 ]];
                         then if [[ $l2nchar == $l2inslen ]]; #test whether there is a base from the second read hanging off the end of the insertion in line 2
                              then m=$(echo "$m" | sed '2s/$/\*/'); #add an asterisk to the end of the second line to act as the unknown read 2 base
                              fi;
                           l2=$(echo "$m" | head -2 | tail -1); #recalculate line 2 input now that it has been modified with an asterisk
                           n=$(paste -d' ' <(echo "$m" | head -1) <(echo "$l2" | rev | cut -c1) | sed 's/ //'); #char states composed of line 1 + last character of line 2
                           for i in $(seq 1 1 $l2inslen);
                             do n=$(echo "$n";echo "$l2" | tr -d 0-9 | cut -c$i)"*"; #combine character n of line 2 with an asterisk to form the read pair
                             done;
                          
                         elif [[ $nl == 3 ]]; #if there are three lines then there are inserts after both reads, which means that the last character of line 2 is always read 2 character 1
                         then n=$(paste -d' ' <(echo "$m" | head -1) <(echo "$l2" | rev | cut -c1) | sed 's/ //'); #char states composed of line 1 + last character of line 2
                           l3=$(echo "$m" | head -3 | tail -1); #get line 3 input
                           l3inslen=$(echo "$l3" | tr -d -c 0-9); #deletes all non-numeric characters from line 3 to read the length of the insertion
                    
                           if [[ $l2inslen == $l3inslen ]]; #case where insertion lengths are the same
                           then for i in $(seq 1 1 $l2inslen);
                             do n=$(echo "$n";(echo "$l2" | tr -d 0-9 | cut -c$i;echo "$l3" | tr -d 0-9 | cut -c$i) | tr -d '\n'); #combine character n of line 2 with character n of line 3 to form the read pair
                             done;
                           elif [[ $l2inslen < $l3inslen ]]; #case where read1 (l2) insertion length is less than read 2 (l3)
                           then for i in $(seq 1 1 $l2inslen);
                             do n=$(echo "$n";(echo "$l2" | tr -d 0-9 | cut -c$i;echo "$l3" | tr -d 0-9 | cut -c$i) | tr -d '\n'); #combine character n of line 2 with character n of line 3 to form the read pair
                             done;
                             #now combine asterisk for line 2 (which is shorter) with character n of line 3 to form the read pair
                             for i in $(seq $(($l2inslen+1)) 1 $l3inslen);
                               do n=$(echo "$n";echo -n "*";echo "$l3" | tr -d 0-9 | cut -c$i); 
                               done;
                           elif [[ $l2inslen > $l3inslen ]]; #case where read1 (l2) insertion length is greater than read 2 (l3)
                           then for i in $(seq 1 1 $l3inslen);
                             do n=$(echo "$n";(echo "$l2" | tr -d 0-9 | cut -c$i;echo "$l3" | tr -d 0-9 | cut -c$i) | tr -d '\n'); #combine character n of line 2 with character n of line 3 to form the read pair
                             done;
                             #now combine character n of line 2 with asterisk for line 3 (which is shorter) to form the read pair
                             for i in $(seq $(($l3inslen+1)) 1 $l2inslen);
                               do n=$(echo "$n";echo "$l2" | tr -d 0-9 | cut -c$i)"*"; 
                               done;
                           fi;    
                         fi;
                      
                    elif [[ $l1len == 2 ]];
                      then n=$(echo "$m" | head -1);
                        for i in $(seq 1 1 $l2inslen);
                          do n=$(echo "$n";echo -n "*";echo "$l2" | tr -d 0-9 | cut -c$i); #place an asterisk for read 1 in combination with character n of line 2 for read 2
                          done;
                    else echo "ERROR: incompatible string: $m";
                    fi;
                    
                    nls=$(( $(echo "$n" | wc -l) - 1 )); #count the number of lines in the processed insertion, decrement by one since first line will not get a decimal place
                    #calculate the line numbering scheme for the processed insertion, will be like: 474 474.1 474.2
                    lz=$(echo "$pos";
                      for k in $(seq 1 1 $nls);
                        do echo "$pos.$k";
                        done;)
                        
                    #paste line numbers to processed insertion base calls and print out.
                    paste -d$'\t' <(echo "$lz") <(echo "$n");
                    
                  done <<< "$j";); #collect processed insertions in variable $lwi
              fi;
              
              #combine lines with insertions with those without, sort, print out without line breaks to be sent to myiupac in a pipe
              echo "$ij"$'\n'"$lwi" | sort -t$'\t' -k1,1n | tr "\n" " ";

}
export -f myinsertion;

#Process alignments of paired reads via mpileup into a single padded consensus sequence representing the haplotype
#Could also parallelize this
for i in $(find . -name "*.2x.sam" | tr "\n" " ");
  do samin="$i";
    echo "$samin";
    rname=$(echo "$samin" | rev | cut -d'/' -f1 | rev | sed 's/\.sam//g'); #get a root name for the read pair/contig alignment
    mp=$(samtools mpileup "$samin"); #form the pileup file
    
    #Extract the pos and base columns, remove read start ^., read end $, remove deletion indicators (e.g. -2AC, they are followed with *), pass to myinsertion subroutine to control insertions, then to myiupac to recode conflicts as ambiguous and produce the consensus sequence
    base1=$(echo "$mp" | cut -d$'\t' -f2,5 | sed 's/\^.//g' | sed 's/\$//g' | sed 's/-.*//g' | tr "\n" " " | myinsertion | myiupac); #first run at haplotype extraction, unpadded as of now
    
    #Deal with contiguity and padding
    csq=$(conseq "$(echo "$mp" | cut -d$'\t' -f2)" | grep -v "^1\-" | sed 's/-/ 1 /g'); #identify regions where reads are not aligned consecutively to the reference, exclude the range from 1-start of overlap, set up to use as in interval for seq command
    
    if [[ "$csq" == "" ]];
    then pads=""; #no pads if no non-contiguous sections
    else pads=$(while read isq;
      do for i in $(seq $isq);
           do echo "$i N";
           done;
      done <<<"$csq"; #possible multi lines of non-contiguous regions that need to be padded with NNNNs
      );
    fi;
    
    #combine to pad haplotype with respect to reference, convert x for deletion to - for deletion
    #base 2 contains the processed haplotype in a pileup like format relative to reference
    base2=$(echo "$base1"$'\n'"$pads" | sort -t$' ' | sed 's/x/-/g' | awk 'NF');
    
    #base3 contains the processed haplotype as a fasta file, deletions relative to reference removed.
    #base3 is as close as we can come to reconstructing the native molecule
    base3=$(echo ">$rname";echo "$base2" | cut -d' ' -f2 | grep -v '-' | tr "\n" " " | sed 's/ //g');
    echo "$base3" > "$rname.fa";
    
    #clean up
    rm "$samin";
  done;










###ALIGNMENT WITH GEM###
#index reference
gem-indexer -t24 -i /share/space/reevesp/patellifolia/ref/Ppanfinal.genome.scf.fasta -o /share/space/reevesp/patellifolia/ref/PpanfinalGEMref;
gem-indexer -t24 -i /share/space/reevesp/patellifolia/ref/PatellifoliaSampleRef.fasta -o /share/space/reevesp/patellifolia/ref/PatellifoliaSampleGEMref;

#perform alignment
gem-mapper --index=/share/space/reevesp/patellifolia/ref/PatellifoliaSampleGEMref.gem \
    --i1=/share/space/reevesp/patellifolia/data/51output_forward_paired.fq.gz \
    --i2=/share/space/reevesp/patellifolia/data/51output_reverse_paired.fq.gz \
    --gzip-output \
    --output=/share/space/reevesp/patellifolia/gem/51on10genes.gem.sam.gz;
    

#This saves all output, intermediate and otherwise, to the local node. Things can get messy if it fails.
mypp() { 
       i=$1;
       f="$wd"/"$i"output_forward_paired.fq.gz; #forward reads
       r="$wd"/"$i"output_reverse_paired.fq.gz; #reverse reads
       tmpd="/state/partition1/tmpGEM$i";
       thr=$(lscpu | grep "^CPU(s):" | awk '{print $2}'); #max threads
       if [ ! -d "$tmpd" ]; then mkdir "$tmpd"; fi; #make a local tmp directory on the main disk
       readgroup="@RG\tID:$i\tSM:$i\tPL:illumina\tLB:na\tPU:na";

       #do GEM alignment, save result to local directory
       /share/apps/gem-mapper --index="$v" \
           --i1="$f" \
           --i2="$r"\
           --sam-read-group-header="$readgroup" \
           --gzip-output \
           --output="$tmpd"/$i"_aligned_reads.gem.sam";

       #sam to bam, then index
       /home/reevesp/bin/samtools sort -O BAM -T "$tmpd"/tmp.$i.bam --threads "$thr" -o "$tmpd"/$i"_aligned_reads.gem.bam" "$tmpd"/$i"_aligned_reads.gem.sam";

       #markdups
       #x=$(ls /state/partition1/tmp/*_aligned_reads.gembam | rev | cut -d_ -f3 | cut -c1-2 | rev); #get the sample number of the aligned reads file on this node (in case nodes have changed due to failure and repeat)
       /home/reevesp/bin/sambamba markdup -t "$thr" --tmpdir="$tmpd" --overflow-list-size 6000000 --remove-duplicates "$tmpd"/$i"_aligned_reads.gem.bam" "$tmpd"/$i"_markdups.gem.bam";
}
export -f mypp;
cd /share/space/reevesp/patellifolia/data;
v="/share/space/reevesp/patellifolia/ref/PpanfinalGEMref.gem"; export v;
#v="/share/space/reevesp/patellifolia/ref/PatellifoliaSampleGEMref.gem"; export v;
wd=$(pwd); export wd;
seq 55 1 55 | parallel --sshloginfile ~/machinesgem --env v --env wd --env mypp mypp;

#retrieve results from tmp directories of remote nodes
for i in {0..9};
  do echo -n "compute-0-$i ";
    ssh compute-0-$i 'ls -l /state/partition1/tmpGEM*/*_aligned_reads.gem.sam 2>/dev/null';
  done;
#results present on compute-0-[5679], rsync them to head node
for i in 3 4 5 7 9;
  do rsync -aP compute-0-$i:"/state/partition1/tmpGEM*/5[012345]_markdups.gem.bam*" /share/space/reevesp/patellifolia/map;
  done;
##then go in and delete the tmp directory contents### 
for i in 5 6 7 9;
  do ssh compute-0-$i 'rm -r /state/partition1/tmpGEM*';
  done;







#Experiment 1/17/20, try an array of different settings meant to force end to end alignment


mypp() { 
       xx=$1; #load in the setting to play with
       i=51;
       f="$wd"/"$i"output_forward_paired.fq.gz; #forward reads
       r="$wd"/"$i"output_reverse_paired.fq.gz; #reverse reads
       tmpd="/state/partition1/tmpGEM$i";
       echo "$i $xx" > "$tmpd"/log.txt; #write the setting played with
       thr=$(lscpu | grep "^CPU(s):" | awk '{print $2}'); #max threads
       if [ ! -d "$tmpd" ]; then mkdir "$tmpd"; fi; #make a local tmp directory on the main disk
       readgroup="@RG\tID:$i\tSM:$i\tPL:illumina\tLB:na\tPU:na";

       #do GEM alignment, save result to local directory
       #include variable $xx as a option
       /share/apps/gem-mapper --index="$v" \
           --i1="$f" \
           --i2="$r" \
           --paired-end-alignment \
           --sam-read-group-header="$readgroup" \
           --gzip-output "$xx" \
           --output="$tmpd"/$i"_aligned_reads.gem.sam";

       #sam to bam, then index
       /home/reevesp/bin/samtools sort -O BAM -T "$tmpd"/tmp.$i.bam --threads "$thr" -o "$tmpd"/$i"_aligned_reads.gem.bam" "$tmpd"/$i"_aligned_reads.gem.sam";

       #markdups
       #x=$(ls /state/partition1/tmp/*_aligned_reads.gembam | rev | cut -d_ -f3 | cut -c1-2 | rev); #get the sample number of the aligned reads file on this node (in case nodes have changed due to failure and repeat)
       /home/reevesp/bin/sambamba markdup -t "$thr" --tmpdir="$tmpd" --overflow-list-size 6000000 --remove-duplicates "$tmpd"/$i"_aligned_reads.gem.bam" "$tmpd"/$i"_markdups.gem.bam";
       
       #mpileup
       /home/reevesp/bin/samtools mpileup "$tmpd"/$i"_markdups.gem.bam" > "$tmpd"/$i".mpl";
}
export -f mypp;
cd /share/space/reevesp/patellifolia/data;
#v="/share/space/reevesp/patellifolia/ref/PpanfinalGEMref.gem"; export v;
v="/share/space/reevesp/patellifolia/ref/PatellifoliaSampleGEMref.gem"; export v;
wd=$(pwd); export wd;
echo "--max-reported-matches=1
--min-template-length=150
--alignment-global-min-identity=90
--alignment-local=never
--discordant-pair-search=never
" | parallel --sshloginfile ~/machinesgem --env v --env wd --env mypp mypp;























/home/reevesp/bin/samtools sort -O BAM -T "$tmpd"/tmp.$i.bam --threads "$thr" -o "$tmpd"/$i"_aligned_reads.gem.bam" "$tmpd"/$i"_aligned_reads.gem.sam";
/home/reevesp/bin/samtools sort -O BAM -T "$tmpd"/tmp.$i.bam -o "$tmpd"/$i"_aligned_reads.bam"



echo 51$'\n'55 | parallel "/home/reevesp/bin/samtools sort -O BAM -T $tmpd/tmp.{}.bam --threads $thr -o $tmpd/{}_aligned_reads.gem.bam $tmpd/{}_aligned_reads.gem.sam";
echo 51$'\n'55 | parallel --dry-run "/home/reevesp/bin/sambamba markdup -t $thr --tmpdir=$tmpd --overflow-list-size 6000000 --remove-duplicates $tmpd/{}_aligned_reads.gem.bam $tmpd/{}_markdups.gem.bam";

echo 51 | parallel "/home/reevesp/bin/samtools sort -O BAM -T $tmpd/tmp.{}.bam --threads $thr -o $tmpd/{}_aligned_reads.gem.bam $tmpd/{}_aligned_reads.gem.sam";
echo 51 | parallel "/home/reevesp/bin/sambamba markdup -t $thr --tmpdir=$tmpd --overflow-list-size 6000000 --remove-duplicates $tmpd/{}_aligned_reads.gem.bam $tmpd/{}_markdups.gem.bam";






#Combine sorted, deduped bam files, index, then calculate depth:
cd /share/space/reevesp/patellifolia/map;
time(sambamba merge -t40 -p AllP.merged.bam 5[012345]_markdups.bam;)
#time(sambamba index -t40 -p AllP.merged.bam AllP.merged.bai;)
time(sambamba depth base -t40 AllP.merged.bam -o AllP.merged.coverage;) #~ 4hrs




./bin/gem-mapper -I hsapiens_v37.gem -1 sample.Illumina.pe.1.fastq -2 sample.Illumina.pe.2.fastq -o sample.Illumina.pe.sam

####DEVELOPMENT####
#figure out some test cases:
#get list of sam files with both read pairs
zz=$((for i in $(find . -name "*.2x.sam");
  do j=$(grep 'RG:Z' "$i" | wc -l);
    echo "$j $i";
  done;) | awk -F' ' '$1==2{print $0}' | cut -d' ' -f2;)
#print regions of alignment with no overlap to reference, if >=2 such regions, there is likely non-overlap btw read pairs.
xx=$(for i in $zz; 
  do echo "$i" >> xxx.tmp;
    mp=$(samtools mpileup "$i" 2>/dev/null);
    csq=$(conseq "$(echo "$mp" | cut -d$'\t' -f2"); #identify regions where reads are not aligned consecutively to the reference, append an 'N' to indicate that the range should be filled with Ns,
    if [[ $(echo "$csq" | wc -l) > 1 ]];
    then echo "$i";
      #echo "$csq";
    fi;
  done;) 
#identify samples with insertions or deletions indicated in the mpileup from the set zz (complete read pair, non overlapping)
for i in $xx;
  do mp=$(samtools mpileup "$i" 2>/dev/null);
    if [[ $(echo "$mp" | awk -F$'\t' '$5~"+"{print $0}' | wc -l) > 0 ]];
    #if [[ $(echo "$mp" | awk -F$'\t' '$5~"-"{print $0}' | wc -l) > 0 ]];
    then echo "$i";
    fi;
  done;

#The alignment ./jcf7180008527965_177_E00558_144_HHGCMCCXY_1_2207_28513_7989_RG_Z_51.2x.sam has non-overlapping read pairs and an insertion relative to reference.
#There are no alignments with two non-overlapping reads (zz) that have deletions.
#Below alignments have overlapping read pairs and deletions
#./jcf7180008527965_177_E00558_144_HHGCMCCXY_4_1206_16853_42429_RG_Z_52.2x.sam
#./jcf7180008395354_294_E00558_144_HHGCMCCXY_4_1115_25286_46244_RG_Z_54.2x.sam
#./jcf7180008527965_177_A00197_23_H5TKNDMXX_1_2131_24957_4789_RG_Z_52.2x.sam
#./jcf7180008527965_177_E00526_165_HHWL7CCXY_2_2220_22191_68271_RG_Z_52.2x.sam
#./jcf7180008395354_294_E00558_144_HHGCMCCXY_2_2224_12966_8657_RG_Z_55.2x.sam


#develop code to address every nucleotide coding possibility in the mpileup.  This might vary by dataset so it would be good to check this each time.
#all positions with -3NNN (or similar, - generally) will also have the following three positions marked as unknown using * or **, so it does not appear to be necessary to deal with '-' other than to strip it plus the appended NNNs from the line


#the problem arises with insertions vis a vis original myiupac subroutine.  here are the classes of insertion possibilities (not sure if this is exhaustive, but should be):
str="T+3CAAC+3CAA";# TC CC AA AA (works)
str="T+3CAAC+2CA";# TC CC AA A* (works)
str="T+2CAC+3CAA";# TC CC AA *A (works)
str="T+3CAA";# T* C* A* A* (works)
str="T+3CAAT";# TT C* A* A* (works)
str="TT+3CAA";# TT *C *A *A (works)





m=$(echo "$str" | tr '+' '\n'); #break into lines on +
l1len=$(echo "$m" | head -1 | awk '{print length}'); #length of line one
nl=$(echo "$m" | wc -l); #number of lines
l2=$(echo "$m" | head -2 | tail -1); #get line 2 input
l2inslen=$(echo "$l2" | tr -d -c 0-9); #deletes all non-numeric characters from line 2 to read the length of the insertion
l2nchar=$(echo "$l2" | tr -d 0-9 | awk '{ print length }'); #number of chars in line 2, counts all nun numeric characters

if [[ $l1len == 1 ]]; #case where there in an insertion in the first read there is one character on line 1. here, you need to get character from first line plus last character of the second line, if present
then if [[ $nl == 2 ]];
     then if [[ $l2nchar == $l2inslen ]]; #test whether there is a base from the second read hanging off the end of the insertion in line 2
          then m=$(echo "$m" | sed '2s/$/\*/'); #add an asterisk to the end of the second line to act as the unknown read 2 base
          fi;
       l2=$(echo "$m" | head -2 | tail -1); #recalculate line 2 input now that it has been modified with an asterisk
       n=$(paste -d' ' <(echo "$m" | head -1) <(echo "$l2" | rev | cut -c1) | sed 's/ //'); #char states composed of line 1 + last character of line 2
       for i in $(seq 1 1 $l2inslen);
         do n=$(echo "$n";echo "$l2" | tr -d 0-9 | cut -c$i)"*"; #combine character n of line 2 with an asterisk to form the read pair
         done;
      
     elif [[ $nl == 3 ]]; #if there are three lines then there are inserts after both reads, which means that the last character of line 2 is always read 2 character 1
     then n=$(paste -d' ' <(echo "$m" | head -1) <(echo "$l2" | rev | cut -c1) | sed 's/ //'); #char states composed of line 1 + last character of line 2
       l3=$(echo "$m" | head -3 | tail -1); #get line 3 input
       l3inslen=$(echo "$l3" | tr -d -c 0-9); #deletes all non-numeric characters from line 3 to read the length of the insertion

       if [[ $l2inslen == $l3inslen ]]; #case where insertion lengths are the same
       then for i in $(seq 1 1 $l2inslen);
         do n=$(echo "$n";(echo "$l2" | tr -d 0-9 | cut -c$i;echo "$l3" | tr -d 0-9 | cut -c$i) | tr -d '\n'); #combine character n of line 2 with character n of line 3 to form the read pair
         done;
       elif [[ $l2inslen < $l3inslen ]]; #case where read1 (l2) insertion length is less than read 2 (l3)
       then for i in $(seq 1 1 $l2inslen);
         do n=$(echo "$n";(echo "$l2" | tr -d 0-9 | cut -c$i;echo "$l3" | tr -d 0-9 | cut -c$i) | tr -d '\n'); #combine character n of line 2 with character n of line 3 to form the read pair
         done;
         #now combine asterisk for line 2 (which is shorter) with character n of line 3 to form the read pair
         for i in $(seq $(($l2inslen+1)) 1 $l3inslen);
           do n=$(echo "$n";echo -n "*";echo "$l3" | tr -d 0-9 | cut -c$i); 
           done;
       elif [[ $l2inslen > $l3inslen ]]; #case where read1 (l2) insertion length is greater than read 2 (l3)
       then for i in $(seq 1 1 $l3inslen);
         do n=$(echo "$n";(echo "$l2" | tr -d 0-9 | cut -c$i;echo "$l3" | tr -d 0-9 | cut -c$i) | tr -d '\n'); #combine character n of line 2 with character n of line 3 to form the read pair
         done;
         #now combine character n of line 2 with asterisk for line 3 (which is shorter) to form the read pair
         for i in $(seq $(($l3inslen+1)) 1 $l2inslen);
           do n=$(echo "$n";echo "$l2" | tr -d 0-9 | cut -c$i)"*"; 
           done;
       fi;    
     fi;
  
elif [[ $l1len == 2 ]];
  then n=$(echo "$m" | head -1);
    for i in $(seq 1 1 $l2inslen);
      do n=$(echo "$n";echo -n "*";echo "$l2" | tr -d 0-9 | cut -c$i); #place an asterisk for read 1 in combination with character n of line 2 for read 2
      done;
else echo "ERROR: incompatible string: $m";
fi;

















e=$(while read j;
do d=$(echo "$j" | cut -d$'\t' -f1,5 | tr "\t" ":");
  echo "$d";
done <<<"$c";)





mp=$(samtools mpileup jcf7180008587925_106_E00558_144_HHGCMCCXY_5_1214_8948_25605_RG_Z_55.bam);
mp=$(samtools mpileup jcf7180008587925_106_A00197_23_H5TKNDMXX_1_1323_12337_8390_RG_Z_53.bam);
if [[ "$mp" != "" ]];
  then echo not equal;


  g echo equal;
  
  
  fi;




samtools mpileup -uf ref.fa aln.bam | bcftools view -cg - | vcfutils.pl vcf2fq > cns.fq

bcftools mpileup --fasta-ref "$h"_ref.txt "$i".bam | bcftools view -cg - | vcfutils.pl vcf2fq > cns.fq


samtools mpileup --reference "$h"_ref.txt "$i".bam | bcftools view -cg - | vcfutils.pl vcf2fq > cns.fq



> xxx.tmp;
for i in $(ls jcf7180008587925_106_*.bam);
  do echo "// $i //" >> xxx.tmp
    #samtools mpileup --reference jcf7180008587925_ref.txt "$i" >> xxx.tmp;
    samtools mpileup --reference jcf7180008587925_ref.txt "$i" >> xxx.tmp;
  done;

samtools mpileup --reference jcf7180008587925_ref.txt jcf7180008587925_106_A00197_23_H5TKNDMXX_1_1323_12337_8390_RG_Z_53.bam; #two sequence with overlap example
samtools mpileup jcf7180008587925_106_A00197_23_H5TKNDMXX_1_1323_12337_8390_RG_Z_53.bam; #two sequence no reference
samtools mpileup --reference jcf7180008587925_ref.txt jcf7180008587925_106_E00558_144_HHGCMCCXY_5_1214_8948_25605_RG_Z_55.bam; #two sequence with no overlap at ref base 150-164
samtools mpileup jcf7180008587925_106_E00558_144_HHGCMCCXY_5_1214_8948_25605_RG_Z_55.bam; #two sequence with no overlap at ref base 150-164














cd /scratch/reevesp/patellifolia;
time(bwa mem -t 40 ref/Ppanfinal.genome.scf.fasta data/50output_forward_paired.fq.gz data/50output_reverse_paired.fq.gz data/51output_forward_paired.fq.gz data/51output_reverse_paired.fq.gz \ 
  data/52output_forward_paired.fq.gz data/52output_reverse_paired.fq.gz \
  data/53output_forward_paired.fq.gz data/53output_reverse_paired.fq.gz \
  data/54output_forward_paired.fq.gz data/54output_reverse_paired.fq.gz \
  data/55output_forward_paired.fq.gz data/55output_reverse_paired.fq.gz \
  -o data/Ppanpe.sam;) #~xhrs





cd /share/space/reevesp/patellifolia/;
time(bwa mem -t 24 ref/Pprofinal.genome.scf.fasta data/17134D-01-55-P1_S19_L001_R1_001.fastq.gz data/17134D-01-55-P1_S19_L001_R2_001.fastq.gz -o data/55pe.sam;) #~9 hrs
^^^

samtools view -q 20 -b map/"$i"se.sam | samtools sort -O BAM - -o map/"$i"filPprose.bam; #~xhrs
sambamba markdup -t40 --tmpdir=/scratch --remove-duplicates map/"$i"filPprose.bam map/"$i"filPproseDedup.bam; #x minutes


#use bazam to extract read pairs at site
time(java -jar /home/reevesp/bin/bazam/build/libs/bazam.jar -bam map/53_markdups.bam --regions jcf7180008368139:50-51;)

