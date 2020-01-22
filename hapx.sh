#!/bin/bash

# Usage: see README.txt

##SUBROUTINES##

#conseq locates gaps in sequence
conseq() {
          echo "$1" | awk '$1!=p+1{print p+1"-"$1-1}{p=$1}';
}
export -f conseq;

#myiupac() replaces character strings with their iupac equivalent or proper insertion and deletion coding.
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

#myinsertion() processes insertions coded in mpileup
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

##END SUBROUTINES##


#acquire command line variables
ref=$1; #path to the reference genome sequence in multi-fasta format
aln=$2; #path to the bam file containing the alignment of reads to ref
sites=$3; #genomic regions to use
e=$(echo "$sites" | tr "," "\n"); #format sites for proper parsing by samtools view

#extract read pairs at site using samtools. Extract proper pairs with mapping quality > 60 excluding unmapped pairs
echo "$e" | parallel 'samtools view -f 2 -F 4 -q 60 '"$aln"' "{}" | sort > {}.tmp; \
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
echo "Number of reads found at sites:"$'\n'"contig:sites pairedreads singlereads"$'\n'"$f"$'\n';

#           sample output:
#           Number of reads found at sites
#           contig:sites pairedreads singlereads
#           jcf7180008454378:303-304 0 8
#           jcf7180008531951:103-495 12 94


g=$(echo "$f" | awk -F' ' '$2+$3!=0 {print $0}'); #remove any contigs from consideration when there are 0 reads (no paired reads:$2, no unpaired reads:$3) that align to them (this could also be a variable that supports a cutoff)

#extract each contig from the reference genome to be used by itself as a reference for realignment with the extracted read pairs, then bwa index
echo "$g" | cut -d: -f1 | parallel 'sed -n -e "/^>{}$/,/>/p" '"$ref"' | sed "$ d" > {}_ref.txt; bwa index {}_ref.txt'; #sed extracts lines from name of contig of interest until next contig name line, second sed deletes the last line

#reconstruct fastq file for reads that aligned to each contig, adding readgroup and f/r designation, r is for reads marked with a negative alignment length in column 9
#then sort reads into individual read pairs, resulting files are like {contigname_readname.rp.fq}
echo "Reconstructing fastq files for reads found at sites:";
for i in $(echo "$g" | cut -d' ' -f1 | tr ':' '_' | tr '\n' ' ');
  do awk -F$'\t' '{print "@"$1,$9,$10,"+",$11,$17}' "$i".tmp | awk -F' ' '{if ($2 > 0) {print $1"_"$6":f",$3,$4,$5} else if ($2 < 0) {print $1"_"$6":r",$3,$4,$5}}' | tr " " "\n" > "$i"_reads.fq;

    m=$(grep ^'@' "$i"_reads.fq | cut -d: -f1-9 | sort -u); #capture list of unique read pair names (includes also unpaired reads with unique names) 
    for j in $m;
      do echo "$i $j";
        fn=$(echo "$j" | sed 's/@//g' | tr ':' '_'); #filename for output
        grep -A3 "$j" "$i"_reads.fq > "$i"_"$fn".rp.fq; #save read pair into its own reconstructed fastq file
      done;
    
    rm "$i"_reads.fq; #clean up
    rm "$i".tmp;
  done;

#           Here is a sample reconstructed *.rp.fq file, note readgroup in col9 of name, and "orientation" or "endedness of read pair" in col10:
#          @E00558:144:HHGCMCCXY:5:2224:6654:40055_RG:Z:55:r
#          CAAAAATATTTTGGTAATTATTCTCAACAAAATGATTTGAAAGGTGTTCATACAACACAAATCGCCTAAGAGACTATGACGGTTTTATCCTCTGATTTGAATTGAGTTTGATCCAAGGGCTTCATATGATTGAAATATAC
#          +
#          FFJJJJAFJ<JJFF7FJAJAA-A7AJJJJJJ<JJJJJJJJJJFJJJJJJJJJJJ<JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFFFAA
#          @E00558:144:HHGCMCCXY:5:2224:6654:40055_RG:Z:55:f
#          ATAAATATGCAAGGAGTCAAAAATATTTGGTAATTAAGCACAAAAACTGATTTGAAATGTGTTCGTACACCACAAATCACCTAAGAGACTATGACGGTTTTACCCTTTGATTTGAATTGAGTTTGATCCAAGGGCTTCATATGCTTTAAA
#          +
#          AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJ

#align *.rp.fq to *_ref.txt, output as bam file, then sambamba index (output files have names like jcf7180008378511_277.bam (contigname_site.bam)
#this step aligns, then verifies mapping quality since it changes from the original values, then realigns
#This step should probably be parallelized
echo "Realigning reads:"
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


#turn on some diagnostic code by un-commenting:
#print mapping quality after 1st realignment (it changes from original mapping quality)
#awk -F$'\t' '{print $5}' *.sam | grep -v ^CL | grep -v ^$ | sort | uniq -c;
#print mapping quality after RE-realignment (all 60)
#awk -F$'\t' '{print $5}' *.2x.sam | grep -v ^CL | grep -v ^$ | sort | uniq -c;

#figure out the possibilities for coding in the nucleotide column ($5)
#> xxx.tmp;
#for i in $(find . -name "jcf*_RG_Z_5[012345].2x.sam");
#  do echo "$i" >> xxx.tmp;
#    samtools mpileup "$i" >> xxx.tmp;
#  done;
#awk -F$'\t' '{print $4,$5}' xxx.tmp | sort | uniq -c > ntcodes.tmp; #print all unique possibilities for nucleotide coding, these have to be accounted for when assembling the consensus





#Process alignments of paired reads via mpileup into a single padded consensus sequence representing the haplotype
if [ ! -d "haplotypes" ]; then mkdir "haplotypes"; fi; #make a directory to hold fasta sequences containing extracted haplotypes

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
    echo "$base3" > haplotypes/"$rname.fa";
    
    #clean up
    rm "$samin";
  done;





#Produce an alignment of the haplotypes to their reference, which might be viewed in something like IGV (Integrative Genomics Viewer)
if [ ! -d "alignments" ]; then mkdir "alignments"; fi; #make a directory to hold fasta sequences containing extracted haplotypes
ams=$(find ./haplotypes -name "*.2x.fa" | cut -d_ -f1-2 | sort -u | rev | cut -d'/' -f1 | rev);
echo "$ams" | parallel 'cat ./haplotypes/{}_*.2x.fa > ./alignments/{}.mfa'; #concatenate all fasta haplotypes belonging to a single contig:range

myalignhaps() {
              i=$1;
              thr=$(lscpu | grep "^CPU(s):" | awk '{print $2}'); #max threads
              rr=$(echo "$i" | cut -d'_' -f1 ); #reference contig name, e.g. jcf7180008454378
              ss=$(echo "$i" | cut -d'.' -f1 ); #refcontig+siterange, e.g. jcf7180008454378_303-304
              #bwa index "$pd"/"$rr"_ref.txt;
              /home/reevesp/bin/bwa mem -t "$thr" "$pd"/"$rr"_ref.txt "$pd"/alignments/"$i"| /home/reevesp/bin/samtools sort -O BAM --threads "$thr" -o "$pd"/alignments/"$ss"_aligned_haps.bam;
              sambamba index "$pd"/alignments/"$ss"_aligned_haps.bam "$pd"/alignments/"$ss"_aligned_haps.bam.bai;
}
export -f myalignhaps;

pd=$(pwd); export pd;
find $(pwd)"/alignments" -name "*.mfa" | rev | cut -d'/' -f1 | rev | parallel --env pd myalignhaps;
#start here, do bwa-mem in parallel with *_ref.txt


#clean up
rm *_ref.txt*;

  
