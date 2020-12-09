#!/bin/bash

# Usage: see README.txt

##SUBROUTINES##

#myconseq locates gaps in numeric sequence
myconseq() {
          echo "$1" | awk '$1!=p+1{print p+1"-"$1-1}{p=$1}';
}
export -f myconseq;

#myiupac() replaces character strings with their iupac equivalent or proper insertion and deletion coding.
myiupac() {
        read i; #input is a piped row of data
        a=$(tr " " "\n" <<<"$i" | cut -d$'\t' -f1); #get the position numbering
        b=$(tr " " "\n" <<<"$i" | cut -d$'\t' -f2 \
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
#echo replaced with here-strings "<<<" for possible speed increase            
              read i;
              i=$(tr ' ' '\n' <<<"$i"); #delimit with line breaks
              j=$(grep "+" <<<"$i"); #collect all lines with insertions
              ij=$(grep -v "+" <<<"$i"); #collect all lines without insertions

              #process lines with insertions ("$j"), save in a variable $lwi
              if [[ "$j" == "" ]];
              then lwi=""; #there are no insertions so create a dummy, empty variable
              else
                lwi=$(while read str;
                  do pos=$(cut -d$'\t' -f1 <<<"$str"); #get the starting position of the insertion
                    m=$(cut -d$'\t' -f2 <<<"$str" | tr '+' '\n'); #break into lines on +
                    l1len=$(head -1 <<<"$m" | awk '{print length}'); #length of line one
                    nl=$(wc -l <<<"$m"); #number of lines
                    l2=$(head -2 <<<"$m" | tail -1); #get line 2 input
                    l2inslen=$(tr -d -c 0-9 <<<"$l2"); #deletes all non-numeric characters from line 2 to read the length of the insertion
                    l2nchar=$(tr -d 0-9 <<<"$l2" | awk '{ print length }'); #number of chars in line 2, counts all non numeric characters
                    
                    if [[ $l1len == 1 ]]; #case where there in an insertion in the first read there is one character on line 1. here, you need to get character from first line plus last character of the second line, if present
                    then if [[ $nl == 2 ]];
                         then if [[ $l2nchar == $l2inslen ]]; #test whether there is a base from the second read hanging off the end of the insertion in line 2
                              then m=$(sed '2s/$/\*/' <<<"$m"); #add an asterisk to the end of the second line to act as the unknown read 2 base
                              fi;
                           l2=$(head -2 <<<"$m" | tail -1); #recalculate line 2 input now that it has been modified with an asterisk
                           n=$(paste -d' ' <(echo "$m" | head -1) <(echo "$l2" | rev | cut -c1) | sed 's/ //'); #char states composed of line 1 + last character of line 2
                           for i in $(seq 1 1 $l2inslen);
                             do n=$(echo "$n";tr -d 0-9 <<<"$l2" | cut -c$i)"*"; #combine character n of line 2 with an asterisk to form the read pair
                             done;
                          
                         elif [[ $nl == 3 ]]; #if there are three lines then there are inserts after both reads, which means that the last character of line 2 is always read 2 character 1
                         then n=$(paste -d' ' <(head -1 <<<"$m") <(rev <<<"$l2" | cut -c1) | sed 's/ //'); #char states composed of line 1 + last character of line 2
                           l3=$(head -3 <<<"$m" | tail -1); #get line 3 input
                           l3inslen=$(tr -d -c 0-9 <<<"$l3"); #deletes all non-numeric characters from line 3 to read the length of the insertion
                    
                           if [[ $l2inslen == $l3inslen ]]; #case where insertion lengths are the same
                           then for i in $(seq 1 1 $l2inslen);
                             do n=$(echo "$n";(tr -d 0-9 <<<"$l2" | cut -c$i;tr -d 0-9 <<<"$l3" | cut -c$i) | tr -d '\n'); #combine character n of line 2 with character n of line 3 to form the read pair
                             done;
                           elif [[ $l2inslen < $l3inslen ]]; #case where read1 (l2) insertion length is less than read 2 (l3)
                           then for i in $(seq 1 1 $l2inslen);
                             do n=$(echo "$n";(tr -d 0-9 <<<"$l2" | cut -c$i;tr -d 0-9 <<<"$l3" | cut -c$i) | tr -d '\n'); #combine character n of line 2 with character n of line 3 to form the read pair
                             done;
                             #now combine asterisk for line 2 (which is shorter) with character n of line 3 to form the read pair
                             for i in $(seq $(($l2inslen+1)) 1 $l3inslen);
                               do n=$(echo "$n";echo -n "*";tr -d 0-9 <<<"$l3" | cut -c$i); 
                               done;
                           elif [[ $l2inslen > $l3inslen ]]; #case where read1 (l2) insertion length is greater than read 2 (l3)
                           then for i in $(seq 1 1 $l3inslen);
                             do n=$(echo "$n";(tr -d 0-9 <<<"$l2" | cut -c$i;tr -d 0-9 <<<"$l3" | cut -c$i) | tr -d '\n'); #combine character n of line 2 with character n of line 3 to form the read pair
                             done;
                             #now combine character n of line 2 with asterisk for line 3 (which is shorter) to form the read pair
                             for i in $(seq $(($l3inslen+1)) 1 $l2inslen);
                               do n=$(echo "$n";tr -d 0-9 <<<"$l2" | cut -c$i)"*"; 
                               done;
                           fi;    
                         fi;
                      
                    elif [[ $l1len == 2 ]];
                      then n=$(echo "$m" | head -1);
                        for i in $(seq 1 1 $l2inslen);
                          do n=$(echo "$n";echo -n "*";tr -d 0-9 <<<"$l2" | cut -c$i); #place an asterisk for read 1 in combination with character n of line 2 for read 2
                          done;
                    else echo "ERROR: incompatible string: $m";
                    fi;
                    
                    nls=$(( $(wc -l <<<"$n") - 1 )); #count the number of lines in the processed insertion, decrement by one since first line will not get a decimal place
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

#mycountqualreadpairs counts how many paired reads and unpaired reads map to the target sites, receives info in $tmpf and $site from mycon1()
mycountqualreadpairs() {
                       p=$(cut -d$'\t' -f1 <<<"$tmpf" | sort | uniq -c | grep ^" \+2 "); #paired reads
                       s=$(cut -d$'\t' -f1 <<<"$tmpf" | sort | uniq -c | grep ^" \+1 "); #single ends of paired reads
                       pp=$(awk 'NF' <<<"$p" | wc -l); #number of paired reads (awk 'NF' removes empty lines)
                       ss=$(awk 'NF' <<<"$s" | wc -l); #number of single ends of paired reads
                       echo "$site $pp $ss";
}
export -f mycountqualreadpairs;

mygetends() {
              i=$1; #contig:site-range
              lesr=$(cut -d: -f2 <<<"$i" | cut -d'-' -f1); #left end site range
              resr=$(cut -d: -f2 <<<"$i" | cut -d'-' -f2); #right end site range
              
              #f=$(samtools view -q 1 Hs1pro1l1.finalaln.bam 51jcf7180007742276:"$i"-"$i" | cut -d$'\t' -f4 | sort -nr | head -1);
              
              #get LE and cigar string for each sequence
              #g=$(/share/apps/samtools view -F 2048 -q "$stq" "$bam" "$i" | cut -d$'\t' -f4,6); #-F 2048 excludes supplementary alignments
              g=$(samtools view -F 2048 -q "$stq" "$bam" "$i" | cut -d$'\t' -f4,6); #-F 2048 excludes supplementary alignments
              if [[ "$g" == "" ]]; then return; fi; #bail out if there are no aligned reads at the position
              
              #calculate closest ends to left side of site range
              le=$(cut -d$'\t' -f1 <<<"$g" | sort -nr | head -1); #position of left read pair end nearest to left end of site range
              led=$(( $lesr - $le )); #distance to closest LE
              
              #parse cigar string for right end positions
              #allres=$(echo "$g" | tr '\t' ':' | parallel --env myparsecigar myparsecigar); #variable to hold all right end positions
              allres="";
              for j in $(echo "$g" | tr '\t' ':' | tr '\n' ' ');
                do lpos=$(cut -d: -f1 <<<"$j"); #position of first aligned base
                  cig=$(cut -d: -f2 <<<"$j"); #cigar string
                  cigops=$(sed 's/[0-9]*//g' <<<"$cig" | grep -o . | sort -u | tr -d '\n'); #gather all cigar string 'operators'
                  
                  #process the cigar string. 1) split string on operations DIMSX=, 2) remove trailing blank line
                  #3) remove any line with I (insertion) op, you are looking for the right end position relative to
                  #the reference so insertions in the read are not counted, 4) remove DIMSX= characters, leaving the number
                  #of base pairs for each op, 5) count all lines, each containing the number of base pairs per op,
                  #6) subtract 1 because the cigar ops start on the lpos so that rpos is the last actual base of the 
                  #read.
                  if [[ "$cigops" == *"P"* ]];
                  then echo "CIGAR string $cig contains invalid character P, skipping";
                    return; #skip the current contig:site-range
                  else lengthcig=$(sed 's:\([DIMSX=]\):\1\n:g' <<<"$cig" | sed '/^$/d' | grep -v "I" | sed 's/[DIMSX=]//' | awk '{s+=$1}END{print s}');
                    rpos=$(( $lpos + $lengthcig - 1 )); #position of read pair right end
                    allres+="$rpos"$'\n';
                   fi;
                done;
              allres=$(sed '/^$/d' <<<"$allres"); #remove trailing blank line
              
              #calculate closest ends to right side of site range
              #distance to all REs from RE of site range
              
              re=$(sort -n <<<"$allres" | head -1); #position of closest RE
              red=$(( $re - $resr )); #distance to closest RE
              
              #report microhaplotype range to calling statement
              echo "$lesr-$resr $le-$re "$(( $re - $le + 1 ));
}
export -f mygetends;

#mycon1 combines old functions into one so that data can be passed in memory instead of using the disk
mycon1() {
       site="$1"; #incoming data is a description of sites to process in contigname:site-range format like jcf7180008531951:276-301
       contigname=$(cut -d: -f1 <<< "$site");

       #extract read pairs at target sites using samtools. Obey include, exclude and quality rules from command line
       #tmpf1=$(/share/apps/samtools view -f "$stf" -F "$stF" -q "$stq" "$bam" "$site" | sort); #get read pairs mapped to contig:site-range in original bwa or gem alignment
       tmpf1=$(samtools view -f "$stf" -F "$stF" -q "$stq" "$bam" "$site" | sort); #get read pairs mapped to contig:site-range in original bwa or gem alignment
       if [[ "$debug" == "YES" ]]; then echo "$tmpf1" > "$pd"/"$site".tmpf1; fi;
       
       #reads extracted above must contain -both- pairs within the range $site. This is undesirable, for example
       #if the site range is 1bp, then reads from pair must overlap in order to recover both.
       #With stf=1, you can recover a single read that exists as a mapped pair, but where both reads don't map to the site-range.
       #From the current contig, include pair of any read that -exists as a mapped pair- (stf=1) and that falls within
       #the site range.  These pairs may be eliminated anyway in later steps that realign and filter by quality,
       #but they should not be excluded by virtue of non-overlap with the targeted site range
       
       #htmpf=$(grep ^'@' <<< "$tmpf1"); #sam header from $tmpf1
       unptmpf=$(cut -d$'\t' -f1 <<< "$tmpf1" | sort -u); #acquire list of all reads that map to the site-range, includes paired and unpaired
       #btmpf=$(/share/apps/samtools view "$bam" "$contigname"); #get all reads associated with contig
       btmpf=$(samtools view "$bam" "$contigname"); #get all reads associated with contig
       tmpf=$(LC_ALL=C grep -w -F -f <(echo "$unptmpf") <<< "$btmpf"); #get read pairs, any one of which mapped to site-range

       #bwa mem, used later, does not use Illumina quality scores for mapping (per Heng Li:https://sourceforge.net/p/bio-bwa/mailman/message/34410817/)
       #bwa mem will not align a fastq sequence if the quality scores are missing (col11 = "*"), so
       #spoof them here, substituting the * in column 11 with column 10 sequence, then converting sequence
       #in column 11 to all JJJs (Q=41), just to get bwa mem to work
       tmpf=$(awk -F$'\t' '{OFS="\t"};$11=="*"{$11=$10;gsub(/./,"J",$11)};{print $0}' <<<"$tmpf"); #list of lengths of sequences 2 lines before asterisks
       
       if [[ "$debug" == "YES" ]]; then echo "$tmpf" > "$pd"/"$site".tmpf; fi;

       #count number of read pairs and single ends that have met the samtools -f/-F/-q rules
       #there will only be single ends if the pair doesn't map to the contig
       if [[ "$tmpf" == "" ]]; then return;
       else f=$(mycountqualreadpairs);
       fi;
       
       g=$(awk -F' ' '$2+$3!=0 {print $0}' <<< "$f"); #remove contig from consideration if there are 0 reads (no paired reads:$2, no unpaired reads:$3) that align (this could also be a variable that supports a cutoff)
       if [[ "$g" == "" ]]; then return;
       fi;



       #myreconfq reconstructs a single fastq file that contains exclusively both reads of a read pair from a sam file
       i="$site"; #input is a description of sites to process in contigname:range format like jcf7180008531951:276-301
               #the sites are accessed from the variable $tmpf below
       h=$(cut -d':' -f1 <<<"$i"); #just the contig name
 
       s=$(awk -F$'\t' -v rgf=$rgf '{print "@"$1"_"$rgf,$10,"+",$11}' <(echo "$tmpf") | tr " " "\n" | sed '/^$/d');
       
       m=$(grep ^'@' <(echo "$s") | cut -d: -f1-9 | sort -u | sed '/^$/d'); #capture list of unique read pair names (includes also unpaired reads with unique names) that map to contig:site $i 

       
       
       
       
       mfa=""; #initialize variable to contain all the read pair haplotypes
       #below j contains a readgroup-annotated read pair name like @E00558:144:HHGCMCCXY:4:2210:23074:34043_RG:Z:55, that aligns to the current contig:site range
       for j in $m;
         do rpfq=$(grep -A3 "$j" <(echo "$s") | sed '0,/\(^@\S\+$\)/ s/\(^@\S\+$\)/\1 1:N:0:AAAAAA/' | sed 's/\(^@\S\+$\)/\1 2:N:0:AAAAAA/'); #save read pair into a variable containing a reconstructed fastq, labeled as a read pair for bwa
           #above, the first sed searches from 0 to the pattern then performs substitution 1 on the fastq defline, which causes the insertion of whitespace.  The second sed then searches the whole file and performs substitution 2 on the second fastq defline. \S means 'not whitespace'
           #above, create single files containing both reads of a read pair
           
           #myrealign takes the reads extracted from the primary alignment ($rpfq) and realigns them independently to the contig as reference
           #mapping quality changes upon realignment, extract only reads with user defined properties (/share/apps/samtools view -fFq), convert those to proper fastq file and re-align again
           #-p, paired-end mode assumes read pairs are consecutive; -M mark split hits as secondary alignments (so they can be ignored in samtools view -F step)
           #t=$(/share/apps/bwa mem -p -M "$pd"/"$h"_ref.txt <(echo "$rpfq") 2>/dev/null | \
           #    /share/apps/samtools sort -O SAM 2>/dev/null | \
           #    /share/apps/samtools view -f "$stf" -F "$stF" -q "$stq" 2>/dev/null | \
           #    awk -F$'\t' '{print "@"$1,$10,"+",$11}' | tr " " "\n" | sed '0,/\(^@\S\+$\)/ s/\(^@\S\+$\)/\1 1:N:0:AAAAAA/' | sed 's/\(^@\S\+$\)/\1 2:N:0:AAAAAA/');
           t=$(bwa mem -p -M "$pd"/"$h"_ref.txt <(echo "$rpfq") 2>/dev/null | \
               samtools sort -O SAM 2>/dev/null | \
               samtools view -f "$stf" -F "$stF" -q "$stq" 2>/dev/null | \
               awk -F$'\t' '{print "@"$1,$10,"+",$11}' | tr " " "\n" | sed '0,/\(^@\S\+$\)/ s/\(^@\S\+$\)/\1 1:N:0:AAAAAA/' | sed 's/\(^@\S\+$\)/\1 2:N:0:AAAAAA/');

           if [[ $(echo "$t") != "" ]]; #only re-realign if there are reads left after latest quality filter
           then x2xsam=$(bwa mem -p -M "$pd"/"$h"_ref.txt <(echo "$t") 2>/dev/null | awk -F$'\t' '$4!="*"{print $0}'); #perform second realignment and manually remove any unmapped reads that may have been found as indicated by a * in the position column $4
             #x2xsam=$(/share/apps/bwa mem -p -M "$pd"/"$h"_ref.txt <(echo "$t") 2>/dev/null | awk -F$'\t' '$4!="*"{print $0}'); #perform second realignment and manually remove any unmapped reads that may have been found as indicated by a * in the position column $4

           else continue; #no reads left to align, continue to next
           fi;
 
           #test whether re-realignment still contains reads after removing unmapped reads, if not skip to next read pair
           lnsam=$(wc -l <<<"$x2xsam");
           if [[ "$lnsam" < 3 ]]; then continue; fi; 

           #mymakehapblocks() processes read pairs into a single contiguous sequence, adds NNNs where opposing read pairs do not overlap to pad relative to the reference, adds IUPAC redundancy codes for conflicts between pairs of a read, removes deletion coding, retains insertions
           rp=$(awk -F$'\t' -v h=$h '$3==h{print $1}' <<<"$x2xsam" | sort -u | tr ":" "_"); #extract the unique read name from the sam formatted data that contains a single mapped read pair
           rname=$(tr ':' '_' <<<"$i")"_$rp"; #get a root name for the read pair/contig alignment

           #mp=$(/share/apps/samtools mpileup -A <(echo "$x2xsam") 2>/dev/null); #form the pileup file, use -A to include orphan reads. samtools mpileup by default excludes improper pairs. in hapx, bwa mem will find no proper pairs because too few reads are used during alignment
           mp=$(samtools mpileup -A <(echo "$x2xsam") 2>/dev/null); #form the pileup file, use -A to include orphan reads. samtools mpileup by default excludes improper pairs. in hapx, bwa mem will find no proper pairs because too few reads are used during alignment

           #Extract the pos and base columns, remove read start ^., read end $, remove deletion indicators (e.g. -2AC, they are followed with *), pass to myinsertion subroutine to control insertions, then to myiupac to recode conflicts as ambiguous and produce the consensus sequence
           base1=$(cut -d$'\t' -f2,5 <<<"$mp"| sed 's/\^.//g' | sed 's/\$//g' | sed 's/-.*//g' | tr "\n" " " | myinsertion | myiupac); #haplotype extraction, unpadded as of now

           #Deal with contiguity and padding
           #csq=$(myconseq "$(cut -d$'\t' -f2)" <<<"$mp" | grep -v "^1\-" | sed 's/-/ 1 /g'); #identify regions where reads are not aligned consecutively to the reference, exclude the range from 1-start of overlap, set up to use as an interval for seq command
           csq=$(myconseq $(cut -d$'\t' -f2 <<<"$mp") | grep -v "^1\-" | sed 's/-/ 1 /g'); #identify regions where reads are not aligned consecutively to the reference, exclude the range from 1-start of overlap, set up to use as an interval for seq command
           
           #calculate the maximum pad size, if any pad is larger than $maxp continue to next contig:read-pair
           maxn=$(awk -F' ' '{print $3 - $1}' <<<"$csq" | sort -nr | head -1); 
           if [[ "$maxn" > "$maxp" ]];
           then continue;
           fi;
           
           #pad non-contiguous sections
           if [[ "$csq" == "" ]];
           then pads=""; #no pads if no non-contiguous sections
           else pads=$(while read isq;
                     do for k in $(seq $isq);
                       do echo "$k N";
                       done;
                     done <<<"$csq"; #possible multi lines of non-contiguous regions that need to be padded with NNNNs
                     );
           fi;
           
           #combine to pad haplotype with respect to reference, convert x for deletion to - for deletion
           #base 2 contains the processed haplotype in a pileup like format relative to reference
           base2=$(echo "$base1"$'\n'"$pads" | sort -t' ' -k1,1n | sed 's/x/-/g' | awk 'NF');
           
           #revise the read pair name so that readgroup appears in front
           nrname1=$(rev <<<"$rname" | cut -d'_' -f1-3 | rev); #readgroup info from end of string
           nrname2=$(rev <<<"$rname" | cut -d'_' -f4- | rev); #rest of string, excluding readgroup info
           nrname="$nrname1"_"$nrname2";
           
           #base3 contains the processed haplotype as a fasta file, deletions relative to reference removed.
           #base3 is as close as we can come to reconstructing the native molecule
           base3=$(echo ">$nrname";echo "$base2" | cut -d' ' -f2 | grep -v '-' | tr "\n" " " | sed 's/ //g');

           mfa+="$base3"$'\n';
           
         done; # for j in $m
         
       #remove terminal line break in major output variable $mfa
       mfa=$(sed '/^$/d' <<< "$mfa");
       if [[ "$debug" == "YES" ]]; then echo "$mfa" > "$pd"/"$site".mfa; fi;

       if [[ "$mfa" == "" ]]; then return; #skip out of this contig:site-range, there is no qualified data
       fi;




       #mydedup() counts and optionally removes duplicate sequences and subsequences that are exactly contained within longer sequences, from the processed multi fasta file
       inf=$(tr ':' '_' <<<"$i"); #make value like jcf7180008587925_40-41 from jcf7180008587925:40-41
       
       #count duplicated sequences and subsequences
       # to this process supply each read group separately, and everything together
       da=$(sed -e '/^>/s/$/@/' -e 's/^>/#/' <(echo "$mfa") | tr -d '\n' | tr "#" "\n" | tr "@" " " | sed '/^$/d' | awk '{print ">rp" NR "_" $0}'); #linearize the data set
       rgs=$(cut -d'_' -f2-4 <<<"$da"| sort -u | tr '\n' ' '); #acquire a list of readgroups
       rgs="$rgs>rp"; #add an item to the list of readgroups that will allow all sequences to be grepped from the file '>rp'
       
       allhash=""; #initialize variable to hold hashes of dna sequences for allele frequency calculation
       basicstats=""; #will contain NumHblocks:NumIdenticalHblocks:NumIdenticalHblockSubsequences:NumDistinctHblocks
       for k in $rgs;
       do 
       
         fon=$(sed 's/>rp/global/' <<<"$k"); #file output name, #replace grep item '>rp' for retrieving all sequences with "global", "global" means all unique haplotypes are counted considering all readgroups simultaneously.
         fonmfa=$(grep "$k" <<<"$da" | sort); #make a readgroup specific mfa file. this has all haploblocks, including duplicated and subsequence
         fonlmfa=$(awk -F' ' '!_[$2]++' <<<"$fonmfa"| sort); #make a readgroup specific, linearized lmfa file. this keeps only one of end-to-end identical sequences via the clever awk clause
         e2eident=$(comm -2 -3 <(echo "$fonmfa") <(echo "$fonlmfa")); #collect sequences that were removed as end-to-end identical, these will be added back later to calculate allele frequencies
       
         #remove identical subsequences as a means of counting them
         rem="";
         while read llu;
           do qu=$(cut -d' ' -f2 <<<"$llu"); #get sequence string
             qt=$(cut -d' ' -f1 <<<"$llu"); #get sample name
             ct=$(grep "$qu" <(echo "$fonlmfa") | wc -l); #count the number of lines that contain the sequence string, if > 1 it is a substring and can be deleted
             if [[ "$ct" > 1 ]];
             then rem+="$llu"$'\n'; #add sequence to list to remove if it is a subsequence of other read pairs
             fi;
           done <<<"$fonlmfa";
         rem=$(sed '/^$/d' <<<"$rem"); #remove trailing blank line
       
         if [[ "$rem" == "" ]];
         then fonfa=$(sort <<<"$fonlmfa"); #there are no identical subsequences so just copy to a final file name
         else
           fonfa=$(grep -v -F -f <(echo "$rem") <(echo "$fonlmfa") | sort); #remove all lines that contain sequences that are subsequences of other lines
         fi;
         
         #count identical sequences and subsequences removed
         ts1=$(wc -l <<<"$fonmfa"); #number of sequences at start
         ts2=$(wc -l <<<"$fonlmfa"); #number of sequences after removing identical
         ts3=$(wc -l <<<"$fonfa");  #number of sequence after removing identical and subsequences
         ni=$(( $ts1 - $ts2 )); #number of identical sequences
         ns=$(( $ts2 - ($ts3) )); #number of identical subsequences



###accumulate counts to report to parallel statement for log.txt###
#echo "#""$inf"."$fon"$'\t'"$ts1":"$ni":"$ns":"$ts3";
basicstats+=$(echo "#""$inf"."$fon"$'\t'"$ts1":"$ni":"$ns":"$ts3")$'\n'; 
###                                               ###



         #Hash sequences from the current readgroup ($k), add to a variable with the readgroup ID.
         #Set up to calculate allele frequencies.  You must add back the end to end identical sequences to $fonfa, but not the identical subsequences.
         #then hash sequences from the current readgroup ($k), add to a variable with the readgroup ID.
       
         fonexcsub=$(echo "$fonfa"$'\n'"$e2eident"); #add back end to end identical sequences to unique sequences contained in $fonfa
         hashrg=$(cut -d' ' -f2 <<<"$fonexcsub" | python -c "exec(\"import hashlib, sys\nfor line in sys.stdin:\n\tprint hashlib.sha224(line).hexdigest()\")" | sed "s/^/$k /"); #hash dna sequences, label with readgroup
         allhash+="$hashrg"$'\n'; #add hashed sequences per readgroup to variable

       done; #for k in $rgs
       allhash=$(sed '/^$/d' <<< "$allhash"); #remove empty line at end


###report accumulated counts to parallel statement for log.txt###
echo "$basicstats" | sed '/^$/d';


       #at this point you have all unique sequences in $fonfa since the last member of $rgs is ">rp", which includes all sequences
       #use this as the basis for searching through each of the readgroup sequence sets to calculate frequencies
       uniqhash=$(cut -d' ' -f2 <<<"$fonfa" | python -c "exec(\"import hashlib, sys\nfor line in sys.stdin:\n\tprint hashlib.sha224(line).hexdigest()\")"); #get hash for unique sequences

       #calculate allele frequency for all unique alleles in each readgroup for current contig:site-range
       reportctsfreqs=""
       for k in $rgs;
       do fon=$(sed 's/>rp/global/' <<<"$k"); #file output name, #replace grep item '>rp' for retrieving all sequences with "global", "global" means all unique haplotypes are counted considering all readgroups simultaneously.
         summ=$(grep ^"$k " <<<"$allhash" | wc -l); #total number of alleles observed
         cts=""; #allele counts
         freqs=""; #allele frequencies
         for u in $uniqhash;
           do nla=$(grep ^"$k $u"$ <<<"$allhash" | wc -l); #number of allele $u
             nlf=$(bc <<<"scale=6;$nla/$summ"); #calculate within pop allele frequency
             cts+="$nla":;
             freqs+="$nlf":; 
           done;
         #echo "@""$inf"."$fon"$'\t'"$cts"$'\t'"$freqs" | sed 's/:$//' | sed "s/:\t/\t/"; #accumulate variable to report allele counts and frequencies
         reportctsfreqs+=$(echo "@""$inf"."$fon"$'\t'"$cts"$'\t'"$freqs" | sed 's/:$//' | sed "s/:\t/\t/")$'\n'; #accumulate variable to report allele counts and frequencies
       done;
###report counts to parallel statement for log.txt###
echo "$reportctsfreqs" | sed '/^$/d'; #report to parallel statement

       #remake the global output file (this is the one for the final alignment, if -mm or -mb)
       #with nothing removed, if -d option not selected. Otherwise, when $dodedup == "YES" (-d switch on),
       #make the global output file with identical (sub)sequences removed
       #only do any of this if the user has agreed to print output (no -x option)
       if [[ "$nooutput" == "NO" ]] && [[ "$dodedup" == "NO" ]]; 
       then sed -e '/^>/s/$/@/' -e 's/^>/#/' <(echo "$mfa") | tr -d '\n' | tr "#" "\n" | tr "@" " " | sed '/^$/d' | awk '{print ">rp" NR "_" $0}' | tr " " "\n" > "$pd"/alignments/"$inf".global.fa; #produce undeduped output for alignment
       elif [[ "$nooutput" == "NO" ]] && [[ "$dodedup" == "YES" ]]
       then echo "$fonfa" | tr ' ' '\n' > "$pd"/alignments/"$inf".global.fa; #produce deduped output for alignment
       #echo "$fonfa" > "$pd"/alignments/"$inf".global.fa; #produce deduped output for alignment
       fi;

}
export -f mycon1;

myalignhaps() {
              i=$1; #i is structured like *jcf7180008848314_28717-32294.global.fa
              thr=$(lscpu | grep "^CPU(s):" | awk '{print $2}'); #max threads
              #rr=$(cut -d'_' -f1 <<<"$i"); #reference contig name, e.g. jcf7180008454378
              #ss=$(cut -d'.' -f1 <<<"$i"); #refcontig+siterange, e.g. jcf7180008454378_303-304
              #tt=$(cut -d'_' -f2 <<<"$ss"); #siterange, e.g. 303-304
              rr=$(rev <<< "$i" | cut -d'_' -f2- | rev ); #reference contig name, e.g. jcf7180008454378
              ss=$(cut -d'.' -f1 <<<"$i"); #refcontig+siterange, e.g. 50_ORF803_jcf7180008454378_303-304
              tt=$(rev <<< "$ss" | cut -d'_' -f1 | rev); #siterange, e.g. 303-304

              #set up for final alignments, include a fragment of the reference contig overlapping the haplotypes
              #add reference sequence to the multi fasta of processed haplotypes
              flankingl=30; #number of bp to extract on each side of the theoretical max and min boundaries of the haplotypes aligned to the reference
              longesth=$(grep -v ^'>' "$pd"/alignments/"$i" | awk '{print length}' | sort -nr | head -1); #find the longest haplotype
              reflength=$(tail -n +2 "$pd"/"$rr"_ref.txt | tr -d '\n' | awk '{print length}');
              le=$(( $(cut -d'-' -f1 <<<"$tt") - $longesth - $flankingl )); #determine the left end of the subsequence to extract from the reference contig
              if (( $le < 0 )); then le=1; fi; #no negative positions allowed
              re=$(( $(cut -d'-' -f2 <<<"$tt") + $longesth + $flankingl )); #determine the right end of the subsequence to extract from the reference contig
              if (( $re > $reflength )); then re=$reflength; fi; #cannot exceed the right end of the reference or you will get an error when you try to extract the 30 bp on each end elongated fragment
              
              refname=$(head -1 "$pd"/"$rr"_ref.txt | sed 's/$/_'$le'_'$re'/'); #name for reference sequence fragment to be included in muscle alignment
              trs=$(grep -v ^'>' "$pd"/"$rr"_ref.txt | tr -d '\n' | cut -c"$le"-"$re"); #add the trimmed reference subsequence to the fasta files for muscle alignment 
              muscin="$refname"$'\n'"$trs"$'\n'$(cat "$pd"/alignments/"$i"); #combine trimmed reference fragment
              echo "$muscin" | tr ' ' '\n' > "$pd"/alignments/"$i"; #overwrite original file containing reads to align with new one that also contains the truncated reference sequence fragment

              if [[ "$debug" == "YES" ]];
              then echo "i = $i
                thr = $thr
                rr = $rr
                ss = $ss
                tt = $tt
                flankingl = $flankingl
                longesth = $longesth
                reflength = $reflength
                le = $le
                re = $re
                refname = $refname
                trs = $trs" >> "$pd"/alignments/aligndb.txt;
              fi;

              #final mapping with bwa
              if [[ $dobwa == "YES" ]];
              then
#              use below for troubleshooting, samtools sort --threads will fail in slurm with n=1 (ncores=1), bwa mem -t "$thr" will not (?)
#                /share/apps/bwa mem -t "$thr" "$pd"/"$rr"_ref.txt "$pd"/alignments/"$i" | \
#                  /share/apps/samtools sort -O BAM --threads "$thr" | \
#                  /share/apps/samtools view -F 2048 -O BAM > "$pd"/alignments/"$ss"_aligned_haps.bam; #-F 2048 excludes supplementary alignments
#                /share/apps/samtools index "$pd"/alignments/"$ss"_aligned_haps.bam;

                #/share/apps/bwa mem -t "$thr" "$pd"/"$rr"_ref.txt "$pd"/alignments/"$i" 2>/dev/null | \
                # /share/apps/samtools sort -O BAM 2>/dev/null | \
                #  /share/apps/samtools view -F 2048 -O BAM 2>/dev/null > "$pd"/alignments/"$ss"_aligned_haps.bam; #-F 2048 excludes supplementary alignments
                #/share/apps/samtools index "$pd"/alignments/"$ss"_aligned_haps.bam 2>/dev/null;

                 bwa mem -t "$thr" "$pd"/"$rr"_ref.txt "$pd"/alignments/"$i" 2>/dev/null | \
                  samtools sort -O BAM 2>/dev/null | \
                  samtools view -F 2048 -O BAM 2>/dev/null > "$pd"/alignments/"$ss"_aligned_haps.bam; #-F 2048 excludes supplementary alignments
                samtools index "$pd"/alignments/"$ss"_aligned_haps.bam 2>/dev/null;
             fi;
              
              #final multiple alignment with muscle
              if [[ $domuscle == "YES" ]];
              then
                #faTMP=$(/share/apps/muscle -quiet -in "$pd"/alignments/"$i"); #capture muscle fasta alignment output in a variable
                faTMP=$(muscle -quiet -in "$pd"/alignments/"$i"); #capture muscle fasta alignment output in a variable

                #put reference sequence at top of muscle output file and sort by read group
                faTMP2=$(sed -e '/^>/s/$/@/' -e 's/^>/#>/' <<<"$faTMP" | tr -d '\n' | tr "#" "\n" | tr "@" " " | sed '/^$/d'); #make output 1 line per sequence
                lrs=$(grep ^"$refname" <<<"$faTMP2" | tr " " "\n"); #get the linearized reference sequence fragment, to add to top of muscle output
                smuscout=$(grep -v ^"$refname" <<<"$faTMP2"| sort -t'_' -k4,4 | tr " " "\n"); #sort muscle aligned haplotypes, to add to revised muscle output
                echo "$lrs"$'\n'"$smuscout" >  "$pd"/alignments/"$ss"_aligned_haps.fa;
              fi;
}
export -f myalignhaps;


##END SUBROUTINES##


#define variables and establish defaults
#ref, path to the reference genome sequence in multi-fasta format
#bam, path to the bam file containing the alignment of reads to ref
#sites, genomic regions to use
stf=1; #samtools view -f option
stF=3852; #samtools view -F option, see https://broadinstitute.github.io/picard/explain-flags.html
stq=60; #samtools view -q option
maxp=1000; #max number of NNNNs used as padding between proper read pairs, implicitly sets an upper limit on insert size
ssh1=""; #default is no --sshloginfile, user may enter a path to machines file with -ssh option
suppar=""; #suppress parallel contig extraction, default is allow GNU parallel --jobs equal to max, -sp switch will set $suppar to --jobs=1 for contig extraction
dodedup=NO; #by default do not remove duplicate (sub)sequences
doalign=NO; #by default do not produce final alignments of extracted haploblocks to their reference
domuscle=NO; #by default do not perform a final muscle alignment of qualified haploblocks
dobwa=NO; #by default do not perform a final bwa mem map of qualified haploblocks
nooutput=NO; #by default do not suppress printing of all output files except the log
debug=NO; #turn debugging off, do not save internal data structures as files

#acquire command line variables
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -r)
    ref="$2"
    shift # past argument
    shift # past value
    ;;
    -b)
    bam="$2"
    shift # past argument
    shift # past value
    ;;
    -o)
    outfol="$2"
    shift # past argument
    shift # past value
    ;;
    -s)
    sites="$2"
    shift # past argument
    shift # past value
    ;;
    -f)
    stf="$2"
    shift # past argument
    shift # past value
    ;;
    -F)
    stF="$2"
    shift # past argument
    shift # past value
    ;;
    -q)
    stq="$2"
    shift # past argument
    shift # past value
    ;;
    -p)
    maxp="$2"
    shift # past argument
    shift # past value
    ;;
    -ssh)
    ssh1="--sshloginfile $2"
    shift # past argument
    shift # past value
    ;;
    -d)
    dodedup=YES
    shift # past argument
    ;;
    -mm)
    doalign=YES
    domuscle=YES
    shift # past argument
    ;;
    -mb)
    doalign=YES
    dobwa=YES
    shift # past argument
    ;;
    -x)
    nooutput=YES
    shift # past argument
    ;;
    -db)
    debug=YES
    shift # past argument
    ;;
    -sp)
    suppar="--jobs 1";
    shift # past argument
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

#determine column number containing read group field
rgf=$(samtools view "$bam" | head -1 | tr "\t" "\n" | sed -n '/RG:Z/=');

#suppress doalign option if -x nooutput option is selected
if [[ "$nooutput" == "YES" ]];
then doalign=NO;
fi;

pd=$(pwd)"/$outfol"; export pd; #path to working directory
if [ ! -d "$pd" ]; then mkdir "$pd"; fi; #make the working directory if not already existing

e=$(cat "$sites"); #content of file $sites
export dodedup;
export nooutput;
export stf;
export stF; 
export stq;
export maxp;
export rgf;
export bam;
export domuscle;
export dobwa;
export debug;

#log
log="$outfol"/log.txt;
date > "$log";
echo >> "$log";
echo "Executable: $0" >> "$log";
echo "Reference sequence: $ref" >> "$log";
echo "Alignment file (bam): $bam" >> "$log";
echo "Target sites: $sites" >> "$log";
echo "samtools view -f $stf" >> "$log";
echo "samtools view -F $stF" >> "$log";
echo "samtools view -q $stq" >> "$log";
echo "Remove duplicate hapblocks: $dodedup" >> "$log";
echo "Map/align hapblocks: $doalign" >> "$log";
echo "Do not generate output files: $nooutput" >> "$log";
echo >> "$log";


#create new reference sequences
#extract each contig listed in $sites from the reference genome to be used by itself as a reference for realignment with the extracted read pairs, then bwa index
echo "Isolating contigs:"; #this step is fast
#Below modify the reference genome to have a '>' as the last line before extracting contig with sed.
#This a hack to take care of the case where the contig of interest is the last one in the reference
echo "$e" | cut -d: -f1 | sort -u | parallel --bar $suppar 'sed -n -e "/^>{}$/,/>/p" '"<(cat "$ref" <(echo \>))"' | sed "$ d" > '"$pd/"'{}_ref.txt'; #sed extracts lines from name of contig of interest until next contig name line, second sed deletes the last line
#echo "$e" | cut -d: -f1 | sort -u | parallel --bar 'sed -n -e "/^>{}$/,/>/p" '"$ref"' | sed "$ d" > '"$pd/"'{}_ref.txt'; #sed extracts lines from name of contig of interest until next contig name line, second sed deletes the last line

echo "Indexing contigs:"; #this step should also be fast
echo "$e" | cut -d: -f1 | sort -u | parallel --bar $suppar 'bwa index '"$pd/"'{}_ref.txt 2>/dev/null'; #index the isolated contigs


#The following single parallel step consolidates independent functions from an earlier prototype of hapx, in order to keep more in memory and less on the drives
#reconstruct fastq file for reads that aligned to each contig, adding readgroup and f/r designation, r is for reads marked with a negative alignment length in column 9
#and
#align *.rp.fq to *_ref.txt, output as bam file, then sambamba index (output files have names like jcf7180008378511_277.bam (contigname_site.bam)
#this step aligns, then verifies mapping quality since it changes from the original values, then realigns
#and
#Process alignments of paired reads via mpileup into a single padded consensus sequence representing the haplotype
if [[ "$nooutput" == "NO" ]];
then
  if [ ! -d "$pd"/"alignments" ]; then mkdir "$pd"/"alignments"; fi; #make a directory to hold alignments if not already existing
fi;

echo "Site.Readgroup"$'\t'"NumHblocks:NumIdenticalHblocks:NumIdenticalHblockSubsequences:NumDistinctHblocks" >> "$log";

echo "Reconstructing haploblocks:";

(echo "$e" | parallel --bar $ssh1 $suppar \
       --env PATH \
       --env pd --env dodedup --env nooutput --env maxp --env rgf --env stf --env stF --env stq --env bam --env debug \
       --env mycon1 --env myiupac --env myinsertion --env myconseq --env mycountqualreadpairs \
       mycon1) >> "$log";

#mycon1 #          Here is a sample reconstructed *.rp.fq file for a single read pair:
#mycon1 #          @E00558:144:HHGCMCCXY:5:2224:6654:40055 1:N:0:AAAAAA
#mycon1 #          CAAAAATATTTTGGTAATTATTCTCAACAAAATGATTTGAAAGGTGTTCATACAACACAAATCGCCTAAGAGACTATGACGGTTTTATCCTCTGATTTGAATTGAGTTTGATCCAAGGGCTTCATATGATTGAAATATAC
#mycon1 #          +
#mycon1 #          FFJJJJAFJ<JJFF7FJAJAA-A7AJJJJJJ<JJJJJJJJJJFJJJJJJJJJJJ<JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFFFAA
#mycon1 #          @E00558:144:HHGCMCCXY:5:2224:6654:40055 2:N:0:AAAAAA
#mycon1 #          ATAAATATGCAAGGAGTCAAAAATATTTGGTAATTAAGCACAAAAACTGATTTGAAATGTGTTCGTACACCACAAATCACCTAAGAGACTATGACGGTTTTACCCTTTGATTTGAATTGAGTTTGATCCAAGGGCTTCATATGCTTTAAA
#mycon1 #          +
#mycon1 #          AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJ

#Produce a pairwise local alignment of the haplotypes to their reference using bwa, which might be viewed in something like IGV (Integrative Genomics Viewer).
#Also produce a multiple alignment of haplotypes to each other using muscle
#final bwa-mem and muscle alignments in parallel
if [[ $doalign == "YES" ]];
then
  echo "Final mapping and alignment:"
  
  #if [ ! -f "$pd"/alignments/*.global.fa ];
  if [[ $(find "$pd"/alignments -name "*.global.fa") == "" ]];
  then echo "0" > "$pd"/alignments/NoReadsSoNoAlignmentPossible; #mark that no reads were found so no alignment is possible
  else find "$pd"/alignments -name "*.global.fa"  | rev | cut -d'/' -f1 | rev | parallel --bar $suppar $ssh1 --env PATH --env pd --env debug --env domuscle --env dobwa --env myalignhaps myalignhaps;
  fi;
fi;

#clean up
if [[ $debug == "NO" ]];
then find "$pd" -name "*_ref.txt*" -print0 | xargs -0 rm;
  #rm "$pd"/*_ref.txt*;
fi;

date >> "$log";

  
