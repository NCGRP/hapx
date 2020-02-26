#!/bin/bash

# Usage: see README.txt

##SUBROUTINES##

#myconseq locates gaps in sequence
myconseq() {
          echo "$1" | awk '$1!=p+1{print p+1"-"$1-1}{p=$1}';
}
export -f myconseq;

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

#mycountqualreadpairs counts how many paired reads and unpaired reads map to the target sites
mycountqualreadpairs() {
                       j=$1;
                       fn=$(echo "$j" | tr ":" "_")".tmp"; #name of input file
                       p=$(cut -d$'\t' -f1 "$pd"/"$fn" | sort | uniq -c | grep ^" \+2 "); #paired reads
                       s=$(cut -d$'\t' -f1 "$pd"/"$fn" | sort | uniq -c | grep ^" \+1 "); #single ends of paired reads
                       pp=$(echo "$p" | awk 'NF' | wc -l); #number of paired reads (awk 'NF' removes empty lines)
                       ss=$(echo "$s" | awk 'NF' | wc -l); #number of single ends of paired reads
                       echo "$j $pp $ss";
}
export -f mycountqualreadpairs;

#mycon1 combines old functions into one so that data can be passed in memory instead of using the disk
mycon1() {
       #myreconfq reconstructs a single fastq file that contains exclusively both reads of a read pair from a sam file
       i="$1"; #input is a description of sites to process in contigname_range format like jcf7180008531951_276-301
               #the sites are accessed from the file "$pd"/"$i".tmp below
       h=$(echo "$i" | cut -d'_' -f1); #just the contig name
 
       s=$(awk -F$'\t' -v rgf=$rgf '{print "@"$1"_"$rgf,$10,"+",$11}' "$pd"/"$i".tmp | tr " " "\n" | sed '/^$/d');
       
       m=$(grep ^'@' <(echo "$s") | cut -d: -f1-9 | sort -u | sed '/^$/d'); #capture list of unique read pair names (includes also unpaired reads with unique names) that map to contig_site $i 

       mfa=""; #initialize variable to contain all the read pair haplotypes
       #below j contains a readgroup-annotated read pair name like @E00558:144:HHGCMCCXY:4:2210:23074:34043_RG:Z:55, that aligns to the current contig:site range
       for j in $m;
         do #echo "j=$j";
           rpfq=$(grep -A3 "$j" <(echo "$s") | sed '0,/\(^@\S\+$\)/ s/\(^@\S\+$\)/\1 1:N:0:AAAAAA/' | sed 's/\(^@\S\+$\)/\1 2:N:0:AAAAAA/'); #save read pair into a variable containing a reconstructed fastq, labeled as a read pair for bwa
           #grep -A3 "$j" <(echo "$s) | sed '0,/\(^@\S\+$\)/ s/\(^@\S\+$\)/\1 1:N:0:AAAAAA/' | sed 's/\(^@\S\+$\)/\1 2:N:0:AAAAAA/' > "$pd"/"$i"_"$fn".rp.fq; #save read pair into its own reconstructed fastq file, labeled as a read pair for bwa
           #above, the first sed searches from 0 to the pattern then performs substitution 1 on the fastq defline, which causes the insertion of whitespace.  The second sed then searches the whole file and performs substitution 2 on the second fastq defline. \S means 'not whitespace'
           #above, create single files containing both reads of a read pair




           #myrealign takes the reads extracted from the primary alignment ($rpfq) and realigns them independently to the contig as reference
           sbm=$(/share/apps/bwa mem -p -M "$pd"/"$h"_ref.txt <(echo "$rpfq") 2>/dev/null | /share/apps/samtools sort -O SAM 2>/dev/null); #capture human readable sam file in a variable, it's small since only involved 2 reads and a reference (-p paired interleaved, -M label split reads as secondary so samtools mpileup excludes them)
 
           #mapping quality changes upon realignment, extract only reads with user defined properties (/share/apps/samtools view -fFq), convert those to proper fastq file and re-align again
           t=$(/share/apps/samtools view -f "$stf" -F "$stF" -q "$stq" <(echo "$sbm")  2>/dev/null | awk -F$'\t' '{print "@"$1,$10,"+",$11}' | tr " " "\n" | sed '0,/\(^@\S\+$\)/ s/\(^@\S\+$\)/\1 1:N:0:AAAAAA/' | sed 's/\(^@\S\+$\)/\1 2:N:0:AAAAAA/');
           
           if [[ $(echo "$t") != "" ]]; #only re-realign if there are reads left after latest quality filter
           then x2xsam=$(/share/apps/bwa mem -p -M "$pd"/"$h"_ref.txt <(echo "$t") 2>/dev/null | awk -F$'\t' '$3!="*"{print $0}'); #perform second realignment and manually remove any unmapped reads that may have been found
           else continue; #no reads left to align, continue to next
           fi;
 
           #test whether re-realignment still contains reads after removing unmapped reads, if not skip to next read pair
           lnsam=$(echo "$x2xsam" | wc -l);
           if [[ "$lnsam" < 3 ]]; then continue; fi; 
  




           #mymakehapblocks() processes read pairs into a single contiguous sequence, adds NNNs where opposing read pairs do not overlap to pad relative to the reference, adds IUPAC redundancy codes for conflicts between pairs of a read, removes deletion coding, retains insertions
           rp=$(echo "$x2xsam" | awk -F$'\t' -v h=$h '$3==h{print $1}' | sort -u | tr ":" "_"); #extract the unique read name from the sam formatted data that contains a single mapped read pair
           rname="$i"_"$rp"; #get a root name for the read pair/contig alignment
           
           mp=$(/share/apps/samtools mpileup -A <(echo "$x2xsam") 2>/dev/null); #form the pileup file, use -A to include orphan reads. samtools mpileup by default excludes improper pairs. in hapx, bwa mem will find no proper pairs because too few reads are used during alignment
           
           #Extract the pos and base columns, remove read start ^., read end $, remove deletion indicators (e.g. -2AC, they are followed with *), pass to myinsertion subroutine to control insertions, then to myiupac to recode conflicts as ambiguous and produce the consensus sequence
           base1=$(echo "$mp" | cut -d$'\t' -f2,5 | sed 's/\^.//g' | sed 's/\$//g' | sed 's/-.*//g' | tr "\n" " " | myinsertion | myiupac); #haplotype extraction, unpadded as of now
           
           #Deal with contiguity and padding
           csq=$(myconseq "$(echo "$mp" | cut -d$'\t' -f2)" | grep -v "^1\-" | sed 's/-/ 1 /g'); #identify regions where reads are not aligned consecutively to the reference, exclude the range from 1-start of overlap, set up to use as an interval for seq command
           
           #calculate the maximum pad size, if any pad is larger than $maxp continue to next contig:read-pair
           maxn=$(echo "$csq" | awk -F' ' '{print $3 - $1}' | sort -nr | head -1); 
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
           nrname1=$(echo "$rname" | rev | cut -d'_' -f1-3 | rev); #readgroup info from end of string
           nrname2=$(echo "$rname" | rev | cut -d'_' -f4- | rev); #rest of string, excluding readgroup info
           nrname="$nrname1"_"$nrname2";
           
           #base3 contains the processed haplotype as a fasta file, deletions relative to reference removed.
           #base3 is as close as we can come to reconstructing the native molecule
           base3=$(echo ">$nrname";echo "$base2" | cut -d' ' -f2 | grep -v '-' | tr "\n" " " | sed 's/ //g');

           mfa+="$base3"$'\n';
           
         done; # for j in $m
         
       #remove terminal line break in major output variable $mfa
       mfa=$(echo "$mfa" | sed '/^$/d');


#clean up
#rm "$pd"/"$i".tmp;
       


        #mydedup() counts and optionally removes duplicate sequences and subsequences that are exactly contained within longer sequences, from the processed multi fasta file
        inf="$i"; #value like jcf7180008587925_40-41
        
        #count duplicated sequences and subsequences
        # to this process supply each read group separately, and everything together
        da=$(sed -e '/^>/s/$/@/' -e 's/^>/#/' <(echo "$mfa") | tr -d '\n' | tr "#" "\n" | tr "@" " " | sed '/^$/d' | awk '{print ">rp" NR "_" $0}'); #linearize the data set
        rgs=$(echo "$da" | cut -d'_' -f2-4 | sort -u | tr '\n' ' '); #acquire a list of readgroups
        rgs="$rgs>rp"; #add an item to the list of readgroups that will allow all sequences to be grepped from the file '>rp'
        
        for k in $rgs;
        do fon=$(echo "$k" | sed 's/>rp/global/'); #file output name, #replace grep item '>rp' for retrieving all sequences with "global", "global" means all unique haplotypes are counted considering all readgroups simultaneously.

          fonmfa=$(echo "$da" | grep "$k"); #make a readgroup specific mfa file. this has all haploblocks, including duplicated and subsequence
          fonlmfa=$(echo "$da" | grep "$k" | awk -F' ' '!_[$2]++'); #make a readgroup specific, linearized lmfa file. this has end-to-end identical sequences removed via the clever awk clause 
        
          #remove identical subsequences as a means of counting them
          rem="";
          while read llu;
            do qu=$(echo "$llu" | cut -d' ' -f2); #get sequence string
              qt=$(echo "$llu" | cut -d' ' -f1); #get sample name
              ct=$(grep "$qu" <(echo "$fonlmfa") | wc -l); #count the number of lines that contain the sequence string, if > 1 it is a substring and can be deleted
              if [[ "$ct" > 1 ]];
              then rem+="$llu"$'\n'; #add sequence to list to remove if it is a subsequence of other read pairs
              fi;
            done <<<"$fonlmfa";
          rem=$(echo "$rem" | sed '/^$/d'); #remove trailing blank line
        
          if [[ "$rem" == "" ]];
          then fonfa=$(echo "$fonlmfa" | tr " " "\n"); #there are no identical subsequences so just copy to a final file name
          else
            fonfa=$(grep -v -F -f <(echo "$rem") <(echo "$fonlmfa") | tr " " "\n"); #remove all lines that contain sequences that are subsequences of other lines
          fi;
        
          #count identical sequences and subsequences removed
          ts1=$(echo "$fonmfa" | wc -l); #number of sequences at start
          ts2=$(echo "$fonlmfa" | wc -l); #number of sequences after removing identical
          ts3=$(echo "$fonfa" | wc -l);  #number of sequence after removing identical and subsequences (2 lines per sequence)
          ni=$(( $ts1 - $ts2 )); #number of identical sequences
          ns=$(( $ts2 - ($ts3/2) )); #number of identical subsequences



          ###report counts to parallel statement for log.txt###
          echo "#""$inf"."$fon" $'\t'"$ts1":"$ni":"$ns":$(( $ts3/2 )); 



        done;
         

        #remake the global output file, with nothing removed, if -d option not selected. Otherwise, when $dodedup == "YES", make the global output file with identical (sub)sequences removed
        if [[ "$dodedup" == "NO" ]]; 
        then sed -e '/^>/s/$/@/' -e 's/^>/#/' <(echo "$mfa") | tr -d '\n' | tr "#" "\n" | tr "@" " " | sed '/^$/d' | awk '{print ">rp" NR "_" $0}' | tr " " "\n" > "$pd"/alignments/"$inf".global.fa; 
        else echo "$fonfa" > "$pd"/alignments/"$inf".global.fa;
        fi;
      
}
export -f mycon1;




#v1 #mydedup() removes duplicate sequences, and subsequences that are exactly contained within longer sequences, from the processed multi fasta file
#v1 mydedup() {
#v1         inf="$1";
#v1         
#v1         #count duplicated sequences and subsequences
#v1         # to this process supply each read group separately, and everything together
#v1         da=$(sed -e '/^>/s/$/@/' -e 's/^>/#/' "$pd"/alignments/"$inf".mfa | tr -d '\n' | tr "#" "\n" | tr "@" " " | sed '/^$/d' | awk '{print ">rp" NR "_" $0}'); #linearize the data set
#v1         rgs=$(echo "$da" | cut -d'_' -f2-4 | sort -u | tr '\n' ' '); #acquire a list of readgroups
#v1         rgs="$rgs>rp"; #add an item to the list of readgroups that will allow all sequences to be grepped from the file '>rp'
#v1         
#v1         for i in $rgs;
#v1         do fon=$(echo "$i" | sed 's/>rp/global/'); #file output name, #replace grep item '>rp' for retrieving all sequences with "global", meaning all unique haplotypes are counted when "global" is reported
#v1 
#v1           fonmfa=$(echo "$da" | grep "$i"); #make a readgroup specific mfa file. this has all haploblocks, including duplicated and subsequence
#v1           fonlmfa=$(echo "$da" | grep "$i" | awk -F' ' '!_[$2]++'); #make a readgroup specific, linearized lmfa file. this has end-to-end identical sequences removed via the clever awk clause 
#v1         
#v1           #remove identical subsequences as a means of counting them
#v1           rem="";
#v1           while read llu;
#v1             do qu=$(echo "$llu" | cut -d' ' -f2); #get sequence string
#v1               qt=$(echo "$llu" | cut -d' ' -f1); #get sample name
#v1               ct=$(grep "$qu" <(echo "$fonlmfa") | wc -l); #count the number of lines that contain the sequence string, if > 1 it is a substring and can be deleted
#v1               if [[ "$ct" > 1 ]];
#v1               then rem+="$llu"$'\n'; #add sequence to list to remove if it is a subsequence of other read pairs
#v1               fi;
#v1             done <<<"$fonlmfa";
#v1           rem=$(echo "$rem" | sed '/^$/d'); #remove trailing blank line
#v1         
#v1           if [[ "$rem" == "" ]];
#v1           then fonfa=$(echo "$fonlmfa" | tr " " "\n"); #there are no identical subsequences so just copy to a final file name
#v1           else
#v1             fonfa=$(grep -v -F -f <(echo "$rem") <(echo "$fonlmfa") | tr " " "\n"); #remove all lines that contain sequences that are subsequences of other lines
#v1           fi;
#v1         
#v1           #count identical sequences and subsequences removed
#v1           ts1=$(echo "$fonmfa" | wc -l); #number of sequences at start
#v1           ts2=$(echo "$fonlmfa" | wc -l); #number of sequences after removing identical
#v1           ts3=$(echo "$fonfa" | wc -l);  #number of sequence after removing identical and subsequences (2 lines per sequence)
#v1           ni=$(( $ts1 - $ts2 )); #number of identical sequences
#v1           ns=$(( $ts2 - ($ts3/2) )); #number of identical subsequences
#v1 
#v1           #report counts to parallel statement
#v1           echo "#""$inf"."$fon" $'\t'"$ts1":"$ni":"$ns":$(( $ts3/2 )); 
#v1 
#v1         done;
#v1         
#v1           #clean up
#v1           #find "$pd"/alignments -name "*$inf.*.mfa" | xargs rm; #remove mfa files for this site
#v1           #find "$pd"/alignments -name "*$inf.*.lmfa" | xargs rm;
#v1           #find "$pd"/alignments -name "*$inf.*.fa" | grep -v global | xargs rm; #retain the global final fa file, does not contain duplicate sequences and subsequences
#v1          
#v1           #remake the global output file, with nothing removed, if -d option not selected
#v1           if [[ "$dodedup" == "NO" ]]; 
#v1           then sed -e '/^>/s/$/@/' -e 's/^>/#/' "$pd"/alignments/"$inf".mfa | tr -d '\n' | tr "#" "\n" | tr "@" " " | sed '/^$/d' | awk '{print ">rp" NR "_" $0}' | tr " " "\n" > "$pd"/alignments/"$inf".global.fa; 
#v1           fi;
#v1 
#v1 }
#v1 export -f mydedup;


#myalignhaps() aligns the extracted read pair haplotypes to their reference contig (bwa mem), or to each other (multiple alignment with muscle) for visualization or further processing by user
myalignhaps() {
              i=$1;
              thr=$(lscpu | grep "^CPU(s):" | awk '{print $2}'); #max threads
              rr=$(echo "$i" | cut -d'_' -f1); #reference contig name, e.g. jcf7180008454378
              ss=$(echo "$i" | cut -d'.' -f1); #refcontig+siterange, e.g. jcf7180008454378_303-304
              tt=$(echo "$ss" | cut -d'_' -f2); #siterange, e.g. 303-304

              #final mapping with bwa
              #echo "Final mapping with bwa: $i";
              /share/apps/bwa mem -t "$thr" "$pd"/"$rr"_ref.txt "$pd"/alignments/"$i" 2>/dev/null | \
                  /share/apps/samtools sort -O BAM --threads "$thr" 2>/dev/null | \
                  /share/apps/samtools view -F 2048 -O BAM 2>/dev/null > "$pd"/alignments/"$ss"_aligned_haps.bam; #-F 2048 excludes supplementary alignments
              /share/apps/samtools index "$pd"/alignments/"$ss"_aligned_haps.bam 2>/dev/null;
              
              #final multiple alignment with muscle, include a fragment of the reference contig overlapping the haplotypes
              #add reference sequence to the multi fasta of processed haplotypes
              flankingl=30; #number of bp to extract on each side of the theoretical max and min boundaries of the haplotypes aligned to the reference
              longesth=$(grep -v ^'>' "$pd"/alignments/"$i" | awk '{print length}' | sort -nr | head -1); #find the longest haplotype
              le=$(( $(echo "$tt" | cut -d'-' -f1) - $longesth - $flankingl )); #determine the left end of the subsequence to extract from the reference contig
              if (( $le < 0 )); then le=1; fi; #no negative positions allowed
              re=$(( $(echo "$tt" | cut -d'-' -f2) + $longesth + $flankingl )); #determine the right end of the subsequence to extract from the reference contig
              
              refname=$(head -1 "$pd"/"$rr"_ref.txt | sed 's/$/_'$le'_'$re'/'); #name for reference sequence fragment to be included in muscle alignment
              #echo "$refname" > "$pd"/alignments/"$i".muscle; #add fasta header line of trimmed reference sequence with range noted for muscle alignment
              trs=$(grep -v ^'>' "$pd"/"$rr"_ref.txt | tr -d '\n' | cut -c"$le"-"$re"); #add the trimmed reference subsequence to the fasta files for muscle alignment 
              #cat "$pd"/alignments/"$i" >> "$pd"/alignments/"$i".muscle; #add processed haplotypes to multi fasta for muscle
              
              muscin="$refname"$'\n'"$trs"$'\n'$(cat "$pd"/alignments/"$i"); #combine trimmed reference fragment
              
              #cat "$pd"/"$rr"_ref.txt "$pd"/alignments/"$i" > "$pd"/alignments/"$i".muscle;
              #echo "Multiple alignment with muscle: $i";
              faTMP=$(/share/apps/muscle -quiet -in <(echo "$muscin")); #capture muscle fasta alignment output in a variable
              
              #put reference sequence at top of muscle output file and sort by read group
              faTMP2=$(echo "$faTMP" | sed -e '/^>/s/$/@/' -e 's/^>/#>/' | tr -d '\n' | tr "#" "\n" | tr "@" " " | sed '/^$/d'); #make output 1 line per sequence
              lrs=$(echo "$faTMP2" | grep ^"$refname" | tr " " "\n"); #get the linearized reference sequence fragment, to add to top of muscle output
              smuscout=$(echo "$faTMP1" | grep -v ^"$refname" | sort -t'_' -k4,4 | tr " " "\n"); #sort muscle aligned haplotypes, to add to revised muscle output
              echo "$lrs"$'\n'"$smuscout" >  "$pd"/alignments/"$ss"_aligned_haps.fa
              
              #clean up
              #rm "$pd"/alignments/"$ss"_aligned_haps.faTMP;
              #rm "$pd"/alignments/"$ss"_aligned_haps.faTMP1;
}
export -f myalignhaps;


##END SUBROUTINES##


#define variables and establish defaults
#ref, path to the reference genome sequence in multi-fasta format
#bam, path to the bam file containing the alignment of reads to ref
#alnr, aligner used to create bam file (gem, bwamem)
#sites, genomic regions to use
stf=1; #samtools view -f option
stF=3852; #samtools view -F option, see https://broadinstitute.github.io/picard/explain-flags.html
stq=60; #samtools view -q option
maxp=1000; #max number of NNNNs used as padding between proper read pairs, implicitly sets an upper limit on insert size
dodedup=NO; #by default do not remove duplicate (sub)sequences
doalign=NO; #by default do not produce final alignments of extracted haploblocks to their reference

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
    -a)
    alnr="$2"
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
    -d)
    dodedup=YES
    shift # past argument
    ;;
    -m)
    doalign=YES
    shift # past argument
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters


#process command line parameters
if [[ "$alnr" == "gem" ]];
then rgf=12; #read group info in field 12 for gem bam output
elif [[ "$alnr" == "bwamem" ]];
then rgf=17; #read group info in field 17 for bwamem bam output
else echo "Unrecognized aligner: $alnr.  Quitting...";
  return;
  exit;
fi;

pd=$(pwd)"/$outfol"; export pd; #path to working directory
mkdir "$pd";

e=$(cat "$sites"); #format sites for proper parsing by samtools view
export dodedup;
export stf;
export stF; 
export stq;
export maxp;
export rgf;

#log
log="$outfol"/log.txt;
date > "$log";
echo >> "$log";
echo "Reference sequence: $ref" >> "$log";
echo "Alignment software: $alnr" >> "$log";
echo "Alignment file (bam): $bam" >> "$log";
echo "Target sites: $sites" >> "$log";
echo "samtools view -f $stf" >> "$log";
echo "samtools view -F $stF" >> "$log";
echo "samtools view -q $stq" >> "$log";
echo "Remove duplicate hapblocks: $dodedup" >> "$log";
echo "Map/align hapblocks: $doalign" >> "$log";
echo >> "$log";





#extract read pairs at target sites using samtools. Obey include, exclude and quality rules from command line
echo "Extracting read pairs:";
echo "$e" | parallel --bar '/share/apps/samtools view -f '"$stf"' -F '"$stF"' -q '"$stq"' '"$bam"' "{}" | sort > '"$pd/"'{}.tmp; \
  mv '"$pd/"'{}.tmp '"$pd/"'$(echo {} | tr ":" "_").tmp';
  
#describe number of read pairs and single ends of a read pair that have met the samtools -f/-F/-q rules, extracted into a tabular format
#this should be parallelized, can operate on head node only
echo "Counting qualified read pairs:"
f=$(echo "$e" | parallel --keep-order --bar mycountqualreadpairs);

g=$(echo "$f" | awk -F' ' '$2+$3!=0 {print $0}'); #remove any contigs from consideration when there are 0 reads (no paired reads:$2, no unpaired reads:$3) that align to them (this could also be a variable that supports a cutoff)
echo "$f" | awk -F' ' '$2+$3==0 {print $0}' | cut -d' ' -f1 | sed 's/$/.tmp/' | sed 's/:/_/' | sed 's:^:'$pd'/:' | xargs rm; #remove .tmp files from contigs with 0 reads
if [[ "$g" == "" ]];
then echo "No sites with reads. Quitting...";
  return;
  exit;
fi;

#extract each contig from the reference genome to be used by itself as a reference for realignment with the extracted read pairs, then bwa index
echo "Isolating contigs:"; #this step is fast
echo "$g" | cut -d: -f1 | sort -u | parallel --bar 'sed -n -e "/^>{}$/,/>/p" '"$ref"' | sed "$ d" > '"$pd/"'{}_ref.txt'; #sed extracts lines from name of contig of interest until next contig name line, second sed deletes the last line
echo "Indexing contigs:"; #this step should also be fast
echo "$g" | cut -d: -f1 | sort -u | parallel --bar 'bwa index '"$pd/"'{}_ref.txt 2>/dev/null'; #index the isolated contigs


#The following single parallel step consolidates independent functions from an earlier prototype of hapx, in order to keep more in memory and less on the drives
#reconstruct fastq file for reads that aligned to each contig, adding readgroup and f/r designation, r is for reads marked with a negative alignment length in column 9
#and
#align *.rp.fq to *_ref.txt, output as bam file, then sambamba index (output files have names like jcf7180008378511_277.bam (contigname_site.bam)
#this step aligns, then verifies mapping quality since it changes from the original values, then realigns
#and
#Process alignments of paired reads via mpileup into a single padded consensus sequence representing the haplotype
#if [ ! -d "haplotypes" ]; then mkdir "$pd"/"haplotypes"; fi; #make a directory to hold fasta sequences containing extracted haplotypes
if [ ! -d "alignments" ]; then mkdir "$pd"/"alignments"; fi; #make a directory to hold alignments

echo "Site.Readgroup"$'\t'"NumHblocksStart:NumIdenticalReads:NumIdenticalSubsequences:NumHblocksFinal" >> "$log";

echo "Evaluating mapping quality:";
mydd=$(echo "$g" | cut -d' ' -f1 | tr ':' '_' | parallel --bar --sshloginfile /home/reevesp/machines --env pd --env dodedup --env maxp --env rgf --env mycon1 --env myiupac --env myinsertion --env myconseq mycon1);

echo "$mydd" >> "$log";





#mycon1 #          Here is a sample reconstructed *.rp.fq file for a single read pair:
#mycon1 #          @E00558:144:HHGCMCCXY:5:2224:6654:40055 1:N:0:AAAAAA
#mycon1 #          CAAAAATATTTTGGTAATTATTCTCAACAAAATGATTTGAAAGGTGTTCATACAACACAAATCGCCTAAGAGACTATGACGGTTTTATCCTCTGATTTGAATTGAGTTTGATCCAAGGGCTTCATATGATTGAAATATAC
#mycon1 #          +
#mycon1 #          FFJJJJAFJ<JJFF7FJAJAA-A7AJJJJJJ<JJJJJJJJJJFJJJJJJJJJJJ<JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFFFAA
#mycon1 #          @E00558:144:HHGCMCCXY:5:2224:6654:40055 2:N:0:AAAAAA
#mycon1 #          ATAAATATGCAAGGAGTCAAAAATATTTGGTAATTAAGCACAAAAACTGATTTGAAATGTGTTCGTACACCACAAATCACCTAAGAGACTATGACGGTTTTACCCTTTGATTTGAATTGAGTTTGATCCAAGGGCTTCATATGCTTTAAA
#mycon1 #          +
#mycon1 #          AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJ


#v1 #Count duplicate sequences and identical subsequences
#v1 (echo -n "Counting ";
#v1 if [[ "$dodedup" == YES ]];
#v1 then echo -n "and removing ";
#v1 fi;
#v1 echo "duplicates:";);
#v1 
#v1 echo "Site.Readgroup"$'\t'"NumHblocksStart:NumIdenticalReads:NumIdenticalSubsequences:NumHblocksFinal" >> "$log";
#v1 mydd=$(find "$pd"/alignments -name "*.mfa" | rev | cut -d'/' -f1 | rev | cut -d. -f1 | parallel --bar --sshloginfile /home/reevesp/machines --env pd --env dodedup --env mydedup mydedup);
#v1 echo "$mydd" >> "$log";



#Produce a pairwise local alignment of the haplotypes to their reference using bwa, which might be viewed in something like IGV (Integrative Genomics Viewer).
#Also produce a multiple alignment of haplotypes to each other using muscle
#final bwa-mem and muscle alignments in parallel
if [[ $doalign == "YES" ]];
then
  echo "Final mapping and alignment:"
  find "$pd"/alignments -name "*.global.fa"  | rev | cut -d'/' -f1 | rev | parallel --bar --sshloginfile /home/reevesp/machines --env pd --env myalignhaps myalignhaps;
fi;

#clean up
#rm *_ref.txt*;

  
