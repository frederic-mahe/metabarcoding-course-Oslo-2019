# Metabarcoding
### easy-peasy


<!-- ----------------------------------------------------------------

                            Prologue

    ----------------------------------------------------------------
-->
## Metabarcoding
#### child of phylogeny
Note: to build a phylogeny, one needs to sequence the same gene from
    many known species. Once such a reference database is available,
    barcoding and metabarcoding become possible (sequence-based
    identification of one or many organisms).


## What you will not find here

* primer design
* experiment design
* DNA extraction and sequencing
* numerical ecology and stats


## What you will find here

* FASTQ merging
* primer clipping
* clustering
* quality filtering
* chimera detection
* table of occurrences



<!-- ----------------------------------------------------------------

                            Code reminder

    ----------------------------------------------------------------
-->
## a primer on code

Command line interface (shell)


#### Tools

* [vsearch](https://github.com/torognes/vsearch),
* [cutadapt](https://github.com/marcelm/cutadapt/),
* [swarm](https://github.com/torognes/swarm),
* bash v4+,
* GNU tools (awk, grep, sed, etc.)
* python 2


#### code block

My pipeline is made of blocks of shell code:

``` bash
# variables
THREADS=4
ENCODING=33

# some comments
vsearch \
    --threads ${THREADS} \
    --fastq_mergepairs R1.fastq.gz \
    --reverse R2.fastq.gz \
    --fastq_ascii ${ENCODING} \
    --fastq_allowmergestagger \
    --quiet \
    --fastqout out.fastq
```


#### redirect

``` bash
# basics
command > output.fastq
command 2> output.log
command 2> /dev/null
command < input.fastq

# but also
>>  2>>  2>&1  <(...)
```

Note: I've been writing shell scripts for 15 years, and I may use
    lesser known aspects of bash. Feel free to ask if my code is
    unclear.


#### wrap

``` bash
# too long to read:
vsearch --threads 4 --fastq_mergepairs R1.fastq.gz --reverse R2.fastq.gz --fastq_ascii 33 --fastq_allowmergestagger --quiet --fastqout out.fastq

# wrapping makes it more readable:
vsearch \
    --threads 4 \
    --fastq_mergepairs R1.fastq.gz \
    --reverse R2.fastq.gz \
    --fastq_ascii 33 \
    --fastq_allowmergestagger \
    --quiet \
    --fastqout out.fastq
```


#### pipe

make data flow

``` bash
# slow
command1 input.fastq > tmp1.fastq
command2 tmp1.fastq > tmp2.fastq
command3 tmp2.fastq > final_output.fastq

# piping avoids temporary files:
command1 input.fastq | \
    command2 | \
    command3 > final_output.fastq
```


#### tee

``` bash
# use a tee to save an intermediary result:
command1 input.fastq | \
    command2 | \
    tee output2.fastq | \
    command3 > final_output.fastq
```


#### test

``` bash
# create toy-examples:
printf ">s_1\nA\n"

# use them to test software behavior:
printf ">s_1\nA\n" | \
    swarm
```

Note: documentation rarely is 100% complete, when you have a doubt
about a tool, create a toy-example to test its behavior



<!-- ----------------------------------------------------------------

                                  FASTQ

    ----------------------------------------------------------------
-->

## FASTQ format


## FASTQ?

go to `./data/`

``` shell
zcat ES1A_S2_L001_R1_001.fastq.gz | less -S
```

``` text
@M05074:97:000000000-BPW9G:1:1101:10203:1383 1:N:0:2
CATAATTTCCTCCGCTTATTGATATGCTTAAGTTCAGCGGGTATCCCTACCTGATCCGAGTTCAACCTAAGAAAGTTGGGGGTTCTGGCGGGTGGACGGCTGAACCCTGTAGCGACAAGTATTACTACGCTTAGAGCCAGACGGCACCGCCACTGCTTTTAAGTGCCGCCGGTACAGCGGGCCCCAAGGCCAAGCAGAGCTTGATTGGTCA
+
@-A-9EFGFFFFD7BFF7FE9,C9F<EFG99,CEF9,@77+@+CCC@F9FCF9,C@C,,+,8C9<CEF,,,,,,,CF,,+++8FEF9,?+++@+++B++@C+,,B?FE8E,,<+++++3C,CF9DF9>>CFE7,,3=@7,,@++@:FC7BC*CC:,7>DF9,,,,7?*=B*5?*:++7***=?EE3***2;***:*0*/;@C8*<C+*<<+
```
note that the quality line starts with a @ ...

Note: Q values are a way to encode on one character numerical values
ranging from 0 to 40 (usually). These values represent the probability
of a wrong base calling for that particular position. Q20 means 1% of
risk, Q30 means 0.1% and Q40 means 0.01%.


## FASTQ format

* [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format#Encoding)
* dominating
* human-readable
* encode quality values (probability of error for each position)
* hard to parse,
* encoding type must be guessed


## Paired-ends

``` shell
zcat ES1A_S2_L001_R1_001.fastq.gz | head -n 1
zcat ES1A_S2_L001_R2_001.fastq.gz | head -n 1
```

``` text
@M05074:97:000000000-BPW9G:1:1101:10203:1383 1:N:0:2
@M05074:97:000000000-BPW9G:1:1101:10203:1383 2:N:0:2
```
each R1 entry has a R2 counterpart.


## check your files

usually provided by the sequencing facility, checks that files were
not modified during transfer:

``` bash
md5sum -c MD5SUM
```

``` text
ES1A_S2_L001_R1_001.fastq.gz: OK
ES1A_S2_L001_R2_001.fastq.gz: OK
```


## in-depth stats on fastq files

[fastqc](https://github.com/s-andrews/FastQC) reports are often provided by sequencing facilities:
``` bash
# run FastQC
for f in *.fastq.gz ; do
    fastqc --quiet ${f}
done

# summarize with multiqc
multiqc --title "Boreal_Forest_Soils"
```
but I prefer to use vsearch


## cumulative expected error (EE)

`$$ EE = \sum_{i=0}^n 10^{-Q/10} $$`

``` bash
# Summarize read quality
for PAIR in 1 2 ; do
    vsearch \
        --quiet \
        --fastq_eestats ES1A_S2_L001_R${PAIR}_001.fastq.gz \
        --output R${PAIR}_eestats.log
done
```
Visualize with some `R` and `tidyverse`

Note: EE is a positive, length-dependent, ever-growing value. EE = 1.0
is often used to discard reads (50% chance of no-error in the read).


## good run

![bad](./images/R1_vs_R2_quality_good.png)


## bad run

![bad](./images/R1_vs_R2_quality_bad.png)


## mixed run

![mix](./images/R1_vs_R2_quality_mix.png)


## guess quality encoding

yes, we need to guess

``` bash
vsearch --fastq_chars ES1A_S2_L001_R1_001.fastq.gz 
```

``` text
vsearch v2.11.1_linux_x86_64, 62.8GB RAM, 8 cores
https://github.com/torognes/vsearch

Reading FASTQ file 100%  
Read 3262836 sequences.
Qmin 40, QMax 71, Range 32
Guess: -fastq_qmin 7 -fastq_qmax 38 -fastq_ascii 33
Guess: Original Sanger format (phred+33)

Letter          N   Freq MaxRun
------ ---------- ------ ------
     A  249004917  25.4%     38
     C  243506540  24.8%    222
...
```
Note: most fastq files use an offset of 33, but beware that you might
encounter older fastq files using a 64 offset.



<!-- ----------------------------------------------------------------

                                Merge

    ----------------------------------------------------------------
-->

## Merging


## Merge R1 and R2 reads

``` bash
FORWARD="ES1A_S2_L001_R1_001.fastq.gz"
REVERSE="${FORWARD/_R1_/_R2_}"
OUTPUT="${FORWARD/_L001*/.fastq}"
LOG="${FORWARD/_L001*/.log}"
ENCODING=33
THREADS=1

vsearch \
    --threads ${THREADS} \
    --fastq_mergepairs ${FORWARD} \
    --reverse ${REVERSE} \
    --fastq_ascii ${ENCODING} \
    --fastq_allowmergestagger \
    --fastqout ${OUTPUT} 2> ${LOG}
```
creates `ES1A_S2.fastq` and `ES1A_S2.log`
Note: allowmergestagger allows to merge reads that are shorter than
the read-length.


## Merging log

``` text
Merging reads
 100%
   3262836  Pairs
   3071725  Merged (94.1%)
    191111  Not merged (5.9%)

Pairs that failed merging due to various reasons:
      3065  too few kmers found on same diagonal
      7936  multiple potential alignments
    115676  too many differences
     64430  alignment score too low, or score drop to high
         4  overlap too short

Statistics of all reads:
    301.00  Mean read length

Statistics of merged reads:
    408.80  Mean fragment length
     53.34  Standard deviation of fragment length
      0.46  Mean expected error in forward sequences
      1.21  Mean expected error in reverse sequences
      0.17  Mean expected error in merged sequences
      0.28  Mean observed errors in merged region of forward sequences
      0.95  Mean observed errors in merged region of reverse sequences
      1.23  Mean observed errors in merged region
```
Note: Torbjørn can help explaining each line.


## Merging theory

![fastq_merging_theory](./images/fastq_merging_theory.png)

Note: from Edgar and Flyvberg, 2015. Merging is an important step that
can change radically the apparent diversity profile of a community
(some popular mergers do not cope well with variable length markers)


## my EE filtering method

* EE / length is constant for read-long markers,
* quality is very high (double-read),
* longer markers accumulate EE, 
* false positives are inevitable
* no early quality-based filtering in my pipeline,
* filtering is done after clustering



<!-- ----------------------------------------------------------------

                                Demultiplex

    ----------------------------------------------------------------
-->

## Demultiplexing


* for each sample,
* add a unique sequence (tag) to all the sequences of a sample,
* pool all the samples into a library,
* sequence that library,
* use the tags to link each sequence to a sample


## Tag list

``` bash
cat sample_barcode.tsv
```

``` text
S001	NAACAAC	NNAACAAC
S002	NNAACCGA	NNNAACCGA
S003	NNNCCGGAA	NCCGGAA
S004	NAGTGTT	NNAGTGTT
S005	NNCCGCTG	NNNCCGCTG
S006	NNNAACGCG	NAACGCG
S007	NGGCTAC	NNGGCTAC
S008	NNTTCTCG	NNNTTCTCG
S009	NNNTCACTC	NTCACTC
S010	NGAACTA	NNGAACTA
S011	NNCCGTCC	NNNCCGTCC
S012	NNNAAGACA	NAAGACA
S013	NCGTGCG	NNCGTGCG
S014	NNGGTAAG	NNNGGTAAG
S015	NNNATAATT	NATAATT
S016	NCGTCAC	NNCGTCAC
S017	NNTTGAGT	NNNTTGAGT
S018	NNNAAGCAG	NAAGCAG
S019	NTTGCAA	NNTTGCAA
S020	NNCACGTA	NNNCACGTA
S021	NNNTAACAT	NTAACAT
S022	NTGCGTG	NNTGCGTG
S023	NNGGTCGA	NNNGGTCGA
S024	NNNCACTCT	NCACTCT
S025	NCTTGGT	NNCTTGGT
S026	NNTCCAGC	NNNTCCAGC
S027	NNNACTTCA	NACTTCA
S028	NGCGAGA	NNGCGAGA
S029	NNTGGAAC	NNNTGGAAC
S030	NNNGTACAC	NGTACAC
S031	NAAGTGT	NNAAGTGT
S032	NNTCTTGG	NNNTCTTGG
S033	NNNAAGGTC	NAAGGTC
S034	NGGCGCA	NNGGCGCA
S035	NNTCGACG	NNNTCGACG
S036	NNNCCTGTC	NCCTGTC
S037	NAGAAGA	NNAGAAGA
S038	NNAATAGG	NNNAATAGG
S039	NNNGGTTCT	NGGTTCT
S040	NTAATGA	NNTAATGA
S041	NNGTAACA	NNNGTAACA
S042	NNNAATCCT	NAATCCT
S043	NAGACCG	NNNAGACCG
S044	NNNTGGCGG	NTGGCGG
S045	NCTATAA	NNCTATAA
S046	NNAATGAA	NNNAATGAA
S047	NNNCGAATC	NCGAATC
S048	NAGAGAC	NNAGAGAC
```
Note: does anyone can explain the Ns?


## Illumina needs variability

and metabarcoding has little

``` text
CCAGCAGCTGCGGTAATTCCAGCTCCAATAGCGTATATTTAAGTTGTTGC...
CCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTTAAGTTGTTGC...
CCAGCACCTGCGGTAATTCCAGCTCCAATAGCGTATATTTAAGTTGTTGC...
CCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTTAAGTTGTTGC...
CCAGCAGCTGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGC...
CCAGCAGCTGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGC...
CCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTTAAGTTGTTGC...
```
* first 5 nucleotides are used to delineate clusters,
* first 25 nucleotides are used to calibrate base-calling

Note: Ns introduce variability in the first flow cycles (4 possible
nucleotides), and variable numbers of N create shifts and even more
variability in the rest of the flow cycles.



<!-- ----------------------------------------------------------------

                                  Trim

    ----------------------------------------------------------------
-->

## Primer clipping

Let's test cutadapt's behavior:

``` bash
# -g, cut after the occurence closest to 5'
printf ">test\naaaACGTggggACGTaaaaa\n" | \
    cutadapt -g ACGT -O 4 - 2> /dev/null

# -a, cut before the occurence most distant to 3'
printf ">test\naaaACGTggggACGTaaaaa\n" | \
    cutadapt -a ACGT -O 4 - 2> /dev/null
```

``` text
>test
ggggACGTaaaaa
>test
aaa
```



## Global dereplication

``` text
vsearch v2.12.0_linux_x86_64, 376.6GB RAM, 48 cores
https://github.com/torognes/vsearch

Dereplicating file - 100%
134902580 nt in 356347 seqs, min 34, max 539, avg 379
Sorting 100%
248072 unique sequences, avg cluster 4.3, median 1, max 221209
Writing output file 100%
```


## Clustering

``` text
Swarm 2.2.2 [Jan 21 2019 22:29:45]
Copyright (C) 2012-2019 Torbjorn Rognes and Frederic Mahe
https://github.com/torognes/swarm

Mahe F, Rognes T, Quince C, de Vargas C, Dunthorn M (2014)
Swarm: robust and fast clustering method for amplicon-based studies
PeerJ 2:e593 https://doi.org/10.7717/peerj.593

Mahe F, Rognes T, Quince C, de Vargas C, Dunthorn M (2015)
Swarm v2: highly-scalable and high-resolution amplicon clustering
PeerJ 3:e1420 https://doi.org/10.7717/peerj.1420

CPU features:      mmx sse sse2 sse3 ssse3 sse4.1 sse4.2 popcnt avx avx2
Database file:     Boreal_forest_soils_18SV4_48_samples.fas
Output file:       Boreal_forest_soils_18SV4_48_samples_1f.swarms
Statistics file:   Boreal_forest_soils_18SV4_48_samples_1f.stats
Resolution (d):    1
Threads:           1
Break OTUs:        Yes
Fastidious:        Yes, with boundary 3

Reading database:  100%  
Indexing database: 100%  
Database info:     93905369 nt in 248072 sequences, longest 539 nt
Hashing sequences: 100%  
Clustering:        100%  

Results before fastidious processing:
Number of swarms:  86193
Largest swarm:     36648

Counting amplicons in heavy and light swarms 100%  
Heavy swarms: 7537, with 167534 amplicons
Light swarms: 78656, with 80538 amplicons
Total length of amplicons in light swarms: 30508515
Bloom filter: bits=16, m=3416953680, k=11, size=407.3MB
Adding light swarm amplicons to Bloom filter 100%  
Generated 204824711 variants from light swarms
Checking heavy swarm amplicons against Bloom filter 100%  
Heavy variants: 425554746
Got 135072 graft candidates
Grafting light swarms on heavy swarms 100%  
Made 31505 grafts

Writing swarms:    100%  
Writing seeds:     100%  
Writing structure: 100%  
Writing stats:     100%  

Number of swarms:  54688
Largest swarm:     43300
Max generations:   15
```


## Chimera detection

``` text
vsearch v2.12.0_linux_x86_64, 376.6GB RAM, 48 cores
https://github.com/torognes/vsearch

Reading file Boreal_forest_soils_18SV4_48_samples_1f_representatives.fas 100%  
20633387 nt in 54688 seqs, min 34, max 539, avg 377
Masking 100% 
Sorting by abundance 100%
Counting k-mers 100% 
Detecting chimeras 100%  
Found 26354 (48.2%) chimeras, 26942 (49.3%) non-chimeras,
and 1392 (2.5%) borderline sequences in 54688 unique sequences.
Taking abundance information into account, this corresponds to
60430 (5.7%) chimeras, 988411 (93.5%) non-chimeras,
and 8076 (0.8%) borderline sequences in 1056917 total sequences.
```




<!-- ----------------------------------------------------------------

                               Limbo

    ----------------------------------------------------------------
-->

## Boundary

Save what can be saved below a certain threshold. Add the weight of
small OTUs to bigger ones.

<!-- two newlines for a new vertical slide  -->
<!-- three newlines for a new horizontal slide  -->
