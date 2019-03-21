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
* sample design
* DNA extraction and sequencing
* numerical ecology and stats


## What you will find here

* FASTQ format
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
* python 2 or 3?,
* biopython?


#### code block

My pipeline is made of blocks of shell code:

``` bash
THREADS=4
ENCODING=33

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
command < input.fastq

# lesser known
>>  2>>  2>&1  <(...)
```

Note: I've been writing shell scripts for 15 years, and I may use
    lesser known aspects of bash. Feel free to ask if my code is
    unclear.


#### wrap

``` bash
# too long to read:
vsearch --threads 4 --fastq_mergepairs R1.fastq.gz --reverse R2.fastq.gz --fastq_ascii 33 --fastq_allowmergestagger --quiet --fastqout out.fastq

# wrapping makes it better:
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
# use a tee to save an intermediary file:
command1 input.fastq | \
    command2 | \
    tee output2.fastq | \
    command3 > final_output.fastq
```



<!-- ----------------------------------------------------------------

                               FASTQ

    ----------------------------------------------------------------
-->
## FASTQ format

* dominating
* human-readable
* not ideal
* 4 different quality encodings
* see [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format#Encoding)



#### check your files

``` bash
md5sum -c MD5SUM
```

``` text
ES1A_S2_L001_R1_001.fastq.gz: OK
ES1A_S2_L001_R2_001.fastq.gz: OK
```


#### guess quality encoding

yes, we have to guess.

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



## Merging

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


## Global dereplication

``` text
Dereplicating file - 100%
27696348 nt in 73299 seqs, min 34, max 539, avg 378
Sorting 100%
51693 unique sequences, avg cluster 3.3, median 1, max 37137
```


## Clustering

``` text
Swarm 2.2.2 [Feb 15 2019 15:24:29]
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
Database info:     19536633 nt in 51693 sequences, longest 539 nt
Hashing sequences: 100%
Clustering:        100%

Results before fastidious processing:
Number of swarms:  21827
Largest swarm:     7352

Counting amplicons in heavy and light swarms 100%
Heavy swarms: 2025, with 31269 amplicons
Light swarms: 19802, with 20424 amplicons
Total length of amplicons in light swarms: 7708077
Bloom filter: bits=16, m=863304624, k=11, size=102.9MB
Adding light swarm amplicons to Bloom filter 100%
Generated 51747582 variants from light swarms
Checking heavy swarm amplicons against Bloom filter 100%
Heavy variants: 79396807
Got 24981 graft candidates
Grafting light swarms on heavy swarms 100%
Made 6433 grafts

Writing swarms:    100%
Writing seeds:     100%
Writing structure: 100%
Writing stats:     100%

Number of swarms:  15394
Largest swarm:     8941
Max generations:   14
```


## Chimera detection

``` text
vsearch v2.11.1_linux_x86_64, 62.8GB RAM, 8 cores
https://github.com/torognes/vsearch

Reading file Boreal_forest_soils_18SV4_48_samples_1f_representatives.fas 100%
5786056 nt in 15394 seqs, min 34, max 539, avg 376
Masking 100% 
Sorting by abundance 100%
Counting k-mers 100% 
Detecting chimeras 100%
Found 6689 (43.5%) chimeras, 8406 (54.6%) non-chimeras,
and 299 (1.9%) borderline sequences in 15394 unique sequences.
Taking abundance information into account, this corresponds to
11100 (6.6%) chimeras, 157061 (92.7%) non-chimeras,
and 1197 (0.7%) borderline sequences in 169358 total sequences.
```




<!-- ----------------------------------------------------------------

                               Limbo

    ----------------------------------------------------------------
-->

## Demonstration

Let's start with redirections and then
[swarm](https://github.com/torognes/swarm).

``` bash
for i in {1..3} ; do
    echo "test"
done
```
`$$ J(\theta_0,\theta_1) = \sum_{i=0} $$`

Note: don't forget to say that.


accents are supported: Frédéric and Torbjørn, 中华人民共和国.
```python
def do_something():
    action = "did something"
    print(action)
    return 0
```

Note: This will only display in the notes window.


Insert an image:

![How to draw an owl](https://typewritermonkeys.files.wordpress.com/2016/03/how-to-draw-an-owl.jpg)


## Boundary

Save what can be saved below a certain threshold. Add the weight of
small OTUs to bigger ones.

<!-- two newlines for a new vertical slide  -->
<!-- three newlines for a new horizontal slide  -->
