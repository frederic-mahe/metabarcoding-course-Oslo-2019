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


#### redirect

``` bash
# basics
command > output.fastq
command 2> output.log
command < input.fastq

# lesser known
>>  2>>  2>&1  <>
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
