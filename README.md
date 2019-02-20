# Metabarcoding course (Oslo 2019)

Slides and notes for the metabarcoding course given at Oslo University
March 25-29. The presented metabarcoding pipeline relies on
[vsearch](https://github.com/torognes/vsearch),
[cutadapt](https://github.com/marcelm/cutadapt/),
[swarm](https://github.com/torognes/swarm) and custom python 2
scripts. The different steps of the pipeline are glued together with
lightweight bash scripting.

To follow the course, please [download this
repository](https://github.com/frederic-mahe/metabarcoding-course-Oslo-2019/archive/master.zip)
and open the file `metabarcoding_course_Oslo_2019.html` in a browser.

Please use the [issue
tracker](https://github.com/frederic-mahe/metabarcoding-course-Oslo-2019/issues)
to report bugs or to ask questions.

## Course's course

(preliminary program)

- primer on input-output, pipes and files,
- basic fastq format checks (vsearch, fastqc),
  - segway to legacy formats (SFF) and new formats (HDF5),
- fastq quality encoding (vsearch),
- fastq merging (vsearch),
- primer trimming (cutadapt or atropos),
- compute expected error (vsearch),
- convert to fasta and dereplicate (vsearch),
- global dereplication (vsearch),
- clustering (swarm),
- OTU refinement (swarm and python),
- chimera detection (vsearch),
- taxonomic assignment (vsearch, against Silva 16S rRNA),
- OTU table (python),
- Statistical ecology (R and PhyloSeq)
