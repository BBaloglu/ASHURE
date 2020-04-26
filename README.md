# ASHURE
A Sure Heuristic under Random events - Python script for analysis of metabarcoding data with Nanopore sequencing.

ASHURE is a metagenome analysis tool for linear consensus or 1d2 like fastq data from nanopore sequencing devices. A library preparation procedure such as rolling circle amplification (RCA) generates ssDNA containing repeats of their originating template. When read by a nanopore sequencing device, these highly errored but repetative reads can be realigned to generate a more error free consensus.

ASHURE works by leveraging priori information such as a reference database of related sequences, primer information, and amplicon size to find, align, and generate a polished consensus read. Reference sequences do not have to match the amplicon. They just have to be similar enough (70-80 percent identity) for fuzzy kmer based aligners to work. If reference sequences are not available, primer information can be used to generate a pseudo reference database for the search.

Generally, the sensitivity and speed of the pipeline will be dependent on the quality of the reference or pseudo reference used. A large library of reference sequences will result in slow search speed. A reference database of wholly unrelated sequences will result in zero hits. For the best results, you should curate a compact and accurate database of reference sequences relevant to your raw fastq data. The tools and steps for generating a good pseudo reference are shown in the jupyter notebooks.

In metagenome analysis, a major challenge is the presence of novel sequences in a sample. Sequencing errors are hard to distinguish from novel sequences because both are unique. For unknown error profiles such as consensus reads from nanopore sequencing, our best guess is through density based clustering. If a novel sequence exists, you would expect a uniform cloud of sequencing errors around the true sequence.

## Dependencies
Runtime dependencies
[minimap2](https://github.com/lh3/minimap2)
[spoa](https://github.com/lh3/minimap2)
[sklearn](https://scikit-learn.org/stable/install.html)
[hdbscan]()
[pandas]()

Numpy, Scipy, Pandas, and Sklearn can be installed via follow commands for pip
```bash
pip install pandas
pip install scikit-learn
pip install hdbscan
```

The relevant install instructions or binaries for spoa and minimap2 can be found that their respective github pages.

Optional dependences
Wrappers for the following aligners are available in bilge_pype.py as long as they are callable from your current path
[mafft]()
[bwa]()
[bowtie2]()

The following libraries are used in jupyter notebooks for plotting the final data
```bash
pip install bokeh
pip install seaborn
pip install matplotlib
```

## Installation
Ashure is written in python and contains wrapper scripts that call a number of useful bioinformatic tools. Pandas dataframes are primarily used to organize the underlying data. These can be easy manipulated, filtered, and export as csv via jupyter notebook or ipython shell.

## General Usage

General usage of ashure is as following:

```bash
./ashure.py run -fq fastq/*.fastq -p primers.csv -o1 cons.csv                   # runs full pipeline with default parameters
./ashure.py run -fq fastq/*.fastq -p primers.csv -o1 cons.csv -c config.json    # runs full pipeline with custom parameters 

./ashure.py run -fq fastq/*.fastq -p primers.csv -o1 cons.csv -r prfg           # runs only pseudo reference generator
./ashure.py run -fq fastq/*.fastq -db ref_db.fa -o1 cons.csv -r fgs             # runs only fragment search with ref_db.fa sequences

# runs only pseudo reference generator, fragment search, and multi-sequence alignment with default parameters
./ashure.py run -fq fastq/*.fastq -p primers.csv -o1 cons.csv -r prfg,fgs,msa
```

Sub options allow more customization over each subset of the pipeline

```bash
usage: ashure.py prfg [-h] [-fs PRFG_FS] [-c PRFG_CONFIG] [-s CONFIG_FILE] [-r] [-fq FASTQ [FASTQ ...]] [-p PRIMER_FILE] [-o DB_FILE] [--low_mem]

optional arguments:
  -h, --help            show this help message and exit
  -fs PRFG_FS           fastq read size ranges to search for primers
                            -fs 100-200             searches for fragments in reads of length 100-200bp
                            -fs 100-200,500-900     searches for fragments in reads of length 100-200bp and 500-900bp
  -c PRFG_CONFIG        config passed to minimap2
  -s CONFIG_FILE        write settings to configuration file
  -r                    generate pseudo reference database now with the current configuration
  -fq FASTQ [FASTQ ...]
                        fastq reads to search
  -p PRIMER_FILE        csv file containing forward and reverse primers used.
                            This must have at least columns [fwd_id, fwd_seq, rev_id, rev_seq]
  -o DB_FILE            output csv file of pseudo reference sequences
  --low_mem             enable optimizations that reduce RAM used
```

## Library
bilge_pype.py contains useful functions for parsing and calling alignment tools. These functions can be import as a python module.

The following code exerpt enables reading and writing fastq files with bilge_pype

```python
import bilge_pype as bpy

# read a fastq file and output a pandas dataframe
out = bpy.read_fastq('mydata.fq')
# print the data
print(out)
out.to_csv('mydata.csv') # writes fastq dataframe to a csv file

# write data to another fastq file
bpy.write_fastq(out.values, 'another.fq')
```

## Contact information
For additional information, help and bug reports please send an email to: ~humorous_bilge@gmail.com
Follow and tweet to: [@bilgeMolEcol](https://twitter.com/bilgeMolEcol)

## Acknowledgement
Thank you thank you me. I am the best. They dont call me the wise for nothing.

