# Getting started
```bash
echo write me write me write me write me write me i am a lazy bastard
```

# ASHURE
A Sure Heuristic under Random events - A Python based toolkit for analysis of metagenomic data from nanopore sequencing.

ASHURE is designed for analysis of linear consensus or 1d2 like fastq data from nanopore sequencing devices. A library preparation procedure such as rolling circle amplification (RCA) generates ssDNA containing repeats of their originating template. When read by a nanopore sequencing device, these highly errored but repetative reads can be realigned to generate a more error free consensus.

ASHURE works by leveraging priori information about the amplicon such as primers used, reference database of related sequences, and amplicon size to find, align, and generate a polished consensus read. Reference sequences do not have to match the amplicon. They just have to be similar enough (70-80 percent identity) for fuzzy kmer based aligners to work.

Generally, the sensitivity and speed of the pipeline will be dependent on the quality of the reference or pseudo reference used. A large library of reference sequences will result in slow search speed. A reference database of wholly unrelated sequences will result in zero hits. For the best results, you should curate a compact but accurate database of reference sequences relevant to your raw fastq data. The tools and steps for generating a good pseudo reference sequences are shown in [pseudo_references.ipynb](https://github.com/bbaloglu/ashure/demos/pseudo_references.ipynb).

For metagenome analysis, de novo sequence clusters are generated in the last step of the pipeline using an adaptive density based clustering approach called OPTICS. Traditional OTU thresholding approaches do not work with nanopore data because error profile is highly variable and subject to sequence identity of the amplicons. Density based clustering is advantageous because it can adaptively call OTU boundaries based on the divergence in local sequence identity. The error profile of the sequencing device does not need to be fully understood for this clustering approach to work. You just need enough read coverage around a true amplicon for you to detect novel clusters. A demo of this approach can be found in [clustering.ipynb](https://github.com/bbaloglu/ashure/demos/clustering.ipynb)

Finally, read polishing can be done with medaka, nanopolish, or racon pipelines. We do not implement these in our code as installation of some of these packages is rather complicated for the end user. We already get good results via clustering find these further steps unneccessary. However, if you wish, our toolkit can prepare the relevant input data for these polishing tools. A demo of this process with racon is shown in [polishing](https://github.com/bbaloglu/ashure/polishing.ipynb).

## Dependencies
Runtime dependencies
[minimap2](https://github.com/lh3/minimap2)
[spoa](https://github.com/lh3/minimap2)
[sklearn](https://github.com/scikit-learn/scikit-learn)
[hdbscan](https://github.com/scikit-learn-contrib/hdbscan)
[pandas](https://github.com/pandas-dev/pandas)

Numpy, Scipy, Pandas, and Sklearn can be installed via follow commands for pip

```bash
pip install pandas
pip install scikit-learn
pip install hdbscan
```

The install instructions for spoa and minimap2 can be found that their respective github pages. If you download or compile these binaries without installing them to /usr/bin/, then do the following to add the binaries to your local path.

In Mac OSX

```bash
mv minimap2 /Users/username/.local/bin
export PATH=$PATH':/Users/username/.local/bin'
```

In linux
```bash
mv minimap2 /home/username/.local/bin
export PATH=$PATH':/home/username/.local/bin'
```

In windows
```bash
echo haha fuck me. I use microsoft.
```

Optional dependences
Wrappers for the following aligners are available in bilge_pype.py as long as they are callable from your current path
[mafft](https://mafft.cbrc.jp/alignment/software/source.html)
[bwa](https://github.com/lh3/bwa)
[bowtie2](https://github.com/BenLangmead/bowtie)

The following libraries are used in demo notebooks for data visualization
```bash
pip install bokeh
pip install seaborn
pip install matplotlib
```

## Installation
Ashure is written in python and is compatible with python3.6 and above. To install, simply clone repository with the following commands.

```bash
git clone https://github.com/bbaloglu/ashure
cd ashure
chmod +x ashure.py      # make it executable
ashure.py run -h        # look at the help commands
```

As long the relevant python libraries are installed, the code should work as is. Pandas dataframes are primarily used to organize the underlying data. These can be easy manipulated, filtered, and export as csv via jupyter notebook or ipython shell.

For some parts of the code, run speed of vector operations can be significantly accelerated if openblas or mkl are installed.
[benchmarks](https://markus-beuckelmann.de/blog/boosting-numpy-blas.html)
[openblas in ubuntu](https://stackoverflow.com/questions/29979539/how-can-i-make-numpy-use-openblas-in-ubuntu#42647590)

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

