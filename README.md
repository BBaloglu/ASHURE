# Getting started
```bash
./ashure.py -h                                                 # prints help
./ashure.py run -h                                             # prints help on a submodule
./ashure.py run -fq folder/*.fq -p primers.csv -o cons.csv -r  # runs everything with default parameters
./ashure.py run -c config.json                                 # writes parameters to config.json
./ashure.py run -fq folder/*.fq -c config.json -r              # loads parameters from config.json and run everything
./ashure.py prfg -fq folder/*.fq -p primers.csv -o seq.csv -r           # runs prfg module
./ashure.py clst -i seq.csv -o cl.csv -iter 10 -r                       # runs clst module for 10 iterations
./ashure.py run -fq folder/*.fq -db cl.csv -o consensus.csv -r          # runs everything with cl.csv as reference
./ashure.py run -fq folder/*.fq -db cl.csv -o cons.csv -r fgs,msa,cons  # runs fgs, msa, and cons modules only 
```

# ASHURE
A Sure Heuristic under Random events - A Python based toolkit for analysis of environmental DNA from nanopore sequencing.

ASHURE is designed for analysis of linear consensus or 1d2 like fastq data from nanopore sequencing devices. Library preparation procedure such as rolling circle amplification (RCA) generate ssDNA containing repeats of their originating template. When these concatemers are read by a nanopore sequencing device, the repetative regions can be realigned to generate a more error free consensus.

ASHURE works by leveraging priori information about the amplicon such as primers used, reference database of related sequences, and amplicon size to find, align, and generate a polished consensus read. Reference sequences do not have to match the amplicon, just similar enough (70-80 percent identity) for fuzzy kmer based alignment.

The sensitivity and speed of the pipeline is dependent on the quality of the reference or pseudo reference used. Searching large libraries will take time and compute effort is sometimes wasted on comparing against unrelated or redundant sequences. A bad reference sequences will result in poor alignments that contaminate multi-alignment. The best results are obtain by searching  against a compact but accurate database of reference sequences relevant to your raw fastq data. The tools and steps for generating good pseudo reference sequences are shown in [pseudo_references.ipynb](https://github.com/bbaloglu/ashure/demos/pseudo_references.ipynb).

Sequence clusters are generated in the last step of the pipeline using a density based clustering approach called OPTICS. Traditional OTU thresholding approaches do not work with nanopore data because the error profile of each read is unpredictable. Density based clustering is more suited for these situations because OTU boundaries are adaptively called based on the divergence in local sequence identity. This method requires sufficient coverage around a true amplicon for clustering to work. An interative demo of how this approach works can be found in [clustering.ipynb](https://github.com/bbaloglu/ashure/demos/clustering.ipynb).

Read polishing with medaka, nanopolish, or racon is supported by our pipeline. Our toolkit only prepares input data for these tools. A demo of this process is shown in [polishing](https://github.com/bbaloglu/ashure/polishing.ipynb).

## Dependencies
### Runtime dependencies
[minimap2](https://github.com/lh3/minimap2) is used for sequence alignment

[spoa](https://github.com/lh3/minimap2) is used for multi-sequence alignment

Numpy, Scipy, Pandas, and Sklearn should also be installed via the following commands for pip
```bash
pip install pandas          # for organizing underlying data
pip install scikit-learn    # for clustering
pip install hdbscan         # for clustering
```

For non-python libraries such as spoa and minimap2, subprocess.run() is used by wrapper functions to make calls to these tools. These binaries should be accessible from your local path. Install instructions for spoa and minimap2 can be found on their github pages. If you download or compile these binaries without installing them to /usr/bin/, you must add them to your local path by running the following commands.
```bash
mv minimap2 /home/username/.local/bin                       # adds minimap2 to your local binary path
export PATH=$PATH':/home/username/.local/bin'               # makes executables in ~/.local/bin accessible in your shell
echo PATH=$PATH':/home/username/.local/bin' >> .bash_login  # updates these settings everytime you login
```

### Optional dependences
Wrappers for the following aligners can also be called from `bilge_pype.py` These tools are not used by `ashure.py`, but if you want to use them in your python scripts you must also make their binaries callable from your local path 

[mafft](https://mafft.cbrc.jp/alignment/software/source.html) for progressive multi-sequence alignment

[bwa](https://github.com/lh3/bwa) and [bowtie2](https://github.com/BenLangmead/bowtie) for aligning miseq fastq

The following libraries are used in demo notebooks for data visualization
```bash
pip install bokeh          # for making interative plots
pip install seaborn        # for making pretty matplotlib plots
```

## Installation
Ashure is written in python and is compatible with python3.6 and above. The following commands will install `ashure.py` to your local path.
```bash
git clone https://github.com/bbaloglu/ashure  # clones this repository
cd ashure                                     # enter the repository folder
chmod +x ashure.py                            # make it executable
./ashure.py run -h                            # look at the help commands
mv ashure.py ~/.local/bin/ashure              # adds ashure to local path
mv bilge_pype.py ~/.local/bin/                # adds bilge_pype module to local path with ashure
ashure -h                                     # call ashure from local path
```

`ashure.py` imports many functions from `bilge_pype.py` These should be in same folder together for the code to work. I will bilge_pype to pypi later if time permits.

If the dependencies are installed, the code should work as is. Pandas dataframes are primarily used to organize the underlying data. These can be easy manipulated, filtered, and export as csv via jupyter notebook or ipython shell.

## Optimization
For some parts of the code, run speed of vector operations can be accelerated if openblas or mkl are installed.
[benchmarks](https://markus-beuckelmann.de/blog/boosting-numpy-blas.html) for mkl, atlas, and blas libraries show the perform gain in certain numpy operations when openblas, mkl, or atlas is installed.

[how to install](https://stackoverflow.com/questions/29979539/how-can-i-make-numpy-use-openblas-in-ubuntu#42647590) shows a quick guide to enabling openblas on ubuntu. Google is your best friend here.

## General Usage
The whole pipeline or subsets of the pipeline can be run in sequence with the `run` module
```bash
./ashure.py run -fq fastq/*.fastq -p primers.csv -o1 cons.csv                   # runs full pipeline with default parameters
./ashure.py run -fq fastq/*.fastq -p primers.csv -o1 cons.csv -c config.json    # runs full pipeline with custom parameters 
./ashure.py run -fq fastq/*.fastq -p primers.csv -o1 cons.csv -r prfg           # runs only pseudo reference generator
./ashure.py run -fq fastq/*.fastq -db ref_db.fa -o1 cons.csv -r fgs             # runs only fragment search with ref_db.fa sequences
# runs only pseudo reference generator, fragment search, and multi-sequence alignment with default parameters
./ashure.py run -fq fastq/*.fastq -p primers.csv -o1 cons.csv -r prfg,fgs,msa
```

The list of submodules can be found with
```bash
./ashure.py -h
```

Example usage from the pseudo reference generator module
```bash
./ashure.py prfg -h                                                         # prints help
./ashure.py prfg -fq folder/*.fq -p primers.csv -o seq.csv -r               # runs the module
./ashure.py prfg -fq folder/*.fq -p primers.csv -o seq.csv -fs 500-3000 -r  # runs the module with fastq filter for 500-3000bp
./ashure.py prfg -fq folder/*.fq -p primers.csv -o seq.csv -fs 500-3000 -c config.json  # updates config.json with custom parameters
```

Example usage from the clustering module
```bash
./ashure.py clst -h                                           # prints help
./ashure.py clst -i input.csv -o clusters.csv -r              # runs clustering
./ashure.py clst -i input.csv -o clusters.csv -c config.json  # updates config.json with custom parameters
```

## Library
`bilge_pype.py` contains several functions you may find useful for parsing and calling commonly used alignment tools. See the [demos](https://github.com/bbaloglu/ashure/demos) folder for how some of these functions are used.

## Citing ASHURE
If you use ASHURE in your work, please cite:

    Baloğlu, B. et al. (2020). “A Workflow for Accurate Metabarcoding Using Nanopore MinION Sequencing.” bioRxiv: 2020.05.21.108852. http://biorxiv.org/content/early/2020/05/25/2020.05.21.108852

## Contact information

Contact via email: bilgenurb@gmail.com or twitter: [@bilgeMolEcol](https://twitter.com/bilgeMolEcol)

## Acknowledgement
This study was supported by funding through the Canada First Research Excellence Fund. This work represents a contribution to the University of Guelph Food From Thought research program.  

