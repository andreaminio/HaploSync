# HaploSync installation

## Prerequisites

### Python libraries

  - Python ver **>2.7**
  - networkx ver **>2.1** (previous versions are incompatible because of different implementation of some functions)
  - Biopython
  - numpy
  - yaml
  - scipy
  - toml
  - datetime
  - pandas
  - multiprocessing

### Tools

The tools below run the pipeline. It is not necessary for these to be in your `$PATH`; see [Setup](#Setup) to set custom paths for these tools.

* Mandatory
  * [Minimap2](https://github.com/lh3/minimap2)
  * [Gmap](http://research-pub.gene.com/gmap/)
    * **Note**: we suggest installing version *2019.09.12*, which is the last version not affected by a bug (@2020.07.16). In version *2020.06.01* and *2020.04.08*, separators between some attributes are missing in the gff3 output. While HaploSync can handle this known issue, the result may be parsed incorrectly if you use the tool independently.
  * [Mummer v.4](https://mummer4.github.io/)
    * **Note 1**: Mummer v3 is not compatible because it does not permit multi-threading. 
    * **Note 2**: Some releases of Mummer v.4 < beta 5 (4.0.0.beta5) should be avoided; a bug allows the tool to report hits between large regions of Ns in the sequences. This does not affect the assembly performance of HaploSync, but dotplots will report unexpected gap-to-gap hits. Installation of the tool through Conda should be avoided. The Bioconda repository only hosts a beta2 version of the tool (@2020.07.16), which is affected by the bug.   
  * [Samtools](https://github.com/samtools)
  * [Bedtools](https://bedtools.readthedocs.io/en/latest/) 
  * [R](https://www.r-project.org/) v4.0 or higher, with the following libraries: 
    * [plotly](https://plotly-r.com/)
    * [optparse](https://github.com/trevorld/r-optparse) 
    * [plyr](https://github.com/hadley/plyr)
    * svgPanZoom
    * knitr
      * **Note**: version 1.29 required 
    * ggplot2
      * **Note**: version 3.3.2 required
* Optional
  * [Blat](http://www.kentinformatics.com/)

## Virtual environment

If the user does not have sudo credentials and/or needs to avoid conflicts with their existing environment, creating a virtual environment in which to use the tool is **strongly recommended**. Here we report how to build a Miniconda environment.

1. Visit the [Miniconda website](https://docs.conda.io/en/latest/miniconda.html) and download the tool appropriate for your system. Miniconda2 is recommended because python 2.7 is necessary.

2. Create a new environment

   ```bash
   bash  Miniconda2-latest-Linux-x86_64.sh
   ```

   * The installation path is `/path/to/virtualenv/miniconda2_HaploSync`.
   * At the end of the installation, Miniconda will ask if the environment should be used as default. This can be done at the discretion of the user.

3. Load the environment

   ```bash
   source /path/to/virtualenv/miniconda2_HaploSync/bin/activate
   ```

4. Set up the environment

   ```bash
   pip install Cython
   pip install networkx==2.1
   pip install biopython==1.76
   pip install pandas
   pip install numpy scipy toml datetime multiprocessing pyyaml matplotlib
   ```

5. Install a local version of the 3rd party tools

   * Use coda to automatise the installation

     ```bash
     conda install libgit2
     conda install -c conda-forge parallel
     conda install -c bioconda minimap2 samtools bedtools blat
     conda install -c bioconda gmap=2019.09.12
     ```

6. Install Mummer v4

   * Download [Mummer version 4.0.0.beta5](https://github.com/mummer4/mummer/releases/tag/v4.0.0.beta5) or later from GitHub and install the tool in the Conda environment 

     ```bash
     tar -xzf ../mummer-4.0.0beta5.tar.gz
     cd mummer-4.0.0beta5/
     ./configure --prefix=$CONDA_PREFIX
     make install
     ```
   
7. Install R and necessary libraries

     ```bash
     conda install -c conda-forge r-base=4.0.2
     conda install -c conda-forge r-curl=4.3
     ```

     * load an R interactive shell checking to be using the conda installation:

       ```bash
       which R
       R
       ```

     * Install necessary R packages:

       ```R
       install.packages("tidyverse")
       install.packages("optparse")
       install.packages("plotly")
       install.packages("plyr")
       install.packages("flexdashboard")
       install.packages("gridExtra")
       install.packages("devtools")
       require(devtools)
       install_version("knitr", version = "1.29", repos = "http://cran.us.r-project.org")
       install_version("ggplot2", version = "3.3.2", repos = "http://cran.us.r-project.org")
       install.packages("svglite")
       install.packages("gridSVG")
       install.packages("svgPanZoom")
       install.packages("ggrepel")
       ```

8. Load the environment (**this should be done before running HaploSync tools**)

   ```bash
   source /path/to/virtualenv/miniconda2_HaploSync/bin/activate
   ```

## Install HaploSync

### From GitHub repo

Install HaploSync by running the following in your preferred working directory:

* Move to the desired directory:

  ```bash
  cd /path/to/installation/directory
  ```

* Download the executables:

  ```bash
  git clone https://github.com/andreaminio/HaploSync
  ```

  This will create a copy of the tool commands.

### From release tarball

Download and unpack the desired release version from [the releases reporsitory](https://github.com/andreaminio/HaploSync/releases).

## Setup custom paths (Optional)

If you want to use compatible tool installations that are not accessible directly in your `$PATH` you can configure the tool for use in custom paths.

If you set up a virtual environment using conda and installed a local version of the tools, you don’t need to customise:

* Generate the file `HaploSync.conf.toml` and then edit the content with desired custom paths:

  ```bash
  cd /path/to/installation/directory/HaploSync
  cp HaploSync.conf.toml_empty HaploSync.conf.toml
  vim HaploSync.conf.toml
  ```

  The original file looks like:

  ```bash 
  title = "HaploBuild configuration file"
  
  [paths]
  minimap2 = ""
  samtools = ""
  bedtools = ""
  nucmer = ""
  show-coords = ""
  gmap = ""
  blat = ""
  ```

* To set custom paths to the tools, **add the absolute path to the folder containing each tool between the quotation marks above**. Leave the space between the quotation marks empty if the tools are in your `$PATH`.

## References

* [Minimap2](https://github.com/lh3/minimap2): Li H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34:3094-3100.
* [Gmap](http://research-pub.gene.com/gmap/): Wu T. D. and Watanabe C. K. (2005). GMAP: a genomic mapping and alignment program for mRNA and EST sequences. *Bioinformatics*, 21:1859-1875
* [Samtools](https://github.com/samtools): Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009). The Sequence alignment/map (SAM) format and SAMtools. *Bioinformatics*, 25:2078-9 
* [Bedtools](https://bedtools.readthedocs.io/en/latest/): Quinlan A. R., Hall I. M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features, *Bioinformatics*, 26:841–842,
* [Plotly](https://plotly-r.com/): Sievert C (2020). *Interactive Web-Based Data Visualization with R, plotly, and shiny*. Chapman and Hall/CRC. ISBN 9781138331457
* [Blat](http://www.kentinformatics.com/): Kent, W.J. (2002). BLAT -- The BLAST-Like Alignment Tool. *Genome Research* 4: 656-664.
* [Mummer v.4](https://mummer4.github.io/): Marçais G., Delcher A.L., Phillippy A.M., Coston R.,Salzberg S.L., Zimin A., (2018). MUMmer4: A fast and versatile genome alignment system. *PLoS computational biology*, 14(1): e1005944.
