# HaploSync installation

We suggest installing the tool within a dedicated virtual environment (even when possessing sudo credentials) to avoid conflicts with the existing global environment. 

```bash
conda create --name haplosync python=2.7
```

## Download HaploSync

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

### From release tarball

Download and unpack the desired release version from [the releases reporsitory](https://github.com/andreaminio/HaploSync/releases).

## Install dependencies
### Using conda

1. Activate the virtual environment:
```bash
conda activate haplosync
```

2. Install R and R packages
```bash
# install r
conda install -c conda-forge r-base

# install r packages
conda install -c conda-forge r-curl r-tidyverse r-optparse r-plotly r-plyr r-flexdashboard r-gridextra r-devtools r-knitr r-svglite r-gridsvg r-svgpanzoom r-ggrepel
```

3. Install python libraries: 
```bash
# get pip
wget https://bootstrap.pypa.io/pip/2.7/get-pip.py
python get-pip.py

# install
pip install Cython networkx==2.1 pandas scipy numpy toml datetime multiprocessing pyyaml matplotlib biopython==1.76
```

3. Install a local version of the 3rd party tools
```bash
conda install libgit2
conda install -c conda-forge parallel
conda install -c bioconda minimap2 samtools bedtools blat gmap=2019.09.12
```

4. Install Mummer v4
Download [Mummer version 4.0.0.beta5](https://github.com/mummer4/mummer/releases/tag/v4.0.0.beta5) or later from GitHub and install the tool in the Conda environment
```bash
wget https://github.com/mummer4/mummer/releases/download/v4.0.0.beta5/mummer-4.0.0beta5.tar.gz
tar -xzf mummer-4.0.0beta5.tar.gz
cd mummer-4.0.0beta5/
./configure --prefix=$CONDA_PREFIX
make install
```

5. Copy the conf.toml, fields can be left empty since tools are already in the path
```bash
cd /path/to/installation/directory/HaploSync
cp HaploSync.conf.toml_empty HaploSync.conf.toml
```

### Manually
To use HaploSync, be sure to have installed the follwing tools, libraries and packages 

#### Python libraries
Make sure to install the following python libraries:
  - Python ver **>2.7**
  - networkx ver **>2.1** (previous versions are incompatible because of different function implementation)
  - Biopython
  - numpy
  - yaml
  - scipy
  - toml
  - datetime
  - pandas
  - multiprocessing

#### Tools

The tools below run the pipeline. It is not necessary for these to be in your `$PATH`; see [Setup](#Setup) to set custom paths for these tools.

* Mandatory
  * [Minimap2](https://github.com/lh3/minimap2)
  * [Gmap](http://research-pub.gene.com/gmap/)
    * **Note**: we suggest installing version *2019.09.12*, which is the last version not affected by a bug (@2020.07.16). In version *2020.06.01* and *2020.04.08*, separators between some attributes are missing in the gff3 output. While HaploSync can handle this known issue, the result may be parsed incorrectly if you use the tool independently.
  * [Mummer v.4](https://mummer4.github.io/) , [Mummer version 4.0.0.beta5 download](https://github.com/mummer4/mummer/releases/tag/v4.0.0.beta5)
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

#### R Packages ####
* In R, install the package using the following command lines:
    ```
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

#### Setup custom paths with custom manual installation of tools

If you want to use compatible tool installations that are not accessible directly in your `$PATH` you can configure the tool for use in custom paths.

If you set up a virtual environment using conda and installed a local version of the tools, you don’t need to customise:

Generate the file `HaploSync.conf.toml` and then edit the content with desired custom paths:

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

To set custom paths to the tools, **add the absolute path to the folder containing each tool between the quotation marks above**. Leave the space between the quotation marks empty if the tools are in your `$PATH`.

## References

* [Minimap2](https://github.com/lh3/minimap2): Li H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34:3094-3100.
* [Gmap](http://research-pub.gene.com/gmap/): Wu T. D. and Watanabe C. K. (2005). GMAP: a genomic mapping and alignment program for mRNA and EST sequences. *Bioinformatics*, 21:1859-1875
* [Samtools](https://github.com/samtools): Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009). The Sequence alignment/map (SAM) format and SAMtools. *Bioinformatics*, 25:2078-9 
* [Bedtools](https://bedtools.readthedocs.io/en/latest/): Quinlan A. R., Hall I. M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features, *Bioinformatics*, 26:841–842,
* [Plotly](https://plotly-r.com/): Sievert C (2020). *Interactive Web-Based Data Visualization with R, plotly, and shiny*. Chapman and Hall/CRC. ISBN 9781138331457
* [Blat](http://www.kentinformatics.com/): Kent, W.J. (2002). BLAT -- The BLAST-Like Alignment Tool. *Genome Research* 4: 656-664.
* [Mummer v.4](https://mummer4.github.io/): Marçais G., Delcher A.L., Phillippy A.M., Coston R.,Salzberg S.L., Zimin A., (2018). MUMmer4: A fast and versatile genome alignment system. *PLoS computational biology*, 14(1): e1005944.
