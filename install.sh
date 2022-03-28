#!/bin/bash -l

echo "#################################################################################################################"
echo "################## Install virtual envornonment for HaploSync - started @ `date` "
echo "#################################################################################################################"
echo

pip install Cython -qqq 1> /dev/null 2> /dev/null
pip install networkx==2.1 -qqq 1> /dev/null 2> /dev/null
pip install pandas -qqq 1> /dev/null 2> /dev/null
pip install scipy -qqq 1> /dev/null 2> /dev/null
pip install numpy toml datetime multiprocessing pyyaml matplotlib -qqq 1> /dev/null 2> /dev/null
pip install biopython==1.76 -qqq 1> /dev/null 2> /dev/null

conda install -y -q libgit2
conda install -y -q -c conda-forge parallel
conda install -y -q -c bioconda minimap2 samtools bedtools blat
conda install -y -q -c bioconda gmap=2019.09.12

mkdir tmp && cd tmp
        wget https://github.com/mummer4/mummer/releases/download/v4.0.0.beta5/mummer-4.0.0beta5.tar.gz
        tar -xzf mummer-4.0.0beta5.tar.gz
        cd mummer-4.0.0beta5/
        ./configure --prefix=$CONDA_PREFIX
        make install
        cd ..
cd ..
rm -rf tmp

conda install -y -q -c conda-forge r-base=4.0.2
conda install -y -q -c conda-forge r-curl=4.3

Rscript -e '
install.packages( "tidyverse", lib = .libPaths()[2] , repos = "http://cran.us.r-project.org" )
install.packages( "optparse", lib = .libPaths()[2] , repos = "http://cran.us.r-project.org" )
install.packages( "plotly", lib = .libPaths()[2] , repos = "http://cran.us.r-project.org" )
install.packages( "plyr", lib = .libPaths()[2] , repos = "http://cran.us.r-project.org" )
install.packages( "flexdashboard", lib = .libPaths()[2] , repos = "http://cran.us.r-project.org" )
install.packages( "gridExtra", lib = .libPaths()[2] , repos = "http://cran.us.r-project.org" )
install.packages( "devtools", lib = .libPaths()[2] , repos = "http://cran.us.r-project.org" )
require(devtools)
install_version( "knitr", lib = .libPaths()[2] , version = "1.29", repos = "http://cran.us.r-project.org")
install_version( "ggplot2", lib = .libPaths()[2] , version = "3.3.2", repos = "http://cran.us.r-project.org")
install.packages( "svglite", lib = .libPaths()[2] , repos = "http://cran.us.r-project.org" )
install.packages( "gridSVG", lib = .libPaths()[2] , repos = "http://cran.us.r-project.org" )
install.packages( "svgPanZoom", lib = .libPaths()[2] , repos = "http://cran.us.r-project.org" )
install.packages( "ggrepel", lib = .libPaths()[2] , repos = "http://cran.us.r-project.org" )
'

echo
echo "#################################################################################################################"
echo "################## Install virtual envornonment for HaploSync - Completed @ `date` "
echo "#################################################################################################################"