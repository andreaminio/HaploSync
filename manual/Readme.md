# HaploSync - Manual

HaploSync was developed to improve the quality of diploid genome assemblies. This tool:

* builds diploid pseudomolecule assemblies from draft sequences using a reference genome as a guide,
* enhances assembly contiguity,
* fills sequence gaps by comparing haplotypes independent of a reference,
* generates maps of the phased haplotypes,
* assesses the quality of assembled sequences, and
* edits sequences. 

## Citations
  **Assembly of complete diploid phased chromosomes from draft genome sequences**  
  Andrea Minio, Noé Cochetel, Amanda M Vondras, Mélanie Massonnet, Dario Cantu  
  _G3 Genes|Genomes|Genetics, 2022, jkac143_; doi: https://doi.org/10.1093/g3journal/jkac143

## Tools

* **HaploSplit**

  HaploSplit reconstructs pseudomolecules from a haploid or diploid draft assembly using collinearity with a reference genome as a guide (mandatory). [[Usage](Usage/HaploSplit_usage.md)]

* **HaploDup**

  HaploDup compares the gene content of two haplotypes, delivers interactive graphs, and identifies chimeric scaffolding where alternative haplotypes are included in the same pseudomolecule. A gene annotation is necessary to perform the analysis. [[Usage](Usage/HaploDup_usage.md)]

* **HaploBreak**

  HaploBreak searches for and corrects breaks in chimeric scaffolds. If HaploBreak is used after HaploSplit or if scaffolding information is provided, breaking scaffolds where gaps occur is prioritized. [[Usage](Usage/HaploBreak_usage.md)]

* **HaploFill**

  HaploFill uses unplaced sequences to fill haplotype-specific gaps produced during the assembly or pseudomolecule reconstruction procedures. When gaps in diploid regions cannot be filled with unplaced sequences, the alternative haplotype is used to fill the unassembled region. [[Usage](Usage/HaploFill_usage.md)]

* **HaploMake**

  HaploMake produces novel sequences and translates their annotation. [[Usage](Usage/HaploMake_usage.md)]

* **HaploMap**

  HaploMap does a pairwise comparison between alternative haplotypes and generates a map of the colinear regions. [[Usage](Usage/HaploMap_usage.md)]

## Installation

 For dependencies and installation instructions, see [the dedicated page](Install.md).

## Workflows

 For examples of how to use HaploSync for different tasks, see the [workflow section](Workflows/Readme.md).

## Browser compatibility

Quality control tools rely on interactive plots to identify errors and correct them. These plots are generated using Rmarkdown, knitr, and [Plotly](https://plotly.com/). We discovered variable compatibility with different platforms and web browsers. See the table below for information about compatability. 

* :white_check_mark: Fully compatible.
* :x: Incompatible.
* :o: Limited compatibility, refer to the exception for more information.
* \- : unavailable/untested

|             |       Chrome       |      Chromium      |      Firefox       |       Safari       |   Edge   |  InternetExplorer  |       Opera        |       Brave        |
| ----------: | :----------------: | :----------------: | :----------------: | :----------------: | :------: | :----------------: | :----------------: | :----------------: |
| **Windows** |       :o:  <sup>1</sup>        |         -          | :white_check_mark: |        :x:         | :o:  <sup>2,3</sup> | :white_check_mark: | :white_check_mark: |       :o:  <sup>1</sup>       |
|     **Mac** | :white_check_mark: |         -          | :white_check_mark: | :white_check_mark: |    -     |         -          | :white_check_mark: | :white_check_mark: |
|   **Linux** | :white_check_mark: | :white_check_mark: | :white_check_mark: |         -          |    -     |         -          | :white_check_mark: | :white_check_mark: |

Exceptions:

1. Need to zoom in to show the dotplot lines.
2. Does not show horizontal grid.
3. Runs the plotting scripts slowly.

