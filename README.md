# HaploSync
Tools for haplotype-wise reconstruction of pseudomolecules

## Aim
The aim of the tool is improving continuity of a given genome draft assembly up to pseudomolecule representation using a closely relate species as guide to sort and direct the query sequences.  

## Description
### General
Despite the recent advance in sequencing technologies allowing less fragmented assemblies, achieving a chromosome scale representation is still challenging. 
By using a closely related species as a guide, the tool aims at use the conserved genomic structures to sort and direct the given sequences, then a representation of the pseudomolecules is given as a non-redundant tiling path of the query assembly     
### Diploid representation
When dealing with heterozygous genomes, picturing the full genomic variability in a diploid representation adds a second layer of complexity. A way to deal with the issue is to haploidize the assembly, but, while complexity is actually reduced, this leads to loss of information on the "alternative" state. 

Moreover, some tools do not deliver phasing information, it is lost during the assembly pipeline or innacurate because of haplotype switches. Even more difficult is the situation of high heterozygosity rates that do trick the assembly tools in overestimating the the "primary" assembly or, better, underestimating the genomic variability, ending in imperfect primary/alternative sequence paring information.

The tool takes care of segregating sequences that do actually represent alternatives to the very same genomic locus. A reconstruction of pseudomolecules is then performed caring of separate usage of the alternative sequences in different tiling paths  

## Prerequisites
- Python libraries  
  - Python ver **>2.7**
  - networkx ver **>2.1** (previous version are incompatible due to missing functions)
  - Biopython
- External tools
  - [Minimap2](https://github.com/lh3/minimap2)

## Output
### Assembly data
- `file.paf`: Query to reference alignment file. Reusable
- Tiling path #1: 
  - `out.1.agp` : AGP file describing output sequences in function of input query used to assemble them  
  - `out.1.fasta` : Output sequences in FASTA format
  - `out.1.list` : Sorted and directed list of input query used 
- Tiling path #2 (if `--N2` is not set)
  - `out.2.agp` : AGP file describing output sequences in function of input query used to assemble them
  - `out.2.fasta` : Output sequences in FASTA format
  - `out.2.list` : Sorted and directed list of input query used
- Unplaced sequences
  - `out.Un.agp` : AGP file describing unplaced sequences in function of unused input query in any tiling path
  - `out.Un.fasta` : Unplaced sequences in FASTA format
  - `out.Un.list` : List of unused input query sequences

List files, like the ones delivered as output, can be used to control the reconstruction procedure through `-1`, `-2`, `--B1`, `--B2`, `--R1`, and `--R2` parameters  
   

### Temporary files
- Intermediate mapping files:
  - tmp.map1.paf.for
  - tmp.map1.paf.multimapping
  - tmp.map1.paf.multimapping.uniq
  - tmp.map1.paf.rev


## Usage
### Command line
```
HaploSplit.py [-h] [-i file.paf [Required]]
                     [-r reference.fasta [Required]]
                     [-q query.fasta [Required]] [-v] [-o NAME] [-p PREFIX]
                     [-a] [-f] [-g N] [-c N] [--distance1 N] [--distance2 N]
                     [-m] [-u] [-t "minimap2 --cs -t 4 -x asm20 -r 1000"]
                     [-n N] [-1 1st.txt] [-2 2nd.txt] [--R1 1st.txt]
                     [--min1 N] [--R2 2nd.txt] [--min2 N] [--B1 2nd.txt]
                     [--B2 2nd.txt] [--N2]
```
### Arguments:
```
  -h, --help            show this help message and exit
  -i file.paf [Required], --input file.paf [Required]
                        File of hits. Used as input if already mapped, output
                        result of mapping procedure if -m/--map
  -r reference.fasta [Required], --reference reference.fasta [Required]
                        FASTA file with reference sequence
  -q query.fasta [Required], --query query.fasta [Required]
                        FASTA file with query sequence
  -v, --dry             Dry run: find best tiling paths and excluded nodes but
                        no FASTA or AGP exported
  -o NAME, --out NAME   Output NAME prefix for files [default: out]
  -p PREFIX, --prefix PREFIX
                        Output NAME prefix for files [default: NEW]
  -a, --agp             Export paths in AGP format to NAME.{1,2,Un}.agp
  -f, --fasta           Export paths sequences in NAME.{1,2,Un}.fasta
  -g N, --gapsize N     Minimum gap size placeholder for AGP and FASTA (in bp)
                        [default: 1,000bp]
  -c N, --concatenate N
                        Set the gap size used to concatenate unplaced
                        sequences (in bp) [default: do not concatenated]
  --distance1 N         Set the maximum distance (bp) between two hits to be
                        considered adjacent [default: 2,000,000bp]
  --distance2 N         Set the maximum distance (bp) between two hits to be
                        considered adjacent at second round [default:
                        4,000,000bp]
  -m, --map             Do run mapping [overwrite input file.paf if existing]
  -u, --uniq            Map only
  -t "minimap2 --cs -t 4 -x asm20 -r 1000", --tool "minimap2 --cs -t 4 -x asm20 -r 1000"
                        Mapping command for minimap2
  -n N, --hitgap N      Allowed gap between hits to be merged [default:
                        500000]
  -1 1st.txt, --path1 1st.txt
                        Use ONLY the following list of (stranded) contigs to
                        create the first haplotype. Tab separated file with
                        Target_ID followed by comma separated list of contigs
                        (contig1,contig2,[...],contigN)
  -2 2nd.txt, --path2 2nd.txt
                        Use ONLY the following list of (stranded) contigs to
                        create the second haplotype. Tab separated file with
                        Target_ID followed by comma separated list of contigs
                        (contig1,contig2,[...],contigN)
  --R1 1st.txt          Require to use the following list of (stranded)
                        contigs in the first haplotype. Tab separated file
                        with Target_ID followed by comma separated list of
                        contigs (contig1,contig2,[...],contigN)
  --min1 N              Minimum length of sequence allowed to be a requirement
                        for the first haplotype (default: 0)
  --R2 2nd.txt          Require to use the following list of (stranded)
                        contigs in the second haplotype. Tab separated file
                        with Target_ID followed by comma separated list of
                        contigs (contig1,contig2,[...],contigN)
  --min2 N              Minimum length of sequence allowed to be a requirement
                        for the second haplotype (default: 0)
  --B1 2nd.txt          Blacklisted (stranded) contigs NOT to be used in the
                        first haplotype. Tab separated file with Target_ID
                        followed by comma separated list of contigs
                        (contig1,contig2,[...],contigN)
  --B2 2nd.txt          Blacklisted (stranded) contigs NOT to be used in the
                        second haplotype. Tab separated file with Target_ID
                        followed by comma separated list of contigs
                        (contig1,contig2,[...],contigN)
  --N2                  Don't run the search for the 2nd path
```
