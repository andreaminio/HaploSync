[<< Back | Workflows ::: Manipulate genomic sequences](edit_sequences.md)

# 3.5 - Upgrade sequences or generate ad hoc loci from assemblies 

HaploMake can build new sequence structures from user-created structure files. This is useful when:

* extracting complex loci from the genome for _ad hoc_ analysis or representation.

* manually editing, correcting or scaffolding genomic sequences.

  

## Example
If we want to produce a new sequence, `new_seq1`, from two old sequences, `old_seq1` and `old_seq2`, and :

* the first part of `new_seq1` corresponds to the reverse complement of `old_seq1` region from 1.5 Mb to 2 Mbp.
* the second part of `new_seq1`  corresponds to `old_seq2`, from the beginning of `old_seq2` to 1 Mbp.
* and there is a 1 Kb gap between the two old sequences.

The procedure differs based on the file format used for the conversion file (`BED`, `AGP` or `BLOCK`).

### 1. With a BED file

1. Generate a `conversion.bed` file:

   ```
   old_seq1	1499999	2000000	.	.	-
   old_seq2	0	1000000	.	.   +   
   ```

   

2. Make the new files:
   ```bash
   export HaploSync_path="/path/to/HaploSync"
   ${HaploSync_path}/HaploMake.py -f sequences.fasta -g annotation.gff3 -s conversion.bed --gap 1000 --format BED -p new_seq
   ```
   * `-f sequences.fasta`: Genomic sequences that contain the `old_seq*` sequences in FASTA format.
   * `-g annotation.gff3`: Genomic annotation that contains the annotation for `old_seq*` sequences in GFF3 format.
   * `--gap 1000`: Add a gap between two consecutive sequences of 1000 bp.
   * `-p new_seq`: use `new_seq` as prefix for the newly generated sequence name.
   * [*Optional*] `—skipoverlap`: Do not check for overlap between adjacent parts of `old_seq`.
   * See [HaploMake usage page](../Usage/HaploMake_usage.md).

### 2. With an AGP file

1. Generate a `conversion.agp` file. See [AGP file format specifications](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/).
   ```
    new_seq1	1	500000	1	W	old_seq1	1500000	2000000	-
    new_seq1	500001	501000	2	N	1000	scaffold	yes	align_genus
    new_seq1	501001	1501000	3	W	old_seq2	1	1000000 +   
   ```

2. Make the new files:
   ```bash
   export HaploSync_path="/path/to/HaploSync"
   ${HaploSync_path}/HaploMake.py -f sequences.fasta -g annotation.gff3 -s conversion.agp --format AGP 
   ```
   * `-f sequences.fasta`: Genomic sequences that contain the `old_seq*` sequences in FASTA format.
   * `-g annotation.gff3`: Genomic annotation that contains the loci annotation for `old_seq*` sequences in GFF3 format.
   * [Optional] `—skipoverlap`: Do not check for overlap between adjacent parts of `old_seq`.
   * See [HaploMake usage page](../Usage/HaploMake_usage.md).

### 3. With a BLOCK file

1. Generate a `conversion.block` file. See [this page for further informations](../block_format.md). 

   ```
   >new_seq1
   old_seq1	1499999	2000000	-
   old_seq2	0	1000000	+
   ```

2. Make the new files:

   ```bash
   export HaploSync_path="/path/to/HaploSync"
   ${HaploSync_path}/HaploMake.py -f sequences.fasta -g annotation.gff3 -s conversion.block --gap 1000 --format BLOCK
   ```
   * `-f sequences.fasta`: Genomic sequences that contain the `old_seq*` sequences in FASTA format.
   * `-g annotation.gff3`: Genomic annotation that contains the annotation for `old_seq*` sequences in GFF3 format.
   * `--gap 1000`: Add a gap between two consecutive sequences of 1000 bp.
   * [Optional] `—skipoverlap`: Do not check for overlap between adjacent parts of `old_seq`.
   * See [HaploMake usage page](../Usage/HaploMake_usage.md).

### Output

* `out.fasta`: FASTA of the selected region with the new given name, `New_seq`.
* `out.gff3`: Annotation of the selected region with coordinates translated for the new sequence.
* `out.agp`: AGP file relating the old and new regions. This file is identical to `conversion.agp` unless `--reverse` is selected.
