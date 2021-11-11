[<< Back | Workflows ::: Manipulate genomic sequences](edit_sequences.md)

# 3.4 - Convert annotation between versions 

HaploMake can bidirectionally convert annotation coordinates and/or re-build updated sequences using coordinates in an AGP file that relate two genome versions.

## Example

To produce the sequence `new_seq*` from the old sequences `old_seq*` and translate the gene annotation with them:
```bash
export HaploSync_path="/path/to/HaploSync"
${HaploSync_path}/HaploMake.py -f sequences.fasta -g annotation.gff3 -s conversion.agp --format AGP --skipoverlap
```
* `-f sequences.fasta`: Genomic sequences. It contains the `old_seq*`  sequences in FASTA format.
* `-g annotation.gff3`: Genomic annotation. It contains the loci annotation for `old_seq*` in GFF3 format.
* `-s conversion.agp `: AGP file of coordinates that relate different sequences [AGP file format specifications](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/).
  
  * `conversion.agp` resembles:
    
    ```
    new_seq1	1	500000	1	W	old_seq1	1500000	2000000	-
    new_seq1	500001	501000	2	N	1000	scaffold	yes	align_genus
    new_seq1	501001	1501000	3	W	old_seq2	1	1000000 +
    ```
  * The file describes a new, 1,501 Kb-long sequence named `new_seq1`:
    * The first region, from 1 to 500,000 bp, corresponds to the reverse complement of `old_seq1` region at 1.5 Mb to 2 Mbp.
    * There is a 1 Kb gap.
    * A third region, from 501,001 bp to 1,501,000 bp, corresponds to `old_seq2` from its beginning to 1 Mbp.
  
* `--reverse` inverts the direction of the conversion, from `new_seq*` based sequences and annotations to `old_seq*` based. `sequences.fasta` contains `new_seq*` sequences in FASTA format and `annotation.gff3` contains the loci annotation for `new_seq*` coordinates.

See [HaploMake usage page](../Usage/HaploMake_usage.md) for a complete description of its usage.

### Output

* `out.fasta` : Sequence in the selected region in FASTA format with the new given name (`New_seq`).
* `out.gff3` : Annotation of the selected region with coordinates translated for the new sequence (`New_seq` based).
* `out.agp` : AGP file associating the selected region to the output. See the [AGP file format specifications](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/). This file is identical to `conversion.agp` unless `--reverse` is selected.

