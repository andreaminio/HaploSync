[<< Back | Workflows ::: Fill gaps using diploid structure](fill_and_make.md)

# 2 - Build sequences

After finding the best fillers for the gaps, filling sequences and annotations must be integrated into the pseudomolecules. To do so we will use HaploMake.

In the workflow we assume to have a variable `${HaploSync_path}` set in the environment pointing to the installation directory fo HaploSync. To do so, run:

```bash
export HaploSync_path="/path/to/HaploSync"
```

## HaploMake

1. The `BLOCK` file describing the new sequence structures is `genome.block`. It should look like this:

   ```
   >new_seq1
   old_seq1	1499999	2000000	.	.	-
   old_seq2	0	1000000	.	.	+
   ```

2. Make the new files:

   ```bash
   ${HaploSync_path}/HaploMake.py -f sequences.fasta -g annotation.gff3 -s genome.block --gap 1000 --format BLOCK
   ```

   with

   * `-f sequences.fasta`: Genomic sequences. It contains the `old_seq*`  sequences in FASTA format.
   * `-g annotation.gff3`: Genomic annotation. It contains the loci annotation over `old_seq*` sequences in GFF3 format.
   * `--gap 1000`: add a gap between two consecutive sequences of 1000bp.
   * `-u` : add unplaced sequences to the output.
   * See [HaploMake usage page](../Usage/HaploMake_usage.md).

### Output

* `out.fasta`: Selected region sequence in FASTA format with the new given name (`New_seq`).
* `out.gff3`: Annotation of the selected region with coordinates translated accordingly to the new sequence (`New_seq` based).
* `out.agp`: AGP file associating the selected region to the output.
  * This file is identical to `conversion.agp` unless `--reverse` is selected.

## Next step

1. A QC of the results is suggested. [HaploDup procedure](generate_pseudomolecules.step_2.md) may be useful to evidence errors that my need correction.
2. [Workflows ::: Phase haplotypes](map_and_phase.md) can be performed to obtain a phasing information between the two haplotypes produced.