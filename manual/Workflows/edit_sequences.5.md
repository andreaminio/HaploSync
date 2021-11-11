[<< Back | Workflows ::: Manipulate genomic sequences](edit_sequences.md)

# 3.2 - Generate new sequences after gap filling with HaploFill

After finding the best fillers for gaps, filling sequences and annotations must be integrated into the pseudomolecules using HaploMake.

In the workflow, we assume `${HaploSync_path}` is set in the environment and points to the installation directory for HaploSync. To do this, run:

```bash
export HaploSync_path="/path/to/HaploSync"
```

## HaploMake procedure

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

   * `-f sequences.fasta`: Genomic sequences. It contains the `old_seq*` sequences in FASTA format.
   * `-g annotation.gff3`: Genomic annotation. It contains the annotation for `old_seq*` sequences in GFF3 format.
   * `--gap 1000`: add a gap between two consecutive sequences of 1000bp.
   * `-u` : add unplaced sequences to the output.
   * [*Optional*] `—skipoverlap`: This flag will prevent HaploMake from reducing overlapping regions between juxtaposed sequences
   * See [HaploMake usage page](../Usage/HaploMake_usage.md).

### Output

* `out.fasta`: Selected region sequence in FASTA format with the new given name (`New_seq`).
* `out.gff3`: Annotation of the selected region with coordinates translated to the new sequence (`New_seq` based).
* `out.agp`: AGP file associating the selected region to the output.
  * This file is identical to `conversion.agp` unless `--reverse` is selected.

## Next step

1. A QC of the results is suggested. [HaploDup procedure](pseudomolecule_QC.md) may be useful to find errors that may need correction.
2. [Phase haplotypes](map_and_phase.md) can be performed to obtain phasing information between the two haplotypes.
