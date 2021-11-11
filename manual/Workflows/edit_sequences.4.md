[<< Back | Workflows ::: Manipulate genomic sequences](edit_sequences.md)

# 3.1 - Break scaffolds with errors

When errors in the sequences are found during QC of pseudomolecules, it is necessary to edit the sequences by breaking erroneous scaffolds. This can be done manually by clearing _ad hoc_ AGP files or automatically with HaploBreak.

The manual procedure can be performed with [HaploMake](), which will build the new sequences and translate the known feature coordinates on them. See [3.4 - Convert annotation between versions](edit_sequences.2.md) for a complete description of the procedure.

To perform the procedure automatically, it is sufficient to define approximate breakpoints on the pseudomolecules. HaploBreak will take care of breaking the most appropriate scaffolding junction, produce new sequences, and translate coordinates on them.

In the workflow, `${HaploSync_path}` is set to the installation directory of HaploSync. Run:

```bash
export HaploSync_path="/path/to/HaploSync"
```

## HaploBreak

* Concatenate sequences files, AGP and annotation files (if necessary):

  ```bash
  cat Pseudomolecules_hap1.fasta Pseudomolecules_hap2.fasta Pseudomolecules_Un.fasta > genome.fasta
  cat Pseudomolecules_hap1.agp Pseudomolecules_hap2.agp Pseudomolecules_Un.agp > previous_to_actual.genome.agp
  cat Pseudomolecules_hap1.gff3 Pseudomolecules_hap2.gff3 Pseudomolecules_Un.gff3 > genome.annotation.gff3
  ```

* Run HaploBreak:

  ```bash
  ${HaploSync_path}/HaploBreak.py -f genome.fasta \
  	-b Breakpoints.txt \
  	-g genome.annotation.gff3 \
  	-a previous_to_actual.genome.agp
  	-o ${broken_outname} \
  	-p ${new_sequence_prefix}
  ```

where

* `genome.fasta`: contains all pseudomolecules sequences, eventually assembled with HaploSplit.
* `Breakpoints.txt`: is the tabular file providing information about the approximate breakpoint position desired.
* `genome.annotation.gff3`: contains pseudomolecule annotation. It is optional, but recommended if available because:
  * HaploBreak will translate it onto the edited sequences and
  * during the selection of the gaps to break, HaploBreak will prevent the selection of gaps bridged by gene loci.
* `${broken_outname}`: name prefix for output files.
* `${new_sequence_prefix}`: prefix for edited sequence IDs. 
* See [HaploBreak usage page](../Usage/HaploBreak_usage.md) for further details.

### What will HaploBreak do?

HaploBreak will search for the nearest available junction between contigs. 

1. First it will search for erroneous junctions in the AGP file to break.

   * It assumes there are and tries to correct HaploSplit errors.

   * This will appear as one or more sequences between two breakpoints corresponding to 2 “gap” components in the AGP file.

     If true, the sequences between breakpoints will be blacklisted. This list can be passed to HaploSplit directly to segregate overlapping sequences in different tiling paths and avoid performing the same error again.

2. If the AGP file search fails, gaps within legacy sequences will be searched. The gap nearest to the suggested breakpoint position will be used to split the scaffold in sub-sequences. To further understand how HaploBreak works, follow the link to HaploBreak [usage page](../Usage/HaploBreak_usage.md) or read [Minio et al. (2020)](ref).

### Note: Breakpoint coordinate precision

Breakpoint coordinates given as input are not expected to be precise. HaploBreak will search for a suitable junction to associate to a breakpoint. The `-d | --distance` parameter sets the maximum possible distance between the breakpoint and its associated junction. By default, it is set to 25Kb; this should be relaxed enough to allow a comfortable detection of breakpoints for the user. However, this may be insufficient if the difference between genotypes is too high. In this case, increase `-d | --distance`. Do so cautiously, because `-d | --distance` also triggers the search for an in-sequence gap to break if there is an assembly error.

###  Output files

* `${broken_outname}.new_legacy.agp`: AGP file reporting the relationship between the processed sequences (`out.new_legacy.fasta`) and legacy sequence structure.

* `${broken_outname}.new_legacy.annotation.gff3`: annotation translated onto the processed sequences (`${broken_outname}.new_legacy.fasta`).

* `${broken_outname}.new_legacy.fasta`: processed sequences.

* `${broken_outname}.blacklist.txt`: sequences in between two AGP-reported gaps used as breakpoints.

* `${broken_outname}.breakable_gaps.list.txt`: list of input in-sequence gaps suitable for breaking.

* `${broken_outname}.broken_gaps.given_sequences.list.txt`: list of breaking point search results on input sequence:

  The format reports:

  * *Header*: “>” + the name of the input sequence.

  * *Regions*: a tabular description of each region for break point search with the first column indicating if the breakpoint falls in a known gap (from the AGP file if given) or within an ungapped region (“sequence”) that requires search for in-sequence gaps. Follow column 2 and 3 with the coordinates of the region.

    Ex: 

    ```text
    >chr01
    sequence        12150000        12857650
    sequence        14110000        20989345
    gap     12957818        12958818
    sequence        14330000        20989345
    >chr02
    sequence        13580000        18271468
    gap     18271468        18272468
    ```

* `${broken_outname}.broken_gaps.legacy_sequences.list.txt`: list of the broken gaps in terms of the legacy sequences IDs (if AGP is given). The format is the same as `${broken_outname}.broken_gaps.given_sequences.list.txt`. However, only “gap” breakpoints are reported because the breakpoints derive only from the in-sequence gap search.

* `${broken_outname}.broken_loci.annotation.legacy_updated.txt`: Loci that overlapped AGP-reported gaps and broken by break point selection. 

* `${broken_outname}.legacy_breakable_gaps.list.txt`: list of legacy in-sequence gaps suitable for breaking.

## Next step

**[Rebuild the pseudomolecules](generate_pseudomolecules.md)** from the corrected scaffolds.
