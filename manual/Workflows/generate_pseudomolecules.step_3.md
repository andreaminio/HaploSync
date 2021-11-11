[<< Back | Workflows ::: Generate pseudomolecules and QC scaffolding errors](generate_pseudomolecules.md)

# Step 3: Sequence editing

Once breakpoints are defined, sequences incorrectly placed by HaploSplit and eventually, erroneous scaffolds, can be broken accordingly to the identified issues using HaploBreak.

In the workflow we assume to have a variable `${HaploSync_path}` set in the environment pointing to the installation directory of HaploSync. To do so, run:

```bash
export HaploSync_path="/path/to/HaploSync"
```

## HaploBreak

To run the editing do

* Concatenate sequences files, AGP and and annotation files (if necessary):

  ```bash
  cat Pseudomolecules_hap1.fasta Pseudomolecules_hap2.fasta Pseudomolecules_Un.fasta > genome.fasta
  cat Pseudomolecules_hap1.agp Pseudomolecules_hap2.agp Pseudomolecules_Un.agp > previous_to_actual.genome.agp
  cat Pseudomolecules_hap1.gff3 Pseudomolecules_hap2.gff3 Pseudomolecules_Un.gff3 > genome.annotation.gff3
  ```

* Runs HaploBreak:

  ```bash
  ${HaploSync_path}/HaploBreak.py -f genome.fasta \
  	-b Breakpoints.txt \
  	-g genome.annotation.gff3 \
  	-a previous_to_actual.genome.agp
  	-o ${broken_outname} \
  	-p ${new_sequence_prefix}
  ```

where

* `genome.fasta`: contains all the pseudomolecules sequences, eventually assembled with HaploSplit.
* `Breakpoints.txt`: is the tabular file providing the information on the approximate breakpoint position desired.
* `genome.annotation.gff3`: contains the annotation onto the pseudomolecules. It is optional, but if available it is suggested to provide it as :
  * HaploBreak will translate it onto the edited sequences.
  * During the selection of the gaps to break, HaploBreak will prevent the selection of gaps that are bridged by gene loci and that would be eventually disrupted.
* `${broken_outname}`: name prefix for output files.
* `${new_sequence_prefix}`: prefix for ids of the edited sequence. 
* See [HaploBreak usage page](../Usage/HaploBreak_usage.md) for further details.

### What will HaploBreak do?

HaploBreak will search the nearest available junction between contigs. 
1. First it will search for erroneous junctions in the AGP file to break.
   * It assumes, and tries to correct, HaploSplit errors in first instance. 
   * This will appear as sequences (one or more) comprised between two breakpoints corresponding to 2 “gap” components of the AGP file.

     If this is the case, the sequences between breakpoints will be classified as blacklisted. This list can be passed to HaploSplit directly in order to segregate overlapping sequences in different tiling paths, thus to avoiding to perform the same error again (as the tool is not able to do).
   
2. If AGP file search fails, gaps within legacy sequences will be searched. The nearest compatible to the suggested breakpoint position will be used to split the scaffold in sub-sequences.
To further understand how HaploBreak works, follow the link to HaploBreak [usage page](../Usage/HaploBreak_usage.md) or read [Minio et al. (2020)](ref).

### Note: Breakpoint coordinate precision

Breakpoint coordinates given as input are not expected to be precise. HaploBreak will search a suitable junction to associate to a breakpoint with allowing some tolerance. The `-d | --distance` parameter controls this tolerance by setting a threshold (in bp) on the maximum distance between the breakpoint and its associated junction. By default it is set to 25Kb, which should be relaxed enough to allow a comfortable detection of breakpoints for the user, however in some cases where difference between genotypes may be too high it may be not sufficient. It is suggest therefore to increase slightly the allowed tolerance.

However, as this threshold trigger also the search for in-sequence gap to break in presence of assembly error, it is suggested to increase it carefully in order to hinder the possibility to correct errors.

###  Output files

* `${broken_outname}.new_legacy.agp`: AGP file reporting the processed sequences (`out.new_legacy.fasta`) relationship to the legacy sequence structure.
* `${broken_outname}.new_legacy.annotation.gff3`: annotation translated onto the processed sequences (`${broken_outname}.new_legacy.fasta`).
* `${broken_outname}.new_legacy.fasta`: processed sequences.
* `${broken_outname}.blacklist.txt`: sequences in between two AGP-reported gaps used as breakpoints.
* `${broken_outname}.breakable_gaps.list.txt`: list of input in-sequence gaps suitable for breaking.
* `${broken_outname}.broken_gaps.given_sequences.list.txt`: list of breaking points search result on input sequence:

  The format reports:
  * *Header*: “>” + the name of the input sequence.

  * *Regions*: a tabular description of each region for break point search with the first column indicating if the breakpoint fall in a known gap (from the AGP file if given) or within an ungapped region (“sequence”) that requires search for in-sequence gaps. Follow column 2 and 3 with the coordinates of the region.
  
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
* `${broken_outname}.broken_gaps.legacy_sequences.list.txt`: list of the broken gaps in terms of the legacy sequences ids (If AGP is given). The format is the same than `${broken_outname}.broken_gaps.given_sequences.list.txt` however,  only “gap” breakpoints are reported (as the breakpoints derive only from in-sequence gap search).
* `${broken_outname}.broken_loci.annotation.legacy_updated.txt`: Loci that were overlapping AGP-reported gaps and resulted broken by break point selection. 
* `${broken_outname}.legacy_breakable_gaps.list.txt`: list of legacy in-sequence gaps suitable for breaking.

## Next step

Move on to **[Step 4 - Build corrected pseudomolecules](generate_pseudomolecules.step_4.md)** to create new pseudomolecule sequences from the results of the breaking and solve the previous issues.