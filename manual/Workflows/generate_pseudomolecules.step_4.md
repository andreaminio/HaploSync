[<< Back | Workflows ::: Generate pseudomolecules and QC scaffolding errors](generate_pseudomolecules.md)

# Step 4: Rebuild pseudomolecules

Once sequences have been edited or classified for blacklisting, pseudomolecules can be rebuilt with HaploSplit. To do so:

1. Separate the blacklisting information delivered by HaploBreak (`${broken_outname}.blacklist.txt`) by haplotype, obtaining two lists `B1.list` and `B2.list` for haplotype 1 and haplotype 2, respectively. 

   If there have been previous runs of breaking and sequences, some legacy sequences may have been already blacklisted for the chromosomes. It may be necessary to merge the two list, at least for haplotype 1, in order to avoid the recurring placement of erroneous sequences.

   File format:

   ```
   chr01   
   chr02   seq660|+,seq677|-,seq683|-,seq1044|+,seq270|-
   chr03   seq1577|-,seq1572|-,seq2331|+,seq225|+
   chr04   seq135|-,
   ```

   * Tabular format with 2 columns containing:
     1. id of the reference sequence whose tiling must be controlled.
     2. comma-separated list of oriented draft sequences, reported as `seq_id|direction`.

2. (optional) move to a new, clean folder to avoid confusion.

3. run HaploSplit using the edited sequences and the blacklists.

In the workflow we assume to have a variable `${HaploSync_path}` set in the environment pointing to the installation directory fo HaploSync. To do so, run:

```bash
export HaploSync_path="/path/to/HaploSync"
```

## HaploSplit

HaploSplit will be used to build new pseudomolecules from the edited sequences and using the the blacklists to control their usage in the tiling path generation.

Run:

```bash
${HaploSync_path}/HaploSplit.py -m -t ' --cs -x asm20 -r 1000' \
	-c ${cores} \
	-f -a \
	-o ${HaploSplit_run2_filename} \
	-p ${HaploSplit_run2_sequence_prefix} \
	-i ${alignment_run2_paf} \
	--dedup --dotplot minimap-wholegenome \
	--GFF3 ${broken_outname}.new_legacy.annotation.gff3 \
	-q ${broken_outname}.new_legacy.fasta \
	-r ${reference_fasta} \
	--B1 B1.list --B2 B2.list
	[ --R1 1st.txt --R2 2nd.txt --distance1 1000000 --distance2 1000000 ]
```
where:
* `${reference_fasta}`: is the `/path/to/reference.fasta`.
  * It is suggested the use reference sequence where repeats have been hard-masked. This will clean the alignment results of any misleading cross-repeat hit.
* `${broken_outname}.new_legacy.fasta`: Query sequences as edited by HaploBreak.
* `${broken_outname}.new_legacy.annotation.gff3`: Gene annotation of query sequences as edited by HaploBreak.
* `${HaploSplit_run2_filename}`: prefix for output files.
* `${HaploSplit_run2_sequence_prefix}`: prefix for assembled sequence.
* `${alignment_run2_paf}`: Minimap2 PAF alignment file of draft sequences on reference.
  * This file can be provided as input or, like in the example above, created by this tool: 
    * If the file provided exists, it will be used to produce the tiling path.
    * If mapping is requested (`-m|--map` flag), HaploSplit will perform the mapping and save the results in `file.paf`, either creating or overwriting it if already existing.
* `${cores}`: number of cores to be used in mapping procedures (Minimap2 and, eventually, GMAP).
* `--B1 B1.list` and `--B2 B2.list`: list of the sequences to avoid for a given tiling path, obtained from splitting  `${broken_outname}.blacklist.txt` by haplotype.
* Optional tuning options:
  * ` --R1 1st.txt` and `--R2 2nd.txt`: Require to use the following list of (stranded) contigs in the first haplotype and the second haplotype respectively (Good list sequence).
  * ` --min1 N` and `--min2 N`: Minimum length of sequence allowed to be a requirement for the first haplotype and the second haplotype respectively (default: 0). Useful when good lists are taken directly from.
* Other flags:
  * `-f`: Export paths sequences in 3 fasta files (`${HaploSplit_run2_sequence_prefix}.1.fasta`, `${HaploSplit_run2_sequence_prefix}.2.fasta`, `${HaploSplit_run2_sequence_prefix}.Un.fasta`).
  * `-a`: Export paths in sequence structure in  3 AGP files (`${HaploSplit_run2_sequence_prefix}.1.agp`, `${HaploSplit_run2_sequence_prefix}.2.agp`, `${HaploSplit_run2_sequence_prefix}.Un.agp`) . See the [AGP file format specifications](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/) for more specifications.
  * `--dedup` : Run HaploDup analisys.
  * `--dotplot minimap-wholegenome` : Generate whole genome comparison dotplots using minimal as aligner.
* See [HaploSplit usage page](../Usage/HaploSplit_usage.md) for further details on HaploSplit command line and functionality.

### Output files

Output files are the same than on [Step 1.1](generate_pseudomolecules.step_1.md#11---haplosplit), except their base name. New tiling path have been built correcting errors  the previous pseudomolecules evidenced. However, smaller errors may be still present that need fixing. 



## Next step

Go back to **[Step 2 - Quality control](generate_pseudomolecules.step_2.md)** and check if there is still errors to resolve.