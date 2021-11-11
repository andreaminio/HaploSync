[<< Back | Workflows ::: Fill gaps using diploid structure](fill_and_make.md)

# 1 - Find the fillings

The first step of the pipeline is searching for appropriate fillers for the gaps present in the input sequences.

In the workflow we assume to have a variable `${HaploSync_path}` set in the environment pointing to the installation directory fo HaploSync. To do so, run:

```bash
export HaploSync_path="/path/to/HaploSync"
```

## HaploFill

HaploFill requires the following data to perform the search of proper gap fillers:

* Pseudomolecule sequences (both haplotypes) and unplaced sequences. No hard masking.

* Sequencing data coverage of all sequences. Here we assume to run the whole pipeline starting from sequences reads, however reads alignment or a per-base coverage profile can be used as input.

* Repetitive regions.

  To run HaploFill use the following command:

```bash
${HaploSync_path}/HaploFill.py \
	-1 pseudomolecules.1.fasta \
	-2 pseudomolecules.2.fasta \
	-U unplaced.fasta \
	-c pseudomolecules.1_to_pseudomolecules.2.correspondance.txt \
	-r repeats.gff3 \
	--repeats_format GFF3 \
	-C \
	-s reads.fasta.gz \
	--sequencing_technology PacBio \
	[ --b1 coverage.run2.1.bam --b2 coverage.run2.2.bam --map_threads ${core} ]
```

where:

* `pseudomolecules.1.fasta`,  `pseudomolecules.2.fasta`, and `unplaced.fasta`: FASTA files with haplotype 1, haplotype 2 and unplaced sequences.

* `pseudomolecules.1_to_pseudomolecules.2.correspondance.txt`: Tabular file reporting the correspondences between the sequences in haplotype 1 (1st column) and the haplotype 2 (2nd column).
* `repeats.gff3`: file of repeats annotation and its format (`--repeats_format GFF3`).
  * Allowed formats: `GFF3`, `BED`, `GFF`, and `GTF`.
* `reads.fasta.gz`: sequenced reads and the sequencing technology (`--sequencing_technology PacBio`).
  * File format: FASTA or FASTQ also compressed in Gzip format.
  * sequencing technologies: `PacBio`for PacBio subreads, `ONT` for Oxford NanoPore reads, `Illumina_pe` and  `Illumina_se` for illumina short reads pair-end or single-end sequencing protocols.
* Optional:
  * `--b1 coverage.pseudomolecules.1.bam` and `--b2 coverage.pseudomolecules.2.bam`: output files for the reads alignment over Haplotype 1 sequences and Haplotype 2.
  * `--map_threads ${core}`: Number of cores for mapping procedures.
* See [HaploFill usage page](../Usage/HaploFill_usage.md#haplofill) for further options and description of the performed procedure.

### Running the tool

The tool may require some time to complete. Some of the process are computationally intensive, like the reads mapping, the coverage extraction, the policy classification or the generation of support sequences. 

However the procedure can be controlled on a by step bases (see [usage page](../Usage/HaploFill_usage.md#haplofill)) and cares of automatic recovery in case of failure avoiding to rerun steps already successfully completed.

### Output

* Standard output: Will log the status of the procedure with time information.
* Standard error: Will detail the processes running and summarise some of the results obtained for QC and debug use.
* `tmp_HaploFill`: This folder will contain all the temporary data of the analysis, comprising the files necessary for recovering the procedure or control the steps to rerun (see [usage page](../Usage/HaploFill_usage.md#haplofill) for more details).
* `out.gap_filling_findings.txt`: Summary of the results obtained for each gap in the sequence.
* `out.structure.block`: [BLOCK](../block_format.md) file format of the expected pseudomolecules structure after gap filling.

## Next step

After finding the fillers, `out.structure.block` file with the new expected structure must be passed to HaploMaker in [Step 2](fill_and_make.step_2.md) to integrate fillers sequences and annotations into the pseudomolecules.