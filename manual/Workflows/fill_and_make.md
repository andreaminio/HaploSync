# 4 - Fill gaps using diploid structure

Gaps present in the pseudomolecules produced by HaploSplit or other scaffolding tools can have two different origins:

* Placeholders between legacy sequences: HaploSplit adds gaps where unknown content is expected between adjacent sequences. These gaps are a user-defined fixed length.

* Gaps due to previous assembly procedures: Gaps in the legacy sequences used as input for HaploSplit are inherited in final pseudomolecules unless the gaps are broken by HaploBreak. Gap lengths may vary based on the approach used for scaffolding.

HaploFill searches for sequences to fill gaps with and considers the particular characteristics of the gap and both haplotypes. For more information about how different types of gaps are treated, see [Minio et al. 2020]() and [HaploFill usage page](../Usage/HaploFill_usage.md).

The pipeline has 2 steps:

1. **[Find the filling](fill_and_make.step_1.md)**: HaploFill uses haplotypes' sequences, raw reads coverage, and repeat information to search for the best sequence to fill gaps.
2. **[Build sequences](edit_sequences.5.md)**: Fillers are compared and integrated into the pseudomolecules sequences.

## 4.1Find the filling

The first step of the pipeline is searching for appropriate fillers for the gaps present in the input sequences.

In this workflow, we assume `${HaploSync_path}` is set in the environment and points to the installation directory fo HaploSync. To do this, run:

```bash
export HaploSync_path="/path/to/HaploSync"
```

## HaploFill

HaploFill requires the following data to perform the search of proper gap filling:

* Pseudomolecule sequences (both haplotypes) and unplaced sequences. No hard masking.

* Sequencing data coverage of all sequences. We will run the whole pipeline starting from sequences reads. However, read alignments or a per-base coverage profile can be used as input.

* Repetitive regions

  To run HaploFill, use the following command:

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

* `pseudomolecules.1_to_pseudomolecules.2.correspondance.txt`: Tabular file reporting the relationship between the sequences in haplotype 1 (1st column) and haplotype 2 (2nd column).
* `repeats.gff3`: repeat annotation and its format (`--repeats_format GFF3`).
  * Allowed formats: `GFF3`, `BED`, `GFF`, and `GTF`.
* `reads.fasta.gz`: sequenced reads and the sequencing technology (`--sequencing_technology PacBio`).
  * File format: gzipped FASTA or FASTQ
  * sequencing technologies: `PacBio`for PacBio subreads, `ONT` for Oxford NanoPore reads, `Illumina_pe` and  `Illumina_se` for illumina short, paired-end or single-end read sequencing protocols.
* Optional:
  * `--b1 coverage.pseudomolecules.1.bam` and `--b2 coverage.pseudomolecules.2.bam`: output files for the read alignment over Haplotype 1 and Haplotype 2.
  * `--map_threads ${core}`: Number of cores for mapping procedures.
* See [HaploFill usage page](../Usage/HaploFill_usage.md#haplofill) for for a full description of the tool and its options.

### Running the tool

The tool may require some time to finish running. Some of the processes are computationally intensive, like the read mapping, the coverage extraction, the policy classification, and generating support sequences. 

The procedure can be run step-by-step (see [usage page](../Usage/HaploFill_usage.md#haplofill)) and permits automatic recovery in the case of failure. This helps avoid needing to rerun successfully completed steps.

### Output

* Standard output: Will log the status of the procedure with time information.
* Standard error: Will detail the processes running and summarise some of the results obtained for QC and debug use.
* `tmp_HaploFill`: This folder will contain all the temporary data of the analysis, including the files necessary for recovering the procedure and controlling the steps to rerun (see [usage page](../Usage/HaploFill_usage.md#haplofill) for more details).
* `out.gap_filling_findings.txt`: Summary of the results obtained for each gap in the sequence.
* `out.structure.block`: [BLOCK](../block_format.md) file format of the expected pseudomolecules structure after gap filling.

## 4.2 - Check the homozygous fillers

Before proceeding with the sequence reconstruction, you should evaluate homozygous filler. 

HaploFill will automatically detect and use regions of the alternative haplotype as filler when the region is diploid based on read coverage but no unplaced sequence has been retrieved to fill the gap. This occurs when a diploid-aware assembler (ex. FalconUnzip) doesn't produce two identical sequences for the two haplotypes in a homozygous region. 

This behaviour may not be desired by the user. Homozygous filler can be removed from the block file by deleting the corresponding block lines in a text editor.

## Next step: Build the sequences

After finding filler, `out.structure.block` with the new expected structure must be passed to [HaploMake](edit_sequences.5.md) to integrate filler sequences and annotations into the pseudomolecules.
