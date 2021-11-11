# 1 - Generate pseudomolecules

When starting from a draft assembly, pseudomolecules can be built with HaploSplit. The tool will identify best tiling path to produce a diploid pseudomolecule representation of each chromosome by scaffolding the draft sequences. 

The procedure can be performed 3 different ways depending on the available information used to anchor tho chromosomes:

1. [HaploSplit with genetic map](): The information provided is a map with sorted genetic markers (ex. Linkage map made on SNPs, rhAmpSeq map, genetic maps, etc.). The position of the markers on the draft assembly drives sorting and orientation of input sequences into phased pairs of pseudomolecules.

   The first step HaploSplit does is a QC of the marker hits, exposing errors in input sequences. It is possible to interrupt the HaploSplit procedure after QC. To do this, use `-v | --dry`.

   See the [HaploSplit](../Usage/HaploSplit_usage.md) usage page for a complete explanation of options and the procedure.

2. [HaploSplit with a guide genome](): The information provided to anchor the sequences to chromosomes is a closely related genome assembled into pseudomolecules. In this situation, HaploSplit will use collinear regions between the sequences to sort, orient and phase the draft sequences.
   If the guide genome is used alone, the information available is not sufficient and QC will not be performed. In this scenario a case we recommend using the pseudomolecule QC to detect sequence errors.

   See the [HaploSplit](../Usage/HaploSplit_usage.md) usage page for a complete explanation of options and the procedure.  

3. [Combined genetic map and guide genome approach](): Both a map of genetic markers and a closely related species guide genome are available. In this situation, HaploSplit will first reconstruct the pseudomolecules based on the markers. Then, results will be compared to the guide genome sequence. Collinearity will then be used to fix ambiguous sequences orientation and improve the pseudomolecule content with sequences that were unplaced after the first step.

   HaploSplit does a QC of the marker hits if a genetic map is used. This exposes errors in input sequences.

   See the [HaploSplit](../Usage/HaploSplit_usage.md) usage page for a complete explanation of options and the procedure.  

In this workflow we assume `${HaploSync_path}` is set and points to the installation directory for HaploSync. To do this, run:

```bash
export HaploSync_path="/path/to/HaploSync"
```

## 1.1 - Command lines

### 1.1.1 - HaploSplit with genetic map

```bash
${HaploSync_path}/HaploSplit.py \
	-o ${output_name} \
	-p ${output_prefix} \
	-i ${draft_fasta} \
	-m ${marker_map_tsv} \
	-n ${marker_position_in_draft_bed} \
	-c ${cores} \
	--haplodup \
	--GFF3 ${draft_annotation_gff3} \
	-a ${draft_agp} \
	-k ${draft_known_groups_tsv} \
	--alternative_groups ${draft_alternative_groups_tsv} \
	--legacy_groups ${contigs_legacy_groups_tsv}
```

* `${draft_fasta}`: is the `/path/to/draft_sequences.fasta`. We refer to this dataset as legacy sequences.
* `${draft_annotation_gff3}`: is the `/path/to/draft_sequences.gene_annotation.gff3`.
* `${cores}`: number of cores available for mapping.
* `${output_name}`: prefix for the output file names.
* `${output_prefix}`: prefix for output sequences.
* `${marker_map_tsv}`: Table of markers with sorted position.
   * Tab separated file with 3 columns: 
     1. Chromosome ID,
     2. Sorted marker position
     3. Marker ID.
* `${marker_position_in_draft_bed}`: Map of markers on input sequences. 
   * BED-like tab-separated file with 3 columns
     1. Sequence ID,
     2. Match start (0-based)
     3. Match end (1-based)
     4. Marker ID 
* `${draft_agp}`: AGP file defining the structure of input sequences.
  * ex. If input sequences are scaffolds, `${draft_agp}` reports the scaffold structure in terms of contig composition.
* `${draft_known_groups_tsv}`: Tab separated file of sequences known to be in the same haplotype. 
  - Two columns reporting:
    1. sequence ID
    2. association group ID
* `${draft_alternative_groups_tsv}` : Tab-separated file with groups of sequences that should be in alternative haplotypes of the same chromosome:
  * Two columns containing comma-separated lists.
  * ex. The primary assembly and the list of the haplotigs associated with it.
* `${contigs_legacy_groups_tsv}`: A file of tab-separated values associating groups of legacy input sequences (must match the components ID reported in `${draft_agp}`).
  * Sequence can be associated to multiple groups with multiple rows.
  * Two columns reporting 
    1. component ID 
    2. association group ID
* `${cores}`: number of cores used during mapping procedures
* `--haplodup`: Perform HaploDup QC on the results automatically

See the [HaploSplit](../Usage/HaploSplit_usage.md) usage page for a complete list of options and a description of how it works.

### 1.1.2 - HaploSplit with a guide genome

Default parameters are a good start:

```bash
${HaploSync_path}/HaploSplit.py \
	-o ${output_name} \
	-p ${output_prefix} \
	-i ${draft_fasta} \
	--GFF3 ${draft_annotation_gff3} \
	-a ${draft_agp} \
	-e ${draft_alternative_groups_tsv} \
	-k ${draft_known_groups_tsv} \
	--alternative_groups ${draft_alternative_groups_tsv} \
	-c ${cores} \
	--haplodup \
	--legacy_groups ${contigs_legacy_groups_tsv} \
	-g ${guide_genome_fasta} \
	--local_alignment ${alignement_paf} \
	--align
```

where:

* `${draft_fasta}`: is the `/path/to/draft_sequences.fasta`. From now on we will refer to this dataset as legacy sequences.
* `${draft_annotation_gff3}`: is the `/path/to/draft_sequences.gene_annotation.gff3`.
* `${cores}`: number of cores available for mapping.
* `${output_name}`: name of prefix for the output files.
* `${output_prefix}`: prefix added to the output sequences.
* `${guide_genome_fasta}`: guide genome FASTA ( `/path/to/guide_genome.fasta`).
* It is suggested to use hard-masked sequences. This will remove misleading cross-repeat hits from the alignment results.
* `${alignement_paf}`: alignment file, in PAF format, of the draft sequences mapped on the reference sequences.
* `--align`: Align draft sequences to the guide genome
* `${draft_agp}`: AGP file of input sequence structure composition.
  * ex. If input sequences are scaffolds, `${draft_agp}` reports the scaffold structure in terms of contig composition
* `${draft_known_groups_tsv}`: Tab separated file of sequences known to be in the same haplotype. 
  - Two columns reporting:
    1. sequence ID
    2. association group ID
* `${draft_alternative_groups_tsv}` : Tab separated file with groups of sequences that should be in alternative haplotypes of the same chromosome:
  * Two columns, with comma separated lists
  * ex. Primary and the list of the haplotigs that it is associated with
* `${contigs_legacy_groups_tsv}`: A file of tab-separated values associating groups of legacy input sequences (must match the components ID reported in `${draft_agp}`).
  * Sequence can be associated with multiple groups with multiple rows.
  * Two columns reporting 
    1. component ID 
    2. association group ID
* `${cores}`: number of cores used during mapping procedures
* `--haplodup`: Perform HaploDup QC on the results automatically

HaploSplit will independently run the alignment, create whole genome interactive dotplots, and analyze duplication for the future QC step. Then proceed to [Step 2]().

This is possible only if a genome annotation is available for the draft assembly. It is suggested to have an annotation, even if raw or made with _ab initio_ predictors. This allows a reference-independent QC of possible in-haplotype duplications.

See the [HaploSplit](../Usage/HaploSplit_usage.md) usage page for a complete list of options and a description of how it works.

### 1.1.3 - Combined genetic map and guide genome approach

```bash
${HaploSync_path}/HaploSplit.py \
	-o ${output_name} \
	-p ${output_prefix} \
	-i ${draft_fasta} \
	-m ${marker_map_tsv} \
	-n ${marker_position_in_draft_bed} \
	-g ${guide_genome_fasta} \
	-c ${cores} \
	--haplodup \
	--GFF3 ${draft_annotation_gff3} \
	-a ${draft_agp} \
	-k ${draft_known_groups_tsv} \
	--alternative_groups ${draft_alternative_groups_tsv} \
	--legacy_groups ${contigs_legacy_groups_tsv}
```


* `${draft_fasta}`: is the `/path/to/draft_sequences.fasta`. From now on we will refer to this dataset as legacy sequences.

* `${draft_annotation_gff3}`: is the `/path/to/draft_sequences.gene_annotation.gff3`.

* `${cores}`: number of cores available for mapping.

* `${output_name}`: name of prefix for the output files.

* `${output_prefix}`: prefix added to the output sequences.

* `${marker_map_tsv}`: Table of markers with map sorting position.

   * Tab separated file with 3 columns: 
     1. Chromosome ID,
     2. Marker sorting position
     3. Marker ID.

* `${marker_position_in_draft_bed}`: Map of markers' positions on input sequences. 

   * BED-like tab separated file with 3 columns
     1. Sequence ID,
     2. Match start (0-based)
     3. Match end (1-based)
     4. Marker ID 

* `${guide_genome_fasta}`: guide genome FASTA ( `/path/to/guide_genome.fasta`).

  * It is suggested to use hard-masked sequences. This will remove misleading cross-repeat hits from the alignment results.

* `${draft_agp}`: AGP file of input sequences structure composition.

  * ex. If input sequences are scaffolds, `${draft_agp}` reports the scaffold structure in terms of contig composition.
  
* `${draft_known_groups_tsv}`: Tab separated file of sequences in the same haplotype. 

  - Two columns reporting:
    1. sequence ID
    2. association group ID

* `${draft_alternative_groups_tsv}` : Tab separated file with groups of sequences that should be in alternative haplotypes of the same chromosome:
  * Two columns, with comma separated lists
  * ex. Primary and the list of the haplotigs that it is associated with

* `${contigs_legacy_groups_tsv}`: A file of tab-separated values associating groups of legacy input sequences (must match the components ID reported in `${draft_agp}`). 
  * Sequences can be associated to multiple groups with multiple rows.
  * Two columns reporting 
    1. component ID 
    2. association group ID
  
* `${cores}`: number of cores to be used during mapping.

* `--haplodup`: Perform HaploDup QC on the results automatically

See the [HaploSplit](../Usage/HaploSplit_usage.md) usage page for a complete list of options and a description of how it works.

## 1.2 - Output files

* `${run1_outname}.1.agp`, `${run1_outname}.2.agp`, `${run1_outname}.Un.agp`: AGP file of the relationships between draft sequences and pseudomolecules.

* `${run1_outname}.1.list`, `${run1_outname}.2.list`, `${run1_outname}.Un.list`: sorted and oriented list of draft sequences in the tiling path used to build the pseudomolecules.

  * Tabular format with 2 columns containing:

    1. ID of the tiled reference sequence.
    2. Comma-separated list of oriented draft sequences whose usage may be controlled reported as `seq_id|direction`.

    * Ex.: 

      ```text
      chr01   
      chr02   seq660|+,seq677|-,seq683|-,seq1044|+,seq270|-
      chr03   seq1577|-,seq1572|-,seq2331|+,seq225|+
      chr04   seq135|-,
      ```

* `${run1_outname}.1.fasta`, `${run1_outname}.2.fasta`, `${run1_outname}.Un.fasta`: Pseudomolecules in FASTA format.

* `${run1_outname.annotation.gff3`: Annotation translated onto pseudomolecules sequences.

* `${run1_outname}.dotplot_dir`: pairwise alignments and dot plots. It contains:

  * Hap1.on.Ref , Hap1.on.Hap1, Hap2.on.Ref, Hap2.on.Hap1, Hap2.on.Hap2, Un.on.Ref, Un.on.Hap1 comparisons where Hap1, Hap2 and Un correspond to `${run1_outname}.1.fasta`, `${run1_outname}.2.fasta`, `${run1_outname}.Un.fasta` respectively and Ref to `${reference.fasta}`.
  * For each comparison:
    * forward-to-forward alignments in PAF format with Minimap2 or delta file plus tabular coords for Nucmer.
    * dotplots in PNG, PDF and interactive HTML format.

* `${run1_outname}.dedup_dir`: Results of HaploDup analysis for the newly produced pseudomolecules. It contains:

  * `new.seq.fasta`: pseudomolecules sequences.
  * `gmap_index/`, `gmap_index.log`, `gmap_index.err`: GMAP index folder, indexing procedure standard log and standard error.
  * `gmap_hits_hap1.txt` and `gmap_hits_hap2.txt`: GMAP mapping hits found for input annotation on each haplotype
  * `new.CDS.fasta`: CDS sequences of the annotated genes.
  * `CDS.on.ref.gmap.gff3`: GMAP alignment results in GFF3 format of CDS sequences on pseudomolecules.
  * `diploid_gene_count_trace.hap1.txt`, `diploid_gene_count_trace.hap2.txt`: Tables with hit counts of each gene annotated on the two haplotypes. The columns contain.
    1. Chromosome id
    2. Locus start
    3. Locus end
    4. Locus id
    5. Hap1 counts
    6. Hap2 counts
  * HTML interactive plots separated by chromosome: 
    * `${run1_outprefix}_Hap1_chr01.html`, […], `${run1_outprefix}_Hap1_chrXX.html`
    * `${run1_outprefix}_Hap2_chr01.html`, […], `${run1_outprefix}_Hap2_chrXX.html`
  * RDA files of the R environment used to create the plots.

## 1.3 - Next Step

QC the results with [HaploDup](pseudomolecule_QC.md). This will identify errors or confirm that all issues have been fixed before proceeding with gap closing or phasing.

## Notes

### 1. Selecting a guide genome

When performing HaploSplit using a guide genome, collinear regions between the species are used to bin the draft sequences and distinguish alternative alleles.

See [Minio et al. 2020]() and the [HaploSplit](../Usage/HaploSplit_usage.md) usage page for more details. 

* **The reference must be in pseudomolecules to obtain pseudomolecules**, otherwise it will not be able to deliver chromosome information to assign the sequences.
* **As close as possible evolutionarily to the query genome**. The more distant it is from the subject, the more prone to bias and errors the results will be.
* **The more closely related are the species, the better**. Collinearity between the species is the basic information extracted. Higher genetic distance between the reference and the query means lower collinearity and reduced possibility of reconstructing specific structural differences between the two.

The more contiguous the draft sequences, the more "divergent structures” from the guide genome can be placed correctly in the pseudomolecules. These are more likely to be connected to flanking collinear regions instead of on their own. When run with highly fragmented assemblies, the procedure is more prone to overfit the guide genome structure. Even in presence of translocation, placement will be defined by the location of the collinear region on the guide genome, while the most divergent regions are more likely not to be shared and remain unplaced.

### 2. Mapping tool and mode selection for dotplot generation

HaploSplit allows the user to choose both the mapping tool (either minimap2 or nucmer) and plotting mode (whole genome or one plot for each chromosome). 

While these approaches may seem interchangeable, choosing one strategy versus another may help the QC procedure. Nucmer alignments, though they often require more computational time, are generally more “granular”; they show high numbers of small hits and, when run whole-genome vs whole-genome, are more prone to report cross-pseudomolecules hits. 

While this can help spot centromeres or translocations, the dotplots generated are far more difficult to use. Even if reported separately for each chromosome, HTML files are too difficult to visualise to pinpoint a hit and obtain coordinates. Therefore, it is suggested to run minimap2 as aligner in the first phase. With these results, the whole genome vs. whole genome HTML dotplot is, in general, small enough to be useful to the user.

### 3. Minimap parameters

The pipeline allows the user to pass specific parameters to Minimap2 for tuning the alignment. As default, the suggested parameters are `--cs -x asm20 -r 1000`. These parameters allow aligning different genomes even they diverge substantially. However, it is suggested that the user try different combinations of parameters to find settings best suited for the genomes they are working with.

For further information about minimap2, visit the [minimap2](https://github.com/lh3/minimap2) guide page and its [relative manual](https://lh3.github.io/minimap2/minimap2.html).