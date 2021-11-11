[<< Back | Workflows ::: Generate pseudomolecules and QC scaffolding errors](generate_pseudomolecules.md)

# Step 1: Find best tiling path, run duplication search and produce dotplots

According to the assembly input data, the procedure may require or not the creation of pseudomolecules from a draft assembly. 

* If starting from draft assembly follow the procedure in [Step 1.1 - Draft Assembly](generate_pseudomolecules.step_1.md#11---draft-assembly) section to create the pseudomolecules first.
* If pseudomolecules are already available, follow the steps in [Step 1.2](generate_pseudomolecules.step_1.md#12---pseudomolecules) section that will skip HaploSplit and produce the interactive plots for the QC procedure.

In the workflow we assume to have a variable `${HaploSync_path}` set in the environment pointing to the installation directory for HaploSync. To do so, run:

```bash
export HaploSync_path="/path/to/HaploSync"
```

## 1.1 - Draft Assembly

When starting from a draft assembly, pseudomolecules can be built with HaploSplit. The tool will also take care of producing all the data necessary for the QC of the sequences.

### HaploSplit

Run HaploSplit to identify best tiling paths of the reference genome and scaffold the draft sequences in pseudomolecules. In the first iteration we aim to build a first, coarse tiling path of the reference. Default parameters are a good start:

```bash
${HaploSync_path}/HaploSplit.py \
	-m \
	-t ' --cs -x asm20 -r 1000 ' \
	-c ${cores} \
	-f -a \
	-o ${run1_outname} \
	-p ${run1_outprefix} \
	-i ${alignement_paf} \
	--dedup \
	--dotplot minimap-wholegenome \
	--GFF3 ${annotation_gff3} \
	-q ${draft_fasta} \
	-r ${reference_fasta}
```
where:
* `${reference_fasta}`: is the `/path/to/reference.fasta`.

  It is suggested to use reference sequence where repeats have been hard-masked. This will clean the alignment results of any misleading cross-repeat hit.

* `${draft_fasta}`: is the `/path/to/draft_sequences.fasta`. From now on we will refer to this dataset as legacy sequences.

* `${annotation_gff3}`: is the `/path/to/draft_sequences.gene_annotation.gff3`.

* `${cores}`: number of cores available for mapping.

* `${run1_outname}`: name prefix for the output files.

* `${run1_outprefix}`: prefix to be added to the output sequences.

* `${alignement_paf}`: draft to reference sequences alignment file in PAF format.

In this way, HaploSplit is ran to perform independently the alignment procedure, the creation of whole genome interactive dotplots and the duplication analysis for the future QC step. If this is the case, you can go directly to [Step 2]().

This is possible only if an annotation of the genome is available for the draft assembly. It is suggested to have an annotation, even if raw or made just with _ab initio_ predictors, as this allows a reference-independent QC of the possible in-haplotype duplications.

See [HaploSplit usage page](../Usage/HaploSplit_usage.md) for further insight.

#### Reference genome

Building pseudomolecules from a draft assembly requires an external source of information. HaploSplit uses for this purpose closely related species, taking advantage of the collinear regions between the species to bin the draft sequences in the chromosomes of origin and, at the same time, distinguish between alternative alleles (see [Minio et al. 2020]() and [HaploSplit](../Usage/HaploSplit_usage.md) usage page for more details).

The reference genome sequence will not be used during the procedure, however:

* **The reference must be in pseudomolecules to obtain pseudomolecules**.
* **The more closely related are the species, the better**, as collinearity between the species is the basic information extracted. Higher the genetic distance between the reference and the query the lower is the collinearity and less probable the possibility to reconstruct specific structural differences between the two (As those region are more likely not to be included in the pseudomolecules).

  This situation double links with the fragmentation of the draft assembly. The more contiguous are the draft sequences the more "divergent structures” can be placed correctly in the pseudomolecules as they are more likely connected to flanking collinear regions instead than by their own. On the contrary, when run with the more fragmented assemblies, the procedure is more prone to overfit the reference structure and classify as Unplaced the "divergent structures”.

#### Output files

* `${run1_outname}.1.agp`, `${run1_outname}.2.agp`, `${run1_outname}.Un.agp`: AGP file of the draft sequences to pseudomolecules relationships.

* `${run1_outname}.1.list`, `${run1_outname}.2.list`, `${run1_outname}.Un.list`: sorted and directed list of draft sequences in the tiling path used to build the pseudomolecules.
  
  * Tabular format with 2 columns containing:
    1. id of the reference sequence tiled.
    2. coma separated listing of oriented draft sequences whose usage may be controlled reported as `seq_id|direction`.
    * Ex.: 
      ```text
      chr01   
      chr02   seq660|+,seq677|-,seq683|-,seq1044|+,seq270|-
      chr03   seq1577|-,seq1572|-,seq2331|+,seq225|+
      chr04   seq135|-,
      ```
  
* `${run1_outname}.1.fasta`, `${run1_outname}.2.fasta`, `${run1_outname}.Un.fasta`: Pseudomolecules sequences in FASTA format.

* `${run1_outname.annotation.gff3`: Annotation translated onto pseudomolecules sequences.

* `${run1_outname}.dotplot_dir`: pairwise alignments and dot plots. It contains:
  
  * Hap1.on.Ref , Hap1.on.Hap1, Hap2.on.Ref, Hap2.on.Hap1, Hap2.on.Hap2, Un.on.Ref, Un.on.Hap1 comparisons where Hap1, Hap2 and Un correspond to `${run1_outname}.1.fasta`, `${run1_outname}.2.fasta`, `${run1_outname}.Un.fasta` respectively and Ref to `${reference.fasta}`.
  * For each comparison:
    * forward-to-forward alignments in PAF format with Minimap2 or delta file plus tabular coords for Nucmer.
    * dotplots in PNG, PDF and interactive HTML format.
  
* `gmap_hits_hap1.txt` and `gmap_hits_hap2.txt`: GMAP mapping hits found for input annotation on each haplotype

* `${run1_outname}.dedup_dir`: Results of HaploDup analysis for the newly produced pseudomolecules. It contains:
  * `new.seq.fasta`: pseudomolecules sequences.
  * `gmap_index/`, `gmap_index.log`, `gmap_index.err`: GMAP index folder, indexing procedure standard log and standard error.
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
  
  ### Next Step
  
  Run **[Step 2 - Quality control](generate_pseudomolecules.step_2.md)**

## 1.2 - Pseudomolecules

As pseudomolecules are already available, we need to perform the duplication analysis (HaploDup) and generate the dotplots that are necessary for the QC of the sequences.

### HaploDup

If duplication analysis must be run separately, it can be run like this:
```bash
${HaploSync_path}/HaploDup.py -f ${pseudomolecules_fasta} -g ${pseudomolecules_annotation_gff3} -i ${Hap1_to_Hap2} -o ${run1_outname} -t ${cores} 
```
where
* `${pseudomolecules_fasta}`: is the `/path/to/pseudomolecules_sequences.fasta` file. 

  * Alternatively, if haplotype sequences are in multiple files (like in the output of HaploSync where the sequences are in `${run1_outname}.1.fasta` and  `${run1_outname}.1.fasta` ), the files can be given as a coma separated list.

    Ex. `-f ${run1_outname}.1.fasta,${run1_outname}.2.fasta`.

* `${pseudomolecules_annotation_gff3}`: is the `/path/to/pseudomolecules_sequences.gene_annotation.gff3` file.
  
  * Alternatively, if annotations are in multiple files, these can be fed as a coma separated list of with multiple `-d|--gff` flags.
  * An annotation is necessary to perform the analysis of gene content duplication. To obtain one it is possible to:
    * Perform a complete annotation pipeline. [(ex. the pipeline that was developed for M. rotunidifolia genome annotation)](https://github.com/andreaminio/AnnotationPipeline-EVM_based-DClab)
    * Run an ab initio identification of the gene content (ex. with Augustus, GeneMark, SNAP, etc.)
    * Map known CDSs sequences of closely related species with GMAP and generate a non-redundant annotation.
  
*  `${Hap1_to_Hap2}`: is the path to one or more tabular files associating each sequence of Haplotype 1 to its own corresponding sequence in Haplotype 2. Multiple files are accepted as coma separated list or by using multiple `-i|—pair` flags.
  
  * Ex. :
    ```
    Hap1_chr01   Hap2_chr01
    Hap1_chr02   Hap2_chr02
        [...]
    Hap1_chrNN   Hap2_chrNN 
    ```
  
* `cores`: number of cores available for mapping procedure.

* `run1_outname`: name prefix for the output files.

* See [HaploDup usage page](../Usage/HaploDup_usage.md) for further details.

#####  Output files

* `gmap_hits_hap1.txt` and `gmap_hits_hap2.txt`: GMAP mapping hits found for input annotation on each haplotype.
* `${run1_outname}.dedup_dir`: Results of HaploDup analysis for the newly produced pseudomolecules. It contains:
  * `new.seq.fasta`: pseudomolecules sequences.
  * `gmap_index/`, `gmap_index.log`, `gmap_index.err`: GMAP index folder, indexing procedure standard log and standard error.
  * `new.CDS.fasta`: CDS sequences of the annotated genes.
  * `CDS.on.ref.gmap.gff3`: GMAP alignment results in GFF3 format of CDS sequences on pseudomolecules.
  * `diploid_gene_count_trace.hap1.txt`, `diploid_gene_count_trace.hap2.txt`: Tables with hit counts of each gene annotated on the two haplotypes. The columns contain:
    1. Chromosome id
    2. Locus start
    3. Locus end
    4. Locus id
    5. Hap1 counts
    6. Hap2 counts
  * HTML interactive plots separated by input sequence: 
    * `seq1.html`, […], `seqXX.html`

### Dotplots

If duplication analysis has be run separately from HaploSplit, it is necessary also to run separately a query-to-reference sequence comparison and generate a dotplot. As we’ll see in the next section, the dotplot and the gene duplication plot  provide  complementary (and easy to spot) information to identify assembly errors. Moreover, an interactive dotplot reporting the coordinates is highly recommended as it makes easier to identify precise position of breakpoints in QC phase.

To do so, it is suggest to use the following script provided with HaploSync, and the adaptation of DotPlotly \([ref](ref.com)\) code to produce an interactive HTML file where, hovering the cursor over the plot, the overlaying label reports start and stop coordinates of the hits in both reference and query axis:

```bash
bash ${HaploSync_path}/support_scripts/doDotPlotly.forward_only.no_chr_names.sh ${reference_fasta} ${pseudomolecules_fasta} "/path/to/minimap2" [$hit_size] [ref_seq1,...,ref_seqN] 
```

where:

* `reference_fasta`: is the `/path/to/reference.fasta`, the reference genome sequences used to compare the query genome structure to.
* `${pseudomolecules_fasta}` : is the `/path/to/pseudomolecules.fasta`, the sequences of the query genome in pseudomolecules for which we are running the QC.
* `"/path/to/minimap2”`: is the absolute path to minimap2 executable.
* Optional: 
  * `match_size`: minimum alignment hit size in bp. If not provide, 7,500bp is used.
  * sorted list of reference sequences ids. If not provided, reference ids are sorted by name.

##### Output files

* files are named `${pseudomolecules_sequences}.on.${reference}`.
* forward-to-forward alignments in PAF format.
* dotplots in PNG, PDF and interactive HTML format.
* RDA files of the R environment used to create the plots.

### Next Step

Run **[Step 2 - Quality control](generate_pseudomolecules.step_2.md)**

## Notes

### 1. Mapping tool and mode selection for dotplot generation

HaploSplit allows to select both the mapper tool (either minimap2 or nucmer) and the plot mode (whole genome or one plot for each chromosome). 

While these approaches may seem interchangeable, there are some subtle difference in the results that may help the QC procedure. Nucmer alignments, except for requiring more computational time, are generally more “granular”, i.e smaller hits in a higher number and, when run whole-genome vs  whole-genome, more prone to report cross-pseudomolecules hits. 

While this result may have its use for spotting centromeres or translocations, the dotplots generated are far more difficult to use. Even if reported separately for each chromosome, HTML files are too complex to visualise up to making it almost impossible to point correctly the hit and obtaining the coordinates. Therefore, in the first phase it is suggested to run minimap2 as aligner. With these results, the whole genome vs. whole genome HTML dotplot is, in general, small enough to allow a good user interaction.

### 2. Minimap parameters

The pipeline allows the user to pass specific parameters to Minimap2 for tuning the alignment. As default the suggested parameters `--cs -x asm20 -r 1000`. These paramours allow to align different genomes even when they show substantial sequence divergence. However, it is suggested to the user to try different parameters combination to find the ones that do best fit the genomes under examination. 

For further information, we suggest to visit [minimap2](https://github.com/lh3/minimap2) page and consult the [relative manual](https://lh3.github.io/minimap2/minimap2.html).

