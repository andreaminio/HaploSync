# 5 - Phase haplotypes

During assembly, phase information between haplotypes may be lost while positioning the sequences. HaploMap generates a map of the collinear regions between alternative haplotypes by doing a pairwise comparison of the sequences.

In this workflow, we assume `${HaploSync_path}` is set in the environment and points to the installation directory for HaploSync. To do this, run:

```bash
export HaploSync_path="/path/to/HaploSync"
```

To run HaploMap, run:

```bash
${HaploSync_path}/HaploMap.py -f sequnces.fasta -c Hap1_to_Hap2.pairs.txt
```

where

* `sequnces.fasta`: Genomic sequences in FASTA format.

  * Multiple files are accepted as a comma-separated list.

* `Hap1_to_Hap2.pairs.txt`: Tabular file(s) that relate pairs of sequences in the following format:

  ```text
  Hap1_chr01   Hap2_chr01
  Hap1_chr02   Hap2_chr02
      [...]
  Hap1_chrNN   Hap2_chrNN 
  ```

  * Multiple files are accepted as a comma-separated list.
  * See [HaploMap usage page](../Usage/HaploMap_usage.md) for further details.

#### Optional

* `--tmpdir`: Path to temporary directory. 
  * Default: `./tmp_HaploMap`.
* `-m | --mapper`: Mapping tool to use.
  * Available tools:
    * `minimap` for minimap2.
    * `nucmer` for nucmer+show-coords.
    * `blat` for blat.
  * Default: `minimap`.
* `-t | --threads`:  Number of cores used in mapping process.
  * Default: `4`.
  * Ignored by `blat`, which will run a single core.
* `-o | --out`: Output files name prefix.
  * Default: `out`.

## Output

* `out.paired_regions.txt`: Tabular file reporting the coordinates relating each collinear region of paired sequences. If the sequence pair requested was: 

  ```Hap1_chr01 Hap2_chr01```

  The format looks like:

  * ```
    Hap1_chr01 4450    367901  Hap2_chr01 0       369778  358964  363451
    Hap1_chr01 405665  756715  Hap2_chr01 402962  757735  344425  351050
    Hap1_chr01 757031  791487  Hap2_chr01 758050  792489  34304   34456
    Hap1_chr01 797613  802669  Hap2_chr01 812258  817361  4691    5056
    [...]
    ```

  * and the columns are:

    1. Hap1 ID
    2. Hap1 region start
    3. Hap1 region stop
    4. Hap2 ID
    5. Hap2 region start
    6. Hap2 region stop
    7. Number of matching based
    8. Alignment region length based on Hap1
