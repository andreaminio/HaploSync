[<< Back | Workflows ::: Manipulate genomic sequences](edit_sequences.md)

# 3.3 - Select regions of interest

HaploMake can be used to extract specific parts of a genomic sequence and its annotation by providing a single region. 

## Example

To extract the reverse complement of a region on chromosome 4 (chr04), from 1,500,000~2,000,000 bp, run:

```bash
export HaploSync_path="/path/to/HaploSync"
${HaploSync_path}/HaploMake.py -f genome.fasta -g annotation.gff3 -s region.file --format ${region.file_format} --skipoverlap
```

with 

* `-f genome.fasta`: Genomic sequences.
* `-g annotation.gff3`: Genomic annotation.
* `-p | —prefix`: Control the output sequence ID.
* `-s region.file ` and `--format ${region.file_format}`: The region to extract and the format of input. The structure can be in any of the following formats: 

    * `BED`:
        ```
        chr04	1499999	2000000	.	.	-
        ```
      * Use `bed6` format to include the orientation of the region to extract. Otherwise, it will be assumed to be the same as the reference sequences (`+`).
      * Different regions can be extracted at the same time by submitting a comma-separated list of `bed` files.
      * Starting coordinate is **0-based**.
      * [BED file format specifications](https://m.ensembl.org/info/website/upload/bed.html).
      
    * `AGP`:
      
        ```
        New_seq	1	500000	1	W	chr04	1500000	2000000	-
      ```
      * The file names a new sequence,`New_seq`, that is 500Kb long, and corresponds to the region to extract.
      * Output sequence ID can be controlled from the AGP file and/or by using the `-p | —prefix` flag.
      * Different regions can be extracted by including them on different lines and different output sequence IDs (like `New_seq1`, `New_seq2`, etc.) and/or submitting a comma-separated list of files.
      * Starting coordinate is **1-based**.
    * [AGP file format](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/).
      
    * `BLOCK`:
        ```a
        >New_seq
        chr04	1499999	2000000	-
        ```
      * The file names a new sequence, `New_seq`, made from the selected region.
      * Output sequence ID can be controlled from the BLOCK file and/or by using the `-p | —prefix` flag.
      * Multiple regions can be extracted at the same time by listing different blocks with different output sequence IDs (like `New_seq1`, `New_seq2`, etc.) in the BLOCK file and/or by submitting a comma-separated list of files.
      * Starting coordinate is **0-based**.
      * [BLOCK file format](../block_format.md).

See [HaploMake usage page](../Usage/HaploMake_usage.md) for a complete description of usage options.

### Output

* `out.fasta`: FASTA of the selected sequence region with the new given name (`New_seq`).
* `out.gff3`: Annotation of the selected region with coordinates translated to the new sequence (`New_seq` based).
* `out.agp`: AGP file associating the selected region to the output. See the [AGP file format specifications](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/) for a complete description of its format.
  * If the input is in AGP format (`region.agp`), the two files are identical.
