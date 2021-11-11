# BLOCK file format

BLOCK files are simple text files that define the structure of a new sequence, without noting its final coordinates, based on the legacy sequences it is made of.

Example:

```
>NEW_seq1
legacy_seq1.1	0	1000000	+
legacy_seq1.2	0	2000000	-
legacy_seq1.3	500	10000	+
[...]
#Note on NEW_seq1

>NEW_seq2 Second sequences of the block file
legacy_seq2.1	1000	20000	+
legacy_seq2.2	0	30000	-
legacy_seq2.3	0	10000000	+
[...]
```

1. Empty lines are allowed.

2. Comment lines begin with `#` (number sign, hash).

3. Each new sequence structure is defined in a block, starting with a header line, and ending when the following block starts.

   1. Header:
      * Begins with `>`.
      * `>` is followed by the name of the sequence without spaces.
      * After the name of the sequence, additional information can be written.
   2. Elements:
      * Each element represents a region on a legacy sequence.
      * Tabular format with `BED`-style coordinates. The columns are:
        1. Legacy sequence ID
        2. Start (0-based)
        3. Stop (1-based)
        4. Orientation

   
