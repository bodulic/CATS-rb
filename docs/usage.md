# Example usage 

CATS-rb workflow consists of three scripts which should be run in succession:

## Genome index generation script: `CATS_rb_index`

CATS-rb maps transcripts to the reference genome using Spaln. The first step is performed by `CATS_rb_index`, which generates the Spaln genome index for the reference genome. The reference genome does not need to be repeat-masked.

Example usage:
```bash
CATS_rb_index [OPTIONS] GENOME
```

## Transcriptome assembly mapping script: `CATS_rb_map`

`CATS_rb_map` maps the analysed transcriptome assembly to the reference genome. The script should be run on all analysed assemblies individually.

Example usage:
```bash
CATS_rb_map [OPTIONS] GENOME_INDEX_DIR TRANSCRIPTOME
```

## Mapping comparison script: `CATS_rb_compare`

`CATS_rb_compare` compares the mapped transcriptome assemblies. Optionally, a reference GTF/GFF3 gene annotation file can be supplied.

While `CATS_rb_compare` is primarily designed to compare multiple transcriptome assemblies, it can also be used with a single assembly. However, the output will not contain relative completeness analysis.

The most important optional parameters which should be set according to the analysed organism are maximum intron length (`-i`), maximum transcript set length for completeness analysis (`-m`), and the number of longest genomic scaffolds for exon set genomic location plot ( `-f`). See the detailed options section for further explanation.

Example usage without reference annotation:
```bash
CATS_rb_compare [OPTIONS] GENOME TRANSCRIPTOME_MAP_DIR ...
```

Example usage with reference annotation:
```bash
CATS_rb_compare -F GTF_FILE [OTHER_OPTIONS] GENOME TRANSCRIPTOME_MAP_DIR ...
```