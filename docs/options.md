# Detailed options

CATS-rb offers a comprehensive list of options which allow users to control the analysis parameters.

## Genome index generation

The following options are available for `CATS_rb_index`:

`-m`: Maximum gene length (in bp), default: estimated from genome size

The value of `m` should be adjusted according to the analysed species.

`-t`: Number of CPU threads, default: 10

Spaln genome index generation is parallelized. Recommended number of threads: 10-20.

`-O`: Overwrite the genome index directory, default: off 

`-h`: Show usage information"

## Transcriptome assembly mapping

The following options are available for `CATS_rb_map`:

`-S`: Enable stranded mapping, default: off

Stranded mapping restricts Spaln to align transcripts in their native 5′ to 3′ orientation.

`-N`: Maximum number of mappings per transcript, default: 5

The value of `N` should be increased when analyzing species with complex genomes that contain a high number of paralogous genes, and decreased for smaller or less complex genomes.

`-i`: Minimum intron length (in bp), default: 20

The value of `i` should be adjusted according to the analysed species.

`-p`: Species-specific preset, default: unset

Spaln provides species-specific presets that control various mapping parameters to suit different genomes. A list of supported species and their input values can be found in the table/ directory within the Spaln installation (path_to_spaln_dir/table/), or on the [Spaln GitHub repository](https://github.com/ogotoh/spaln/blob/master/table/gnm2tab).

`-s`: Splice site characterization option, default: 2

The value of `s` controls how Spaln treats splice sites, following a predefined set of rules [adopted from Spaln github repository](https://github.com/ogotoh/spaln/blob/master):

0: accept only the canonical pairs (GT..AG,GC..AG,AT..AC)

1: accept also AT..AN

2: allow up to one mismatch from GT..AG

3: accept any pairs

`-P`: Relative contribution of coding potential to mapping score, default: 1

`-T`: Relative contribution of translation initiation signal to mapping score, default: 1

Values of `P` and `T` should be adjusted according to the leverage that should be given to protein-coding transcripts (higher values -> more leverage).

`-t`: Number of CPU threads, default: 10

Transcriptome assembly mapping by Spaln is parallelized. Recommended number of threads: 10-20.

`-D`: Mapping output directory name, default: TRANSCRIPTOME_CATS_rb_map

`-O`: Overwrite the mapping output directory, default: off

`-h`: Show usage information

## Transcriptome assembly mapping comparison

The following options are available for `CATS_rb_compare`:

`-S`: Enable stranded analysis, default: off

Stranded analysis ensures that CATS-rb only examines transcripts mapping in their native 5' to 3' orientation. Furthermore, element set coordinates will be dependent on the genomic strand to which the transcript maps. Stranded analysis should only be enabled if all analysed transcriptome assemblies were mapped in stranded mode.

`-p`: Minimum exon identity proportion, default: 0.98

More complex genomes should be assigned a higher value of `p` (e.g. 0.99 or 0.995) to minimze off-target mapping. On the other hand, the value of `p` should be reduced if working with reference genome from a related species.

`-e`: Minimum exon length (in bp), default: 20

`-i`: Maximum intron length (in bp), default: 100000

Values of `e` and `i` should be adjusted according to the analysed species.

`-M`: Alignment proportion threshold for structural inconsistency detection, default: 0.9

`-C`: Maximum proportion of allowed transcript segment overlap for identification of segments mapping to disjunct genomic regions, default: 0.3

A transcript is classified as structurally inconsistent if its alignment proportion falls below the threshold `M` or if it contains regions that map to disjunct genomic loci. The latter is assessed by identifying transcript regions overlapping by less than `C` and mapping either to different scaffolds, to opposite strands, or beyond the intron length threshold.

`-l`: Minimum exon set length for completeness analysis (in bp), default: 0

`-L`: Minimum transcript set length for completeness analysis (in bp), default: 100

`-m`: Maximum transcript set length for completeness analysis (in bp), default: 1000000

Thresholds for element set length `l`, `L`, and `m` should be adjusted according to the analysed species. Complex genomes should be assigned higher thresholds. Maximum transcript set length should be adjusted according to the expected maximum gene size.

`-j`: Minimum overlap between exon sets for edge specification (in bp), default: 1

`-J`: Minimum overlap between transcript sets for edge specification (in bp), default: 1

Values of `j` and `J` control the required overlap length for exon and transcript sets to be connected by edges when constructing inter-assembly exon/transcript set graphs.

`-o`: Minimum overlap between transcript set and transcript for isoform specification (in bp), default: 1

Isoforms are defined as transcripts overlapping the associated transcript set with a minimum of `o` bases.

`-P`: Transcript set proximity region length for unique exon set analysis (in bp), default: 5000

Genomic coordinates of unique exon sets from each transcriptome assembly are analyzed across non-origin assemblies. Each unique exon set is classified as either: (1) located within a transcript set, (2) proximal to a transcript set based on a defined threshold `P`, or (3) distant from any transcript set (in non-origin assemblies).

`-x`: Figure extension, default: png

`-d`: Figure DPI, default: 600

Extension (device) and DPI of each plotted figure are controlled with `x` and `d`, respectively.

`-r`: Raincloud plot colors (quoted hexadecimal codes or R color names, specified with x,y,z...), default: adjusted Set1 palette from RColorBrewer package

All color sets (parameters `r`, `b`, `n`, `u`, `v`, `y`, and `c`) should be supplied as R color names or hexadecimal codes separated with commas and enclosed in quotes (e.g. "#FDAF4A,#DC151D"). R color cheatsheet is available [here](https://sites.stat.columbia.edu/tzheng/files/Rcolor.pdf).

`-b`: Barplot colors (quoted hexadecimal codes or R color names, specified with x,y,z...), default: adjusted YlOrRd palette from RColorBrewer package

`-n`: Exon set genomic location plot colors (quoted hexadecimal codes or R color names, specified with x,y,z...), default: adjusted Set1 palette from RColorBrewer package

`-u`: UpSet plot bar and matrix colors (quoted hexadecimal codes or R color names, specified with x,y), default: "#FDAF4A,#DC151D"

`-v`: Venn diagram colors (quoted hexadecimal codes or R color names, specified with x,y,z...), default: adjusted Reds palette from RColorBrewer package

`-y`: Pairwise similarity tileplot colors (quoted hexadecimal codes or R color names, specified with x,y,z...), default: adjusted YlOrRd palette from RColorBrewer package

`-c`: Hierarchical clustering heatmap colors (quoted hexadecimal codes or R color names, specified with x,y,z...), default: adjusted YlOrRd palette from RColorBrewer package

`-q`: Maximum right-tail distribution quantile for raincloud plots, default: 0.995

Raincloud plots omit right-tail extreme values for visualization purposes. The x-axis in all raincloud plots is logarithmically scaled. Raincloud plot densities are normalized for each transcriptome assembly. Boxplots within raincloud plots mark the distribution median, Q<sub>1</sub>, and Q<sub>3</sub>, with whiskers extending from Q<sub>1</sub> - 1.5 * IQR to Q<sub>3</sub> + 1.5 * IQR of the distribution.

`-f`: Number of longest genomic scaffolds for exon set genomic location plot, default: all scaffolds

The value of `f` should be adjusted according to the number of relevant genomic scaffolds.

`-B`: Number of genomic bins for exon set genomic location plot, default: 25000

Higher values of `B` allow for a higher resolution of exon set genomic location plots.

`-V`: Minimum completeness threshold for assigning an element set to a Venn set, default: 0.35.

In Venn diagrams, element set completeness is used to define the plotted Venn sets. If element set completeness exceeds `V`, the set is considered shared with the reference element set. If both compared element sets are shared with the reference set, these element sets are considered common between the compared transcriptome assemblies.

`-H`: Number of longest element sets used in hierarchical clustering, default: 5000

Higher values of `H` will result in more detailed heatmaps, but significantly increase runtime and RAM usage. The value of `H` is capped at 65000.

`-E`: Use raster for heatmap plotting, default: off

Rasterization can be used to improve heatmap quality.

`-A`: Proportion of aligned transcript distribution breakpoints (specified with x,y,z...), default: "0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,1"

Proportion of aligned transcript is split into intervals defined by `A` (e.g. [0-0.2>, [0.2-0.4>...). This category variable is used for plotting. 

All category variable breaks (`A`, `N`, `R`, `I`, and `s`) should be supplied as strings separated with commas and enclosed in quotes (e.g. "0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,1").

`-N`: Number of exons per transcript distribution breakpoints (specified with x,y,z...), default: "1,2,4,6,8,10,15,20"

Number of exons per transcript is split into intervals defined by `N` (e.g. [0-2>, [2-4>...). This category variable is used for plotting.

`-R`: Common element set relative length distribution breakpoints (specified with x,y,z...), default: "0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,1"

Common element set relative length is split into intervals defined by `R` (e.g. [0-0.2>, [0.2-0.4>...). This category variable is used for plotting.

`-I`: Number of isoforms per transcript set distribution breakpoints (specified with x,y,z...), default: "1,2,4,6,8,10,15,20"

Number of isoforms per transcript set is split into intervals defined by `I` (e.g. [0-2>, [2-4>...). This category variable is used for plotting.

`-F`: GTF/GFF3 file for the annotation-based analysis

If a GTF/GFF3 file is supplied, CATS-rb will also perform the annotation-based analysis.

`-g`: Minimum proportion of an exon set that must be covered to be considered a match to a GTF exon set (and vice versa); default: 0.35

`-G`: Minimum proportion of a transcript set that must be covered to be considered a match to a GTF transcript set (and vice versa); default: 0.35

An assembly element set is considered matched to a GTF element set if their overlap exceeds a proportion `g` (for exon sets) or `G` (for transcript sets) of the assembly element set’s length. Conversely, a GTF element set is considered matched to an assembly element set if their overlap exceeds the same proportion of the GTF element set length.

`-s`: Proportion of element sets covered by a GTF set distribution breakpoints (specified with x,y,z...), default: "0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,1"

Proportion of element sets covered by a GTF set is split into intervals defined by `s` (e.g. [0-0.2>, [0.2-0.4>...). This category variable is used for plotting.

`-t`: Number of CPU threads, default: 10

Several steps of CATS-rb transcriptome assembly comparison are parallelized. This mainly includes operations performed by the data.table package. Recommended number of threads: 8-12.

`-D`: Comparison output directory name, default: CATS_rb_comparison

`-O`: Overwrite the comparison output directory, default: off

`-h`: Show usage information