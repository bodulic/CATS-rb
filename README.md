# CATS-rb

<img src="cats_rb_logo.png" alt="Logo" width="750" height="160"/>

# Introduction 

CATS-rb is the reference-based module of the CATS (Comprehensive Assessment of Transcript Sequences) framework. It evaluates transcriptome assembly quality using the reference genome of the corresponding or a closely related species. The pipeline maps transcripts to the reference genome and examines several mapping and completeness metrics.

The main contribution of CATS-rb is the transcriptome assembly completeness analysis, which can be performed in two settings: 

- Relative completeness analysis: requires multiple (two or more) transcriptome assemblies
- Annotation-based completenes analysis: requires one or more transcriptome assemblies and a reference gene annotation

Completeness analysis introduces exon and transcript sets as units for assembly comparison, together denoted as element sets. Precisely, CATS-rb collapses overlapping exon and transcript genomic coordinates of a given assembly into non-redundant exon and transcript sets, respectively. Completeness of exon/transcript sets is compared between the analysed assemblies by constructing an undirected graph in which vertices represent exon/transcript sets and edges indicate overlaps between the corresponding sets of the compared assemblies. Overlapping exon/transcript sets are grouped into connected components, with the longest set designated as the group representative.

Element set completeness is quantified by its relative length compared to the representative set. Relative exon and transcript scores for each transcriptome assembly are computed as the mean of exon and transcript set completeness. Alongside completeness scores, CATS-rb also provides an in-depth analysis of missing, common, and unique element sets.

Additionally, CATS-rb can perform an annotation-based analysis using element sets derived from a GTF file. This workflow follows the same principles as relative completeness analysis, while grouping element sets based on overlaps with GTF sets. As such, GTF sets are considered the representative for each set group. Annotation-based exon and transcript scores are calculated analogously to relative exon and transcript scores, offering an absolute measure of assembly completeness.

CATS-rb exon and transcript scores exhibit a strong correlation with transcriptome assembly quality. Furthermore, relative and annotation-based scores are strongly correlated when applied to assembly sets with varying quality, enabling precise assembly quality assessment without strictly requiring reference annotation. 

For detailed benchmarks and methodology, please refer to the CATS [preprint](test)

## Use cases

A typical CATS-rb analysis should follow one of the following use cases:

- Assessing accuracy and completeness of one or more transcriptome assemblies relative to the reference transcriptome assembly and/or reference gene annotation
- Comparing relative accuracy and completeness of transcriptome assemblies generated from the same RNA-seq library
- Comparing transcriptomic content between different RNA-seq libraries

# Installation 

## Compatibility

### Linux and Windows

For the best compatibility and performance, we recommend running CATS-rb on:
- Any modern Linux distribution (e.g. Ubuntu, Debian, Fedora, etc.)
- WSL (i.e. Ubuntu on Windows)

### MacOS
If you are using MacOS, Bash (version >= 4.0) and GNU versions of core utilities are required. In this case, `PATH` variable should be adjusted so that CATS-rb uses GNU versions of core utilities:

- Install Bash ≥ 4.0 via [Homebrew](https://formulae.brew.sh/formula/bash):

```bash
brew install bash
```

- Install GNU core utilities:

```bash
brew install coreutils findutils gnu-sed gawk grep
```

- Add Bash and GNU tools to your `PATH` (adjust path depending on your architecture):

For Apple Sillicon:
```bash
export PATH="/opt/homebrew/bin:$PATH"
export PATH="/opt/homebrew/opt/coreutils/libexec/gnubin:$PATH"
export PATH="/opt/homebrew/opt/findutils/libexec/gnubin:$PATH"
export PATH="/opt/homebrew/opt/gnu-sed/libexec/gnubin:$PATH"
export PATH="/opt/homebrew/opt/grep/libexec/gnubin:$PATH"
```

For Intel-based configurations
```bash
export PATH="/usr/local/bin:$PATH"
export PATH="/usr/local/opt/coreutils/libexec/gnubin:$PATH"
export PATH="/usr/local/opt/findutils/libexec/gnubin:$PATH"
export PATH="/usr/local/opt/gnu-sed/libexec/gnubin:$PATH"
export PATH="/usr/local/opt/grep/libexec/gnubin:$PATH"
```

- Run the script using the installed Bash version:
  
```bash
bash CATS_rb
```

## Installing CATS-rb via conda

CATS-rb and its dependencies can be directly installed via [Bioconda](https://bioconda.github.io/):

```bash
conda install -c bioconda cats_rb
```

## Installing CATS-rb from source

CATS-rb consists of Bash and R scripts located in the `scripts` directory of this repository. After cloning the repository, all CATS-rb scripts must be included in the `PATH` environment variable. 

The following dependencies are required:

| **Dependency**           | **Tested Version** | **Homepage**                                                                    | **Conda Installation**                                        | **R installation**                             |
|--------------------------|--------------------|---------------------------------------------------------------------------------|---------------------------------------------------------------|------------------------------------------------|
| spaln                    | 3.0.1              | https://github.com/ogotoh/spaln                                                 | `conda install -c bioconda spaln`                             | /                                              | 
| R                        | 4.4.2              | https://www.r-project.org/                                                      | `conda install conda-forge::r-base`                           | /                                              |
| knitr (R)                | 1.49               | https://cran.r-project.org/web/packages/knitr/                                  | `conda install conda-forge::r-knitr`                          | `install.packages("knitr")`                    |
| data.table (R)           | 1.16.4             | https://cran.r-project.org/package=data.table                                   | `conda install conda-forge::r-data.table`                     | `install.packages("data.table")`               |
| ggplot2 (R)              | 3.5.1              | https://cran.r-project.org/web/packages/ggplot2/                                | `conda install conda-forge::r-ggplot2`                        | `install.packages("ggplot2")`                  |
| ggdist (R)               | 3.3.2              | https://cran.r-project.org/web/packages/ggdist/                                 | `conda install conda-forge::r-ggdist`                         | `install.packages("ggdist")`                   | 
| GenomicRanges (R)        | 1.56.2             | https://bioconductor.org/packages/devel/bioc/html/GenomicRanges.html            | `conda install -c bioconda bioconductor-genomicranges`        | `BiocManager::install("GenomicRanges")`        | 
| Matrix (R)               | 1.7.1              | https://cran.r-project.org/web/packages/Matrix/                                 | `conda install conda-forge::r-matrix`                         | `install.packages("Matrix")`                   |
| igraph (R)               | 4.4.2              | https://cran.r-project.org/web/packages/igraph/                                 | `conda install conda-forge::igraph`                           | `install.packages("igraph")`                   |
| UpSetR (R)               | 2.1.3              | https://cran.r-project.org/web/packages/UpSetR/                                 | `conda install conda-forge::r-upsetr`                         | `install.packages("UpSetR")`                   |
| ggVennDiagram (R)        | 1.5.2              | https://cran.r-project.org/web/packages/ggVennDiagram/                          | /                                                             | `install.packages("ggVennDiagram")`            |
| egg (R)                  | 0.4.5              | https://cran.r-project.org/web/packages/egg/index.html                          | `conda install conda-forge::r-egg`                            | `install.packages("egg")`                      |
| ComplexHeatmap (R)       | 2.20.0             | https://www.bioconductor.org/packages/devel/bioc/html/ComplexHeatmap.html       | `conda install -c bioconda bioconductor-complexheatmap`       | `BiocManager::install("ComplexHeatmap")`       |
| GenomicDistributions (R) | 1.12.0             | https://www.bioconductor.org/packages/devel/bioc/html/GenomicDistributions.html | `conda install -c bioconda bioconductor-genomicdistributions` | `BiocManager::install("GenomicDistributions")` |

R (Rscript) and Spaln executables must be included in `PATH`. Tools denoted with (R) correspond to R packages and can be installed via conda or directly in R with the supplied commands. R package BiocManager is required when installing Bioconductor packages (GenomicRanges, ComplexHeatmap, and GenomicDistributions) im R.

# Test data

CATS-rb installation can be tested using instructions and files located in `test_data` directory.

# Example usage 

CATS-rb workflow consists of three scripts which should be run in succession:

## Genome index generation script: `CATS_rb_index`

`CATS_rb_index` generates the Spaln genome index for the reference genome.

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

`CATS_rb_compare` compares the mapped transcriptome assemblies. Optionally, a reference gene annotation file can be supplied.

While `CATS_rb_compare` is primarily designed to compare multiple transcriptome assemblies, it can also be used with a single assembly to analyse general assembly length, composition, and CATS-rb mapping metrics.

Example usage without reference annotation:
```bash
CATS_rb_compare [OPTIONS] GENOME TRANSCRIPTOME_MAP_DIR ...
```

Example usage with reference annotation:
```bash
CATS_rb_compare -F GTF_FILE [OTHER_OPTIONS] GENOME TRANSCRIPTOME_MAP_DIR ...
```

# Detailed options

CATS-rb offers a comprehensive list of options which allow users to control the analysis parameters.

## Genome index generation

The following options are available for `CATS_rb_index`:

`-m`: Maximum gene length (in bp), default: estimated from genome size

Value of `m` should be adjusted according to the analysed species.

`-t`: Number of CPU threads, default: 10

Spaln genome index generation is parallelized. Recommended number of threads: 10-20.

`-O`: Overwrite the genome index directory, default: off" 

`-h`: Show usage information"

## Transcriptome assembly mapping

The following options are available for `CATS_rb_map`:

`-S`: Enable stranded mapping, default: off

Stranded mapping restricts Spaln to align transcripts in their native 5′ to 3′ orientation.

`-N`: Maximum number of mappings per transcript, default: 5

Value of `N` should be increased when analyzing species with complex genomes that contain many paralogous genes, and decreased for smaller or less complex genomes.

`-i`: Minimum intron length (in bp), default: 20

Value of `i` should be adjusted according to the analysed species.

`-p`: Species-specific preset, default: unset

Spaln provides species-specific presets that control various mapping parameters to suit different genomes. A list of supported species and their input values can be found in the table/ directory within the Spaln installation (path_to_spaln_dir/table/), or on the [Spaln GitHub repository](https://github.com/ogotoh/spaln/blob/master/table/gnm2tab).

`-s`: Splice site characterization option, default: 2

Value of `s` controls how Spaln treats splice sites, following a predefined set of rules [adopted from Spaln github](https://github.com/ogotoh/spaln/blob/master):

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

Strandness analysis ensures that CATS-rb only examines transcripts mapping in their native 5' to 3' orientation. Furthermore, element set coordinates are depdendent on the genomic strand to which the transcript maps.

`-p`: Minimum exon identity proportion, default: 0.98

Parameter `p` controls the minimum identity proportion for an exon to be analysed. More complex genomes should be assigned a higher value of `p` (e.g. 0.99 or 0.995) to minimze off-target mapping.

`-e`: Minimum exon length (in bp), default: 20

`-i`: Maximum intron length (in bp), default: 100000

Parameters `e` and `i` should be adjusted according to the analysed species.

`-M`: Alignment proportion threshold for structural inconsistency detection, default: 0.9

`-C`: Maximum proportion of allowed transcript segment overlap for identification of segments mapping to disjunct genomic regions, default: 0.3

A transcript is classified as structurally inconsistent if its alignment proportion falls below the threshold `M` or if it contains regions that map to disjunct genomic loci. The latter is assessed by identifying transcript regions overlapping by less than `C` and mapping either to different scaffolds, to opposite strands, or beyond the intron length threshold.

`-l`: Minimum exon set length for completeness analysis (in bp), default: 0

`-L`: Minimum transcript set length for completeness analysis (in bp), default: 100

`-m`: Maximum transcript set length for completeness analysis (in bp), default: 1000000

Thresholds for element set length should be adjusted according to the analysed species. Complex genomes should be assigned higher thresholds. Maximum transcript set length should be adjusted according to the expected maximum gene size.

`-j`: Minimum overlap between exon sets for edge specification (in bp), default: 1

`-J`: Minimum overlap between transcript sets for edge specification (in bp), default: 1

Parameters `j` and `J` control the required overlap length for exon and transcript sets to be connected by edges when constructing inter-assembly exon/transcript set graphs.

`-o`: Minimum overlap between transcript set and transcript for isoform specification (in bp), default: 1

Isoforms are defined as transcripts overlapping the associated transcript set with a minimum of `o` bases.

`-P`: Transcript set proximity region length for unique exon set analysis (in bp), default: 5000

Genomic coordinates of unique exon sets from each transcriptome assembly are analyzed across non-origin assemblies. Each unique exon set is classified as either: (1) located within a transcript set, (2) proximal to a transcript set based on a defined threshold `P`, or (3) distant from any transcript set (in non-origin transcriptome assemblies).

`-x`: Figure extension, default: png

`-d`: Figure DPI, default: 300

Parameters `x` and `d` control the extension and DPI of each plotted figure. Note that UpSet plots and heatmaps will always be saved as a pdf file.

`-r`: Raincloud plot colors (quoted hexadecimal codes or R color names, specified with x,y,z...), default: adjusted Set1 palette from RColorBrewer package

All color name options (parameters `r`, `b`, `n`, `u`, `v`, `y`, and `c`) should be given as R color names or hexadecimal codes separated with commas and enclosed in quotes (e.g. "#FDAF4A,#DC151D"). R color cheatsheet is available [here](https://sites.stat.columbia.edu/tzheng/files/Rcolor.pdf)

`-b`: Barplot colors (quoted hexadecimal codes or R color names, specified with x,y,z...), default: adjusted YlOrRd palette from RColorBrewer package

`-n`: Exon set genomic location plot colors (quoted hexadecimal codes or R color names, specified with x,y,z...), default: adjusted Set1 palette from RColorBrewer package

`-u`: UpSet plot bar and matrix colors (quoted hexadecimal codes or R color names, specified with x,y,z...), default: "#FDAF4A,#DC151D"

`-v`: Venn diagram colors (quoted hexadecimal codes or R color names, specified with x,y,z...), default: adjusted Reds palette from RColorBrewer package

`-y`: Pairwise similarity tileplot colors (quoted hexadecimal codes or R color names, specified with x,y,z...), default: adjusted YlOrRd palette from RColorBrewer package

`-c`: Hierarchical clustering heatmap colors (quoted hexadecimal codes or R color names, specified with x,y,z...), default: adjusted YlOrRd palette from RColorBrewer package

`-q`: Maximum right-tail distribution quantile for raincloud plots, default: 0.995

Raincloud plots omit right-tail extreme values for visualization purposes. The x-axis in all raincloud plots is logarithmically scaled.

`-f`: Number of longest genomic scaffolds for exon set genomic location plot, default: all scaffolds

Parameter `f` should be adjusted according to the number of relevant genomic scaffolds.

`-B`: Number of genomic bins for exon set genomic location plot, default: 50000

Increased values of `B` allow for a higher resolution of exon set genomic location plots.

`-V`: Minimum completeness threshold for assigning an element set to a Venn set, default: 0.35.

In Venn diagrams, element set completeness is used to define the plotted Venn sets. If element set completeness exceedes `V`, the set is considered shared with the reference element set.

`-H`: Number of longest element sets used in hierarchical clustering, default: 5000

Higher values of parameter `H` will result in more detailed heatmaps, but significantly increase runtime and RAM usage. The parameter is capped at 65000.

`-E`: Use raster for heatmap plotting, default: off

Rasterization can be used to improve the plotting quality of resulting heatmaps.

`-A`: Proportion of aligned transcript distribution breakpoints (specified with x,y,z...), default: "0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,1"

Proportion of aligned transcript is split into intervals defined by `A` (e.g. [0-0.2>, [0.2-0.4>...). This category variable is used for plotting. 

All category variable breaks (options `A`, `N`, `R`, `I`, and `s`) should be given as strings separated with commas and enclosed in quotes (e.g. "0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,1").

`-N`: Number of exons per transcript distribution breakpoints (specified with x,y,z...), default: "0,2,4,6,8,10,15,20"

Number of exons per transcript is split into intervals defined by `N` (e.g. [0-2>, [2-4>...). This category variable is used for plotting.

`-R`: Common element set relative length distribution breakpoints (specified with x,y,z...), default: "0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,1"

Common element set relative length is split into intervals defined by `R` (e.g. [0-0.2>, [0.2-0.4>...). This category variable is used for plotting.

`-I`: Number of isoforms per transcript set distribution breakpoints (specified with x,y,z...), default: "0,2,4,6,8,10,15,20"

Number of isoforms per transcript set is split into intervals defined by `I` (e.g. [0-2>, [2-4>...). This category variable is used for plotting.

`-F`: GTF/GFF3 file for the annotation-based analysis

If a GTF/GFF3 file is supplied, CATS-rb will perform the annotation-based analysis.

`-g`: Minimum proportion of an exon set that must be covered to be considered a match to a GTF exon set (and vice versa); default: 0.35

`-G`: Minimum proportion of a transcript set that must be covered to be considered a match to a GTF transcript set (and vice versa); default: 0.35

An assembly element set is considered matched to a GTF element set if their overlap exceeds a proportion `g` or `G` of the assembly element set’s length. Conversely, a GTF element set is considered matched to an assembly element set if their overlap exceeds the same proportion of the GTF element set length.

`-s`: Proportion of element sets covered by a GTF set distribution breakpoints (specified with x,y,z...), default: "0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,1"

Proportion of element sets covered by a GTF set is split into intervals defined by `s` (e.g. [0-0.2>, [0.2-0.4>...). This category variable is used for plotting.

`-t`: Number of CPU threads, default: 10

Several steps of CATS-rb are paralellized. This mainly includes operations performed by the data.table package. Recommended number of threads: 8-12.

`-D`: Comparison output directory name, default: CATS_rb_comparison

`-O`: Overwrite the comparison output directory, default: off

`-h`: Show usage information

# Output explanation

The analysis is summarized in the `CATS_rb_comparison.html` HTML file. 
An example of the HTML output is provided [here](output_example.html)

## Summary tables

CATS-rb produces several summary files encompassing transcriptome assembly length statistics, various mapping metrics, and completeness analysis results:

`CATS_rb_general_statistics.tsv`: contains descriptive statistics of transcript length (mean, median, interquartile range, range, N50, L50, N90, L90), GC content, and read mapping rate.

`CATS_rb_main_comparison_results.tsv`: contains transcriptome assembly mapping metrics and the results of relative completeness analysis

`CATS_rb_annotation_based_analysis_results.tsv`: contains the results of annotation-based completeness analysis. This table is provided only if the annotation-based analysis is enabled.

## Figures

CATS-rb produces several figures, providing a detailed visualization of quality metrics. 

`transcript_length` visualizes the distribution of transcript length.

`transcript_alignment_proportion` visualizes the distribution of transcript alignment proportion.

`number_of_exons_per_transcript` visualizes the distribution of exon number per transcript.

`exon_length` visualizes the distribution of exon length.

`exon_set_genomic_distribution` visualizes the positional distribution of exon sets in the analysed genome.

`number_of_isoforms_per_transcript_set` visualizes the distribution of isoform number per transcript set. 

`exon_set_length` and `transcript_set_length` visualize the distribution of exon/transcript set length.

`common_exon_set_length` and `common_transcript_set_length` visualize the distribution of common exon/transcript set length. Common sets correspond to set groups found in all analysed trasncriptome assemblies.

`common_exon_set_relative_length` and `common_transcript_set_relative_length` visualize the distribution of relative common exon/transcript set length. Relative length is calculated with respect to the longest set within each group of common sets.

`unique_exon_set_length` and `unique_transcript_set_length` visualize the distribution of unique exon/transcript set length. Unique sets correspond to sets found in only one of the analysed transcriptome asseblies.

`exon_set_upset_plot.pdf` and `transcript_set_upset_plot.pdf` visualize the UpSet plot for exon/transcript sets. The UpSet plot is accompanied by two boxplots: the upper boxplot illustrates the length distribution of exon/transcript sets within each subset, while the lower boxplot displays the distribution of the ratio between the minimum and maximum exon/transcript set length within each subset.

`exon_set_pairwise_comp_similarity_tileplot` and `transcript_set_pairwise_comp_similarity_tileplot` visualize the exon/transcript set pairwise completeness similarity tileplot between the analysed transcriptome assemblies. Completeness similarity is defined as the mean completeness ratio of each corresponding exon/transcript set between each transcriptome assembly pair.

`pairwise_exon_set_venn_diagrams` and `pairwise_transcript_set_venn_diagrams` visualize the exon/transcript set Venn diagrams for each pair of the analysed assemblies. These plots are generated only when the comparison involves ten or fewer transcriptome assemblies.

`exon_set_heatmap.pdf` and `transcript_set_heatmap.pdf` visualize hierarchical clustering heatmaps of exon/transcript sets. Transcriptome assemblies (columns) are clustered based on the relative completeness of exon/transcript sets (rows). Clustering is performed using complete linkage and Euclidean distance.

`unique_exon_set_position_in_non_origin_transcriptomes` visualizes the positional analysis of unique exon sets in non-origin transcriptome assemblies. Each unique exon set is classified as either: (1) located within a transcript set, (2) proximal to a transcript set, or (3) distant from any transcript set (in non-origin transcriptome assemblies).

`missing_exon_set_position` visualizes the positional analysis of missing exon sets. Each missing exon set is classified as either: (1) located in a transcript set, or (2) located outside any transcript set.

`prop_of_exon_set_covered_by_a_gtf_set` and `prop_of_transcript_set_covered_by_a_gtf_set` visualize the distribution of the proportion of exon/transcript sets covered by a GTF set.

`annotation_based_exon_set_upset_plot.pdf` and `annotation_based_transcript_set_upset_plot.pdf` visualize the annotation-based UpSet plot for exon/transcript sets. The UpSet plot is accompanied by two boxplots: the upper boxplot illustrates the length distribution of exon/transcript sets within each subset, while the lower boxplot displays the distribution of the ratio between the minimum and maximum exon/transcript set length within each subset.

`annotation_based_pairwise_exon_set_venn_diagrams` and `annotation_based_pairwise_transcript_set_venn_diagrams` visualize the annotation-based exon/transcript set Venn diagrams for each pair of the analysed assemblies. These plots are generated only when the comparison involves ten or fewer transcriptome assemblies.

`annotation_based_exon_set_heatmap.pdf` and `annotation_based_transcript_set_heatmap.pdf` visualize annotation-based hierarchical clustering heatmaps of exon/transcript sets. Transcriptome assemblies (columns) are clustered based on the relative completeness of exon/transcript sets (rows). Clustering is performed using complete linkage and Euclidean distance.

## Detailed tables

CATS-rb also produces several .tsv files containing detailed per-transcript and element set metrics:

`unmapped_transcripts.tsv` lists unmapped transcripts.

`transcript_aln_prop.tsv` contains the proportion of aligned transcript length for each transcript.

`transcripts_low_aln_prop.tsv` lists transcripts classified as structurally inconsistent due to low alignment rate.

`transcript_mapped_N.tsv` contains the number of mappings for each transcript.

`transcripts_disjunct_genomic_region.tsv` lists transcripts classified as structurally inconsistent, with different transcript regions mapping to disjunct genomic locations.

`str_inconsistent_transcripts.tsv` lists all structurally inconsistent transcripts (unmapped + low alignment rate + segments mapping to disjunct genomic locations).

`per_transcript_exon_n.tsv` contains the exon number per transcript.

`exon_sets.tsv` and `transcript_sets.tsv` contain exon/transcript set coordinates.

`unique_exon_sets.tsv` and `unique_transcript_sets.tsv` contain unique eoxn/transcript set coordinates.

`missing_exon_set_ranges.tsv` contains coordinate ranges of missing exon sets identified in other transcriptome assemblies. Range coordinates are defined by taking the range form minimum to maximum coordinate of the exon set group in all assemblies in which the set was identified.

`exon_set_pairwise_completeness_similarity_matrix.tsv` and `transcript_set_pairwise_completeness_similarity_matrix.tsv` contain the exon/transcript set pairwise completeness similarity tileplot between the analysed transcriptome assemblies. Completeness similarity is defined as the mean completeness ratio of each corresponding exon/transcript set between each transcriptome assembly pair.

# Citation

CATS is an academic software distributed under the MIT license. 

Copyright © 2025 Kristian Bodulić

if you use CATS, please cite the CATS preprint:

(add reference)

# Troubleshooting

Please report all potential bugs in the Issues tracker.

# Changelog

Version 1.0.0. Initial commit, June 5, 2025.
