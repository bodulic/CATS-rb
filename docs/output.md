# Output explanation

The analysis is summarized in the `CATS_rb_comparison.html` HTML file. 
An example of the HTML output is provided [here](CATS_rb_output_example.html)

Note on transcriptome assembly order and names: Assemblies will appear in the order they were provided on the command line when running the tool. For visualization purposes, assembly names are limited to a maximum of 20 characters; names exceeding this limit will be truncated. If multiple assemblies share the same name, a numeric suffix (e.g., .1, .2, etc.) will be appended to distinguish these assemblies.

## Summary tables

CATS-rb produces several summary files encompassing transcriptome assembly length statistics, various mapping metrics, and completeness analysis results:

`CATS_rb_general_statistics.tsv`: contains descriptive statistics of transcript length (mean, median, interquartile range, range, N50, L50, N90, L90( and GC content.

`CATS_rb_main_comparison_results.tsv`: contains transcriptome assembly mapping metrics and the results of relative completeness analysis.

`CATS_rb_annotation_based_analysis_results.tsv`: contains the results of annotation-based completeness analysis. This table is provided only if annotation-based analysis is enabled.

## Figures

CATS-rb produces several figures, providing a detailed visualization of CATS-rb quality metrics. 

`transcript_length` visualizes the distribution of transcript length.

`transcript_alignment_proportion` visualizes the distribution of transcript alignment proportion.

`number_of_exons_per_transcript` visualizes the distribution of exon number per transcript.

`exon_length` visualizes the distribution of exon length.

`exon_set_genomic_distribution` visualizes the positional distribution of exon sets in the analysed genome.

`number_of_isoforms_per_transcript_set` visualizes the distribution of isoform number per transcript set. 

`exon_set_length` and `transcript_set_length` visualize the distribution of exon/transcript set length.

`common_exon_set_length` and `common_transcript_set_length` visualize the distribution of common exon/transcript set length. Common sets correspond to set groups found in all analysed transcriptome assemblies.

`common_exon_set_relative_length` and `common_transcript_set_relative_length` visualize the distribution of common exon/transcript set relative length. Relative length is calculated with respect to the longest set within each group of common sets.

`unique_exon_set_length` and `unique_transcript_set_length` visualize the distribution of unique exon/transcript set length. Unique sets correspond to sets found in only one of the analysed transcriptome assemblies.

`exon_set_upset_plot.pdf` and `transcript_set_upset_plot.pdf` visualize UpSet plots for exon/transcript sets. Each UpSet plot is accompanied by two boxplots: the upper boxplot illustrates the length distribution of exon/transcript sets within each subset, while the lower boxplot displays the distribution of the ratio between the minimum and maximum exon/transcript set length within each subset.

`exon_set_pairwise_comp_similarity_tileplot` and `transcript_set_pairwise_comp_similarity_tileplot` visualize the exon/transcript set pairwise completeness similarity tileplot between the analysed transcriptome assemblies. Completeness similarity is defined as the mean completeness ratio of each corresponding exon/transcript set between each assembly pair.

`pairwise_exon_set_venn_diagrams` and `pairwise_transcript_set_venn_diagrams` visualize the exon/transcript set Venn diagrams for each transcriptome assembly pair. These plots are generated only if the comparison involves ten or fewer assemblies.

`exon_set_heatmap.pdf` and `transcript_set_heatmap.pdf` visualize hierarchical clustering heatmaps of transcriptome assemblies and exon/transcript sets. Assemblies (columns) are clustered based on the relative completeness of clustered exon/transcript sets (rows). Clustering is performed using complete linkage and Euclidean distance.

`unique_exon_set_position_in_non_origin_transcriptomes` visualizes the positional analysis of unique exon sets in non-origin transcriptome assemblies. Each unique exon set is classified as either: (1) located within a transcript set, (2) proximal to a transcript set, or (3) distant from any transcript set (in non-origin assemblies).

`missing_exon_set_position` visualizes the positional analysis of missing exon sets. Each missing exon set is classified as either: (1) located within a transcript set, or (2) located outside any transcript set (in the transcriptome assembly from which the exon set is missing).

`prop_of_exon_set_covered_by_a_gtf_set` and `prop_of_transcript_set_covered_by_a_gtf_set` visualize the distribution of the proportion of exon/transcript sets covered by a GTF set.

`annotation_based_exon_set_upset_plot.pdf` and `annotation_based_transcript_set_upset_plot.pdf` visualize annotation-based UpSet plots for exon/transcript sets. Each UpSet plot is accompanied by two boxplots: the upper boxplot illustrates the length distribution of exon/transcript sets within each subset, while the lower boxplot displays the distribution of the ratio between the minimum and maximum exon/transcript set length within each subset.

`annotation_based_pairwise_exon_set_venn_diagrams` and `annotation_based_pairwise_transcript_set_venn_diagrams` visualize annotation-based exon/transcript set Venn diagrams for each transcriptome assembly pair. These plots are generated only if the comparison involves ten or fewer assemblies.

`annotation_based_exon_set_heatmap.pdf` and `annotation_based_transcript_set_heatmap.pdf` visualize annotation-based hierarchical clustering heatmaps of transcriptome assemblies and exon/transcript sets. Assemblies (columns) are clustered based on the relative completeness of clustered exon/transcript sets (rows). Clustering is performed using complete linkage and Euclidean distance.

## Detailed tables

CATS-rb also produces several .tsv files containing detailed per-transcript and element set metrics:

`unmapped_transcripts.tsv` lists unmapped transcripts.

`transcript_aln_prop.tsv` contains the proportion of aligned transcript length for each transcript.

`transcripts_low_aln_prop.tsv` lists transcripts classified as structurally inconsistent due to low alignment rate.

`transcript_mapped_N.tsv` contains the number of mappings for each transcript.

`transcripts_disjunct_genomic_regions.tsv` lists transcripts classified as structurally inconsistent due to different transcript segments mapping to disjunct genomic regions.

`str_inconsistent_transcripts.tsv` lists all structurally inconsistent transcripts (unmapped + low alignment rate + segments mapping to disjunct genomic regions).

`per_transcript_exon_N.tsv` contains the exon number per transcript.

`exon_sets.tsv` and `transcript_sets.tsv` contain exon/transcript set coordinates.

`unique_exon_sets.tsv` and `unique_transcript_sets.tsv` contain unique exon/transcript set coordinates.

`missing_exon_set_ranges.tsv` contains genomic coordinate ranges of missing exon sets identified in other transcriptome assemblies. Range coordinates are defined by taking the range from minimum to maximum genomic coordinate of the exon set group in all assemblies in which the set was found.

`exon_set_pairwise_completeness_similarity_matrix.tsv` and `transcript_set_pairwise_completeness_similarity_matrix.tsv` contain exon/transcript set pairwise completeness similarity between the analysed transcriptome assemblies. Completeness similarity is defined as the mean completeness ratio of each corresponding exon/transcript set between each assembly pair.

`annotation_based_exon_set_coordinates.tsv` amd `annotation_based_transcript_set_coordinates.tsv` contain annotation-based exon/transcript set coordinates.

`exon_set_annotation_based_comp_matrix.tsv` and `transcript_set_annotation_based_comp_matrix.tsv` contain the completeness of exon/transcript sets relative to the GTF reference set.
