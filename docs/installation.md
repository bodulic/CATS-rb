# Installation 

## Running CATS-rb via Docker

Starting from version 1.0.2, CATS-rb can be run using Docker:

Build the Docker image with:

```bash
docker pull bodulic/cats-rb
```

Run the container with:

```bash
docker run --rm -v "$PWD":/data -w /data bodulic/cats-rb CATS_rb_index

docker run --rm -v "$PWD":/data -w /data bodulic/cats-rb CATS_rb_map

docker run --rm -v "$PWD":/data -w /data bodulic/cats-rb CATS_rb_compare
```

## Running CATS-rb via Singularity

Build the Singularity image with:

```bash
singularity build cats_rb.sif docker://bodulic/cats-rb:latest
```

Run the container with:

```bash
singularity run cats_rb.sif CATS_rb_index

singularity run cats_rb.sif CATS_rb_map

singularity run cats_rb.sif CATS_rb_compare
```

## Installing CATS-rb via conda

CATS-rb and its dependencies can be directly installed via [Bioconda](https://anaconda.org/bioconda/cats-rb):

```bash
conda install -c bioconda cats-rb
```

In case of dependency conflicts, please see the Troubleshooting section.

## Installing CATS-rb from source

CATS-rb consists of Bash and R scripts located in the `scripts` directory of this repository. After cloning the repository, all CATS-rb scripts must be included in the `PATH` environment variable.

The following dependencies are required:

| **Dependency**           | **Tested Version** | **Homepage**                                                                    | **Conda Installation**                                        | **R installation**                             |
|--------------------------|--------------------|---------------------------------------------------------------------------------|---------------------------------------------------------------|------------------------------------------------|
| spaln                    | 3.0.1              | https://github.com/ogotoh/spaln                                                 | `conda install -c bioconda spaln`                             | /                                              | 
| R                        | 4.3.0.-4.4.3       | https://www.r-project.org                                                       | `conda install conda-forge::r-base`                           | /                                              |
| data.table (R)           | 1.16.4             | https://cran.r-project.org/package=data.table                                   | `conda install conda-forge::r-data.table`                     | `install.packages("data.table")`               |
| pandoc                   | 2.19.2             | https://pandoc.org/                                                             | `conda install conda-forge::pandoc`                           | /                                              |
| rmarkdown (R)            | 2.29               | https://cran.r-project.org/package=rmarkdown                                    | `conda install conda-forge::r-rmarkdown`                      | `install.packages("rnarkdown)`                 |
| ggplot2 (R)              | 3.5.1              | https://cran.r-project.org/web/packages/ggplot2                                 | `conda install conda-forge::r-ggplot2`                        | `install.packages("ggplot2")`                  |
| ggdist (R)               | 3.3.2              | https://cran.r-project.org/web/packages/ggdist                                  | `conda install conda-forge::r-ggdist`                         | `install.packages("ggdist")`                   | 
| GenomicRanges (R)        | 1.56.2             | https://www.bioconductor.org/packages/devel/bioc/html/GenomicRanges.html        | `conda install -c bioconda bioconductor-genomicranges`        | `BiocManager::install("GenomicRanges")`        | 
| Matrix (R)               | 1.7.1              | https://cran.r-project.org/web/packages/Matrix                                  | `conda install conda-forge::r-matrix`                         | `install.packages("Matrix")`                   |
| igraph (R)               | 4.4.2              | https://cran.r-project.org/web/packages/igraph                                  | `conda install conda-forge::igraph`                           | `install.packages("igraph")`                   |
| UpSetR (R)               | 2.1.3              | https://cran.r-project.org/web/packages/UpSetR                                  | `conda install conda-forge::r-upsetr`                         | `install.packages("UpSetR")`                   |
| ggVennDiagram (R)        | 1.5.2              | https://cran.r-project.org/web/packages/ggVennDiagram                           | /                                                             | `install.packages("ggVennDiagram")`            |
| egg (R)                  | 0.4.5              | https://cran.r-project.org/web/packages/egg/index.html                          | `conda install conda-forge::r-egg`                            | `install.packages("egg")`                      |
| ComplexHeatmap (R)       | 2.20.0             | https://www.bioconductor.org/packages/devel/bioc/html/ComplexHeatmap.html       | `conda install -c bioconda bioconductor-complexheatmap`       | `BiocManager::install("ComplexHeatmap")`       |
| GenomeInfoDb (R)         | 1.40.1             | https://www.bioconductor.org/packages/devel/bioc/html/GenomeInfoDb.html         | `conda install -c bioconda bioconductor-genomeinfodb`         | `BiocManager::install("GenomeInfoDb")`         |
| GenomicDistributions (R) | 1.12.0             | https://www.bioconductor.org/packages/devel/bioc/html/GenomicDistributions.html | `conda install -c bioconda bioconductor-genomicdistributions` | `BiocManager::install("GenomicDistributions")` |

R (Rscript), Spaln, and pandoc executables must be included in `PATH`. Tools denoted with (R) correspond to R packages and can be installed via conda or directly in R with the supplied commands. R package BiocManager is required when installing Bioconductor packages (GenomicRanges, ComplexHeatmap, GenomeInfoDb, and GenomicDistributions) in R.

### MacOS
If you are using MacOS, Bash (version >= 4.0) and GNU versions of core utilities are required. In this case, `PATH` should be adjusted so that CATS-rb uses GNU versions of core utilities:

- Install Bash ≥ 4.0 via [Homebrew](https://formulae.brew.sh/formula/bash):

```bash
brew install bash
```

- Install GNU utilities:

```bash
brew install coreutils gnu-sed gawk
```

- Add Bash and GNU utilities to your `PATH` (adjust path depending on your architecture):

For Apple Silicon:
```bash
export PATH="/opt/homebrew/bin:$PATH"
export PATH="/opt/homebrew/opt/coreutils/libexec/gnubin:$PATH"
export PATH="/opt/homebrew/opt/gnu-sed/libexec/gnubin:$PATH"
```

For Intel-based configurations
```bash
export PATH="/usr/local/bin:$PATH"
export PATH="/usr/local/opt/coreutils/libexec/gnubin:$PATH"
export PATH="/usr/local/opt/gnu-sed/libexec/gnubin:$PATH"
```

- Run CATS-rb using the installed Bash version:
  
```bash
bash CATS_rb
```

The stated changes can be made permanent by modifying the appropriate .rc file. 
