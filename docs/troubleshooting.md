# Troubleshooting

Please report all potential bugs in the Issues tracker.

## Singularity

If you run into Singularity errors involving the default `TMPDIR` environment variable, set `TMPDIR` to your current working directory so the container uses a writable location for temporary files:

```bash
TMPDR=$(pwd)
```

## Conda installation

In case of dependency conflicts, try setting the following channel priority:

```bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

If this does not work, try setting channel priority to flexible:

```bash
conda config --set channel_priority flexible
```

If installation fails on ARM64 Mac architecture (Apple Silicon), try creating a new environment emulating the x86_64 architecture:

```bash
CONDA_SUBDIR=osx-64 conda create -n cats_rb -c conda-forge -c bioconda -c defaults r-base=4.3 cats-rb
```