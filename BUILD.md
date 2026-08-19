# eccsplorer-recipe

Conda recipe for [ECCsplorer](https://github.com/crimBubble/ECCsplorer), a pipeline
for detecting extrachromosomal circular DNA (eccDNA) from paired-end sequencing data.

This recipe builds a self-contained conda package because ECCsplorer is not available
on bioconda or biocontainers.

## Package details

| Field    | Value            |
|----------|------------------|
| Name     | `eccsplorer`     |
| Version  | `2022.01.1.1`    |
| Build    | 0                |
| Channel  | `siyangming`     |
| License  | GPL-3.0          |
| Upstream | https://github.com/crimBubble/ECCsplorer |

## Dependencies

The recipe pulls in the following runtime dependencies from `conda-forge` and
`bioconda` (channel order: `conda-forge`, `bioconda`, `defaults`):

- Python 3.7, numpy, biopython, scipy, pyRserve
- R (`r-base`) with ggplot2, ggrepel, gridExtra, dplyr
- blast, segemehl, samtools (>=1.9), bedtools (>=2.28.0)
- repeatexplorer2, trimmomatic, seqtk

> **Availability notes**
>
> - `pyrserve` is provided by `conda-forge`.
> - `segemehl`, `blast`, `samtools`, `bedtools`, `trimmomatic`, `seqtk` are
>   available on `bioconda`.
> - `repeatexplorer2` may not be present on every platform/channel combination.
>   If a dependency is missing during build or install, install it manually or
>   adjust the `requirements.run` list in `meta.yaml`.

## Build

Configure the channel order (one-time setup):

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Then build the package:

```bash
conda build eccsplorer-recipe/
```

## Install

```bash
conda install -c siyangming eccsplorer=2022.01.1.1
```

## Upload to Anaconda

```bash
anaconda login
anaconda upload <path-to-built-package>
# e.g. anaconda upload /opt/conda/conda-bld/noarch/eccsplorer-2022.01.1.1-0.tar.bz2
```

After upload the package is available at:
<https://anaconda.org/siyangming/eccsplorer>

## Post-install configuration

`build.sh` attempts to rewrite third-party tool paths in `lib/config.py` to point
at the conda environment's `$PREFIX/bin`. If a tool cannot be found at runtime,
edit `$CONDA_PREFIX/lib/eccsplorer/lib/config.py` manually and set the correct
absolute path (e.g. `$CONDA_PREFIX/bin/blastn`).

## Notes

- The package is declared `noarch: python` because ECCsplorer itself is pure
  Python (scripts only). Its compiled dependencies (blast, segemehl, samtools,
  bedtools, R) are resolved per-platform from bioconda/conda-forge.
- The source is cloned from the `2022.01.1.1` tag. If the tag is unavailable on
  the upstream repository, update `git_rev` in `meta.yaml` to the corresponding
  commit hash.
