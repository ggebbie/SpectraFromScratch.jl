# SpectraFS

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ggebbie.github.io/SpectraFS.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ggebbie.github.io/SpectraFS.jl/dev)
[![Build Status](https://github.com/ggebbie/SpectraFS.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ggebbie/SpectraFS.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ggebbie/SpectraFS.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ggebbie/SpectraFS.jl)

Spectral analysis can get complicated, but the basic concepts are simple. Here I follow Tom Farrar's approach of building up a spectral analysis toolbox from scratch. It's not really from scratch as I rely on Steven Johnson's Fastest Fourier Transform of the West (FFTW.jl). The goal here is not to make the best operational spectral analysis, but instead to facilitate my learning and to make useful tools at the same time.

This package originated as a Julia Colab notebook in ipynb format. The goal is to transform it to the standard Julia package format. 

For examples on how to use this toolbox, see `test/runtests.jl`. 

# How this Julia package was started

PkgTemplates.tl was used for the original template. See `scripts/start_package.jl` for the details. Julia 1.0 is not supported. For integration with GitHub, I used the following steps:

1. Check that remote is set correctly (yes, was done automatically).
2. Check that default branch is set to main not master (done by git configuration on my machine).
3. `upgrade_manifest()`
4. Push output of PkgTemplates into a bare GitHub repository.
5. Check status of GitHub Actions. It failed due to Documenter.jl not being set up. (Requires a documenter key)
