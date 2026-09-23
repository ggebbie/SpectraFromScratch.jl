# SpectraFromScratch

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ggebbie.github.io/SpectraFromScratch.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ggebbie.github.io/SpectraFromScratch.jl/dev)
[![Build Status](https://github.com/ggebbie/SpectraFromScratch.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ggebbie/SpectraFromScratch.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ggebbie/SpectraFromScratch.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ggebbie/SpectraFromScratch.jl)

Spectral analysis can get complicated, but the basic concepts are simple. Here I follow Tom Farrar's approach of building up a spectral analysis toolbox from scratch. It's not really from scratch as I rely on Steven Johnson's Fastest Fourier Transform of the West (FFTW.jl). The goal here is not to make the best operational spectral analysis, but instead to facilitate my learning and to make useful tools at the same time.

This package originated as a Julia Colab notebook in ipynb format written by Tom Farrar <jfarrar@whoi.edu>. Here that notebook is transformed to the standard Julia package format. 

# Usage guide

For examples on how to use this toolbox, see `test/runtests.jl`. 

Taking a Fourier transform by the Fast Fourier Transform method requires having a uniformly spaced timeseries. Put such a timeseries into the expected format first:
```julia
y = RegularTimeseries(yy, t)
```
where the constructor `RegularTimeseries` is defined by this package. Access the timeseries with `y.x` and the times with the private field `y.time`. 

You are now ready to take the Fourier Transform by FFT method:
```julia
ŷ = FourierTransform(y, alg=:centered_fft)
```
or simply by removing the optional keyword argument:
```julia
ŷ = FourierTransform(y) 
```

How much does the FFT speed up the calculation?
Compare to this manual expansion:
```julia
 ŷ = FourierTransform(y, alg=:manual)
```
and you will likely find that the FFT is a huge cost savings.

As the name of this package implies, it is also convenient to find the frequency spectrum. Currently, the raw spectrum (or periodogram) is calculated:
```julia
Ψ = periodogram(y)
```

Utilities for band averaging and computing confidence limits can also be found in this package. 

---

*This package was generated using PkgTemplates.jl.*

