---
layout: page
title: Documentation
---

The SSCHA code is a suite of 3 packages: the *CellConstructor* package, the *python-sscha* package, and the *tdscha* package. Each of them has its own documentation, which is available online.
- `cellconstructor` is a python library that can load, manipulate, and save crystal structures, and phonon dynamical matrices as well as perform various operations on them. It is the core of the SSCHA code, and it is used by both the *python-sscha* and the *tdscha* packages. We provide a documnetation of its API [here](CellConstructor.pdf).
- `python-sscha` is the main package of the SSCHA code, which implements the stochastic self-consistent harmonic approximation. It is used to perform the actual SSCHA calculations, and it provides a high-level interface to the functionalities of the code. The documentation of the *python-sscha* code can be found instead [here](python-sscha.pdf), where the complete API to public methods are discussed.
- `tdscha` is the python library that performs the dynamical linear response calculations, the complete documentation of this package is available in the [documentation pages](https://sscha.eu/tdscha).
