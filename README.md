# hippogryph

[![PyPI - Version](https://img.shields.io/pypi/v/hippogryph.svg)](https://pypi.org/project/hippogryph)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/hippogryph.svg)](https://pypi.org/project/hippogryph)

-----

**Table of Contents**

- [Overview](#overview)
- [Installation](#installation)
- [About the Name](#about-the-name)
- [License](#license)

## Overview

`hippogryph` is a Python package for the generation of mapped, multiblock structured grids. A structured grid is a grid for which there is a relationship between the physical space coordinates $x$, $y$, and $z$ and the three indices $i$, $j$, and $k$ and a multiblock grid is constructed from a set of subgrids ("blocks") with potentially different connections between the coordinates and the indices. A "mapped" grid is a grid for which there is a transformation $f$ such that the nodal locations are determined by applying this mapping to the indices

$$\vec{x}(i,j,k) = f(i,j,k),$$

where the $\vec{x} = (x, y, z)$ are the nodal coordinates. While it is possible to write many (if not all) structured grid generation schemes in this form, the focus here is on mappings that are known beforehand and not computed using the numerical calculations that are commonly used within the field of numerical grid generation. The grid is first constructed in index space using primitive shapes, and then the physical space grid is constructed from a specified mapping. Currently, the focus is on tensor-product grids where $x=f(i)$, $y=f(j)$, and $z=f(k)$, though more mappings may be added in the future.

## Installation

```console
pip install hippogryph
```

## About the Name

Many years ago, the ancestor of this code generated overset grids, also called chimera grids, with all the hole cutting, trilinear interpolation, and so on. Very little of that code remains, with the exception of the grid stretching calculations. The chimera is a fire-breathing creature from Greek mythology with the head of a lion, the body of goat, and the tail of a serpent so it made sense to call the original code something along those lines. The hippogryph (or hippogriff) is half eagle and half horse, so that's what was chosen. 

## License

`hippogryph` is distributed under the terms of the [BSD-3-Clause](https://spdx.org/licenses/BSD-3-Clause.html) license.
