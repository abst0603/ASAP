# ASAP
ASAP - A Sub-sampling Approach forPreserving Topological Structures

## Matlab
For few ten thousand points epsilon_dense.m is fast. For a larger number of points please use the C++ code.

data: datafile as a matrix of NumberOfPoints*Dimensions

ep: the radius for subsampling

samples: the output matrix which contains the samples

## C++
To build the run file please download the nanoflann library. It's a header-only library, therefore you don't need to install it, however, you should add the address of "utils.h" and "nanoflann.hpp" to your search address. This code can only process 3D data points, and it is not implemented for a higher dimension yet.

## Citation
@inproceedings{Taghribi2020,
title = {ASAP - A Sub-sampling Approach for Preserving Topological Structures},
author = {Abolfazl Taghribi and Kerstin Bunte and Michele Mastropietro and Sven De Rijcke and Peter Tino},
booktitle = {Proc. of the  28th "European Symposium on Artificial Neural Networks (ESANN)},
pages = {},
year = {2020},
publisher = {Ciaco - i6doc.com},
publisher2 = {D-facto Publications},
editor = {M. Verleysen},
abstract = {Topological data analysis tools enjoy increasing popularity in a wide range of applications. However, due to computational complexity,    processing large number of samples of higher dimensionality quickly becomes infeasible.    We propose a novel sub-sampling strategy inspired by Coulomb’s law to decrease the number of data points in d-dimensional point clouds while preserving its Homology.    The method is not only capable of reducing the memory and computation time needed for the construction of different types of simplicial complexes but also preserves the size of the voids in d-dimensions, which is crucial e.g. for astronomical applications.  We demonstrate and compare the strategy in several synthetic scenarios and an astronomical particle simulation of a Jellyfish galaxy for the detection of superbubbles  (supernova signatures).},
}

## Manuscript
https://www.esann.org/sites/default/files/proceedings/2020/ES2020-147.pdf

## Contact
If you have a question you can directly raise an issue in the repository. Otherwise, feel free to send me an email (abolfazl.taghribi@gmail.com)
