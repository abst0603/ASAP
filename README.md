# ASAP
ASAP - A Sub-sampling Approach forPreserving Topological Structures

## Matlab
For few ten thousand points epsilon_dense.m is fast. For a larger number of points please use the C++ code.

data: datafile as a matrix of NumberOfPoints*Dimensions

ep: the radius for subsampling

samples: the output matrix which contains the samples

## install
1- In this project we also use two other packages, namely "Nanoflann" and "pybind11". As a result to build the project completely you sould clone the ripository with the external submodules using following code:
```shell
git clone --recurse-submodules https://github.com/abst0603/ASAP.git
```
2- Run `pip install -r requirements.txt` or `pip3 install -r requirements.txt`
3- Build the project using CMake:
```shell
mkdir build
cd build
cmake ..
make
```
If you gave another path as your new PYTHONPATH since you are using anaconda, etc, Then you should use the following lines instead:
```shell
mkdir build
cd build
cmake .. -DPYTHON_INCLUDE_DIR=$(python -c "from distutils.sysconfig import get_python_inc; print(get_python_inc())")  -DPYTHON_LIBRARY=$(python -c "import distutils.sysconfig as sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
make
```
Check the link bellow for more information about this problem:
[cmake is not able to find python libraries](https://stackoverflow.com/questions/24174394/cmake-is-not-able-to-find-python-libraries)
4- When this is done, you should see a shared library is build in python folder. You can insert this library to subsample the data as follows:
```python
import asap
```
5- There is a full example in python folder on how to apply ASAP subsampling, then extracting cycles in every iteration of subsamples, and finally, apply voting procedure as described in our paper.


## Citation
Please consider citing our papers if you use this tool.
```
@article{taghribi_asap_2021,
	title = {{ASAP} - {A} {Sub}-sampling {Approach} for {Preserving} {Topological} {Structures} {Modeledwith} {Geodesic} {Topographic} {Mapping}},
	abstract = {Topological data analysis tools enjoy increasing popularity in a wide range of applications, such as Computer graphics, image analysis, Machine Learning and Astronomy for information extraction that is independent of a specific metric. 
However, due to computational complexity, processing large numbers of samples of higher dimensionality quickly becomes infeasible. 
This contribution is two-fold:
We present an efficient novel sub-sampling strategy inspired by Coulomb's law to decrease the number of data points in \$d\$-dimensional point clouds while preserving its homology.
The method is not only capable of reducing the memory and computation time needed for the construction of different types of simplicial complexes but also preserves the size of the voids in \$d\$-dimensions, which is crucial e.g. for astronomical applications. 
Furthermore, we propose a technique to construct a probabilistic description of the border of significant cycles and cavities inside the point cloud.
We demonstrate and empirically compare the strategy in several synthetic scenarios and an astronomical particle simulation of a dwarf galaxy for the detection of superbubbles (supernova signatures).},
	journal = {Neurocomputing},
	author = {Taghribi, Abolfazl and Canducci, Marco and Mastropietro, Michele and Rijcke, Sven De and Bunte, Kerstin and Tino, Peter},
	year = {2021}
}

@inproceedings{Taghribi2020,
title = {ASAP - A Sub-sampling Approach for Preserving Topological Structures},
author = {Abolfazl Taghribi and Kerstin Bunte and Michele Mastropietro and Sven De Rijcke and Peter Tino},
booktitle = {Proc. of the  28th "European Symposium on Artificial Neural Networks (ESANN)},
pages = {67-72},
year = {2020},
publisher = {Ciaco - i6doc.com},
publisher2 = {D-facto Publications},
editor = {M. Verleysen},
abstract = {Topological data analysis tools enjoy increasing popularity in a wide range of applications. However, due to computational complexity,    processing large number of samples of higher dimensionality quickly becomes infeasible.    We propose a novel sub-sampling strategy inspired by Coulomb’s law to decrease the number of data points in d-dimensional point clouds while preserving its Homology.    The method is not only capable of reducing the memory and computation time needed for the construction of different types of simplicial complexes but also preserves the size of the voids in d-dimensions, which is crucial e.g. for astronomical applications.  We demonstrate and compare the strategy in several synthetic scenarios and an astronomical particle simulation of a Jellyfish galaxy for the detection of superbubbles  (supernova signatures).},
}
```


## Manuscript
https://www.esann.org/sites/default/files/proceedings/2020/ES2020-147.pdf

## Contact
If you have a question you can directly raise an issue in the repository. Otherwise, feel free to send me an email (abolfazl.taghribi@gmail.com)
