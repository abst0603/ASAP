# ASAP
ASAP - A Sub-sampling Approach for Preserving Topological Structures Modeledwith Geodesic Topographic Mapping

## Matlab
For few ten thousand points epsilon_dense.m is fast. For a larger number of points please use the C++ code.

data: datafile as a matrix of NumberOfPoints*Dimensions

ep: the radius for subsampling

samples: the output matrix which contains the samples

## install
1- In this project we also use two other packages, namely "Nanoflann" and "pybind11". As a result, to build the project completely you should clone the repository with the external submodules using the following code:
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
Check the link below for more information about this problem:
[cmake is not able to find python libraries](https://stackoverflow.com/questions/24174394/cmake-is-not-able-to-find-python-libraries)

4- When this is done, you should see a shared library is built in the python folder. You can insert this library to subsample the data as follows:
```python
import asap
```
5- There is a full example in python folder on how to apply ASAP subsampling, then extracting cycles in every iteration of subsamples, and finally, apply the voting procedure as described in our paper.


## Citation
Please consider citing our papers if you use this tool.
```
@article{TAGHRIBI2021,
title = {ASAP - A Sub-sampling Approach for Preserving Topological Structures Modeled with Geodesic Topographic Mapping},
journal = {Neurocomputing},
year = {2021},
issn = {0925-2312},
doi = {https://doi.org/10.1016/j.neucom.2021.05.108},
url = {https://www.sciencedirect.com/science/article/pii/S0925231221011139},
author = {Abolfazl Taghribi and Marco Canducci and Michele Mastropietro and Sven De Rijcke and Kerstin Bunte and Peter Tiňo},
keywords = {Topological Data Analysis, Persistent Homology, Sub-sampling, Generative Topographic Mapping, probabilistic modeling, Particle Simulation, supernova shells},
abstract = {Topological data analysis tools enjoy increasing popularity in a wide range of applications, such as Computer graphics, Image analysis, Machine learning, and Astronomy for extracting information. However, due to computational complexity, processing large numbers of samples of higher dimensionality quickly becomes infeasible. This contribution is two-fold: We present an efficient novel sub-sampling strategy inspired by Coulomb’s law to decrease the number of data points in d-dimensional point clouds while preserving its homology. The method is not only capable of reducing the memory and computation time needed for the construction of different types of simplicial complexes but also preserves the size of the voids in d-dimensions, which is crucial e.g. for astronomical applications. Furthermore, we propose a technique to construct a probabilistic description of the border of significant cycles and cavities inside the point cloud. We demonstrate and empirically compare the strategy in several synthetic scenarios and an astronomical particle simulation of a dwarf galaxy for the detection of superbubbles (supernova signatures).}
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


## Manuscripts
[ASAP - A Sub-sampling Approach for Preserving Topological Structures Modeled with Geodesic Topographic Mapping](https://www.sciencedirect.com/science/article/pii/S0925231221011139)

[ASAP - A Sub-sampling Approach for Preserving Topological Structures](https://www.esann.org/sites/default/files/proceedings/2020/ES2020-147.pdf)

## Contact
If you have a question you can directly raise an issue in the repository. Otherwise, feel free to send me an email (abolfazl.taghribi@gmail.com)
