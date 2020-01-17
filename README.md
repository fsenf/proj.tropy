# Tropy

A collection of python routines developed to read and analyze data at TROPOS

## Requirements
The `tropy`package depends on standard and typically scientific python modules like

* `numpy`
* `scipy`, `sklearn`, `skimage`
* `cv2` : needed for optical flow calculations
* `mahotas` & `SimpleITK`: for segmentation and some cool filtering
* ...

The installation routine needs to run `gfortran`compiler to build a fast interface routine for reading MSG SEVIRI data.


## Installation
Download the code from github

```
git clone https://github.com/fsenf/proj.tropy.git
```

Try to install it locally via `pip`

```
pip install --upgrade .
```

## Problems
The `tropy` installation sometimes and somehow breaks on server machines that have several `gfortran` versions installed (typically via modules). 

Sometimes it help to `unload` all modules going back to the system-wide `gfortran` installation. Unfortunately, I am unaware of a clever solution to that problem. Sorry!  


## Documentation and Tutorial
This two aspects have been started. The status is rather poor and needed to be / will be improved in future.

Basic documenation can found here: [https://tropy.readthedocs.io/en/latest/](https://tropy.readthedocs.io/en/latest/)

A set of jupyter notebooks are collected in a separate project [https://github.com/fsenf/proj.tropy-tutorials](https://github.com/fsenf/proj.tropy-tutorials) which is parsed here: [https://tropy-tutorials.readthedocs.io/en/latest/](https://tropy-tutorials.readthedocs.io/en/latest/)
