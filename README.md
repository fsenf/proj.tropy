# Tropy

A collection of python routines developed to read and analyze data at TROPOS

## Requirements
The `tropy`package depends on standard and typically scientific python modules like

* `numpy`
* `scipy`, `skit-learn`, `skit-
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
pip install .
```

## Problems
The `tropy` installation sometimes and somehow breaks on server machines that have several `gfortran` versions installed (typically via modules). 

Sometimes it help to `unload` all modules going back to the system-wide `gfortran` installation. Unfortunately, I am unaware of a clever solution to that problem. Sorry!  
