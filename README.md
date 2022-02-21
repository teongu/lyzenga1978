# Lyzenga1978

Water column correction from Lyzenga 1978.

This repository implements the water column correction described in Lyzenga 1978, "Passive remote sensing techniques for mapping water depth and bottom features" (https://doi.org/10.1364/AO.17.000379)

That has been strongly inspired by the work of @jkibele, but some modifications have been performed. 

## Introduction

The main idea of Lyzenga's paper is to be able to perform a water column correction of a multibands image, i.e. to remove the effects due to the water column on the reflectance perceived by the sensor. 

To be able to do so, we must have a multibands image at our disposal, and a few bathymetry points to be able to estimate a linear relationship between depth and the logarithm of reflectance. The strength of the algorithm is that it is not required to know the bathymetry in each point of the image.

A corrected image is called a depth invariant image, because it does not depend on the depth anymore. Such an image usually increases the accuracy of any classification performed on it (see Nguyen et al. 2021, "Mapping of Coral Reefs with Multispectral Satellites: A Review of Recent Papers", https://doi.org/10.3390/rs13214470), even though it can have some other advantages depending on the work you want to achieve.

## How to use

This module only requires the use of `numpy`, `matplotlib` and `sklearn`. 

The notebook `Lyzenga 1978 demo.ipynb` easily shows how the function should be used to perform the water column correction. 

Let's assume that you already know the slope of the linear regression between the depth and the logarithm of reflectance (the function `linear_regression` helps you with that), and you have stored them in a list `slopes`. 

If you have, for each band, computed the mean reflectance on deep water zonse and stored them in a list `dw`, and if your image is of the shape `[B_0,B_1,...,B_N]`, then the computation of the depth invariant image can basically be summed up in one call:

```
from water_column_correction import compute_Xi, Aij, depth_invariant, reshape_di

di_image = reshape_di(depth_invariant(Aij(slopes),compute_Xi(image,dw)),image)
```

But, to see more details, I strongly recommend you to have a look at the notebook `Lyzenga 1978 demo.ipynb`.


## Contact

Should you have any questions of remarks regarding this work, you can directly contact me via email: teo.nguyen@univ-pau.fr
