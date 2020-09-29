# phase-contrast-phantom

Code to upsample X-ray phantoms for use in phase-contrast imaging.

## Overview 

Simulations of propagation-based phase-contrast X-ray imaging (PB-PCXI) requires micrometer level sampling. Most available clinical phantoms have sampling in the range 0.1 to 1 mm. Using these phantoms result in simulation artifacts. To solve this the phantoms can be upsampled, i.e. voxels are divided into smaller ones. The code in this repository does this in a careful way to minimize distortion of features and to not introduce new artifacts in the process. The concept is shown on a mammography phantom (see below).


## Mammography Phantom

Original phantom taken from [VICTRE project](https://github.com/DIDSR/VICTRE). This is one of the sample phantoms provided. The sampling is 50 µm. An upsampled version (32x -> 1.5625 µm sampling) is available. In the upsampled version only the projected thickness of each material is saved to reduce data size. 


## Publication

__In Silico Phase-Contrast X-Ray Imaging of Anthropomorphic Voxel-Based Phantoms__
Ilian Häggmark, Kian Shaker, and Hans M. Hertz. under review.



