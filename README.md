# jfreg: Joint fMRI Registration

`jfreg` is a Python script that is intended to facilitate analysis of
functional magnetic resonance imaging (fMRI) datasets. `jfreg` provides a
pipeline for motion correction, functional-structural coregistration, and
spatial normalization. The goal of this pipeline is to process the data so that
a spatially normalized version of the functional dataset is available for
further analysis. The pipeline minimizes interpolation of the functional
dataset. `jfreg` uses the FMRI Software Library (FSL). `jfreg` supports the
NIfTI-1.1 data storage format.

# License

Copyright (c) 2019, Jeffrey M. Engelmann

`jfreg` is released under the revised (3-clause) BSD license.
For details, see [LICENSE.txt](LICENSE.txt).

The [FMRI Software Library](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki) is a
comprehensive library of analysis tools for brain imaging data that is
developed by the Analysis Group at the Wellcome Centre for Integrative
Neuroimaging at the University of Oxford. FSL is released for non-commercial
use under an open-source
[license](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Licence).

The NIfTI file format specification is available from the
[Neuroimaging Informatics Technology Initiative](https://nifti.nimh.nih.gov).
