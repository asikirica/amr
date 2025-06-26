# Adaptive Mesh Refinement

This repository hosts the source code for a collection of tools and libraries for Adaptive Mesh Refinement (AMR) in OpenFOAM 10. Current implementation supports both 2D and 3D problems.

This implementation builds upon significant prior work in the field of Adaptive Mesh Refinement for OpenFOAM. The core methodologies and algorithms are based on the following research and implementations:
- https://doi.org/10.1016/j.softx.2019.100317
- https://github.com/ElsevierSoftwareX/SOFTX_2018_143
- https://github.com/HenningScheufler/multiDimAMR


## Key Features

- Support for both 2D and 3D problems
- Dynamic mesh refinement driven by multi-field criteria with geometric constraints
- Dynamic load balancing based on MPI rank loads and exponential moving average


## Main Libraries and Utilities

The following libraries and utilities provide the core functionality for mesh manipulation and load balancing:

- #### libfvMeshTopoChangers2D
  - Handles mesh refinement in 2D
  - Provides **refiner2D** topoChanger
  - Provides **multiFieldRefiner2D** topoChanger

- #### libfvMeshTopoChangers3D
  - Handles mesh refinement in 3D
  - Provides **multiFieldRefiner3D** topoChanger

- #### decomposeParMesh2D
  - Ensures proper **refinementHistory** decomposition in 2D

- #### decomposeParMesh3D
  - Ensures proper **refinementHistory** decomposition in 3D

- #### reconstructParMesh2D
  - Ensures proper **refinementHistory** reconstruction in 2D

- #### reconstructParMesh3D
  - Ensures proper **refinementHistory** reconstruction in 3D

- #### libfvMeshDistributorsMPI
  - Provides **distributorMPI** distributor
  - Handles runtime 2D and 3D mesh redistribution based on MPI-linked rank computational time and memory usage

- #### libfvMeshDistributorsRollingMPI
  - Provides **distributorsRollingMPI** distributor
  - Handles runtime 2D and 3D mesh redistribution based on MPI-linked rank computational time and memory usage
  - Relies on load history to determine whether to refine


## Validation

Validation directory contains 2D and 3D simulation setup data using RANS and LES models.

- #### cylinderBenchmark
  - Validation case based on https://doi.org/10.1007/978-3-322-89849-4_39
  - Results obtained using laminar model for 2D and 3D variant
  - Case is using pimpleFoam solver

- #### hysingBenchmark
  - Validation case based on https://doi.org/10.1002/fld.1934
  - Results obtained using laminar model for 2D and 3D variant
  - Case is using interFoam solver

- #### damBreak
  - Validation case based on https://doi.org/10.1016/j.jcp.2004.12.007
  - Results obtained using kOmegaSST model for 3D variant
  - Case is using interFoam solver

- #### channelFlow
  - Validation case based on https://doi.org/10.1063/1.869966
  - Results obtained using Smagorinsky and dynamic Smagorinsky models
  - Utility **perturbUChannel** by De Villiers, E. (2006) ported to OpenFOAM 10
  - LES model **dynamicSmagorinsky** by Passalacqua, A. (2014) ported to OpenFOAM 10
  - Case is using pimpleFoam solver

- #### squareCylinder
  - Validation case based on https://doi.org/10.1017/S0022112095004435
  - Results obtained using kOmegaSST model for 2D and 3D variant
  - Results obtained using WALE model for conventional and Wall-Modelled LES (WMLES)
  - LES library **libWallModelledLES** by Mukha, T. et al. (2019) ported to OpenFOAM 10
  - Case is using pimpleFoam solver

- #### scalarTransport
  - Validation case based on https://doi.org/10.1115/GT2010-22709
  - Results obtained using kOmegaSST model for 3D variant
  - Results obtained using WALE model for conventional and Wall-Modelled LES (WMLES)
  - Case is using custom **scalarPimpleFoam** solver


## Documentation

Contains doxygen documentation to be deployed locally.
