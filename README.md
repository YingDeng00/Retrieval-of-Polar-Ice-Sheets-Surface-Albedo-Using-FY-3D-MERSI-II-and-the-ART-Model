## Overview

This repository contains the Python code for retrieving shortwave broadband planar albedo (BBA) over the Greenland Ice Sheet (GrIS) and Antarctic Ice Sheet (AIS) using observations from the Fengyun-3D (FY-3D) Medium Resolution Spectral Imager II (MERSI-II) and the Asymptotic Radiative Transfer (ART) model.

The study adapts the ART-based snow and ice albedo retrieval approach to FY-3D MERSI-II observations by selecting six spectral bands according to their spectral response characteristics. The resulting retrievals are evaluated using independent ground-based observations, the MODIS MCD43A3 albedo product, and Sentinel-3 Ocean and Land Colour Instrument (OLCI) observations.

The repository provides code for MERSI-II data preprocessing, ice optical-property interpolation, ART-based albedo retrieval, and the calculation of shortwave broadband planar albedo.

## Workflow

The retrieval workflow consists of the following main steps:

1. **Data preprocessing:** Extract MERSI-II top-of-atmosphere (TOA) reflectance and observation geometry information from the original HDF files and generate georeferenced raster datasets.
2. **Spectral adaptation:** Select six MERSI-II bands to represent the spectral inputs required by the ART-based retrieval method.
3. **Ice optical properties:** Interpolate the real and imaginary parts of the ice refractive index and calculate the absorption coefficient at the relevant wavelengths.
4. **ART-based retrieval:** Apply the ART model to retrieve  spectral albedo from the preprocessed satellite observations.
5. **Broadband conversion:** Integrate the retrieved spectral information to estimate shortwave broadband planar albedo.
6. **Validation:** Compare the retrieval results with ground-based observations, MODIS MCD43A3, and Sentinel-3 OLCI-derived albedo.

## Repository Structure

```text
Retrieval-of-Polar-Ice-Sheets-Surface-Albedo-Using-FY-3D-MERSI-II-and-the-ART-Model/
│
├── DataPreprocessing.py
│   └── Preprocess FY-3D MERSI-II HDF data and generate georeferenced GeoTIFF files.
│
├── AlbedoRetrieval.py
│   └── Retrieve surface albedo using the ART model and calculate broadband planar albedo.
│
├── ice_refractive_interpolation.py
│   └── Interpolate ice refractive indices and calculate absorption coefficients.
│
└── Figure&Data/
    └── Contains plotting scripts and research-related data used in this study.
```

## Data Requirements

### FY-3D MERSI-II

The primary input consists of FY-3D MERSI-II observations, including TOA reflectance and geolocation information.

The preprocessing workflow uses the MERSI-II multispectral and geolocation HDF files. The relevant inputs include spectral reflectance, solar zenith angle (SZA), sensor viewing zenith angle (VZA), solar azimuth angle (SAA), and sensor azimuth angle (VAA).

Data can be obtained from the [National Satellite Meteorological Center (NSMC)](https://satellite.nsmc.org.cn/).

The retrieval code requires appropriately preprocessed and cloud-screened observations. 

### Validation Data

The study uses the following independent datasets for validation and comparison:

- **PROMICE automatic weather stations (AWS):** Ground-based surface albedo observations over Greenland. Data are available from [PROMICE](https://promice.org/).
- **MODIS MCD43A3:** The MODIS albedo product used for spatial comparison and validation. Data are available through [NASA Earthdata](https://www.earthdata.nasa.gov/).
- **Sentinel-3 OLCI:** Satellite observations used to evaluate the adaptation of the ART-based retrieval approach to MERSI-II.

Additional ground observations over Antarctica, including measurements from Dome C and Dome A, are used in the validation analysis described in the associated study.

## Scripts

### 1. `DataPreprocessing.py`

This script preprocesses FY-3D MERSI-II observations.

Its main functions include:

- Reading MERSI-II multispectral and geolocation HDF files.
- Extracting the required TOA reflectance bands.
- Resampling reflectance and observation geometry variables to a common grid.
- Generating georeferenced GeoTIFF files for subsequent albedo retrieval.


### 2. `ice_refractive_interpolation.py`

This module provides the optical properties of ice required by the ART retrieval.

It includes:

- Spectral data for the real part of the refractive index, \(n\).
- Spectral data for the imaginary part of the refractive index, \(k\).
- Interpolation of \(n\) and \(k\) at MERSI-II specified wavelengths.
- Calculation of the ice absorption coefficient.


### 3. `AlbedoRetrieval.py`

This script implements the ART-based albedo retrieval and broadband conversion.

Its main functions include:

- Reading the preprocessed reflectance and observation geometry data.
- Applying the ART-based retrieval equations.
- Constructing the spectral albedo vector for the selected MERSI-II bands.
- Calculating shortwave broadband planar albedo.
- Extracting and processing observations at selected validation sites.



## MERSI-II Spectral Band Configuration

Six MERSI-II bands were selected to adapt the ART-based retrieval method to FY-3D observations.

| MERSI-II band | Central wavelength (μm) | Corresponding OLCI band |
|---|---:|---|
| B8 | 0.412 | OLCI B1 |
| B11 | 0.555 | OLCI B6 |
| B13 | 0.709 | OLCI B11 |
| B14 | 0.746 | OLCI B12 |
| B15 | 0.865 | OLCI B17 |
| B19 | 1.030 | OLCI B21 |

This configuration is specific to the MERSI-II implementation described in the associated study.

## Retrieval Configuration

The retrieval uses region-specific observation geometry filtering and atmospheric assumptions.

- **Greenland:** Solar zenith angle below 75°.
- **Antarctica:** Solar zenith angle below 80°.
- **Aerosol optical thickness at 500 nm:** 0.07.
- **Ångström exponent:** 1.3.

Other model settings, screening thresholds, and site-specific parameters are defined in the Python scripts. 
