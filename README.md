# socco-slocum-pipeline

Processing pipeline for Slocum glider data 

## Overview
This repository provides a reproducible, configuration-driven workflow to:
- ingest raw Slocum glider files (L0),
- apply transparent QC and calibration steps (flag-first, no silent deletions),
- produce **OG1.0-aligned** outputs for downstream interoperability,
- prepare data for processing of microstructure turbulence data

## Scope
- L0 ingestion (dbd/ebd, tbd/sbd, caches)
- variable mapping + metadata for [OG1.0](https://oceangliderscommunity.github.io/OG-format-user-manual/OG_Format.html)
 *based on [pyglider](https://pyglider.readthedocs.io/en/latest/)*
- align pressure
- thermal lag corrections
- apply flight model
- align sensors (e.g optode, suna)
- QC flagging and comparison with ship and glider data
- writing standard products (NetCDF/xarray)

*Scope is subject to change*


This repository uses the GEOMAR glider toolbox as part of the processing pipeline.

The idea is initially to create a baseline dataset following the
  [GEOMAR Slocum Processing Toolbox](https://git.geomar.de/open-source/geomar_glider_toolbox/-/blob/main/glider/README.md),
  with the longer-term aim is to develop a single, Python-based processing pipeline.


## Data levels (working definitions)
- **L0**: Raw Slocum glider data as recorded by the platform/instrument files.
- **OG1.0**: OceanGliders Level-1 format target (harmonised variables, metadata, and QC flags).

