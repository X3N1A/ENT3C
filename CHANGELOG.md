# Changelog

## [2.1.2] - 2025-08-04
### Added
- sorted legend and colormap reverse
- promt user for recomputing the entropy output table
- append output table after each resolution/chr/file instead of saving at end
- print version ```ENT3C --version```
### Changed 
- version to match git hub


## [2.0.7] - 2025-08-01
### Added
- color schemes for BRs of same cell type
### Changed 
- linux executable command now: ENT3C_exe to avoid confusion with python CLI
- minor fixes: 
 - for API: return ENT3C_OUT and Similarity when calling ENT3C.run_get_functions
 - in case double underscore in config file: str.split("_BR").str[0] in get_similarity_table.py

## [2.0.6] - 2025-07-22
### Changed
- could not get a conda build to work due to dependency clashes. 
- Support for Python>=3.11
- Current minimum version requirements:
dependencies = [
  "cooler>=0.10.3",
  "matplotlib>=3.10.0",
  "numpy>=2.3.0",
  "pandas>=2.3.0"
    ]

## [2.0.5] - 2025-07-22
### Added
- Support for Python 3.10 added to enable conda build.

### Changed
- Relaxed minimum version requirements further for to allow broader compatibility.


## [2.0.4] - 2025-07-22
### Added
- Support for Python 3.11 and 3.12

### Changed
- Relaxed minimum version requirements for to allow broader compatibility.

## [2.0.3] - 2025-07-21
### Added
- Initial stable release


