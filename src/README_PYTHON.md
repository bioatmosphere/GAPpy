# UVAFME Python Translation

This directory contains the Python translation of the UVAFME (University of Virginia Forest Model Enhanced) vegetation model, originally written in Fortran.

## Overview

The Python translation preserves the original model structure while making it more accessible for integration with Python-based scientific computing workflows. All core modules have been translated:

- **constants.py**: Global constants and parameters
- **species.py**: Species data structures and response functions
- **site.py**: Site data structures and climate adjustment functions
- **soil.py**: Soil biogeochemical processes and water balance
- **climate.py**: Climate data processing and conversions
- **uvafme.py**: Main model execution loop and orchestration

## Key Features

### Translated Modules

1. **SpeciesData Class**: Complete translation of species attributes and response functions
   - Temperature response (`temp_rsp`)
   - Drought response (`drought_rsp`) 
   - Light response (`light_rsp`)
   - Fire response (`fire_rsp`)
   - Flood response (`flood_rsp`)

2. **SiteData Class**: Site-level data management
   - Climate data attachment
   - Altitude adjustments
   - Species attachment and filtering

3. **SoilData Class**: Soil processes
   - Decomposition model (`soil_decomp`)
   - Water cycle model (`soil_water`)
   - Three-layer soil structure (A0, A, Base layers)

4. **Climate Functions**: 
   - Monthly to daily conversions (`cov365`, `cov365a`)
   - Extraterrestrial radiation (`ex_rad`)
   - Hargreaves evaporation (`hargrea`)

5. **UVAFMEModel Class**: Main model orchestration
   - Site loop processing
   - Annual simulation cycles
   - Output management

## Usage

### Basic Usage

```python
from vegetation import UVAFMEModel

# Create model instance
model = UVAFMEModel()

# Run with default parameters
model.run()

# Or run with input file list
model.run("input_files.txt")
```

### Integration Example

```python
from vegetation import SpeciesData, SiteData, SoilData

# Create species instance
species = SpeciesData()
species.initialize_species(
    species_id=1,
    genus_name="Quercus",
    taxonomic_name="Quercus alba",
    # ... other parameters
)

# Use species response functions
light_response = species.light_rsp(available_light=0.5)
species.temp_rsp(degree_days=1200)
```

## Implementation Status

### ✅ Completed Components

- Core data structures (Species, Site, Soil)
- Species response functions
- Climate processing functions
- Soil biogeochemical processes
- Basic model execution framework

### 🚧 Partial Implementation

- Tree data structures and growth functions
- Plot-level processes
- Input/output file handling
- Parameter management
- Random number generation

### ❌ Not Yet Implemented

- Complete tree growth and mortality models
- Canopy light interactions
- Detailed biogeochemical climate model
- Full input/output system
- Visualization and analysis tools

## Data Structures

### Species Response Functions

The species module implements all original response functions:

```python
# Temperature response (parabolic)
species.temp_rsp(degree_days)

# Drought response (exponential decay)
species.drought_rsp(dry_days, dry_days_surface)

# Light response (exponential saturation)
light_factor = species.light_rsp(available_light)
```

### Soil Biogeochemistry

Three-layer soil model with carbon and nitrogen cycling:

```python
# Decomposition with temperature and moisture controls
avail_N, C_resp = soil.soil_decomp(
    litter_c1, litter_c2, litter_n1, litter_n2,
    temperature, precipitation, moisture_factors...
)

# Water balance with evapotranspiration
water_results = soil.soil_water(
    slope, lai, lai_water, sigma, freeze,
    rain, potential_evap
)
```

## Key Differences from Fortran

1. **Object-Oriented Design**: Classes replace Fortran derived types
2. **NumPy Arrays**: Replace Fortran arrays for numerical operations
3. **Exception Handling**: Python error handling vs Fortran error codes
4. **Dynamic Memory**: Python lists vs Fortran allocatable arrays
5. **Module System**: Python imports vs Fortran USE statements

## Development Notes

### Preservation of Original Logic

The translation preserves:
- Original mathematical formulations
- Parameter values and constants
- Algorithm flow and decision logic
- Data structure relationships

### Python Enhancements

- Type hints for better code documentation
- Docstrings for all functions and classes
- More readable variable names where appropriate
- Modular design for easier testing and extension

## Testing and Validation

To ensure accuracy, the Python translation should be validated against:
1. Original Fortran model outputs
2. Known test cases
3. Literature benchmarks
4. Cross-validation with other forest models

## Future Development

Priority areas for completion:
1. Complete tree growth and mortality algorithms
2. Implement full I/O system for standard UVAFME file formats
3. Add parameter configuration management
4. Develop testing framework
5. Create visualization and analysis tools
6. Add parallel processing capabilities

## Dependencies

- Python 3.7+
- NumPy
- Additional scientific computing libraries as needed

## License

This translation maintains the same license as the original UVAFME model.