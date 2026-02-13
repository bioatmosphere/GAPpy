# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**GAPpy** (formerly UVAFME-VOC) is a Python translation of the UVAFME (University of Virginia Forest Model Enhanced) individual-based forest gap model. This is a sophisticated ecosystem model that simulates forest composition dynamics, tree growth, mortality, and biogeochemical cycles at the individual tree level.

The project is transitioning from Fortran (src_fortran/) to Python (src/), maintaining the same scientific algorithms while making the model more accessible for modern scientific computing workflows.

## Running the Model

### Basic execution
```bash
python main.py
```

### With custom input file list
```bash
python main.py input_data/file_list.txt
```

### Expected directory structure
- `input_data/`: All model input files (species parameters, climate, site data)
- `output_data/`: Where simulation outputs are written
- `src/`: Python source modules
- `src_fortran/`: Original Fortran implementation (reference)

## Architecture

### Core Hierarchical Structure

The model follows a nested hierarchy from ecosystem to individual:

```
UVAFMEModel (main orchestrator)
  └── SiteData (geographic location with climate)
       └── PlotData (spatial unit, typically ~0.1 ha)
            └── TreeData (individual tree)
```

Each level maintains its own state and delegates to the next level down.

### Main Components

**uvafme.py**: Top-level orchestrator
- Initializes input files, parameters, and species data
- Manages the main site loop and annual simulation cycles
- Coordinates output writing
- Entry point: `UVAFMEModel.run()`

**model.py**: Core forest dynamics engine (`ForestModel` class)
- `bio_geo_climate()`: Processes daily climate variables, calculates degree days, and runs soil biogeochemistry with daily timesteps
- `canopy()`: Calculates light distribution using Beer-Lambert law with LAI
- `growth()`: Three-loop tree growth algorithm (calculate stress → apply growth → check mortality)
- `mortality()`: Individual tree and plot-level disturbance (fire/wind)
- `renewal()`: Seedbank dynamics and tree recruitment based on environmental conditions
- `run_annual_cycle()`: Orchestrates the full annual sequence

**Key data structures**:
- `species.py`: Species traits and response functions (temperature, drought, light, fire, flood)
- `site.py`: Site-level data with climate and multiple plots
- `plot.py`: Plot-level container for trees, light arrays, seedbank/seedling tracking
- `tree.py`: Individual tree extending SpeciesData with growth state
- `soil.py`: Three-layer soil model (`soil_water()`, `soil_decomp()`)

**Supporting modules**:
- `climate.py`: Climate data conversions (monthly→daily, evapotranspiration)
- `parameters.py`: Global simulation parameters (singleton `params` object)
- `constants.py`: Physical/ecological constants
- `io_utils.py`: Input/output file management
- `random_utils.py`: Random number generation for stochastic processes

### Critical Implementation Details

**Species Response Functions** (`species.py`):
All environmental responses are multiplicative factors (0-1) calculated per species:
- `temp_rsp()`: Parabolic function of growing degree days
- `drought_rsp()`: Exponential decay with dry days
- `flood_rsp()`: Flood tolerance response
- `light_rsp()`: Exponential saturation curve
- `poor_soil_rsp()`: Nutrient limitation response
- `fire_rsp()`: Fire tolerance multiplier

These factors are combined multiplicatively in growth calculations to determine realized growth rates.

**Tree Growth Algorithm** (model.py:386-648):
The growth() function uses three sequential loops over all trees:
1. Calculate environmental stress factors and potential growth
2. Apply nitrogen limitation and determine final diameter increment
3. Adjust canopy height and calculate litter inputs

This structure is critical - the loops cannot be combined because nitrogen demand must be calculated for ALL trees before determining nitrogen limitation factors.

**Light Competition** (model.py:298-384):
Light is distributed across height layers using cumulative LAI:
- Deciduous trees contribute fully to `dec_light` array, 80% to `con_light`
- Coniferous trees contribute equally to both arrays
- Light availability calculated top-down using Beer's Law: `exp(-0.4 * cumulative_LAI / plotsize)`

**Recruitment Process** (model.py:792-982):
Complex seedbank → seedling → tree pipeline:
- Seedbank updated annually with seed production + sprouting
- Seedlings established when environmental conditions exceed threshold
- Trees recruited from seedling pool based on probability weighted by growth capacity and light
- First recruitment cycle (seedling_number == 0) has special initialization logic

**Daily Biogeochemistry Loop** (model.py:230-290):
Critical that soil processes run with daily timestep even though model is annual:
- Monthly climate converted to 365 daily values
- Each day: calculate PET → soil water balance → decomposition
- Annual totals accumulated for nitrogen availability

## Development Patterns

### Adding New Species Parameters

1. Add parameter to `species.py`: Update `__init__()` and `initialize_species()`
2. Update input reader in `input_module.py`: Modify CSV parsing
3. If affects growth: Update relevant response function or add to stress calculation in `model.py`

### Modifying Growth Algorithms

Always preserve the three-loop structure in `growth()`. If modifying:
- Loop 1 (lines 428-503): Environmental stress and potential growth
- Loop 2 (lines 522-594): Nitrogen-limited actual growth and mortality
- Loop 3 (lines 596-637): Canopy adjustment and litter

### Testing Strategy

Compare outputs against Fortran reference:
1. Use identical input files from `input_data/`
2. Set `fixed_seed=True` in parameters for reproducibility
3. Compare annual biomass, NPP, soil C/N, and species composition
4. Key validation outputs: BiomC, available N, degree days, light arrays

### Common Pitfalls

- **Tree/Species data synchronization**: Trees inherit species data but maintain independent state. Always call `tree.update_tree(species)` before using species-specific factors in loops.
- **Array indexing**: Height arrays are 0-indexed but ecological height values are 1-based. Use `min(int(height), maxheight-1)` consistently.
- **Unit conversions**: Biomass calculations often require `HEC_TO_M2 / plotsize / numplots` conversion between tree-level (kg) and ecosystem-level (t/ha) units.
- **Mortality markers**: Growth loop marks mortality but doesn't remove trees. Actual removal happens in `mortality()`.
- **Seedling scaling**: Seedlings are scaled by plotsize during recruitment, then scaled back afterward. Don't modify this unless you understand the full recruitment algorithm.

## Model Configuration

Key parameters in `parameters.py` and `input_data/uvafme_config.json`:
- `numyears`: Simulation length (typically 100-500 years)
- `numplots`: Spatial replication (usually 5-20)
- `plotsize`: Area in m² (typically 1000-2000)
- `maxtrees`: Maximum trees per plot (typically 500-2000)
- `maxheight`: Maximum tree height in meters (typically 60-80)
- `fixed_seed`: Set True for reproducible runs
- `with_clim_change`: Enable climate change scenarios

## Scientific Background

This is a gap model in the JABOWA/FORET family:
- Individual tree resolution with size-dependent growth
- Annual timestep for demographics, daily for biogeochemistry
- Gap-scale spatial resolution (forest canopy gaps ~0.1 ha)
- Multiplicative environmental limitations on potential growth

Key publications:
- Wang et al. 2017, Ecological Modelling: Original VOC model
- Shugart et al. 2018, Environmental Research Letters: Gap model review
- Zhang et al. 2021, Forests: Functional group aggregation

## Input Data Requirements

Required CSV files in `input_data/`:
- `UVAFME2012_specieslist.csv`: Species trait parameters
- `UVAFME2012_site.csv`: Site characteristics (slope, elevation, soil)
- `UVAFME2012_climate.csv`: Monthly temperature and precipitation
- `UVAFME2012_climate_stddev.csv`: Climate variability
- `UVAFME2012_sitelist.txt`: List of sites to simulate
- `uvafme_config.json`: Runtime parameters

## Dependencies

Python 3.9+ required. Key packages (see pyproject.toml):
- numpy: Array operations and numerical calculations
- pandas: Input data parsing
- matplotlib, seaborn: Visualization (for analysis scripts)
- scipy: Statistical functions

Note: The project uses `uv` for dependency management based on pyproject.toml.
