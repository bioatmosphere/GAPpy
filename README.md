# GAPpy

**Individual-based forest gap model in Python**

GAPpy is a Python implementation of the UVAFME (University of Virginia Forest Model Enhanced) individual-based forest gap model, originally developed in Fortran. This translation makes the model accessible for integration with Python-based scientific computing workflows, while preserving the original model's scientific algorithms.

## Overview

The model simulates:
- Individual tree growth, mortality, and recruitment
- Species competition for light, water, and nutrients
- Soil biogeochemistry (carbon and nitrogen cycling)
- Response to climate, disturbance (fire, wind), and environmental stress
- Forest succession dynamics over decades to centuries

## Quick Start

### Installation

**Requirements**: Python 3.9 or higher

**Using uv (recommended)**:

```bash
# Install uv if you don't have it
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clone the repository
git clone https://github.com/yourusername/GAPpy
cd GAPpy

# Install dependencies
uv sync
```

**Using pip**:

```bash
pip install -r requirements.txt  # or install from pyproject.toml
```

### Running the Model

**Basic run with default parameters**:

```bash
uv run python main.py
```

**With custom input file list**:

```bash
uv run python main.py input_data/file_list.txt
```

**Expected output**: The model will simulate forest dynamics, displaying progress every 10 years with statistics on live/dead trees and biomass carbon.

## Directory Structure

```
GAPpy/
├── src/                    # Python source code
│   ├── __init__.py
│   ├── uvafme.py          # Main model orchestrator
│   ├── model.py           # Core forest dynamics engine
│   ├── species.py         # Species traits and response functions
│   ├── tree.py            # Individual tree data structure
│   ├── plot.py            # Plot-level container
│   ├── site.py            # Site data with climate
│   ├── soil.py            # Soil biogeochemistry
│   ├── climate.py         # Climate processing
│   └── ...                # Additional support modules
├── src_fortran/           # Original Fortran implementation (reference)
├── input_data/            # Model input files
│   ├── UVAFME2012_specieslist.csv
│   ├── UVAFME2012_site.csv
│   ├── UVAFME2012_climate.csv
│   └── uvafme_config.json
├── output_data/           # Simulation outputs
├── main.py                # Entry point
├── plot_outputs.py        # Output visualization
├── pyproject.toml         # Python dependencies
└── README.md
```

## Model Components

1. **Individual tree resolution**: Each tree tracked with size, age, species identity
2. **Species response functions**: Temperature, drought, flood, light, fire, nutrient responses
3. **Daily biogeochemistry**: Soil water balance and decomposition on daily timestep
4. **Annual demographics**: Growth, mortality, and recruitment on annual timestep
5. **Spatial structure**: Multiple plots per site for landscape heterogeneity

### Configuration

Key parameters (set in `input_data/uvafme_config.json`):
- `numyears`: Simulation length (default: 1000)
- `numplots`: Number of spatial replicates (default: 200)
- `plotsize`: Plot area in m2 (default: 500)
- `maxtrees`: Maximum trees per plot (default: 1000)
- `fixed_seed`: Reproducible random numbers (default: true)

## Running the Fortran Version

The original Fortran implementation is preserved in `src_fortran/` for reference and validation:

```bash
cd src_fortran
make UVAFME.exe
mv UVAFME.exe ..
cd ..
./UVAFME.exe file_list.txt
```

**Note**: Requires a Fortran compiler (e.g., gfortran, ifort)

## Development

### Testing Against Fortran Reference

To validate the Python translation:

1. Run both versions with identical input files
2. Set `fixed_seed=True` for reproducibility
3. Compare annual outputs: biomass C/N, NPP, soil pools, species composition

### Adding New Features

See `CLAUDE.md` for detailed development guidance, including:
- Model architecture and data flow
- Critical implementation details (3-loop growth algorithm, light competition, recruitment)
- Common pitfalls and unit conversion patterns

## Model Outputs

Simulation outputs are written to `output_data/` and include:

- **Site-level data**: Annual climate, soil carbon/nitrogen pools, NPP
- **Species-level data**: Biomass, basal area, and abundance by species
- **Genus-level data**: Aggregated metrics by genus
- **Tree-level data** (optional): Individual tree attributes for detailed analysis

Output format and frequency can be configured via parameters (`year_print_interval`, `tree_level_data`).

## Publications

**Wang, B.**, Shugart, H. H., & Lerdau, M. T. (2017). [An individual-based model of forest volatile organic compound emissions--UVAFME-VOC v1.0](https://doi.org/10.1016/j.ecolmodel.2017.02.006). **Ecological Modelling**, 350, 69-78.

**Wang, B.**, Shugart, H. H., Shuman, J. K., & Lerdau, M. T. (2016). [Forests and ozone: productivity, carbon storage, and feedbacks](https://www.nature.com/articles/srep22133). **Scientific Reports**, 6, 22133.

**Wang, B.**, Shuman, J., Shugart, H. H., & Lerdau, M. T. (2018). [Biodiversity matters in feedbacks between climate change and air quality: a study using an individual-based model](https://doi.org/10.1002/eap.1721). **Ecological Applications**.

Shugart, H. H., **Wang, B.**, Fischer, R., Ma, J., Fang, J., Yan, X., ... & Armstrong, A. H. (2018). [Gap models and their individual-based relatives in the assessment of the consequences of global change](https://doi.org/10.1088/1748-9326/aaaacc). **Environmental Research Letters**, 13, 033001.

**Wang, B.**, Shugart, H.H. & Lerdau, M.T. (2019). [Complexities between plants and the atmosphere](https://doi.org/10.1038/s41561-019-0413-8). **Nature Geoscience** 12, 693-694

H. H. Shugart, Adrianna Foster, **Bin Wang**, Dan Druckenbrod, Jianyong Ma, Manuel Lerdau, Sassan Saatchi, Xi Yang & Xiaodong Yan. (2020). [Gap models across micro- to mega-scales of time and space: examples of Tansley's ecosystem concept](https://doi.org/10.1186/s40663-020-00225-4). **Forest Ecosystems** 7, 14.

Zhang, H.; Shugart, H.H.; **Wang, B.**; Lerdau, M. [The Significance of Aggregation Methods in Functional Group Modeling](https://doi.org/10.3390/f12111560). **Forests** 2021, 12, 1560.

## Citation

If you use GAPpy in your research, please cite:

```bibtex
@article{wang2017individual,
  title={An individual-based model of forest volatile organic compound emissions---UVAFME-VOC v1.0},
  author={Wang, Bin and Shugart, Herman H and Lerdau, Manuel T},
  journal={Ecological Modelling},
  volume={350},
  pages={69--78},
  year={2017},
  publisher={Elsevier},
  doi={10.1016/j.ecolmodel.2017.02.006}
}
```
