# PyGen

PyGen is a dynamic cellular simulation tool designed for modelling cellular interactions within diverse three-dimensional tissue architectures. This tool leverages various mechanical, geometric, and life-event forces to mirror the adaptability of cellular spatial configurations across varied tissue landscapes.

## Requirements

Ensure you have the following packages installed on a python virtual environment:

- numpy
- scipy
- h5py
- vispy
- torch
- numba

You can typically install these using pip, for instance:

`pip install numpy scipy h5py vispy etc`


## Codebase

The project comprises two main files:

1. `pygen.py`: The primary script for executing simulations.
2. `visuals.py`: Designed for executing visualisations either in real time or offline.

## Usage

To use PyGen, follow these steps:

1. Clone or download this repository to your local computer.
2. Navigate to the project directory.
3. Run the simulation using the following command:

```python
python pygen.py
```

By default, the code runs a simulation with predefined parameters, for example `n_steps=5000` and `geometry='Drosophila'` as shown in the `if __name__ == "__main__":` block.

### Parameters:

```python
main(
    N_cells=200,      # Change the initial number of cells
    n_steps=10000,    # Change the number of simulation steps
    dt=1e-3,          # Adjust the time step
    geometry='sphere',# Choose a different background geometry
    division=100,     # Set how many cell division events to occur
    death=50          # Set how many cell death events to occur
)
```

- **geometric_background**: Specifies the target geometry of the tissue architecture. Options include:
    - `sphere`
    - `torus`
    - `wfoil`
    - `wknot`
    - `cylinder`
    - `beltrami`

Alternatively, you can import or generate your own point clouds utilising the notebook `generate_points.ipynb` within tools.

## Example

To run the simulation with a new configuration, simply edit the `main()` call as described above and execute:

```python
python pygen.py
```

For example, changing `N_cells=500`, `division=100`, and `death=50` would simulate a scenario where 500 cells evolve over time, with 100 division events and 50 death events.

### Maintenance

For any further queries or issues, feel free to open an issue on this repository or contact the maintainer.

jesus [dot] fuentes [at] uni [dot] lu



