# PyGen

PyGen is a dynamic cellular simulation tool designed for modelling cellular interactions within diverse three-dimensional tissue architectures. This tool leverages various mechanical, geometric, and life-event forces to mirror the adaptability of cellular spatial configurations across varied tissue landscapes.

## Requirements

Ensure you have the following packages installed:

- Numpy
- Scipy
- h5py
- imageio
- vispy
- python 3

You can typically install these using pip:

`pip install numpy scipy h5py imageio vispy`


## Codebase

The project comprises four main files:


1. `runme.py`: The primary script for executing simulations.
2. `pygen.py`: The central engine fo the simulation.
3. `simulation.py`: Designed for executing visualisations in real time.
4. `playback.py`: It contains a class for visualising evolution through files.

## Usage

To use PyGen, follow these steps:

1. Clone or download this repository to your local computer.
2. Navigate to the project directory.
3. Run the simulation using the following command:

> python runme.py --type [simulation_type] --C [initial_number_of_cells] --S [time_steps] --B [geometric_background]


### Parameters:

- **simulation_type**: Determines the mode of the simulation. 
    - `realtime`: Visualize the evolution in real time.
    - `offline`: For larger simulations. Records the evolution in a file before visualization.
    - `file`: To visualize pre-recorded simulations.
  
- **initial_number_of_cells**: The initial number of cells in the ensemble.
  
- **time_steps**: Define the number of simulation steps. (Integer value)

- **geometric_background**: Specifies the target geometry of the tissue architecture. Options include:
    - `sphere`
    - `torus`
    - `wfoil`
    - `wknot`
    - `cylinder`
    - `beltrami`
    - `cylinder`

### Tweaking Additional Parameters:

Within `runme.py`, users can adjust several parameters for more refined simulations:

```python
    """ 
        Adjustable parameters for the simulation
    """
    l0 = 1e-2                         # Cell-cell equilibrium distance
    k = 1e-2 * np.ones((N_cells))     # Cell-cell adhesion
    a = 1e-2 * np.ones((N_cells))     # Cell-cell repulsion
    gamma = 1 * np.ones((N_cells))    # Cell-tissue surface tension
    division, death = 100, 50         # Integer values for cell divisions & deaths
```

Adjust these parameters as needed to fit your specific requirements.

For any further queries or issues, feel free to open an issue on this repository or contact the maintainer.

jesus.fuentes [at] uni [dot] Luxembourg ccTLD



