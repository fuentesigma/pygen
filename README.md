# PyGen

PyGen is a dynamic cellular simulation tool designed for modelling cellular interactions within diverse three-dimensional tissue architectures. This tool leverages various mechanical, geometric, and life-event forces to mirror the adaptability of cellular spatial configurations across varied tissue landscapes.

## Requirements

Ensure you have the following packages installed:

- numpy
- scipy
- h5py
- imageio
- vispy
- python 3
- torch
- numba

You can typically install these using pip, for instance:

`pip install numpy scipy h5py imageio vispy etc`


## Codebase

The project comprises two main files:

1. `pygen.py`: The primary script for executing simulations.
2. `visuals.py`: Designed for executing visualisations either in real time or offline.

## Usage

To use PyGen, follow these steps:

1. Clone or download this repository to your local computer.
2. Navigate to the project directory.
3. Run the simulation using the following command:

> python pygen.py

By default, the code runs a simulation with predefined parameters, for example n_steps=5000 and geometry='Drosophila' as shown in the if `__name__ == "__main__":` block.

### Parameters:

All parameters are defined within the `main()` function. To adjust them, open the `pygen.py` file and modify the call to `main()` at the bottom of the script:

- **simulation_type**: Determines the mode of the simulation. 
    - `realtime`: Visualize the evolution in real time.
    - `offline`: For larger simulations. Records the evolution in a file before visualization.
    - `file`: To visualize pre-recorded simulations.
  
- **initial_number_of_cells**: The initial number of cells in the ensemble. (Integer value.)
  
- **time_steps**: Define the number of simulation steps. (Integer value.)

- **geometric_background**: Specifies the target geometry of the tissue architecture. Options include:
    - `sphere`
    - `torus`
    - `wfoil`
    - `wknot`
    - `cylinder`
    - `beltrami`
    - `cylinder`

Alternatively, you can import or generate your own point clouds utilising the notebook `generate_points.ipynb` within tools.

### Tweaking Additional Parameters:

Within `runme.py`, users can adjust several parameters for more refined simulations:

```python
    """ 
        Adjustable parameters for the simulation
    """
    a = 1e-2 * np.ones((N_cells))     # Cell-cell repulsion
    k = 1e-2 * np.ones((N_cells))     # Cell-cell adhesion
    l0 = 1e-2                         # Cell-cell equilibrium distance
    gamma = 1.0 * np.ones((N_cells))  # Cell-tissue surface tension
    division = 100                    # Integer values for cell divisions
    death = 50                        # Integer values for cell deaths
```

Adjust these parameters as needed to fit your specific requirements.

## Example

Edit the `runme.py` file, then tweak these parameters:

```python
    a = 1e-2 * np.ones((N_cells))
    k = 1e-2 * np.ones((N_cells))
    l0 = 1e-2
    gamma = 1.0 * np.ones((N_cells))
    division = 100
    death = 50
```

Open a terminal window, go to the folder's project and execute the following command line:

>  python runme.py

This will execute a simulation file containing 500 cells evolving over a foil background. In the simulation, 100 cells will be divided, while 50 of them will be killed.

### Maintenance

For any further queries or issues, feel free to open an issue on this repository or contact the maintainer.

jesus [dot] fuentes [at] uni [dot] lu



