#!/usr/bin/env python

__author__ = "Jesus Fuentes"
__version__ = "PyGenometry v2"

import argparse
import numpy as np
import pygen as pg
import playback as pb
import simulation as sim

def welcome_title():
    import os
    os.system('cls' if os.name == 'nt' else 'clear')
    print(""
    "==================================================\n"
    "    _______  __   __  _______  _______  __    _   \n"
    "   |       ||  | |  ||       ||       ||  |  | |  \n"
    "   |    _  ||  |_|  ||    ___||    ___||   |_| |  \n"
    "   |   |_| ||       ||   | __ |   |___ |       |  \n"
    "   |    ___||_     _||   ||  ||    ___||  _    |  \n"
    "   |   |      |   |  |   |_| ||   |___ | | |   |  \n"
    "   |___|      |___|  |_______||_______||_|  |__|  \n"
    "                                                  \n"
    "==================================================\n"
    "   PyGenometry v2.0.21                            \n"
    "   github.com/fuentesigma                         \n"
    "   Let's start the simulation...                  \n"
    "==================================================\n"
    "")

def main(args):
    # /////////////////////////// #
    #np.random.seed(42)
    # /////////////////////////// #
    """
        Command-line arguments/parameters
        C: Initial number of cells
        S: Number of simulation steps
        B: Background option for the simulation
    """
    # Initial number of cells
    N_cells = args.C
    # Number of simulation steps
    n_steps = args.S
    # Background option
    geometry = args.B
    # Load background geometry
    S = np.loadtxt("data/" + str(geometry) + ".txt")

    """ 
        Adjustable parameters for the simulation
        l0: Cell-cell equilibrium distance
        k: Cell-cell adhesion
        a: Cell-cell repulsion
        gamma: Cell-tissue surface tension
        division: Integer number of cell divisions
        death: Integer number of cell deaths
    """
    l0 = 1e-2
    k = 1e-2 * np.ones((N_cells))
    #k = 5e-2 * np.random.normal(0, 1, N_cells)
    a = 1e-2 * np.ones((N_cells))
    gamma = 1 * np.ones((N_cells))
    division, death = 100, 50
    
    """
        Fixed parameters for the simulation
        dt: Time step
        eta: Cell-fluid friction
        D: Cell diffusion coefficient
    """
    # Time step
    dt = 1e-3
    # Cell-fluid friction
    eta = 1e-6
    # Cell diffusion coefficient
    D = 1e-3

    """
        Cell events and simulation steps
    """
    START, END = int(1/4 * n_steps), int(3/4 * n_steps)
    # /////////////////////////// #
    cell_events = pg.create_cell_events(START, END, division, death, n_steps)
    # /////////////////////////// #
    # Filename for the simulation
    f = f"N_cells_{N_cells}_div_death_{division}_{death}_steps_{n_steps}_bg_{geometry}"
    
    # Welcome title
    welcome_title()

    """  Initial conditions from file  """
    import h5py

    filename = "tissue_regeneration.hdf5"

    with h5py.File("data/" + filename, 'r') as fname:
        X0 = fname['cell_positions'][:]
    
    """  //////////////////////////////  """
    
    # Initial cell positions
    #X0 = 0.5 * pg.sphere(N_cells)

    # Check which type of simulation to run based on command-line arguments
    if args.type == 'realtime':
        # Simulation in real time
        morphogenesis = sim.CellSimulation(X0, k, l0, D, eta, gamma, a, S, dt, n_steps, cell_events)
        morphogenesis.run()

    elif args.type == 'offline':
        # Simulation via file
        pg.run_simulation(X0, k, l0, D, eta, gamma, a, S, dt, n_steps, cell_events, "data/" + f)
        print("\n")
        simulation = pb.Evolution(n_steps, f)
        simulation.run()

    elif args.type == 'file':
        # Replay simulation from file
        simulation = pb.Evolution(n_steps, f)
        simulation.run()

    else:
        print("Invalid simulation type. Please choose: 'realtime' or 'offline' or 'file'.")

# Run main function
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run simulation')
    parser.add_argument('-t', '--type', type=str, help='Type of simulation to run: "realtime" or "offline" or "file"')
    parser.add_argument('--C', type=int, help='Initial number of cells')
    parser.add_argument('--S', type=int, help='Number of simulation steps')
    parser.add_argument('--B', type=str, help='Background option for the simulation')
    args = parser.parse_args()
    main(args)
