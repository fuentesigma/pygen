#!/usr/bin/env python

__author__ = "Jesus Fuentes"

# ////////////////////////////////////////////////
import sys
import h5py
import time
import torch
import numpy as np
import visuals as vi
# ////////////////////////////////////////////////
from numba import njit, prange
from scipy.spatial import KDTree
# ////////////////////////////////////////////////

# Generate initial particle positions on a sphere
def initial_positions(N):
    indices = np.arange(0, N, dtype=float) + 0.5
    phi = np.arccos(1 - 2 * indices / N)
    theta = np.pi * (1 + 5 ** 0.5) * indices
    x = np.cos(theta) * np.sin(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(phi)
    return np.column_stack((x, y, z))

# Compute elastic force between two particles
def elastic_force(xi, xj, k_i, l0):
    delta = xi - xj
    distance = np.sqrt((delta ** 2).sum())
    if distance == 0.0:
        return np.zeros_like(delta)
    direction = delta / distance
    force = -k_i * (distance - l0) * direction
    return force

# Compute control term to drive cells towards S
def control_term(xi, S_tree, S_points, control_strength_i):
    distance, idx = S_tree.query(xi)
    if idx == len(S_points) or np.isinf(distance):
        return np.zeros_like(xi)
    closest_point = S_points[idx]
    control_force = control_strength_i * (closest_point - xi)
    return control_force

# Compute the normal vector at a point xi on surface S
def compute_normal_vector(xi, S_tree, S_points):
    # Find the indices of nearest neighbors on S
    distances, idxs = S_tree.query(xi, k=10)
    neighbors = S_points[idxs]
    # Compute covariance matrix of the neighbors
    cov = np.cov(neighbors - xi, rowvar=False)
    # Compute eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    # Normal vector is the eigenvector corresponding to the smallest eigenvalue
    normal_vector = eigenvectors[:, np.argmin(eigenvalues)]
    return normal_vector

# Project a point onto the target surface using the closest point
def geometric_force(xi, S_tree, S_points, gamma_i):
    # Find the closest point on the target surface
    distance, idx = S_tree.query(xi)
    closest_point = S_points[idx]
    # Compute the force towards the closest point
    force = gamma_i * (closest_point - xi)
    return force

# Compute geometric forces and normals
def compute_geometric_forces(X, S_tree, S_points, gamma):
    N = len(X)
    forces = np.zeros_like(X)
    normals = np.zeros_like(X)
    for i in range(N):
        xi = X[i]
        gamma_i = gamma[i]
        # Compute geometric attraction force
        force = geometric_force(xi, S_tree, S_points, gamma_i)
        forces[i] = force
        # Compute normal vector at xi
        n_i = compute_normal_vector(xi, S_tree, S_points)
        normals[i] = n_i
    return forces, normals

# Parallelised elastic force computation function
def compute_elastic_forces_parallel(X, k, l0, neighbors_list):
    N = X.shape[0]
    forces = np.zeros_like(X)
    pairwise_forces = []
    pair_indices = []
    for i in prange(N):
        xi = X[i]
        k_i = k[i]
        force_i = np.zeros(3)
        neighbors = neighbors_list[i]
        for idx in neighbors:
            if i != idx:
                xj = X[idx]
                fij = elastic_force(xi, xj, k_i, l0)
                force_i += fij
                # Store pairwise force and indices
                pairwise_forces.append(fij)
                pair_indices.append((i, idx))
        forces[i] = force_i
    return forces, pairwise_forces, pair_indices

# Generate the cell events, if any, or leave the list empty
def create_cell_events(start, end, n_divisions, n_deaths, n_steps):
    """ 
    ---------------------------
    Define the cell events.
    ---------------------------
    Examples:
    1.  This simulation will start with N cells, 
        add a new cell at time step 10, remove a cell at time step 20 and so on
        cell_events = [(10, "division"), (20, "death"), (30, "division")]

    2. This simulation will keep the number of cells constant
        cell_events = []
    """
    if n_divisions + n_deaths > n_steps:
        raise ValueError("The number of cell events is greater than the number of time steps")
    if n_divisions == 0 and n_deaths == 0:
        cell_events = []
    else:
        cell_events = [(i, "division") for i in np.random.randint(start, end, n_divisions)]
        cell_events += [(i, "death") for i in np.random.randint(start, end, n_deaths)]
        cell_events = sorted(cell_events, key=lambda x: x[0])
    return cell_events

# Handle particle addition and removal events
def handle_cell_events(X, params, cell_events, current_step):
    # Unpack parameters
    k, gamma = params['k'], params['gamma']
    # Process events scheduled for the current step
    events = [event for event in cell_events if event[0] == current_step]
    for _, event_type in events:
        if event_type == "division":
            # Choose a random particle to divide
            idx = np.random.choice(len(X))
            xi = X[idx]
            # Create a new particle with a small random offset
            offset = np.random.normal(scale=0.1, size=xi.shape)
            new_particle = xi + offset
            X = np.vstack([X, new_particle])
            # Assign parameters to the new particle
            k = np.append(k, k[idx])
            gamma = np.append(gamma, gamma[idx])
        elif event_type == "death":
            # Choose a random particle to remove
            idx = np.random.choice(len(X))
            X = np.delete(X, idx, axis=0)
            k = np.delete(k, idx)
            gamma = np.delete(gamma, idx)
    # Update parameters
    params['k'], params['gamma'] = k, gamma
    return X, params

# Adaptive parameter optimisation using PyTorch
def optimise_parameters(X, S_tree, S_points, params, learning_rate=0.01):
    N = len(X)
    device = torch.device('cpu')  # Or 'cuda'

    # Convert parameters to PyTorch tensors with gradients
    k = torch.tensor(params['k'], dtype=torch.float32, requires_grad=True, device=device)
    gamma = torch.tensor(params['gamma'], dtype=torch.float32, requires_grad=True, device=device)
    control_strength = torch.tensor(params['control_strength'], dtype=torch.float32, requires_grad=True, device=device)

    # Convert positions to tensors
    X_tensor = torch.tensor(X, dtype=torch.float32, device=device)

    # Find closest points on S
    distances, idxs = S_tree.query(X)
    closest_points = S_points[idxs]
    closest_points_tensor = torch.tensor(closest_points, dtype=torch.float32, device=device)

    # Compute adherence loss
    adherence_loss = torch.mean((X_tensor - closest_points_tensor).pow(2).sum(dim=1))

    # Regularisation to prevent parameters from growing too large
    reg_loss = torch.mean(k.pow(2) + gamma.pow(2) + control_strength.pow(2))

    # Total loss
    total_loss = adherence_loss + 0.01 * reg_loss

    # Backpropagation
    total_loss.backward()

    # Update parameters using gradient descent
    with torch.no_grad():
        k -= learning_rate * k.grad
        gamma -= learning_rate * gamma.grad
        control_strength -= learning_rate * control_strength.grad

    # Zero gradients for the next iteration
    k.grad.zero_()
    gamma.grad.zero_()
    control_strength.grad.zero_()

    # Update params dictionary by detaching tensors before converting to NumPy arrays
    params['k'] = k.detach().cpu().numpy()
    params['gamma'] = gamma.detach().cpu().numpy()
    params['control_strength'] = control_strength.detach().cpu().numpy()
    return params

# Main simulation function with parallelisation
def run_simulation(X0, S, params, dt, n_steps, cell_events, filename):
    # Initialise positions
    X = X0.copy()
    # Build KDTree for the target surface (does not change over time)
    S_tree = KDTree(S)
    # Create the HDF5 file
    f = h5py.File(filename + '.hdf5', 'w')
    # Create an empty resizable dataset to store the cell positions
    dset = f.create_dataset("cell_positions", (1, X.shape[0], 3), maxshape=(None, None, 3), dtype='f')
    # ////////////////////////////////////////////////////////////
    start_time = time.time()
    # Time integration loop
    for step in range(n_steps):
        # Handle cell events (division and death)
        X, params = handle_cell_events(X, params, cell_events, step)
        # Unpack updated parameters
        k, gamma, control_strength = params['k'], params['gamma'], params['control_strength']
        # Here l0 and D remain constants
        l0, D = params['l0'], params['D']
        # Number of cells at this time step
        N = len(X)
        # Resize the datasets to accommodate the new data
        dset.resize((step + 1, N, 3))
        # Build KDTree for particle interactions
        tree = KDTree(X)
        # Initialise forces
        forces = np.zeros_like(X)
        # Precompute the neighbour indices for all particles
        neighbors_list = tree.query_ball_point(X, r=l0 * 1.5)
        # Compute elastic forces in parallel
        elastic_forces, pairwise_forces, pair_indices = compute_elastic_forces_parallel(X, k, l0, neighbors_list)
        forces += elastic_forces
        # Compute geometric forces and normals
        geometric_forces, normals = compute_geometric_forces(X, S_tree, S, gamma)
        forces += geometric_forces
        # Compute control forces
        control_forces = np.zeros_like(X)
        for i in range(N):
            xi = X[i]
            control_strength_i = control_strength[i]
            control_force = control_term(xi, S_tree, S, control_strength_i)
            control_forces[i] = control_force
        forces += control_forces
        # Add stochastic term for Brownian motion
        stochastic_term = np.random.normal(0, np.sqrt(2 * D * dt), size=X.shape)
        # Update positions
        X += (forces * dt + stochastic_term)
        # Adaptive parameter optimisation using PyTorch
        params = optimise_parameters(X, S_tree, S, params)
        # Store the current positions in the HDF5 file
        dset[step] = X
        # Compute stress tensors
        stress_tensors = np.zeros((N, 3, 3))
        for (i, j), fij in zip(pair_indices, pairwise_forces):
            rij = X[j] - X[i]
            stress_contribution = np.outer(fij, rij)
            stress_tensors[i] += stress_contribution
        # Divide by effective volume (set V_i = 1.0 for simplicity)
        stress_tensors /= 1.0
        # Compute pressures
        pressures = -np.trace(stress_tensors, axis1=1, axis2=2) / 3.0
        # Resize stress and pressure datasets if necessary
        if step == 0:
            stress_dset = f.create_dataset("stress_tensors", (1, N, 3, 3), maxshape=(None, None, 3, 3), dtype='f')
            pressure_dset = f.create_dataset("pressures", (1, N), maxshape=(None, None), dtype='f')
        else:
            stress_dset.resize((step + 1, N, 3, 3))
            pressure_dset.resize((step + 1, N))
        # Store stress tensors and pressures
        stress_dset[step] = stress_tensors
        pressure_dset[step] = pressures
        # ////////////////////////////////////////////////////////////
        # Compute elapsed time and estimated remaining time
        elapsed_time = time.time() - start_time
        progress = (step + 1) / n_steps
        remaining_time = (elapsed_time / progress) - elapsed_time
        # Show progress
        sys.stdout.write('\r')
        # Print the progress bar, elapsed time, and estimated remaining time
        sys.stdout.write(
            "[%-20s] %d%% | Elapsed: %s | Remaining: %s" % (
                '=' * int(20 * progress),
                100 * progress,
                format_time(elapsed_time),
                format_time(remaining_time)
            )
        )
        # Apply a flush to the screen
        sys.stdout.flush()
    # Save the final parameters to the HDF5 file
    param_group = f.create_group('final_parameters')
    param_group.create_dataset('k', data=params['k'])
    param_group.create_dataset('gamma', data=params['gamma'])
    param_group.create_dataset('control_strength', data=params['control_strength'])
    # Close the file
    f.close()
    return

def format_time(seconds):
    hrs, secs = divmod(seconds, 3600)
    mins, secs = divmod(secs, 60)
    return "%02d:%02d:%02d" % (hrs, mins, secs)

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
    "   PyGen                                          \n"
    "   github.com/fuentesigma/pygen                   \n"
    "   Let's start the simulation...                  \n"
    "==================================================\n"
    "")

def main(N_cells=100, n_steps=1000, dt=1e-2, geometry='sphere', division=0, death=0):
    # Define seed for reproducibility
    np.random.seed(42)

    if geometry == 'Drosophila':
        # Load initial cell positions from CSV file
        X0 = np.loadtxt("drosophila/XYZ_Coordinates_T20.csv", delimiter=",")

        # Number of cells
        N_cells = X0.shape[0]

        # Load background geometry from CSV file
        S = np.loadtxt("drosophila/XYZ_Coordinates_T80.csv", delimiter=",")
    else:
        # Initial cell positions
        X0 = 0.5 * initial_positions(N_cells)

        # Load background geometry
        S = np.loadtxt("data/" + str(geometry) + ".txt")

    # Initialise cell parameters with some variability
    k = np.random.uniform(1e-2, 2, size=N_cells)
    gamma = np.random.uniform(1e-2, 2, size=N_cells)
    control_strength = np.random.uniform(1e-2, 2, size=N_cells)
    l0 = 1e-2
    D = 1e-4

    # Pack parameters into a dictionary
    params = {'k': k, 'l0': l0, 'gamma': gamma, 'D': D, 'control_strength': control_strength}

    # Cell events and simulation steps
    cell_events = create_cell_events(0, n_steps, division, death, n_steps)

    # Filename for the simulation
    filename = f"N_cells_{N_cells}_div_death_{division}_{death}_steps_{n_steps}_bg_{geometry}"

    # Welcome title
    welcome_title()

    # Execute simulation
    run_simulation(X0, S, params, dt, n_steps, cell_events, filename="data/" + filename)
    print("\n")

    # Display trajectories
    # Comment these lines out if run on HPC
    simulation = vi.Evolution(n_steps, filename)
    simulation.run()

# Run main function
if __name__ == "__main__":
    main(n_steps=5000, geometry='Drosophila')
    #   //////////////
    #      _____ 
    #     /     \
    #    | () () |
    #     \  ^  /
    #      |||||
    #      |||||
    #
    #   //////////////