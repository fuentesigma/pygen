#!/usr/bin/env python

__author__ = "Jesus Fuentes"
__version__ = "PyGenometry v2"

import numpy as np
from scipy.stats import gaussian_kde
from scipy.spatial import distance as scidist

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# Generate a sphere of points
def sphere(N):
    i = np.arange(0, N, dtype=float) + 0.5
    phi = np.arccos(1 - 2*i/N)
    theta = np.pi*(1 + 5**0.5)*i
    x, y, z = np.cos(theta)*np.sin(phi), np.sin(theta)*np.sin(phi), np.cos(phi)
    return np.array((x, y, z)).T 

# Compute distances between cells
def distance(xi, xj):
    return np.linalg.norm(xi - xj)

# Neighbours are the nearest cells except the cell itself
def find_nearest_neighbors(X, threshold=None, all_neighbors=False):
    N = len(X)
    # Estimate an appropriate threshold
    if threshold is None:
        threshold = np.mean([distance(X[i], X[j]) for i in range(N) for j in range(i + 1, N)])
    # Initialize neighbors
    neighbors = [set() for _ in range(N)]
    # Find nearest neighbors
    if all_neighbors:
        for i in range(N):
            neighbors[i] = set(range(max(0, i - 1), min(N, i + 2))) - {i}
    else:
        for i in range(N):
            for j in range(i + 1, N):
                if distance(X[i], X[j]) <= threshold:
                    neighbors[i].add(j)
                    neighbors[j].add(i)
    return neighbors

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

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# Elastic force between two cells
def spring_force(xi, xj, k, l0):
    """
    This term serves to model cell-cell adhesion. The formula for the elastic force 
    is given by force = k * (d - l0) * direction, where:
    
    ** k is the spring constant, which determines the strength of the interaction 
    between the cells. A higher k value represents a stronger interaction.
    
    ** d is the actual distance between the cells, calculated by d = distance(xi, xj).
    l0 is the equilibrium distance between the cells. 
    The term d - l0 is the displacement of the "spring" from its equilibrium position.
    
    ** direction is a unit vector pointing from cell j to cell i, 
    calculated by direction = (xi - xj) / d. 
    This vector gives the direction of the force. The force is a vector quantity, 
    so it has both magnitude and direction. 
    
    In this case, the force acts along the line connecting the two cells.
    So, the force exerted on cell i due to cell j is given by k * (d - l0) * direction. 
    This force tends to restore the distance between the cells to the equilibrium 
    distance l0, acting to compress the cells if they're too far apart and to extend 
    them if they're too close together.
    
    One thing to note is that this force is pairwise and acts equally and oppositely 
    on the two cells due to Newton's third law (for every action, there is an equal 
    and opposite reaction). When computing the total force on a cell, we sum over these 
    pairwise forces for all neighboring cells.
    """
    d = distance(xi, xj)
    direction = (xi - xj) / d
    force = -k * (d - l0) * direction
    return force

# Compute damping force
def damping_force(xi, xj, xi_prev, xj_prev, eta, dt):
    """
    The damping force is proportional to the relative velocity of the cells, which can 
    be approximated by the finite difference between the current and previous positions 
    divided by the time step. This is a standard first order finite difference approximation 
    for the derivative. When calculating the damping force, we need to consider the relative
    velocity of the cells, not the absolute. So, the damping force term should take the 
    difference in the displacement of the two cells from the previous time step to the 
    current one. 
    """
    relative_velocity = (xi - xi_prev) - (xj - xj_prev)
    return -eta * relative_velocity / dt

# Compute repulsion forces
def repulsion_force(xi, xj, a):
    """
    This force would act in addition to the current forces acting on the cells.
    a: is a constant that determines the strength of the repulsion 
    direction: is the unit vector pointing from cell j to cell i

    This force acts along the line joining the two cells and is always repulsive, 
    tending to push the cells apart.
    """
    d = distance(xi, xj)
    direction = (xi - xj) / d
    force = a / (d ** 2) * direction
    return force

# Geometric force density
def geometric_force_density(Xi, S, gamma, sigma):
    """
    Compute the geometric force for a point Xi with respect to a set of points S, 
    taking into account the density of the points on S.
    Xi: np.array of shape (3,) representing the position of the cell
    S: np.array of shape (N, 3) representing the positions of the points in the surface
    gamma: float representing the strength of the geometrical force
    sigma: float representing the scale of the Gaussian kernel in the KDE
    """
    # Compute the distances from Xi to all points in S
    distances = scidist.cdist([Xi], S)[0]
    # Find the index of the closest point
    closest_point_idx = np.argmin(distances)
    # Compute the geometric force
    force_closest = gamma * (S[closest_point_idx] - Xi)
    # Compute a density estimate for the points on the surface S
    kde = gaussian_kde(S.T, bw_method=sigma)
    # Compute the density at the position of Xi
    density_Xi = kde(Xi)
    # Compute a force that pushes Xi away from crowded areas
    force_density = -gamma * density_Xi * (Xi - S[closest_point_idx])
    # Combine the two forces
    return force_closest + force_density

# Compute the mean curvature force
def mean_curvature_force(Xi, S, gamma, neighborhood_size):
    """
    Compute a force based on the mean curvature of the surface at a point Xi.
    Xi: np.array of shape (3,) representing the position of the point
    S: np.array of shape (N, 3) representing the positions of the points in the surface
    gamma: float representing the strength of the geometrical force
    neighborhood_size: int representing the number of nearest neighbors to consider
    """
    # Compute the distances from Xi to all points in S
    distances = scidist.cdist([Xi], S)[0]
    # Find the indices of the closest points
    closest_point_idxs = np.argpartition(distances, neighborhood_size)[:neighborhood_size]
    # Define a local coordinate system using the mean of the closest points as the origin
    origin = np.mean(S[closest_point_idxs], axis=0)
    local_coords = S[closest_point_idxs] - origin
    # Compute the covariance matrix of the local coordinates
    C = np.cov(local_coords, rowvar=False)
    # Compute the eigenvalues and eigenvectors of the covariance matrix
    eigenvalues, eigenvectors = np.linalg.eig(C)
    # Find the normal vector, which is the eigenvector associated with the smallest eigenvalue
    normal_vector = eigenvectors[:,np.argmin(eigenvalues)]
    # Compute the mean curvature
    mean_curvature = np.sum((local_coords @ normal_vector) ** 2) / neighborhood_size
    # Compute the force based on the mean curvature
    force_curvature = -gamma * mean_curvature * normal_vector
    return force_curvature

# Compute the force that pushes Xi towards the surface S
def force_to_surface(Xi, S, gamma):
    """
    Compute a force that pushes Xi towards a point on the surface S 
    that is in the same direction from the center of the sphere.
    Xi: np.array of shape (3,) representing the position of the point
    S: np.array of shape (N, 3) representing the positions of the points in the surface
    gamma: float representing the strength of the geometrical force
    """
    # Compute the unit vector from the center of the sphere to Xi
    direction_Xi = Xi / np.linalg.norm(Xi)
    # Compute the unit vectors from the center of the sphere to each point in S
    directions_S = np.divide(S, np.linalg.norm(S, axis=1, keepdims=True), 
                             out=np.zeros_like(S), where=np.linalg.norm(S, axis=1, keepdims=True)!=0, 
                             dtype=np.float64)
    #directions_S = S / np.linalg.norm(S, axis=1, keepdims=True)
    # Compute the dot product of direction_Xi with each direction in directions_S
    dot_products = directions_S @ direction_Xi
    # Find the index of the point in S that has the highest dot product with direction_Xi
    closest_point_idx = np.argmax(dot_products)
    # Compute the force that pushes Xi towards the closest point in the same direction
    force = gamma * (S[closest_point_idx] - Xi)
    return force

# Compute neighbouring force
def neighbouring_force(X, X_prev, neighbors, A_opt=5, alpha=1e-4):
    """
    Compute the topological force on each cell.
    * X is a 2D array of cell positions at time t.
    * X_prev is a 2D array of cell positions at time t-1.
    * neighbors is a list of lists, where neighbors[i] is a list of the indices of cell i's neighbors.
    * A_opt is the optimal number of neighbors for each cell.
    * alpha is a constant force that determines the strength of the penalty for deviating from the optimal 
    number of neighbours. It increases the force whenever a cell has too few or too many 
    neighbours. In the model, cells "want" to minimise their energy, so they will move 
    in a way that brings them closer to the optimal number of neighbours. If a cell has 
    too few neighbours, it is "encouraged" to move closer to other cells. Conversely, 
    if a cell has too many neighbours, it is "encouraged" to move away from some of them.
    """
    N = len(X)
    # Prepare an array to store the forces
    F = np.zeros_like(X)
    # For each cell
    for i in range(N):
        # Compute the difference between the optimal number of neighbors and the actual number of neighbors
        delta_A = A_opt - len(neighbors[i])
        # Compute the vector from cell i's position at time t-1 to its position at time t
        r_i = X[i] - X_prev[i]
        # Add the contribution of cell i's change in position to the force on cell i
        F[i] = -alpha * np.sign(delta_A) * r_i
    # Return the forces
    return F

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# Compute forces on all cells
def compute_forces(X, X_prev, k, l0, eta, gamma, a, S, dt, neighbors):
    """
    The model includes the following forces:
    * Elastic force
    * Damping force
    * Repulsive force
    * Geometrical forces
    """
    N = len(X)
    forces = np.zeros_like(X)
    # Bandwith for the density estimate
    # Check the number of points in the surface S
    M = S.shape[0]
    # Check the dimensionality of the space
    d = S.shape[1]
    # Compute the bandwidth using Scott's rule
    sigma = M**(-1. / (d + 4))
    # Adjust surface tension
    gamma = 10 * gamma
    # Pairwise elastic forces and constraints
    for i in range(N):
        for j in neighbors[i]:
            # Elastic force
            elastic = spring_force(X[i], X[j], k[i], l0)
            # Damping force
            damping = damping_force(X[i], X[j], X_prev[i], X_prev[j], eta, dt)
            # Repulsive force
            repulsive = repulsion_force(X[i], X[j], a[i])
            # Add the forces to the total force
            forces[i] += elastic + damping + repulsive
        """
        Constraints
        """
        # 1. Background constraints
        geometry = geometric_force_density(X[i], S, gamma[i], sigma)
        surface = force_to_surface(X[i], S, 0.5 * gamma[i])
        curvature = mean_curvature_force(X[i], S, 5.0 * gamma[i], 6)
        # 2. Connectivity corrections
        # connectivity = neighbouring_force(X, X_prev, neighbors)
        # Compute the total force
        forces[i] += geometry + surface + curvature
    # Denoise the forces
    forces = denoise_fields(forces)
    return forces

def denoise_fields(points):
    # Calculate the distances between cells and the mean distance
    distances = scidist.squareform(scidist.pdist(points))
    mean_distance = np.mean(distances, axis=1)
    # Calculate the mean distance
    m_mean_distance = np.max(mean_distance)
    # Correct the cell positions
    points_corrected = points * m_mean_distance / mean_distance[:, np.newaxis]
    # Return corrected cell positions
    return points_corrected

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# Brownian dynamics
def cell_evolution(X, X_prev, k, l0, D, eta, gamma, a, S, dt, step, cell_event_idx, cell_events):
    """
    The stochastic force term is modeled as a Gaussian white noise with mean 0 and a
    variance proportional to the diffusion coefficient (D) and the time step size (dt).
    The factor of 2 in the variance comes from the fluctuation-dissipation theorem,
    which relates the strength of the random forces to the dissipation of energy in
    the system. In this case, it ensures that the simulated system reaches the correct
    equilibrium distribution at a given temperature.
    """
    # Set the division offset. It is estimated based on empirical data or 
    # through trial-and-error to match the expected biological behaviour
    division_offset = 0.1
    # Set an initial seed for reproducibility
    #np.random.seed(42)
    # Check size of the input
    N, d = X.shape
    # Add or remove cells from the simulation
    while cell_event_idx < len(cell_events) and cell_events[cell_event_idx][0] == step:
        _, event_type = cell_events[cell_event_idx]
        if event_type == "division":
            # Pick a random cell to divide
            idx = np.random.choice(X.shape[0])
            # Generate a small random offset
            offset = np.random.randn(1, d) * division_offset
            # Add a new cell at the position of the existing cell plus the offset
            X = np.vstack([X, X[idx] + offset])
            # Add a new phantom cell at the previous step
            X_prev = np.vstack([X_prev, np.zeros((1, d))])
            # Do the same for the rest of the parameters
            a = np.hstack([a, a[idx]])
            k = np.hstack([k, k[idx]])
            gamma = np.hstack([gamma, gamma[idx]])
        elif event_type == "death":
            # Remove a random cell
            idx = np.random.choice(X.shape[0])
            # Remove the cell from the simulation
            X = np.delete(X, idx, axis=0)
            # Remove the phantom cell from the previous step
            X_prev = np.delete(X_prev, idx, axis=0)
            # Do the same for the rest of the parameters
            a = np.delete(a, idx, axis=0)
            k = np.delete(k, idx, axis=0)
            gamma = np.delete(gamma, idx, axis=0)
        # Update the number of cells
        N = X.shape[0]
        # And increment the event index
        cell_event_idx += 1
    # Find the nearest neighbors based on the current positions
    neighbors = find_nearest_neighbors(X)
    # Deterministic forces
    F = compute_forces(X, X_prev, k, l0, eta, gamma, a, S, dt, neighbors)
    # Stochastic forces
    R = np.random.normal(0, np.sqrt(2*D*dt), size=(N, d))
    # Update positions
    X_prev = np.copy(X)
    X += (F + R) * dt
    # Return the updated positions and the cell event index
    return X, X_prev, cell_event_idx, a, k, gamma

# Define a function to iterate cell_evolution() and save the output in an HDF5 file
def run_simulation(X0, k, l0, D, eta, gamma, a, S, dt, n_steps, cell_events, filename):
    import h5py, sys
    # Initial conditions
    # Sphere
    X = X0
    X_prev = X - 0.001
    # Number of cells
    N_cells = X.shape[0]
    # Initialise the cell event index
    cell_event_idx = 0
    step = 0
    # Create the HDF5 file
    f = h5py.File(filename + '.hdf5', 'w')
    # Create an empty resizable dataset to store the cell positions
    maxshape = (None, None, 3)  # Maximum shape: (unlimited steps, unlimited cells, 3 coordinates)
    dset = f.create_dataset("cell_positions", (1, N_cells, 3), maxshape=maxshape, dtype='f')
    for step in range(n_steps):
        # Compute cell positions
        X, X_prev, cell_event_idx, a, k, gamma  = cell_evolution(
                                                    X, 
                                                    X_prev, 
                                                    k, 
                                                    l0, 
                                                    D, 
                                                    eta, 
                                                    gamma, 
                                                    a,  
                                                    S, 
                                                    dt, 
                                                    step, 
                                                    cell_event_idx, 
                                                    cell_events)
        # Resize the dataset to accommodate the new data
        dset.resize((step + 1, X.shape[0], 3))
        # Store the current positions in the HDF5 file
        dset[step] = X
        # Show progress of the computation with a bar
        sys.stdout.write('\r')
        # Print the progress bar
        sys.stdout.write("[%-20s] %d%%" % ('='*int(20*(step+1)/n_steps), 100*(step + 1)/n_steps))
        # Flush the output
        sys.stdout.flush()
    # Close the file
    f.close()