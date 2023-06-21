#!/usr/bin/env python

__author__ = "Jesus Fuentes"
__version__ = "PyGenometry v2"

import sys
import numpy as np
import pygen as pg

# Vispy ------------------------
import vispy
from vispy.scene import visuals
from vispy import app
# ------------------------------

# Constants
CELL_COLOR  = "#89CFF0"
ARROW_COLOR = "#0098DB"
CANVAS_SIZE = (800, 600)
SAVE_STEP = 100

class CellSimulation:
    """
    Class to run the cell simulation
    """
    def __init__(self, X0, k, l0, D, eta, gamma, a, S, dt, n_steps, cell_events):
        self.init_simulation_parameters(X0, k, l0, D, eta, gamma, a, S, dt, n_steps, cell_events)
        self.init_visualization()

    def init_simulation_parameters(self, X0, k, l0, D, eta, gamma, a, S, dt, n_steps, cell_events):
        """
        Initialise simulation parameters
        """
        self.X = X0
        self.X_prev = self.X - 0.001
        self.cell_event_idx = 0
        self.step = 0
        self.k = k
        self.l0 = l0
        self.D = D
        self.eta = eta
        self.gamma = gamma
        self.a = a
        self.S = S
        self.dt = dt
        self.n_steps = n_steps
        self.cell_events = cell_events
    
    def init_visualization(self):
        """
        Initialise visualisation components
        """
        self.frame = vispy.scene.SceneCanvas(title='PyGenometry Window', keys='interactive', bgcolor='#EFEFEF', show=True)
        self.frame.size = CANVAS_SIZE
        self.pov = self.frame.central_widget.add_view()
        self.pov.camera = 'turntable'

        self.cells = visuals.Markers(scaling=True, spherical=False, symbol='disc', alpha=0.8)
        self.orientation = visuals.Arrow(arrow_type='triangle_30', antialias=True, connect="segments", width=1)
        self.pov.add(self.cells)
        self.pov.add(self.orientation)

        self.pov.camera.fov = 30
        self.pov.camera.distance = 5

        self.axis = visuals.XYZAxis(parent=self.pov.scene)

        self.timer = app.Timer()
        self.timer.connect(self.update)
        self.timer.start(interval=0, iterations=self.n_steps)
        self.frame.events.close.connect(self.stop_simulation)
        
    def stop_simulation(self, event):
        """
        Stop the simulation
        """
        self.timer.stop()
        app.quit()

    def update(self, event):
        """
        Update the simulation state
        """
        self.X, self.X_prev, self.cell_event_idx, self.a, self.k, self.gamma = pg.cell_evolution(
                                                                        self.X, 
                                                                        self.X_prev, 
                                                                        self.k, 
                                                                        self.l0, 
                                                                        self.D, 
                                                                        self.eta, 
                                                                        self.gamma, 
                                                                        self.a,
                                                                        self.S, 
                                                                        self.dt, 
                                                                        self.step,
                                                                        self.cell_event_idx, 
                                                                        self.cell_events)
        self.cells.set_data(self.X, size=400, edge_width=10, face_color=CELL_COLOR, edge_color=CELL_COLOR)
        # Add orientation arrow
        U = self.X / np.linalg.norm(self.X)
        self.orientation.set_data(np.stack((self.X, self.X + U), axis=1), color=ARROW_COLOR)

        # Print progress
        progress = (self.step + 1) / self.n_steps
        sys.stdout.write(f'\r[{"=" * int(20 * progress): <20}] {100 * progress:.0f}%')
        sys.stdout.flush()

        self.step += 1
            
    def run(self):
        """
        Run the simulation
        """
        self.frame.show()
        app.run()