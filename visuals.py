#!/usr/bin/env python

__author__ = "Jesus Fuentes"

# Standard libraries
import sys
import numpy as np

# Third-party libraries
import h5py
from vispy.scene import visuals
from vispy import app, scene

# Constants
COLOR_CELL = "#b38eed"
CANVAS_SIZE = (1200, 960)
CANVAS_BGCOLOR = '#EFEFEF'

def progressbar(j, count, prefix="", size=30, out=sys.stdout):
    x = int(size * j / count)
    if j < count:
        print(f"{prefix}[{'=' * (x + 1)}{'.' * (size - x - 1)}] {j}/{count}", end='\r', file=out, flush=True)
    else:
        print(f"{prefix}[{'=' * size}] {j}/{count}", end='\r', file=out, flush=True)

class Canvas:
    def __init__(self):
        self.canvas = scene.SceneCanvas(
            title='PyGenometry Window',
            keys='interactive',
            show=True,
            bgcolor=CANVAS_BGCOLOR,
            size=CANVAS_SIZE
        )
        self.view = self.canvas.central_widget.add_view()
        self.view.camera = 'turntable'
        self.view.camera.fov = 45
        self.view.camera.distance = 100

        self.scatter = visuals.Markers(
            scaling=True,
            spherical=False,
            symbol='disc',
            alpha=0.8
        )

        self.view.add(self.scatter)
        self.axis = visuals.XYZAxis(parent=self.view.scene)

class Evolution:
    def __init__(self, n_steps, datafile='file'):
        self.canvas = Canvas()
        self.t = 0
        self.n_steps = n_steps
        self.df = 'data/' + datafile + '.hdf5'

        with h5py.File(self.df, 'r') as f:
            self.data = f['cell_positions'][:]

    def simulation(self, t):
        cells = self.data[t]
        idx = np.where(np.linalg.norm(cells, axis=1) < 1e-6)[0]
        cells = np.delete(cells, idx, axis=0)
        return cells

    def update(self, event):
        r = self.simulation(self.t)
        self.canvas.scatter.set_data(
            r,
            face_color=COLOR_CELL,
            size=100,
            edge_color=COLOR_CELL,
            edge_width=5
        )

        progressbar(self.t + 1, self.n_steps, "Simulating...")
        self.t += 1

    def run(self):
        timer = app.Timer(
            interval=0.01,
            connect=self.update,
            start=True,
            iterations=self.n_steps
        )
        self.canvas.canvas.app.run(1)