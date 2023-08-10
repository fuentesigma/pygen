#!/usr/bin/env python

__author__ = "Jesus Fuentes"
__version__ = "PyGenometry v2"

# Standard libraries
import sys
import numpy as np

# Third-party libraries
import h5py
import imageio
import vispy
from vispy.scene import visuals
from vispy import app

# Constants
COLOR_CELL = "#89CFF0"
ARROW_COLOR = "#0098DB"
CANVAS_SIZE = (800, 600)
CANVAS_BGCOLOR = '#EFEFEF'

class Canvas:
    def __init__(self):
        self.canvas = vispy.scene.SceneCanvas(title='PyGenometry Window', 
                                              keys='interactive', 
                                              show=True, 
                                              bgcolor=CANVAS_BGCOLOR, 
                                              size=CANVAS_SIZE)
        self.view = self.canvas.central_widget.add_view()
        self.view.camera = 'turntable'
        self.view.camera.fov = 45
        self.view.camera.distance = 5
        self.scatter = visuals.Markers(scaling=True, spherical=False, symbol='disc', alpha=0.8)
        self.orientation = visuals.Arrow(arrow_type='triangle_30', antialias=True, connect="segments", width=1)
        self.view.add(self.scatter)
        self.view.add(self.orientation)
        self.axis = visuals.XYZAxis(parent=self.view.scene)

class VideoWriter:
    def __init__(self, filename):
        self.filename = 'data/' + filename + '.mp4'
        self.writer = imageio.get_writer(self.filename)

    def append(self, frame):
        self.writer.append_data(frame)

    def close(self):
        self.writer.close()

class Evolution:
    def __init__(self, n_steps, datafile='file', video=False, filename='video.mp4'):
        self.canvas = Canvas()
        self.t = 0
        self.n_steps = n_steps
        self.df = 'data/' + datafile + '.hdf5'
        
        if video:
            self.writer = VideoWriter(filename)
        
        with h5py.File(self.df, 'r') as f:
            self.data = f['cell_positions'][:]
        
    def simulation(self, t):
        cells = self.data[t]
        # Identify and remove cells with coordinates (0, 0, 0)
        idx = np.where(np.linalg.norm(cells, axis=1) < 1e-6)[0]
        cells = np.delete(cells, idx, axis=0)
        return cells
    
    def update(self, event):
        r = self.simulation(self.t)
        self.canvas.scatter.set_data(r, face_color=COLOR_CELL, size=100, edge_color=COLOR_CELL, edge_width=5)
        # Unitary vector
        u = r / np.linalg.norm(r)
        self.canvas.orientation.set_data(np.stack((r, r + u), axis=1), color=ARROW_COLOR)
        progressbar(self.t+1, self.n_steps, "Simulating...")
        self.t += 1
    
    def movie(self, event):
        self.update(event)
        # Record frames
        #image = self.canvas.canvas.render()
        #self.writer.append(image)
        # Update canvas
        self.writer.append(self.canvas.render(alpha=True))
        # Step forward
        if self.t+1 == self.n_steps:
            self.writer.close()

    def run(self, video=False):
        if video:
            timer = app.Timer(interval=0.01, connect=self.movie, start=True, iterations=self.n_steps)
        else:
            timer = app.Timer(interval=0.01, connect=self.update, start=True, iterations=self.n_steps)
        timer.events.stop.connect(lambda x: app.quit())
        self.canvas.canvas.app.run()

def progressbar(j, count, prefix="", size=30, out=sys.stdout):
    x = int(size*j/count)
    if j < count:
        print(f"{prefix}[{'='*(x+1)}.{'.'*(size-x-1)}] {j}/{count}", end='\r', file=out, flush=True)
    else:
        print(f"{prefix}[{'='*size}] {j}/{count}", end='\r', file=out, flush=True)