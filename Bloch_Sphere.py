import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio

# Create the sphere
phi, theta = np.mgrid[0:2*np.pi:180j, 0:np.pi:90j]
x = np.sin(theta) * np.cos(phi)
y = np.sin(theta) * np.sin(phi)
z = np.cos(theta)

# Initial state parameters
theta_init = np.pi / 4  # Initial polar angle
phi_init = 0  # Initial azimuthal angle

# Calculate initial coordinates
x_init = np.sin(theta_init) * np.cos(phi_init)
y_init = np.sin(theta_init) * np.sin(phi_init)
z_init = np.cos(theta_init)

# Time parameters for animation
frames = 100
theta_vals = np.linspace(theta_init, np.pi, frames)
phi_vals = np.linspace(phi_init, 2 * np.pi, frames)

# Create a list to store frames
data_frames = []

# Generate frames for the animation
for i in range(frames):
    x_val = np.sin(theta_vals[i]) * np.cos(phi_vals[i])
    y_val = np.sin(theta_vals[i]) * np.sin(phi_vals[i])
    z_val = np.cos(theta_vals[i])
    
    frame = go.Scatter3d(x=[0, x_val], y=[0, y_val], z=[0, z_val],
                         mode='lines+markers', marker=dict(size=[0, 5]),
                         line=dict(width=5, color='red'))
    data_frames.append(frame)

# Create the figure
fig = go.Figure(
    data=[go.Surface(x=x, y=y, z=z, opacity=0.5, colorscale='Viridis')],
    layout=go.Layout(
        title="Bloch Sphere",
        scene=dict(
            xaxis=dict(nticks=5, range=[-1, 1]),
            yaxis=dict(nticks=5, range=[-1, 1]),
            zaxis=dict(nticks=5, range=[-1, 1]),
        ),
        updatemenus=[{
            "buttons": [
                {
                    "args": [None, {"frame": {"duration": 50, "redraw": True},
                                    "fromcurrent": True, "mode": "immediate"}],
                    "label": "Play",
                    "method": "animate"
                },
                {
                    "args": [[None], {"frame": {"duration": 0, "redraw": True},
                                      "mode": "immediate"}],
                    "label": "Pause",
                    "method": "animate"
                }
            ],
            "direction": "left",
            "pad": {"r": 10, "t": 87},
            "showactive": False,
            "type": "buttons",
            "x": 0.1,
            "xanchor": "right",
            "y": 0,
            "yanchor": "top"
        }]
    ),
    frames=[go.Frame(data=[frame]) for frame in data_frames]
)

# Show the figure
pio.show(fig)
