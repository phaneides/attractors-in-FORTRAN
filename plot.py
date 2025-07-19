import numpy as np
from vispy import scene, app

# === Load data ===
data = np.loadtxt('solution.dat')  # Columns: t, x, y, z
positions = data[:, 1:4]
n_points = positions.shape[0]

# === Create canvas ===
canvas = scene.SceneCanvas(keys='interactive', bgcolor='black', size=(800, 600), show=True)
view = canvas.central_widget.add_view()
view.camera = scene.cameras.TurntableCamera(fov=60, azimuth=60, elevation=30, distance=50)

# === White trail line ===
line = scene.Line(pos=positions[:1], color='white', width=2, parent=view.scene)

# === Red head marker ===
head = scene.Markers(pos=positions[:1], face_color='red', size=8, parent=view.scene)

# === Animation state ===
index = 1

# === Update function ===
def update(event):
    global index
    if index < n_points:
        line.set_data(positions[:index])
        head.set_data(positions[index-1:index])  # just the latest point
        index += 1
    else:
        timer.stop()

# === Timer ===
timer = app.Timer(interval=1/1000.0, connect=update, start=True)

# === Run ===
if __name__ == '__main__':
    app.run()

