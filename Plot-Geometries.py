import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

output_dat = r'C:\Users\chem-chem1612\Documents\Max\output.dat'
save_gif = r'C:\Users\chem-chem1612\Documents\Max\traj.gif'

with open(output_dat,'r') as f:
    data = f.readlines()[32:] # Skip the first 32 header lines of the file

geometries = []
write_data = False
for line in data:
    line = line.strip()
    if line == '! 11 Geometry in a.u.': 
        write_data = True
        continue
    if line == '! 12 Velocities in a.u.': 
        write_data = False
        continue
    if not write_data: continue
    geometries.append(line.split())
geometries = np.array(geometries,dtype=float)

N = geometries[0::3]
O = geometries[1::3]
Ar = geometries[2::3]

# init the figure
fig, ax = plt.subplots(figsize=(5,5))

delimiter = 1
x_axis = 1
y_axis = 2

def update(index):
    frame = index * delimiter
    # clear the axis each frame
    ax.clear()
    # replot things
    ax.scatter(N[:,x_axis][frame], N[:,y_axis][frame], c='blue', label='N')
    ax.scatter(O[:,x_axis][frame], O[:,y_axis][frame], c='red', label='O')
    ax.scatter(Ar[:,x_axis][frame], Ar[:,y_axis][frame], c='green', label='Ar')
    ax.set_title(f'{frame * 10 :.1f} fs')

    # reformat things
    ax.set_xlim(-8,8)
    ax.set_ylim(-8,8)
    ax.set_xlabel('Bohr')
    ax.set_ylabel('Bohr')
    ax.legend()

ani = animation.FuncAnimation(fig, update, frames=int(N.shape[0]), interval=1)
ani.save(save_gif, writer='pillow')