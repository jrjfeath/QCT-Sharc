import os
import queue
import select
import shutil
import subprocess
import threading
import time

import h5py
import numpy as np
import PESV2
import periodictable

from tqdm import tqdm

file_path = os.path.realpath(__file__)
file_path = os.path.dirname(file_path)
iscatter_input = '/mnt/scratch/users/jrjfeath/QCT_Sharc/Example/NO+Ar/input'
trajectories = 100
numpy_seed = 177013
cores = 8
# If debug=True it saves all geometry/velocity data for each trajectory
debug = False
A2b = 1.8897259886

# Import the iScatter functions required for generating the inputs
from iScatter.functions import iscattering
from iScatter.generate_conditions import generate_conditions

q_in = queue.Queue()
q_out = queue.Queue()

# load masses of atoms
massN = periodictable.N.mass
massO = periodictable.O.mass
massAr = periodictable.Ar.mass

#Load potential energy surface data
coeffs = np.loadtxt(f'{file_path}/potcoeff.dat', usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9))
Rs = coeffs[:10001, 0]
coeffs = coeffs[:10001, 1:]/219474.6

def format_floats(value : float) -> str:
    float_to_exp = f'{value:.9E}'.split("E")
    # Needs to be E+000 or E-000
    format_exp = f'{int(float_to_exp[1][1:]):03}'
    formatted_string = f'{float_to_exp[0]}E{float_to_exp[1][0]}{format_exp}'
    if formatted_string[0] != '-':
        formatted_string = ' ' + formatted_string
    return formatted_string

def make_QM_out(pot, deriv) -> str:
    string = ''.join(['! 1 Hamiltonian Matrix (1x1, complex)\n',
        '1 1\n',
        f'{format_floats(pot)} {format_floats(0.0)}\n',
        '! 2 Dipole Moment Vectors (3x1x1, complex)\n',
        '1 1\n',
        f'{format_floats(0.0)} {format_floats(0.0)}\n',
        '1 1\n',
        f'{format_floats(0.0)} {format_floats(0.0)}\n',
        '1 1\n',
        f'{format_floats(0.0)} {format_floats(0.0)}\n',
        '! 3 Gradient Vectors (1x3x3, real)\n',
        '3 3\n',
        f'{format_floats(deriv[0, 0])} {format_floats(deriv[0, 1])} {format_floats(deriv[0, 2])}\n',
        f'{format_floats(deriv[1, 0])} {format_floats(deriv[1, 1])} {format_floats(deriv[1, 2])}\n',
        f'{format_floats(deriv[2, 0])} {format_floats(deriv[2, 1])} {format_floats(deriv[2, 2])}\n',
        '! 6 Overlap matrix (1x1, complex)\n',
        '1 1\n',
        f'{format_floats(0.0)} {format_floats(0.0)}\n',
        '! 8 Runtime\n',
        f'{format_floats(0.0)}\n'])
    return string

def make_qct_input(gi, xyz, vxyz, atoms) -> None:
    '''Create all the files needed for running.'''
    dir_name = f'{file_path}/Traj/{gi}'
    # Make or empty trajectory folders
    if os.path.isdir(dir_name): shutil.rmtree(dir_name)
    os.mkdir(dir_name)

    velocities = []
    for velocity in vxyz:
        velocities.append(velocity.split()[1:])
    np.savetxt(f'{dir_name}/veloc',np.array(velocities,dtype=float),fmt='%.6f')

    geoms = ''
    for id in range(len(xyz)):
        geom = xyz[id].split()
        geoms += f"{atoms[id]['symbol']} {atoms[id]['number']:.1f} {geom[1]} {geom[2]} {geom[3]} {atoms[id]['mass']}\n"
    
    with open(f'{dir_name}/geom','w') as opf:
        opf.write(geoms)

    with open(f'{dir_name}/rattle','w') as opf:
        opf.write('1 2')

    input_str = 'printlevel 0\ngeomfile "geom"\nveloc external\nvelocfile "veloc"\nrattle\n'
    input_str += 'rattlefile "rattle"\nnstates 1\nstate 1 mch\ncoeff auto\nstepsize 10.0\n'
    input_str += 'nsubsteps 1\ntmax 20000000\nnorestart'
    
    with open(f'{dir_name}/input','w') as opf:
        opf.write(input_str)

def worker(worker_id) -> None:
    while True:
        gi, xyz, vxyz, natoms, atoms = q_in.get()
        if gi is None:
            return
        runtime = time.time()
        
        # Make the input files for sharc
        make_qct_input(gi, xyz, vxyz, atoms)
        # Create the working directory path for subprocess
        dir_name = f'{file_path}/Traj/{gi}'

        command = [f"/programs/sharc-3.0.2/bin/sharc.x", "input"]

        process = subprocess.Popen(
            command, 
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=f'{dir_name}'
        )

         # Set up file descriptors for stdout and stderr
        stdout_fd = process.stdout.fileno()
        stderr_fd = process.stderr.fileno()

        # Empty lists to hold the outputs from sharc
        geometries = []
        velocities = []
        terminate = 0
        start = time.time() # Kill the job if it gets stuck for some reason
        while (terminate == 0) and (start-time.time() < 5.0):
            # Use select to wait for any output from stdout or stderr
            readable_fds, _, _ = select.select([stdout_fd, stderr_fd], [], [])
            for fd in readable_fds:
                # Check to see if Sharc output an error or text
                if fd == stdout_fd:
                    stdout = process.stdout.readline().strip()
                    if len(stdout) == 0: continue
                    # If Sharc wants us to Run the QM Calculation
                    if stdout == 'QM Please':
                        data = []
                        # Loop over data until we have the geometry
                        while len(data) < natoms:
                            stdout = process.stdout.readline().strip()
                            if len(stdout) == 0: continue
                            data.append(stdout.split()[1:4])
                        # After we get data format it for PESV2
                        N  = data[0]
                        O  = data[1]
                        Ar = data[2]
                        pot, deriv, terminate = PESV2.main(coeffs, Rs, N, O, Ar, massN, massO, massAr)
                        # If the geometry is bad do not make a QM file, terminate trajectory
                        if np.isnan(pot):
                            terminate = 1
                        else:
                            output = make_QM_out(pot, deriv)
                            # Tell fortran we are done processing the PES
                            process.stdin.write('Done\n')
                            process.stdin.write(output)
                            process.stdin.write("EOF\n")  # Send custom end marker
                            process.stdin.flush()
                    # Get geometry from Sharc
                    if stdout == 'GEOM':
                        geom = []
                        while len(geom) < natoms:
                            stdout = process.stdout.readline().strip()
                            if len(stdout) == 0: continue
                            geom.append(stdout.split())
                        geometries.append(np.array(geom, dtype=float))
                        # Tell Sharc it got all the data
                        process.stdin.write('Done\n')
                        process.stdin.flush()
                    # Get velocity from Sharc
                    if stdout == 'VELOC':
                        velocity = []
                        while len(velocity) < natoms:
                            stdout = process.stdout.readline().strip()
                            if len(stdout) == 0: continue
                            velocity.append(stdout.split())
                        velocities.append(np.array(velocity, dtype=float))
                        # Tell Sharc it got all the data
                        process.stdin.write('Done\n')
                        process.stdin.flush()
                if fd == stderr_fd:
                    line = process.stderr.readline().strip()
                    if not line: continue
                    print(f'{worker_id} running ({gi}) encountered issues: {line}')
                    terminate = 1
        
        # It is possible to reach the limits on the PES before sharc completes (PESV2), if so terminate sharc
        process.terminate()
        process.wait()

        # Delete sharc files after run
        shutil.rmtree(dir_name)

        # Write everything to files if debug selected
        if debug:
            with open(f'{dir_name}_geom.npy', 'wb') as f:
                np.save(f, geometries)
            with open(f'{dir_name}_veloc.npy', 'wb') as f:
                np.save(f, velocities)
        
        # Some geometries will fail immediately, so pass empty lists
        if len(geometries) != 0:
            # Put the first and last sets of data into the output queue
            q_out.put([
                gi, 
                time.time() - runtime, 
                [geometries[0], geometries[-1]], 
                [velocities[0], velocities[-1]]
            ])
        else:
            q_out.put([
                gi,
                time.time() - runtime,
                [],
                []
            ])

def run_parallel_qct_batch():
    threads = [threading.Thread(target=worker, args=(id,)) for id in range(cores)]
    for thread in threads:
        thread.start()

    #Make directory for storing Trajectory data
    if os.path.isdir(f'{file_path}/Traj'): shutil.rmtree(f'{file_path}/Traj')
    os.mkdir(f'{file_path}/Traj')

    # If you run too many trajectories at a time iScatter slows down dramatically
    limit = 50000
    iterations = int(np.ceil(trajectories / limit))

    atoms = []
    natoms = 0
    for current_step in range(iterations):
        print(f'Beginning step {current_step+1} of {iterations}')
        samples = limit
        if current_step == (iterations-1): samples = trajectories % limit
        sc = iscattering.iscatter(samples = samples, seed = numpy_seed+current_step)
        sc.ReadInput(iscatter_input)
        sc = generate_conditions(sc, samples)
        print(f'Beginning Sharc-QCT process')

        # When it is first initialized grab the atomic numbers, mass, and number of atoms
        if current_step == 0:
            for geom in sc.xyz[0]:
                geom = geom.split()
                if hasattr(periodictable,geom[0]):
                    atom = getattr(periodictable,geom[0])
                else:
                    print(f'{geom[0]} is not a valid atom!')
                atoms.append({'symbol':atom, 'mass': atom.mass, 'number':atom.number})
                natoms += 1

        for id in range(len(sc.xyz)):
            # Add the file ids to the queue
            q_in.put([id, sc.xyz[id], sc.vxyz[id], natoms, atoms])

        # The workers will return the job id, runtime, geometry, and velocity
        for _ in tqdm(range(len(sc.xyz))):
            id, runtime, geometries, velocities = q_out.get()
            # print(id, runtime)
            with open(f'{file_path}/Traj/geometries.npy', 'ab') as f:
                np.save(f, geometries)
            with open(f'{file_path}/Traj/velocities.npy', 'ab') as f:
                np.save(f, velocities)

    # When we run out of data kill the threads
    for thread in threads:
        q_in.put([None, None, None, None, None])

if __name__ == '__main__':
    run_parallel_qct_batch()