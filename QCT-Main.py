import os
import queue
import select
import shutil
import subprocess
import threading
import time

import numpy as np
import PESV2
import periodictable

A2b = 1.8897259886
file_path = os.path.realpath(__file__)
file_path = os.path.dirname(file_path)
iscatter_input = '/mnt/scratch/users/jrjfeath/QCT_Sharc/Example/NO+Ar/input'
trajectories = 100
numpy_seed = 1
cores = 4

# Import the iScatter functions required for generating the inputs
from iScatter.functions import iscattering
from iScatter.generate_conditions import generate_conditions

q = queue.Queue()

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

def make_qct_input(gi, xyz, vxyz) -> None:
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
    for geom in xyz:
        geom = geom.split()
        if hasattr(periodictable,geom[0]):
            atom = getattr(periodictable,geom[0])
        else:
            print(f'{geom[0]} is not a valid atom!')
        atomic_number = atom.number
        atomic_mass = atom.mass
        geoms += f'{atom} {atomic_number:.1f} {geom[1]} {geom[2]} {geom[3]} {atomic_mass}\n'
    
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
    start =time.time()
    while True:
        gi, xyz, vxyz = q.get()
        if gi is None:
            print(f'{worker_id} terminated after: {time.time() - start}')
            return
        
        # Make the input files for sharc
        make_qct_input(gi, xyz, vxyz)
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

        count = 0
        terminate = 0
        while terminate == 0:
            # Use select to wait for any output from stdout or stderr
            readable_fds, _, _ = select.select([stdout_fd, stderr_fd], [], [])
            for fd in readable_fds:
                # Loop until QM Please is passed by fortran
                if fd == stdout_fd:
                    line = process.stdout.readline()
                    if not line: continue
                                        
                    # Read the stdout from fortran written by write(6)
                    # data is coord data in the order of the atoms in the input
                    data = []
                    # Loop through atom data until Done! is reached
                    for line in iter(process.stdout.readline, b''):
                        # Decode the byte string and remove \n
                        line_str = line.strip()
                        if line_str == 'Done!': break
                        # Ignore lines that are not array data
                        if len(line_str.split()) != 7: continue
                        # Sharc returns geom and velocity, we only need geom
                        data.append(line_str.split()[1:4])
                    N  = data[0]
                    O  = data[1]
                    Ar = data[2]
                    pot, deriv, terminate = PESV2.main(coeffs, Rs, N, O, Ar, massN, massO, massAr)
                    # If the geometry is bad do not make a QM file, kills Sharc
                    if np.isnan(pot):
                        pass
                    else:
                        output = make_QM_out(pot, deriv)
                    # Tell fortran we are done processing the PES
                    process.stdin.write('Done\n')
                    process.stdin.write(output)
                    process.stdin.write("EOF\n")  # Send custom end marker
                    process.stdin.flush()
                    count+=1

                if fd == stderr_fd:
                    line = process.stderr.readline()
                    if not line: continue
                    print(f'{worker_id} running ({gi}) encountered issues: {line}')
                    terminate = 1
        
        # It is possible to reach the limits on the PES before sharc completes (PESV2), if so terminate sharc
        process.terminate()
        process.wait()

def run_parallel_qct_batch(batch_size=100):
    threads = [threading.Thread(target=worker, args=(id,)) for id in range(cores)]
    for thread in threads:
        thread.start()

    #Make directory for storing Trajectory data
    if os.path.isdir(f'{file_path}/Traj'): shutil.rmtree(f'{file_path}/Traj')
    os.mkdir(f'{file_path}/Traj')

    # If you run too many trajectories at a time iScatter slows down dramatically
    limit = 50000
    iterations = int(np.ceil(trajectories / limit))

    for current_step in range(iterations):
        print(f'Beginning step {current_step+1} of {iterations}')
        samples = limit
        if current_step == (iterations-1): samples = trajectories % limit
        sc = iscattering.iscatter(samples = samples, seed = numpy_seed+current_step)
        sc.ReadInput(iscatter_input)
        sc = generate_conditions(sc, samples)
        print(f'Beginning Sharc-QCT process')

        for id in range(len(sc.xyz)):
            # Add the file ids to the queue
            q.put([id, sc.xyz[id], sc.vxyz[id]])

    # After all the files have been created, start the QCT
    for thread in threads:
        q.put([None, None, None])

if __name__ == '__main__':
    run_parallel_qct_batch()