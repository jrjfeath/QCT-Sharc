import os
import subprocess
import time
import numpy as np
import periodictable
import PESV2 # Made with numpy.f2py, see README!

sharc_directory = r'/home/hartree2/mb/chem1612/programs/sharc-3.0.2/bin/sharc.x'

A2b = 1.8897259886
file_path = os.path.realpath(__file__)
file_path = os.path.dirname(file_path)

def format_floats(value):
    float_to_exp = f'{value:.9E}'.split("E")
    # Needs to be E+000 or E-000
    format_exp = f'{int(float_to_exp[1][1:]):03}'
    formatted_string = f'{float_to_exp[0]}E{float_to_exp[1][0]}{format_exp}'
    if formatted_string[0] != '-':
        formatted_string = ' ' + formatted_string
    return formatted_string

def make_QM_out(pot, deriv):
    with open(f'{dir_name}/QM/QM.out', 'w') as f:
        f.write(''.join(['! 1 Hamiltonian Matrix (1x1, complex)\n',
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
        f'{format_floats(0.0)}\n']))
    return 

files = [x for x in os.listdir(f'{file_path}/outputs')]
#Make directory for storing Trajectory data
if os.path.isdir(f'{file_path}/Traj'):
    os.system(f'rm -rf {file_path}/Traj')
os.mkdir(f'{file_path}/Traj')
print('Creating Sharc inputs')
for index in range(int(len(files)/3)):
    # Make or empty trajectory folders
    if os.path.isdir(f'{file_path}/Traj/{index}'):
        os.system(f'rm -rf {file_path}/Traj/{index}')        
    os.system(f'mkdir {file_path}/Traj/{index}')
    os.system(f'mkdir {file_path}/Traj/{index}/QM')
    
    #Make veloc file
    with open(f'{file_path}/outputs/out_{index}.vel','r') as opf:
        data = opf.read().split('\n')
    na = int(data[0])
    velocities = []
    for velocity in data[2:2+na]:
        velocities.append(velocity.split()[1:])
    np.savetxt(f'{file_path}/Traj/{index}/veloc',np.array(velocities,dtype=float),fmt='%.6f')

    #Make geom file
    with open(f'{file_path}/outputs/out_{index}.xyz','r') as opf:
        data = opf.read().split('\n')
    na = int(data[0])
    geoms = ''
    for geom in data[2:2+na]:
        geom = geom.split()
        if hasattr(periodictable,geom[0]):
            atom = getattr(periodictable,geom[0])
        else:
            print(f'{geom[0]} is not a valid atom!')
        atomic_number = atom.number
        atomic_mass = atom.mass
        geoms += f'{atom} {atomic_number:.1f} {geom[1]} {geom[2]} {geom[3]} {atomic_mass}\n'
    with open(f'{file_path}/Traj/{index}/geom','w') as opf:
        opf.write(geoms)

    with open(f'{file_path}/Traj/{index}/rattle','w') as opf:
        opf.write('1 2')

    input_str = 'printlevel 0\ngeomfile "geom"\nveloc external\nvelocfile "veloc"\nrattle\n'
    input_str += 'rattlefile "rattle"\nnstates 1\nstate 1 mch\ncoeff auto\nstepsize 10.0\n'
    input_str += 'nsubsteps 1\ntmax 20000000\nnorestart'
    with open(f'{file_path}/Traj/{index}/input','w') as opf:
        opf.write(input_str)

coeffs = np.loadtxt(f'{file_path}/potcoeff.dat', usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9))
Rs = coeffs[:10001, 0]
coeffs = coeffs[:10001, 1:]/219474.6

massN = periodictable.N.mass
massO = periodictable.O.mass
massAr = periodictable.Ar.mass

start_i = time.time()
for index in range(int(len(files)/3)):
    if index % 50 == 0: print(f'Running trajectory {index} of {int(len(files)/3)}')
    start_traj = time.time()
    dir_name = f'{file_path}/Traj/{index}'
    # Remove files before executing code
    main_files = ['output.dat','output.lis','output.log','output.xyz','STOP']
    for file in main_files:
        if os.path.exists(f'{dir_name}/{file}'): os.remove(f'{dir_name}/{file}')
    qm_files = ['QM.out']
    for file in qm_files:
        if os.path.exists(f'{dir_name}/QM/{file}'): os.remove(f'{dir_name}/QM/{file}')

    os.chdir(f'{dir_name}')
    command = [f"{sharc_directory} input"]

    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True, preexec_fn=os.setsid)

    while not os.path.exists(f'{dir_name}/STOP'):
        line = process.stdout.readline()
        # Loop until QM Please is passed by fortran
        if len(line) > 0:
            # Read the stdout from fortran written by write(6)
            # data is coord data in the order of the atoms in the input
            data = []
            # Loop through atom data until Done! is reached
            for line in iter(process.stdout.readline, b''):
                # Decode the byte string and remove \n
                line_str = line.decode().strip()
                if line_str == 'Done!': break
                # Sharc returns geom and velocity, we only need geom
                data.append(line_str.split()[1:4])
            N  = data[0]
            O  = data[1]
            Ar = data[2]
            start = time.time()
            pot, deriv = PESV2.main(coeffs, Rs, N, O, Ar, massN, massO, massAr)
            make_QM_out(pot, deriv)
            # Tell fortran we are done processing the PES
            process.stdin.write('Done\n'.encode())
            process.stdin.flush()
print(f'All trajectories finished in {time.time()-start_i:.2f} seconds')



