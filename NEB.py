import os
import numpy as np
from ase import Atoms
from ase.calculators.aims import Aims, AimsCube
from ase.optimize import BFGS
from ase.io import read, write
from ase.dyneb import DyNEB
from ase.neb import NEB

os.environ['ASE_AIMS_COMMAND'] = 'srun --distribution=block:block --hint=nomultithread aims.mpi.x > aims_epcc.out'
os.environ['AIMS_SPECIES_DIR'] = '/work/e05/e05/mat92/Codes/Aims/defaults_2020/intermediate/'

#spin_cube = AimsCube(plots=('total_density',
#                             'spin_density'))
initial = read('ini.in')
final = read('fin.in')


#images = read('neb.traj@-10:')

#for image in images:
#    calc = Aims(xc=('hse06',0.25),
#            hybrid_xc_coeff=0.5,
#            hse_unit='B',
#            mixer='linear',
#            charge_mix_param=0.08,
#            sc_accuracy_etot=1e-4,
#            sc_accuracy_eev=1e-2,
#            sc_accuracy_rho=1e-5,
#            sc_accuracy_forces=1e-2,
#            spin='collinear',
#            charge=-1,
#            fixed_spin_moment=1.0,
#            load_balancing='.true.',
#            collect_eigenvectors='.false.',
##            distributed_spline_storage='.true.',
#            n_max_pulay=6,
##            default_initial_moment=1.0,
#            relativistic=('atomic_zora','scalar'),
#            k_grid=(1,1,1))
#    image.set_calculator(calc)

images = []
for i in range(0,9):
    image = initial.copy()
    images.append(image)
images.append(final)

#Create initial moment arrays with one electron on atom 236 and 237 respectively

moments = np.zeros(len(initial.get_atomic_numbers()))

#
moments1 = initial.get_initial_magnetic_moments()
#moments1[0] = 1.0
#
moments2 = final.get_initial_magnetic_moments()
#moments2[1] = 1.0

for i in range(5):
    images[i].set_initial_magnetic_moments(moments1)


for i in range(5,10):
    images[i].set_initial_magnetic_moments(moments2)

#neb = NEB(images)
neb = DyNEB(images, k=0.6, fmax=0.05, dynamic_relaxation=True)
neb.interpolate(method='idpp')

#for i in range(len(images)):
#    write('0{i}.in'.format(i=i), images[i], format='aims')
#    write('POSCAR_0{i}'.format(i=i), images[i], format='vasp')


for image in images:
    image.calc = Aims(xc=('hse06',0.25),
        hybrid_xc_coeff=0.5,
        hse_unit='B',
        sc_accuracy_etot=1e-4,
        sc_accuracy_eev=1e-2,
        sc_accuracy_rho=1e-5,
        sc_accuracy_forces=1e-2,
        calculate_fock_matrix_version=5,
        charge_mix_param=0.01,
        spin='collinear',
        charge=-1,
        fixed_spin_moment=1.0,
        relativistic=('atomic_zora','scalar'),
        k_grid=(1,1,1))

qn = BFGS(neb, trajectory='neb.traj')
qn.run(fmax=0.05)
