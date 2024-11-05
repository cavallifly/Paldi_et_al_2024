from tadphys.modelling.lammps_modelling import generate_chromosome_rosettes_conformation
from math import pi
import sys

chrlength  = XXXchrlengthXXX
resolution = XXXresolutionXXX
nparticles = XXXsizeXXX
nchrs = XXXncopiesXXX
b=XXXbXXX
chromosome_particle_numbers=[nparticles]*nchrs
print(chrlength,nparticles, nchrs, chromosome_particle_numbers)

radius = XXXradiusXXX

print("%d Chromosomes of length %d in a sphere of radius %s at density %f bp/nm³" % (nchrs, nparticles, radius, chrlength/(4./3.*pi*radius*radius*radius*b*b*b)))
sys.stdout.flush()


r = XXXreplicaXXX
seed = XXXseedXXX

sigma=XXXsigmaXXX

#At the new coarse-graining, 14nm~1kb: What should be the size of a single bead?
for replica in range(r, r+1,1):
    initial_conformation = "Initial_rosette_conformation_sphere_replica_%s.dat" % (replica)
    generate_chromosome_rosettes_conformation ( chromosome_particle_numbers, fractional_radial_positions=None,
                                                confining_environment=['sphere',radius], rosette_radius=12.0 , particle_radius=sigma/2 ,
                                                seed_of_the_random_number_generator=seed ,
                                                number_of_conformations=1,
                                                outfile = initial_conformation,
                                                atom_types = 3)

