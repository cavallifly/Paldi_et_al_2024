from os import rename

from tadphys.modelling.lammps_modelling import run_lammps

# From 52795 to 58989 both included
nparticles = XXXnparticlesXXX
chrlength  = int(nparticles*5000)

nchrs = XXXncopiesXXX
chromosome_particle_numbers=[nparticles]*nchrs
#print(chrlength,nparticles, nchrs, chromosome_particle_numbers)

timesteps_relaxation = 0

lk = XXXlkXXX
b  = XXXbXXX
run_duration=0
timesteps_per_loop_extrusion_step = 0

replica = XXXreplicaXXX
radius=XXXradiusXXX

initial_conformation = "XXXinitial_conformationXXX"
reset_timestep=0

wdir = "XXXwdirXXX/"
runtime=XXXruntimeXXX
extrusionTime=XXXextrusionTimeXXX
sigma=XXXsigmaXXX

run_lammps(initial_conformation=initial_conformation,
           minimize  = False, 
           tethering = False,
           kseed = int(replica),
           to_dump = extrusionTime,
           persistence_length = lk / b / 2,
           timestep=0.006,
           connectivity   = {1:["FENE", 30.0, 1.5*sigma, 1.0, 1.0*sigma]},
           excludedVolume = {(1,1):["LJ",1.0, 1.0*sigma, 1.12246152962189*sigma]},
           confining_environment = ['sphere',radius,0.0,0.0,0.0,1.0,sigma],
           run_time = runtime,
           lammps_folder = wdir,
           chromosome_particle_numbers = chromosome_particle_numbers,
           reset_timestep=reset_timestep
)
