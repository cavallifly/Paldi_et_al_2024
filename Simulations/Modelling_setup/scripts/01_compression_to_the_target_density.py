from os import rename

from tadphys.modelling.lammps_modelling import run_lammps

# From 52795 to 58989 both included
nparticles = XXXnparticlesXXX

nchrs = 1
chromosome_particle_numbers=[nparticles]*nchrs

timesteps_relaxation = 0

lk = XXXlkXXX
b  = XXXbXXX
run_duration=0
timesteps_per_loop_extrusion_step = 0

replica = XXXreplicaXXX
radius=XXXradiusXXX
sigma=XXXsigmaXXX

initial_conformation = "./final_conformation.txt"
reset_timestep=XXXreset_timestepXXX

wdir = "XXXwdirXXX/"
runtime=XXXruntimeXXX

run_lammps(initial_conformation=initial_conformation,
           minimize   = True,
           tethering  = False, 
           kseed = int(replica),
           to_dump = runtime,
           persistence_length = lk / b / 2,
           timestep=0.006,
           connectivity   = {1:["FENE", 30.0, 1.5*sigma, 1.0, 1.0*sigma]},
           excludedVolume = {(1,1):["LJ",1.0, 1.0*sigma, 1.12246152962189*sigma]},
           run_time = 0,
           lammps_folder = wdir,
           compress_without_pbc = [0.0,0.0,0.0,XXXIradiusXXX,XXXTradiusXXX,runtime],
           confining_environment = ['sphere',radius,0.0,0.0,0.0,1.0,sigma],
           chromosome_particle_numbers = chromosome_particle_numbers,
           hide_log=False,
           reset_timestep=reset_timestep
)
