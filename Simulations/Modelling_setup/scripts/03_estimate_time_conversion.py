from os import rename

from tadphys.modelling.lammps_modelling import run_lammps

nparticles = XXXnparticlesXXX

lk = XXXlkXXX
b  = XXXbXXX
replica = XXXreplicaXXX
radius  = XXXradiusXXX

sigma=XXXsigmaXXX


initial_conformation = "final_conformation.txt"
initial_relaxation = { "relaxation_time" : XXXrunXXX,
                       "MSD"             : 100,
                      }

run_lammps(initial_conformation=initial_conformation,
           minimize   = False, 
           initial_relaxation = initial_relaxation,
           connectivity   = {1:["FENE", 30.0, 1.5*sigma, 1.0, 1.0*sigma]},
           excludedVolume = {(1,1):["LJ",1.0, 1.0*sigma, 1.12246152962189*sigma]},
           confining_environment = ['sphere',radius,0.0,0.0,0.0,1.0,sigma],
           kseed = int(replica),
           to_dump = 5000,
           timestep = XXXtimestepXXX,
           persistence_length = lk / b / 2,
           reset_timestep = 0,
           run_time = 0,
           lammps_folder = "./",
           chromosome_particle_numbers = [nparticles]*2,
           hide_log=False
           )
