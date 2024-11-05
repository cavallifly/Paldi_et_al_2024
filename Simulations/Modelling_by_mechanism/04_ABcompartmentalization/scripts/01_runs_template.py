from tadphys.modelling.lammps_modelling import run_lammps
from numpy import zeros
import itertools

# From 52795 to 58989 both included
nparticles = 4000
resolution = 5000
nchrs = 1
chromosome_particle_numbers=[nparticles]*nchrs
print("Nparticles",nparticles,"Ncopies",nchrs,"Particles per copy",chromosome_particle_numbers)

timestep=0.006

phi=0.10 #XXXphiXXX
if phi == 0.03:
    secondInTauLJ = 16 # tau_LJ for 1s in real time
if phi == 0.10:
    secondInTauLJ = 9  # tau_LJ for 1s in real time    
if phi == 0.20:
    secondInTauLJ = 9  # tau_LJ for 1s in real time
if phi == 0.30:
    secondInTauLJ = 13 # tau_LJ for 1s in real time    

second        = int(secondInTauLJ / 0.012 * int(0.012/timestep)) # tau_LJ for 1s in real time
print("1 Second in simulation time = ",secondInTauLJ,"TauLJ")
print("1 Second in simulation time = ",second,"dt")
print("")

runtime = int(4*60*60*second) # TSA treatment lasts 4 hours
print("Run time =",4*60*60,"s")
print("Run time =",runtime,"dt")
dumpTime = int(10*60*second)
print("Dump time =",dumpTime,"dt")

# Compartmentalization
attrEnergy=0.0
Asize=1500000 # Average size of A compartment
Bsize=1500000 # Average size of B compartment
Ssize=Asize+Bsize
offset1 = int(  1000000/resolution)
offset2 = int( 21000000/resolution)
offset3 = int( 41000000/resolution)
offset4 = int( 61000000/resolution)
offset5 = int( 81000000/resolution)
offset6 = int(101000000/resolution)
offset7 = int(121000000/resolution)
offset8 = int(141000000/resolution)
print("First part from 0 to 1000000bp is neutral. It spans %d beads" % offset1)
print("Segment with A-B partition spans %d bp or %d beads" % (Ssize,int(Ssize/resolution)))
# 1- 3 First chromosome

ncopies = XXXncopiesXXX

if ncopies == 1: 
    compartmentalization = {
        "runtime"   : runtime,
        "partition" : {
            1 : itertools.chain(range(offset1+0*int(Ssize/resolution)+1,offset1+0*int(Ssize/resolution)+int(Asize/resolution)+1,1), # First Chain
                                range(offset1+1*int(Ssize/resolution)+1,offset1+1*int(Ssize/resolution)+int(Asize/resolution)+1,1),
                                range(offset1+2*int(Ssize/resolution)+1,offset1+2*int(Ssize/resolution)+int(Asize/resolution)+1,1),
                                range(offset1+3*int(Ssize/resolution)+1,offset1+3*int(Ssize/resolution)+int(Asize/resolution)+1,1),
                                range(offset1+4*int(Ssize/resolution)+1,offset1+4*int(Ssize/resolution)+int(Asize/resolution)+1,1),
                                range(offset1+5*int(Ssize/resolution)+1,offset1+5*int(Ssize/resolution)+int(Asize/resolution)+1,1)),
            2 : itertools.chain(range(offset1+0*int(Ssize/resolution)+int(Asize/resolution)+1,offset1+0*int(Ssize/resolution)+int(Asize/resolution)+int(Bsize/resolution)+1,1), # First chain
                                range(offset1+1*int(Ssize/resolution)+int(Asize/resolution)+1,offset1+1*int(Ssize/resolution)+int(Asize/resolution)+int(Bsize/resolution)+1,1),
                                range(offset1+2*int(Ssize/resolution)+int(Asize/resolution)+1,offset1+2*int(Ssize/resolution)+int(Asize/resolution)+int(Bsize/resolution)+1,1),
                                range(offset1+3*int(Ssize/resolution)+int(Asize/resolution)+1,offset1+3*int(Ssize/resolution)+int(Asize/resolution)+int(Bsize/resolution)+1,1),
                                range(offset1+4*int(Ssize/resolution)+int(Asize/resolution)+1,offset1+4*int(Ssize/resolution)+int(Asize/resolution)+int(Bsize/resolution)+1,1),
                                range(offset1+5*int(Ssize/resolution)+int(Asize/resolution)+1,offset1+5*int(Ssize/resolution)+int(Asize/resolution)+int(Bsize/resolution)+1,1))},
        "radii" : { 1: 0.5,  2: 0.5},
        "interactions" : {( 1, 1) : ["attraction",AAattrEnergy],( 2, 2) : ["attraction",BBattrEnergy]}
    }

for chromType in compartmentalization["partition"]:
    print(chromType,compartmentalization["partition"][chromType])

b=XXXbXXX
lk=XXXlkXXX
radius=XXXradiusXXX

r = XXXreplicaXXX

#At the new coarse-graining, 14nm~1kb: What should be the size of a single bead?
for replica in range(r, r+1,1):
    initial_conformation = "initial_conformation.txt"
    
    run_lammps(initial_conformation=initial_conformation,
               minimize  = False, 
               tethering = False,
               compartmentalization    = compartmentalization,
               kseed = int(replica),
               run_time = runtime,
               to_dump = dumpTime,
               pbc=True,
               persistence_length = lk / b / 2,               
               lammps_folder = "./",
               chromosome_particle_numbers = [nparticles]*nchrs,
               hide_log = False,
               confining_environment = ['sphere',radius],
               timestep=timestep,
               )
