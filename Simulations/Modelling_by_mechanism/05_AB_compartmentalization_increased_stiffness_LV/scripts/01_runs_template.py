from tadphys.modelling.lammps_modelling import run_lammps
from numpy import zeros, array_equal
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

offset = [0]*XXXncopiesXXX
for o in range(len(offset)):
    offset[o] = int((1000000+o*20000000)/resolution)
    print(o+1,offset[o])
    
ncopies = XXXncopiesXXX

partition = {}
partition[1] = []
partition[2] = []    
for i in range(ncopies):
    partition[1] += itertools.chain(range(offset[i]+0*int(Ssize/resolution)+1,offset[i]+0*int(Ssize/resolution)+int(Asize/resolution)+1,1), # i-th Chain
                                    range(offset[i]+1*int(Ssize/resolution)+1,offset[i]+1*int(Ssize/resolution)+int(Asize/resolution)+1,1),
                                    range(offset[i]+2*int(Ssize/resolution)+1,offset[i]+2*int(Ssize/resolution)+int(Asize/resolution)+1,1),
                                    range(offset[i]+3*int(Ssize/resolution)+1,offset[i]+3*int(Ssize/resolution)+int(Asize/resolution)+1,1),
                                    range(offset[i]+4*int(Ssize/resolution)+1,offset[i]+4*int(Ssize/resolution)+int(Asize/resolution)+1,1),
                                    range(offset[i]+5*int(Ssize/resolution)+1,offset[i]+5*int(Ssize/resolution)+int(Asize/resolution)+1,1))
    partition[2] += itertools.chain(range(offset[i]+0*int(Ssize/resolution)+int(Asize/resolution)+1,offset[i]+0*int(Ssize/resolution)+int(Asize/resolution)+int(Bsize/resolution)+1,1), # i-th chain
                                    range(offset[i]+1*int(Ssize/resolution)+int(Asize/resolution)+1,offset[i]+1*int(Ssize/resolution)+int(Asize/resolution)+int(Bsize/resolution)+1,1),
                                    range(offset[i]+2*int(Ssize/resolution)+int(Asize/resolution)+1,offset[i]+2*int(Ssize/resolution)+int(Asize/resolution)+int(Bsize/resolution)+1,1),
                                    range(offset[i]+3*int(Ssize/resolution)+int(Asize/resolution)+1,offset[i]+3*int(Ssize/resolution)+int(Asize/resolution)+int(Bsize/resolution)+1,1),
                                    range(offset[i]+4*int(Ssize/resolution)+int(Asize/resolution)+1,offset[i]+4*int(Ssize/resolution)+int(Asize/resolution)+int(Bsize/resolution)+1,1),
                                    range(offset[i]+5*int(Ssize/resolution)+int(Asize/resolution)+1,offset[i]+5*int(Ssize/resolution)+int(Asize/resolution)+int(Bsize/resolution)+1,1))         

compartmentalization = {
    "runtime"   : runtime,
    "fraction"  : { 1 : fractionA , 2 : fractionB},
    "partition" : {
        1 : partition[1],                    
    2 : partition[2]},
    "radii" : { 1: 0.5,  2: 0.5},
    "interactions" : {( 1, 1) : ["attraction",AAattrEnergy],( 2, 2) : ["attraction",BBattrEnergy]}
}
#print(compartmentalization)

#for p in compartmentalization["partition"].keys():
#    print(p)
#    b = []
#    for j in partition[int(p)]:
#        print(j)
#    b = b.sort()
#    print(b)
#    continue
#    for i in range(len(b)):
#        print(b[i])
#exit(1)
     
for chromType in compartmentalization["partition"]:
    print(chromType,compartmentalization["partition"][chromType])

b=XXXbXXX
lk=XXXlkXXX
lkFA=XXXlkFAXXX
lkFB=XXXlkFBXXX
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
               persistence_length = [lk / b / 2 * lkFB, lk / b / 2 * lkFA],
               #persistence_length = lk / b / 2 * lkFA,               
               lammps_folder = "./",
               chromosome_particle_numbers = [nparticles]*nchrs,
               hide_log = False,
               confining_environment = ['sphere',radius],
               timestep=timestep,
               )
