# Simulations

To perform the simulations, we used the TADphys package, which you can download from [here](https://github.com/cavallifly/TADphys) and install using the instructions therein.

To prepare the simulations with 1 and 20 chains, you can launch the following command from the Simulation folder:


```bash
cd Modelling_setup
bash ./scripts/Modelling_setup.sh
```
<!---
<div style="text-align: justify">
Each chain is initially organized in a rod-like folding featuring rosettes along the main axis and placed in random positions inside a confining sphere of radius R<sup>*</sup> to set up the volume density of the system to 3% (R<sup>*</sup>= 25.5σ for the 1-chain system and R<sup>*</sup>=69.3σ for the 20-chain system), avoiding clashes with other chains (step 00). After energy minimization, each polymeric system is compressed to reach the DNA density of 10% (step 01). These conditions were achieved by minimization (LAMMPS command: minimize 1.0e-4 1.0e-6 100000 100000) followed by molecular dynamics simulations of 600 τ<sub>LJ</sub> (100,000 ∆t with ∆t = 0.006 τ<sub>LJ</sub>) during which the radius of confining spheres is reduced from the minimum radius to include all the particles of the chains at time 0 to the target radius R (R= 17.1σ for the 1-chain system and R=46.4σ for the 20-chain system). At the target volume density of 10%, the polymer chains have parameters σ~54.2nm and 𝐾<sub>𝜃</sub> ~92.3nm. These estimates were done by considering a fine-scale chromatin model with 𝜈<sub>FS</sub>=100 𝑏𝑝, 𝜎<sub>FS</sub>=20 𝑛𝑚 and 𝐾<sub>𝜃</sub><sup>FS</sup>=50 𝑛𝑚. Finally, each polymeric system is relaxed with a molecular dynamics run of 30,600 τ<sub>LJ</sub> (5,100,000 ∆t with ∆t = 0.006 τ<sub>LJ</sub>). By comparing the average monomer Mean-Squared Displacement (MSD) in these relaxation runs and the MSD of non-transcribed genes measured 96 by live-cell imaging, we obtained an approximated estimate of the simulated time (in τ<sub>LJ</sub>) corresponding to 1s ~ 9 τ<sub>LJ</sub>. These conformations are next used as the initial conformations for the downstream simulations (Steps 02 and 03). </div>
--->
