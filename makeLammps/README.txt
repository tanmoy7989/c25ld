All Lammps data are to be kept under data/lammpstraj. 
The directory hierarchy format that we'll be following is:

Target_dir--> constrained_lj, unconstrained_lj, constrained_wca, unconstrained_wca


the steps are as:

### To generate a Lammps data file for the system anywhere

1) Enter data at a python prompt:
   
   datafileprefix = c25   (Idea: The datafile prefix will usually refer to the chain length)
   n_water = 1700
   temp = 298.0
   Target_dir = your desired directory

2) Generate data file inside required directory:
	
python genLammpsin.py temp n_water datafileprefix Target_dir



### To generate complete directory hierarchy, each complete with its own set of python and job scripts and so on.

1) Please read and modify the file writelammpsjob.py to set the temperature, n_water etc correctly before you run it. This is to make sure that the final lammps input scripts are in sync with the datafiles they'll call. Discrepancies in the number of water molecules in the two cases can lead to major computational overload and ultimate program crash.
2) Run it from the python prompt as:

python writelammpsjob.py target_dir biasfile_lj biasfile_wca

target_dir should preferably be named after the month or week and year. So, something like jan_1_2015_runs and so on.
biasfile_lj and biasfile_wca are the files where a list of umbrella centers and corresponding force constants are to be maintained for that particular set of umbrella sampling runs. Frequently these will be the same files. A quick and dirty way to create umbrellas that are equispaced with the same force-constant is encapsulated in makebiasfile.py present in the c25ld parent directory. Please look it up.

