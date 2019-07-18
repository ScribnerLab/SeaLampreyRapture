This is the file structure that I used to run the colony runs in parallel on the MSU cluster. The dat.files folder holds the individual
simulation files 1-100 (There are 8 (of 9) other folders like this for the sea lamprey paper. The jobs are submitted in slurm. When COLONY
is done with identifying a BestConfiguration pedigree that file is post processed. See the scripts in the Source directory for more information.
#output/ is where all the outputs are pointed and std.eo is where standard error files go.