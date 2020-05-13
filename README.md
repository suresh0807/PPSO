# PPSO
A generalized parallel PSO code

1. Edit input-pso fine to describe your problem
2. Edit the subroutine fitness to include your objective function
3. Compile with mpif90 or mpiifort : mpif90 PSO.f90 -o PSO.x | mpiifort PSO.f90 -o PSO.x
4. Run with mpirun -np X PSO.x

If the objective function is an openMP program:
For SLURM manager use:
  #SBATCH --cpus-per-task=X*Y - X is the numper of cores for PSO program and Y is the number of openmp cores for the objective function

