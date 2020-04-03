# PPSO
A generalized parallel PSO code

1. Edit input-pso fine to describe your problem
2. Edit the subroutine fitness to include your objective function
3. Compile with mpif90 or mpiifort : mpif90 PSO.f90 -o PSO.x | mpiifort PSO.f90 -o PSO.x
4. Run with mpirun -np X PSO.x
