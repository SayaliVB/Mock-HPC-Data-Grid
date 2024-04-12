OpenMP_pollutant_average : ./omp "path to data folder"
MPI_min_max_AQI : mpirun -n 8 mpi "path to data folder" no_of_recods order_asc_desc
MPI_total_average : mpirun -n 8 mpi "path to data folder"
