program main
    use mpi
    use PaScaL_TDMA
    implicit none
    integer :: ierr, nprocs, myrank

    integer :: i, j, k

    call MPI_Init(ierr)
    call MPI_Comm_size( MPI_COMM_WORLD, nprocs, ierr)
    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr)



    call MPI_Finalize(ierr)

end program