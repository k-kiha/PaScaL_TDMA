program main
    use mpi
    use PaScaL_TDMA
    implicit none

    integer, parameter :: n1 = 1024 ,n2 = 1024 ,n3 = 1025
    integer, parameter :: np_dim(1:3) = (/1, 1, 4/)
    integer :: ierr, nprocs, myrank
    integer :: n1sub,n2sub,n3sub
    real*8, allocatable, dimension(:,:,:) :: a,b,c,d
    integer :: mpiutil_para

    integer :: i, j, k
    integer :: indx_tmpa, indx_tmpb

    call MPI_Init(ierr)
    call MPI_Comm_size( MPI_COMM_WORLD, nprocs, ierr)
    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr)

    n1sub = mpiutil_para(1, n1, myrank, np_dim(1), indx_tmpa, indx_tmpb)
    n2sub = mpiutil_para(1, n2, myrank, np_dim(2), indx_tmpa, indx_tmpb)
    n3sub = mpiutil_para(1, n3, myrank, np_dim(3), indx_tmpa, indx_tmpb)

    write(*,*) "myrank=", myrank, " n1sub=", n1sub, " n2sub=", n2sub, " n3sub=", n3sub

    allocate(a(1:n1sub,1:n2sub,1:n3sub), b(1:n1sub,1:n2sub,1:n3sub))
    allocate(c(1:n1sub,1:n2sub,1:n3sub), d(1:n1sub,1:n2sub,1:n3sub))



    deallocate(a, b, c, d)
    call MPI_Finalize(ierr)

end program



function mpiutil_para(sta_g, end_g, myrank, nprocs, indx_a, indx_b)result(nsub)
    integer :: sta_g, end_g, myrank, nprocs, indx_a, indx_b
    integer :: nsub

    integer :: n, tmp1, tmp2, aa,bb

    n    = end_g-sta_g+1
    tmp1 = int(n/nprocs)
    tmp2 = mod(n,nprocs)
    
    if ( myrank<tmp2 ) then
        aa = myrank
        bb = 1
    else
        aa = tmp2
        bb = 0
    endif 

    indx_a = myrank*tmp1 + aa + sta_g;
    indx_b = indx_a + tmp1 + bb -1;
    nsub   = tmp1 + bb;
end function mpiutil_para