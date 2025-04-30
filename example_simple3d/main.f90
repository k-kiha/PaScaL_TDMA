program main
    use mpi
    use PaScaL_TDMA
    implicit none

    integer, parameter :: n1 = 256 ,n2 = 256 ,n3 = 1024
    integer, parameter :: np_dim(1:3) = (/1, 1, 8/)
    integer :: ierr, nprocs, myrank
    integer :: n1sub,n2sub,n3sub
    real*8, allocatable, dimension(:,:,:) :: a,b,c,d
    integer :: mpiutil_para
    real*8 :: timeA,timeB

    integer :: i, j, k
    integer :: indx_tmpa, indx_tmpb

    type(ptdma_plan_many) :: pz_many

    call MPI_Init(ierr)
    call MPI_Comm_size( MPI_COMM_WORLD, nprocs, ierr)
    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr)

    n1sub = mpiutil_para(1, n1, myrank, np_dim(1), indx_tmpa, indx_tmpb)
    n2sub = mpiutil_para(1, n2, myrank, np_dim(2), indx_tmpa, indx_tmpb)
    n3sub = mpiutil_para(1, n3, myrank, np_dim(3), indx_tmpa, indx_tmpb)

    write(*,*) "myrank=", myrank, " n1sub=", n1sub, " n2sub=", n2sub, " n3sub=", n3sub

    ![[ ==PaScaL_TDMA ====
        allocate(a(1:n1sub,1:n2sub,1:n3sub), b(1:n1sub,1:n2sub,1:n3sub))
        allocate(c(1:n1sub,1:n2sub,1:n3sub), d(1:n1sub,1:n2sub,1:n3sub))

        a(1:n1sub,1:n2sub,1:n3sub) = 1.d0
        b(1:n1sub,1:n2sub,1:n3sub) =-2.d0
        c(1:n1sub,1:n2sub,1:n3sub) = 1.d0
        d(1:n1sub,1:n2sub,1:n3sub) = 1.d0

        call PaScaL_TDMA_plan_many_create(pz_many, (n1sub*n2sub), myrank, nprocs, MPI_COMM_WORLD)
        timeA = MPI_Wtime()
        call PaScaL_TDMA_many_solve(pz_many, a,b,c,d,(n1sub*n2sub),n3sub)
        timeB = MPI_Wtime()
        call PaScaL_TDMA_plan_many_destroy(pz_many,nprocs)

        deallocate(a, b, c, d)

        write(*,*) "myrank=", myrank, " time=", timeB-timeA
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
    !==PaScaL_TDMA ==== ]]
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