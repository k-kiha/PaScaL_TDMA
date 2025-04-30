program main
    use mpi
    use PaScaL_TDMA
    implicit none

    integer, parameter :: n1 = 3 ,n2 = 3 ,n3 = 15
    ! integer, parameter :: np_dim(1:3) = (/1, 1, 8/)
    integer :: ierr, nprocs, myrank
    integer :: n1sub,n2sub,n3sub,n12ssub
    real*8, allocatable, dimension(:,:) :: a,b,c,d
    real*8, allocatable, dimension(:,:) :: d_center
    integer :: mpiutil_para
    real*8 :: timeA,timeB

    integer :: i, j, k, index, count
    integer :: indx_tmpa, indx_tmpb

    type(ptdma_plan_many) :: pz_many

    integer, allocatable, dimension(:,:) :: ssscount_c2z, rrrcount_c2z, sssdist_c2z, rrrdist_c2z
    integer, allocatable, dimension(:,:) :: ssscount_z2c, rrrcount_z2c, sssdist_z2c, rrrdist_z2c

    integer, allocatable, dimension(:) :: sendcount_c2z, recvcount_c2z, senddist_c2z, recvdist_c2z
    integer, allocatable, dimension(:) :: sendcount_z2c, recvcount_z2c, senddist_z2c, recvdist_z2c

    real*8, allocatable, dimension(:) :: packbuf_c2z, unpackbuf_c2z
    real*8, allocatable, dimension(:) :: packbuf_z2c, unpackbuf_z2c

    call MPI_Init(ierr)
    call MPI_Comm_size( MPI_COMM_WORLD, nprocs, ierr)
    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr)

    n1sub = n1
    n2sub = n2
    n3sub = mpiutil_para(1, n3, myrank, nprocs, indx_tmpa, indx_tmpb)

    write(*,*) "myrank=", myrank, " n1sub=", n1sub, " n2sub=", n2sub, " n3sub=", n3sub

    !=====================
    !=====================
    ![[ ==PaScaL_TDMA ====
        allocate(a(1:n1sub*n2sub,1:n3sub), b(1:n1sub*n2sub,1:n3sub))
        allocate(c(1:n1sub*n2sub,1:n3sub), d(1:n1sub*n2sub,1:n3sub))

        a(1:n1sub*n2sub,1:n3sub) = 1.d0
        b(1:n1sub*n2sub,1:n3sub) =-2.d0
        c(1:n1sub*n2sub,1:n3sub) = 1.d0
        d(1:n1sub*n2sub,1:n3sub) = 1.d0

        call PaScaL_TDMA_plan_many_create(pz_many, (n1sub*n2sub), myrank, nprocs, MPI_COMM_WORLD)

        timeA = MPI_Wtime()
        call PaScaL_TDMA_many_solve(pz_many, a,b,c,d,(n1sub*n2sub),n3sub)
        timeB = MPI_Wtime()
        
        call PaScaL_TDMA_plan_many_destroy(pz_many,nprocs)

        deallocate(a, b, c, d)

        write(*,*) "myrank=", myrank, " time=", timeB-timeA
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
    !==PaScaL_TDMA ==== ]]
    !=====================
    !=====================


    !=======================
    !=======================
    ![[ ==Oridinal TDMA ====
        allocate(d_center(1:n1sub*n2sub,1:n3sub))

        n12ssub = mpiutil_para(1, n1sub*n2sub, myrank, nprocs, indx_tmpa, indx_tmpb)
        allocate(a(1:n12ssub,1:n3sub), b(1:n12ssub,1:n3sub))
        allocate(c(1:n12ssub,1:n3sub), d(1:n12ssub,1:n3sub))

        allocate(ssscount_c2z(1:2,0:nprocs-1), rrrcount_c2z(1:2,0:nprocs-1))
        allocate(ssscount_z2c(1:2,0:nprocs-1), rrrcount_z2c(1:2,0:nprocs-1))
        allocate(sssdist_c2z(1:2,0:nprocs-1), rrrdist_c2z(1:2,0:nprocs-1))
        allocate(sssdist_z2c(1:2,0:nprocs-1), rrrdist_z2c(1:2,0:nprocs-1))
        allocate(sendcount_c2z(0:nprocs-1), recvcount_c2z(0:nprocs-1))
        allocate(sendcount_z2c(0:nprocs-1), recvcount_z2c(0:nprocs-1))
        allocate(senddist_c2z(0:nprocs-1), recvdist_c2z(0:nprocs-1))
        allocate(senddist_z2c(0:nprocs-1), recvdist_z2c(0:nprocs-1))

        do i = 0, nprocs-1
            ssscount_c2z(1,i) = mpiutil_para(1, (n1sub*n2sub), i     , nprocs, indx_tmpa, indx_tmpb)
            ssscount_c2z(2,i) = n3sub
            rrrcount_c2z(1,i) = mpiutil_para(1, (n1sub*n2sub), myrank, nprocs, indx_tmpa, indx_tmpb)
            rrrcount_c2z(2,i) = mpiutil_para(1, n3           , i     , nprocs, indx_tmpa, indx_tmpb)
        end do

        do i = 0, nprocs-1
            sssdist_c2z(1,i) = sum(ssscount_c2z(1,0:i)) - ssscount_c2z(1,i)
            sssdist_c2z(2,i) = 0
            rrrdist_c2z(1,i) = 0
            rrrdist_c2z(2,i) = sum(rrrcount_c2z(2,0:i)) - rrrcount_c2z(2,i)
        end do

        ssscount_z2c(:,:) = rrrcount_c2z(:,:)
        rrrcount_z2c(:,:) = ssscount_c2z(:,:)

        do i = 0, nprocs-1
            sendcount_c2z(i) = ssscount_c2z(1,i)*ssscount_c2z(2,i)
            recvcount_c2z(i) = rrrcount_c2z(1,i)*rrrcount_c2z(2,i)
        end do
        do i = 0, nprocs-1
            senddist_c2z(i) = sum(sendcount_c2z(0:i)) - sendcount_c2z(i)
            recvdist_c2z(i) = sum(recvcount_c2z(0:i)) - recvcount_c2z(i)
        end do

        sendcount_z2c(:) = recvcount_c2z(:)
        recvcount_z2c(:) = sendcount_c2z(:)
        senddist_z2c(:) = recvdist_c2z(:)
        recvdist_z2c(:) = senddist_c2z(:)

        write(*,'(1A,1I3,1A,4I5)') "myrank=", myrank, " sendcount_c2z=", sendcount_c2z
        write(*,'(1A,1I3,1A,4I5)') "myrank=", myrank, " recvcount_c2z=", recvcount_c2z
        write(*,'(1A,1I3,1A,4I5)') "myrank=", myrank, " senddist_c2z=", senddist_c2z
        write(*,'(1A,1I3,1A,4I5)') "myrank=", myrank, " recvdist_c2z=", recvdist_c2z

        allocate(packbuf_c2z(0:sum(sendcount_c2z(:))-1), unpackbuf_c2z(0:sum(recvcount_c2z(:))-1))
        allocate(packbuf_z2c(0:sum(sendcount_z2c(:))-1), unpackbuf_z2c(0:sum(recvcount_z2c(:))-1))


        a(1:n12ssub,1:n3sub) = 1.d0
        b(1:n12ssub,1:n3sub) =-2.d0
        c(1:n12ssub,1:n3sub) = 1.d0
        d(1:n12ssub,1:n3sub) = 0.d0

        do k = 1, n3sub
        do i = 1, n1sub*n2sub
            d_center(i,k) = i*(100) + k
        end do
        end do
        

        ! alltoall pack
        count = 0
        do index = 0, nprocs-1
            do k = 1, ssscount_c2z(2,index)
            do i = 1, ssscount_c2z(1,index)
                packbuf_c2z(count) = d_center(sssdist_c2z(1,index)+ i, sssdist_c2z(2,index)+ k)
                count = count+1
            end do
            end do
        end do


        do index = 0, nprocs-1
            if(myrank==index) write(*,*) "myrank=", myrank, packbuf_c2z
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
        end do
        
        
        ! ! alltoall c to z

        ! ! alltoall unpack
        ! count = 0
        ! do index = 0, nprocs-1
        !     do k = 1, rrrcount_c2z(2,index)
        !     do i = 1, rrrcount_c2z(1,index)
        !         d(rrrdist_c2z(1,index)+ i, rrrdist_c2z(2,index)+ k) = unpackbuf_c2z(count) 
        !         count = count+1
        !     end do
        !     end do
        ! end do

        ! ! tdma many
        
        ! ! alltoall pack
        
        ! ! alltoall c to z
        
        ! ! alltoall unpack
        
        deallocate(ssscount_c2z, rrrcount_c2z)
        deallocate(ssscount_z2c, rrrcount_z2c)
        deallocate(sssdist_c2z, rrrdist_c2z)
        deallocate(sssdist_z2c, rrrdist_z2c)
        deallocate(packbuf_c2z, unpackbuf_c2z)
        deallocate(packbuf_z2c, unpackbuf_z2c)
        deallocate(senddist_c2z, recvdist_c2z)
        deallocate(senddist_z2c, recvdist_z2c)
        deallocate(sendcount_c2z, recvcount_c2z)
        deallocate(sendcount_z2c, recvcount_z2c)

        deallocate(a, b, c, d)

        deallocate(d_center)
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
    !==Oridinal TDMA ==== ]]
    !=======================
    !=======================
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