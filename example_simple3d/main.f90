program main
    use mpi
    use PaScaL_TDMA
    implicit none

    ! integer, parameter :: n1 = 256 ,n2 = 256 ,n3 = 1024
    integer :: n1 ,n2 ,n3

    integer :: ierr, nprocs, myrank
    integer :: n1sub,n2sub,n3sub,n12ssub
    real*8, allocatable, dimension(:,:) :: a,b,c,d
    real*8, allocatable, dimension(:,:) :: d_center
    integer :: mpiutil_para
    real*8 :: timeA,timeB
    real*8 :: timeA0,timeB0
    real*8 :: timeA_pack1,timeB_pack1
    real*8 :: timeA_unpack1,timeB_unpack1
    real*8 :: timeA_comm1,timeB_comm1
    real*8 :: timeA_tdma,timeB_tdma
    real*8 :: timeA_pack2,timeB_pack2
    real*8 :: timeA_unpack2,timeB_unpack2
    real*8 :: timeA_comm2,timeB_comm2

    integer :: i, j, k, index, count
    integer :: indx_tmpa, indx_tmpb

    type(ptdma_plan_many) :: pz_many
    
    integer :: ios
    character*8 :: tmpchar
    integer, allocatable, dimension(:,:) :: ssscount_c2z, rrrcount_c2z, sssdist_c2z, rrrdist_c2z
    integer, allocatable, dimension(:,:) :: ssscount_z2c, rrrcount_z2c, sssdist_z2c, rrrdist_z2c

    integer, allocatable, dimension(:) :: sendcount_c2z, recvcount_c2z, senddist_c2z, recvdist_c2z
    integer, allocatable, dimension(:) :: sendcount_z2c, recvcount_z2c, senddist_z2c, recvdist_z2c

    real*8, allocatable, dimension(:) :: packbuf_c2z, unpackbuf_c2z
    real*8, allocatable, dimension(:) :: packbuf_z2c, unpackbuf_z2c

    call MPI_Init(ierr)
    call MPI_Comm_size( MPI_COMM_WORLD, nprocs, ierr)
    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr)
    ! Read values from the input file "in.in"
    
    open(unit=10, file="in.in", status="old", action="read", iostat=ios)
    read(10,*) tmpchar, n1
    read(10,*) tmpchar, n2
    read(10,*) tmpchar, n3
    close(10)

    n1sub = n1
    n2sub = n2
    n3sub = mpiutil_para(1, n3, myrank, nprocs, indx_tmpa, indx_tmpb)

    ! write(*,*) "myrank=", myrank, " n1sub=", n1sub, " n2sub=", n2sub, " n3sub=", n3sub

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

        timeA0 = MPI_Wtime()
        call PaScaL_TDMA_many_solve(pz_many, a,b,c,d,(n1sub*n2sub),n3sub)
        timeB0 = MPI_Wtime()
        
        call PaScaL_TDMA_plan_many_destroy(pz_many,nprocs)

        deallocate(a, b, c, d)

        if ( myrank == 0 ) write(*,*) "myrank=", myrank, " time=", timeB0-timeA0
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
    !==PaScaL_TDMA ==== ]]
    !=====================
    !=====================

    if ( myrank == 0 ) write(*,*) "~~~~~~~~~~~"

    !=======================
    !=======================
    ![[ ==Oridinal TDMA ====
        allocate(d_center(1:n1sub*n2sub,1:n3sub))

        n12ssub = mpiutil_para(1, n1sub*n2sub, myrank, nprocs, indx_tmpa, indx_tmpb)
        allocate(a(1:n12ssub,1:n3), b(1:n12ssub,1:n3))
        allocate(c(1:n12ssub,1:n3), d(1:n12ssub,1:n3))

        !~~~ alltoall info ~~~~
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
            sssdist_z2c(:,:) = rrrdist_c2z(:,:)
            rrrdist_z2c(:,:) = sssdist_c2z(:,:)

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

            allocate(packbuf_c2z(0:sum(sendcount_c2z(:))-1), unpackbuf_c2z(0:sum(recvcount_c2z(:))-1))
            allocate(packbuf_z2c(0:sum(sendcount_z2c(:))-1), unpackbuf_z2c(0:sum(recvcount_z2c(:))-1))

        !~~~ alltoall info ~~~~

        a(1:n12ssub,1:n3sub) = 1.d0
        b(1:n12ssub,1:n3sub) =-2.d0
        c(1:n12ssub,1:n3sub) = 1.d0
        d(1:n12ssub,1:n3sub) = 0.d0

        do k = 1, n3sub
        do i = 1, n1sub*n2sub
            d_center(i,k) = dble(i*(100) + k)
        end do
        end do
        
        timeA = MPI_Wtime()
        ! alltoall pack
        timeA_pack1     = MPI_Wtime()
        count = 0
        do index = 0, nprocs-1
            do k = 1, ssscount_c2z(2,index)
            do i = 1, ssscount_c2z(1,index)
                packbuf_c2z(count) = d_center(sssdist_c2z(1,index)+ i, sssdist_c2z(2,index)+ k)
                count = count+1
            end do
            end do
        end do
        timeB_pack1     = MPI_Wtime()
        
        ! alltoall c to z
        timeA_comm1     = MPI_Wtime()
        call MPI_Alltoallv(  packbuf_c2z, sendcount_c2z, senddist_c2z, MPI_DOUBLE, &
                           unpackbuf_c2z, recvcount_c2z, recvdist_c2z, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
        timeB_comm1     = MPI_Wtime()
        ! alltoall unpack
        timeA_unpack1   = MPI_Wtime()
        count = 0
        do index = 0, nprocs-1
            do k = 1, rrrcount_c2z(2,index)
            do i = 1, rrrcount_c2z(1,index)
                d(rrrdist_c2z(1,index)+ i, rrrdist_c2z(2,index)+ k) = unpackbuf_c2z(count) 
                count = count+1
            end do
            end do
        end do
        timeB_unpack1   = MPI_Wtime()

        ! tdma many
        timeA_tdma      = MPI_Wtime()
        call tdma_many(A,B,C,D, n12ssub, n3)
        timeB_tdma      = MPI_Wtime()

        ! alltoall pack
        timeA_pack2     = MPI_Wtime()
        count = 0
        do index = 0, nprocs-1
            do k = 1, ssscount_z2c(2,index)
            do i = 1, ssscount_z2c(1,index)
                packbuf_z2c(count) = d(sssdist_z2c(1,index)+ i, sssdist_z2c(2,index)+ k)
                count = count+1
            end do
            end do
        end do
        timeB_pack2     = MPI_Wtime()
        
        ! alltoall c to z
        timeA_comm2     = MPI_Wtime()
        call MPI_Alltoallv(  packbuf_z2c, sendcount_z2c, senddist_z2c, MPI_DOUBLE, &
                           unpackbuf_z2c, recvcount_z2c, recvdist_z2c, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
        timeB_comm2     = MPI_Wtime()
        
        ! alltoall unpack
        timeA_unpack2   = MPI_Wtime()
        count = 0
        do index = 0, nprocs-1
            do k = 1, rrrcount_z2c(2,index)
            do i = 1, rrrcount_z2c(1,index)
                d_center(rrrdist_z2c(1,index)+ i, rrrdist_z2c(2,index)+ k) = unpackbuf_z2c(count) 
                count = count+1
            end do
            end do
        end do
        timeB_unpack2   = MPI_Wtime()
        
        timeB = MPI_Wtime()
        if ( myrank == 0 ) then
            write(*,*) "myrank=", myrank, " time=", timeB-timeA
            write(*,'(1A13,1A1,F30.20)') "time",":", (timeB_pack1 - timeA_pack1)    &
                                                    +(timeB_unpack1 - timeA_unpack1)&
                                                    +(timeB_comm1 - timeA_comm1)    &
                                                    +(timeB_tdma - timeA_tdma)      &
                                                    +(timeB_pack2 - timeA_pack2)    &
                                                    +(timeB_unpack2 - timeA_unpack2)&
                                                    +(timeB_comm2 - timeA_comm2)
            write(*,'(1A13,1A1,F30.20)') "time_pack1",":",timeB_pack1 - timeA_pack1
            write(*,'(1A13,1A1,F30.20)') "time_unpack1",":",timeB_unpack1 - timeA_unpack1
            write(*,'(1A13,1A1,F30.20)') "time_comm1",":",timeB_comm1 - timeA_comm1
            write(*,'(1A13,1A1,F30.20)') "time_tdma",":",timeB_tdma - timeA_tdma
            write(*,'(1A13,1A1,F30.20)') "time_pack2",":",timeB_pack2 - timeA_pack2
            write(*,'(1A13,1A1,F30.20)') "time_unpack2",":",timeB_unpack2 - timeA_unpack2
            write(*,'(1A13,1A1,F30.20)') "time_comm2",":",timeB_comm2 - timeA_comm2
        end if
        call MPI_Barrier(MPI_COMM_WORLD, ierr)

        if ( myrank == 0 ) then
            open(unit=20, file="out.txt", status="unknown", position="append", action="write")
            write(20,'(5I15,9E30.20)') nprocs,n1sub,n2sub,n3,n3sub   &
            ,timeB0-timeA0                  &
            ,timeB-timeA                    &
            ,timeB_pack1 - timeA_pack1      &
            ,timeB_unpack1 - timeA_unpack1  &
            ,timeB_comm1 - timeA_comm1      &
            ,timeB_tdma - timeA_tdma        &
            ,timeB_pack2 - timeA_pack2      &
            ,timeB_unpack2 - timeA_unpack2  &
            ,timeB_comm2 - timeA_comm2
            close(20)
        end if
        
        !~~~ alltoall info ~~~~
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
        !~~~ alltoall info ~~~~

        deallocate(a, b, c, d)

        deallocate(d_center)
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
    !==Oridinal TDMA ==== ]]
    !=======================
    !=======================
    call MPI_Finalize(ierr)

end program

! subroutine tdma_many(a, b, c, d, n1, n2)

!     implicit none

!     integer, intent(in) :: n1,n2
!     double precision, intent(inout) :: a(n1,n2), b(n1,n2), c(n1,n2), d(n1,n2)
    
!     integer :: i,j
!     double precision, allocatable, dimension(:) :: r

!     allocate(r(1:n1))

!     do i=1,n1
!         d(i,1)=d(i,1)/b(i,1)
!         c(i,1)=c(i,1)/b(i,1)
!     enddo

!     do j=2,n2
!         do i=1,n1
!             r(i)=1.d0/(b(i,j)-a(i,j)*c(i,j-1))
!             d(i,j)=r(i)*(d(i,j)-a(i,j)*d(i,j-1))
!             c(i,j)=r(i)*c(i,j)
!         enddo
!     enddo

!     do j=n2-1,1,-1
!         do i=1,n1
!             d(i,j)=d(i,j)-c(i,j)*d(i,j+1)
!         enddo
!     enddo

!     deallocate(r)

! end subroutine tdma_many

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