
!==================================================================
        module mod_mpi_io
!==================================================================

        use shympi

        real, allocatable, dimension(:,:) :: inTempv
        real, allocatable, dimension(:,:) :: inSaltv
        integer, allocatable, dimension(:) :: inIlhkv
        real, allocatable, dimension(:) :: inHkv
        real, allocatable, dimension(:) :: inHev

        real, allocatable, dimension(:) :: outZnv,outV1v,outRdist
        real, allocatable, dimension(:) :: outHev
        real, allocatable, dimension(:,:) :: outZenv,outSaux
        integer, allocatable, dimension(:) :: outIlhv,outIlhkv
        real, allocatable, dimension(:,:) :: outUtlnv,outVtlnv
        real, allocatable, dimension(:,:) :: outSaltv,outTempv

        INTERFACE rebuild_2d_nodes
        	MODULE PROCEDURE  &
     			  rebuild_2d_nodes_i  &
     			  ,rebuild_2d_nodes_r
        END INTERFACE

        INTERFACE rebuild_3d_elems3
        	MODULE PROCEDURE  & 
     			  rebuild_3d_elems3_r
        END INTERFACE

        INTERFACE rebuild_2d_elems
        	MODULE PROCEDURE &
     			  rebuild_2d_elems_i &
     			  ,rebuild_2d_elems_r
        END INTERFACE

        INTERFACE rebuild_3d_nodes
        	MODULE PROCEDURE & 
     			  rebuild_3d_nodes_r
        END INTERFACE

        INTERFACE rebuild_3d_elems
        	MODULE PROCEDURE & 
     			  rebuild_3d_elems_r
        END INTERFACE

        contains


!******************************************************************

        subroutine rebuild_2d_nodes_i(array, reb_array)

        use basin
        use mpi_common_struct

        implicit none

        integer array(nkn)
        integer, allocatable, optional, dimension(:) :: reb_array

        if(my_id .eq. 0) then
          if(.not. allocated(reb_array)) then
            allocate(reb_array(nkndi))
          end if
          call rebuild_integ_2d(array, nkndi, numberNodes, 2,reb_array)
        else
          call rebuild_integ_2d(array, nkndi, numberNodes, 2)
        end if

        end subroutine rebuild_2d_nodes_i

!******************************************************************

        subroutine rebuild_2d_nodes_r(array, reb_array)

        use basin
        use mpi_common_struct

        implicit none

        real array(nkn)
        real, allocatable, optional, dimension(:) :: reb_array

        if(my_id .eq. 0) then
          if(.not. allocated(reb_array)) then
            allocate(reb_array(nkndi))
          end if
          call rebuild_real_2d(array, nkndi, numberNodes, 2,reb_array)
        else
          call rebuild_real_2d(array, nkndi, numberNodes, 2)
        end if


        end subroutine rebuild_2d_nodes_r

!******************************************************************

        subroutine rebuild_3d_elems3_r(array, reb_array)

        use basin
        use mpi_common_struct
        implicit none

        real array(3,nel)
        real, allocatable, optional, dimension(:,:) :: reb_array

        if(my_id .eq. 0) then
          if(.not. allocated(reb_array)) then
            allocate(reb_array(3,neldi))
          end if
          call rebuild_real_3d(array,3, neldi, numberElements,1,reb_array)
        else
          call rebuild_real_3d(array,3, neldi, numberElements,1)
        end if

        end subroutine rebuild_3d_elems3_r

!******************************************************************

        subroutine rebuild_2d_elems_i(array,reb_array)

        use basin
        use mpi_common_struct
        implicit none

        integer array(nel)
        integer, allocatable, optional, dimension(:) :: reb_array

        if(my_id .eq. 0) then
          if(.not. allocated(reb_array)) then
            allocate(reb_array(neldi))
          end if
          call rebuild_integ_2d(array, neldi, numberElements, 1,reb_array)
        else
          call rebuild_integ_2d(array, neldi, numberElements, 1)
        end if

        end subroutine rebuild_2d_elems_i

!******************************************************************

        subroutine rebuild_2d_elems_r(array,reb_array)

        use basin
        use mpi_common_struct
        implicit none

        real array(nel)
        real, allocatable, optional, dimension(:) :: reb_array

        if(my_id .eq. 0) then
          if(.not. allocated(reb_array)) then
            allocate(reb_array(neldi))
          end if
          call rebuild_real_2d(array, neldi, numberElements, 1,reb_array)
        else
          call rebuild_real_2d(array, neldi, numberElements, 1)
        end if

        end subroutine rebuild_2d_elems_r

!******************************************************************

        subroutine rebuild_3d_nodes_r(array, reb_array)

        use basin
        use levels
        use mpi_common_struct

        implicit none

        real array(nlvdi,nkn)
        real, allocatable, optional, dimension(:,:) :: reb_array

        if(my_id .eq. 0) then
          if(.not. allocated(reb_array)) then
            allocate(reb_array(nlvdi,nkndi))
          end if
          call rebuild_real_3d(array,nlvdi,nkndi,numberNodes, 2,reb_array)
        else
          call rebuild_real_3d(array,nlvdi,nkndi,numberNodes, 2)
        end if

        end subroutine rebuild_3d_nodes_r


!******************************************************************

        subroutine rebuild_3d_elems_r(array, reb_array)

        use basin
        use levels
        use mpi_common_struct

        implicit none

        real array(nlvdi,nel)
        real, allocatable, optional, dimension(:,:) :: reb_array

        if(my_id .eq. 0) then
          if(.not. allocated(reb_array)) then
            allocate(reb_array(nlvdi,neldi))
          end if
          call rebuild_real_3d(array,nlvdi, neldi, numberElements, 1,reb_array)
        else
          call rebuild_real_3d(array,nlvdi, neldi, numberElements, 1)
        end if

        end subroutine rebuild_3d_elems_r

!******************************************************************

        subroutine rebuild_structures

        use basin
        use mod_hydro
        use levels

        implicit none

        if(my_id .eq. 0) then
          call rebuild_2d_nodes(znv,outZnv)
          call rebuild_3d_elems3(zenv,outZenv)
          call rebuild_2d_elems(ilhv,outIlhv)
          call rebuild_3d_elems(utlnv,outUtlnv)
          call rebuild_3d_elems(vtlnv,outVtlnv)
        else
          call rebuild_2d_nodes(znv)
          call rebuild_3d_elems3(zenv)
          call rebuild_2d_elems(ilhv)
          call rebuild_3d_elems(utlnv)
          call rebuild_3d_elems(vtlnv)
        end if

        end subroutine rebuild_structures

!******************************************************************

        subroutine rebuild_ous_header

        use basin
        use levels
        use mod_depth

        implicit none
        
        if(my_id .eq. 0) then
          call rebuild_2d_elems(ilhv,outIlhv)
          call rebuild_2d_elems(hev,outHev)
        else
          call rebuild_2d_elems(ilhv)
          call rebuild_2d_elems(hev)
        end if

        end subroutine rebuild_ous_header

!******************************************************************

        subroutine rebuild_nos_header

        use basin
        use levels
        use mod_depth

        implicit none
        
        if(my_id .eq. 0) then
          call rebuild_2d_nodes(ilhkv, outIlhkv)
          call rebuild_2d_elems(hev,outHev)
        else
          call rebuild_2d_nodes(ilhkv)
          call rebuild_2d_elems(hev)
        end if

        

        end subroutine rebuild_nos_header

!******************************************************************

        subroutine rebuild_scalar(what,array)

        use basin
        use levels
        !use mod_ts

        implicit none

        real array(nlvdi,nkn)
        
        character*(5) what

        if(what .eq. 'saltv') then
          if(my_id .eq. 0) then
            call rebuild_2d_nodes(ilhkv, outIlhkv)
            call rebuild_3d_nodes(array,outSaltv)
          else
            call rebuild_2d_nodes(ilhkv)
            call rebuild_3d_nodes(array)
          end if
        else if(what .eq. 'tempv') then
          if(my_id .eq. 0) then
            call rebuild_2d_nodes(ilhkv, outIlhkv)
            call rebuild_3d_nodes(array,outTempv)
          else
            call rebuild_2d_nodes(ilhkv)
            call rebuild_3d_nodes(array)
          end if
        end if
        

        end subroutine rebuild_scalar

!******************************************************************

        subroutine compute_my_tracer(fullStruct,myStruct)

        use basin
        use levels
        use mpi_communication_struct

        implicit none

        integer i
        real fullStruct(nlvdi,nkndi)
        real myStruct(nlvdi,nkndi)

        do i=1,nkn
          myStruct(:,i) = fullStruct(:,mynodes%globalID(i))
        end do

        end subroutine compute_my_tracer

!******************************************************************

        subroutine shympi_ts_init(ilhkv_fem,hk_fem,he_fem,hlv_fem, &
                                        nkn_fem,nel_fem,nlv_fem)

        use levels, only: ilhkv,hlv,ilhv,nlvdi
        use mod_depth
        use basin, only: neldi,nkndi

        implicit none
        !argument
        integer nkn_fem,nel_fem,nlv_fem
        integer, allocatable, dimension(:) :: ilhkv_fem
        real, allocatable, dimension (:) :: hk_fem, he_fem, hlv_fem
        real, allocatable, dimension (:) :: new_hlv_fem

        integer ierr, hlv_dim

          allocate(inTempv(nlvdi,nkndi))
          allocate(inSaltv(nlvdi,nkndi))

          if(my_id .eq. 0) then
            call rebuild_2d_nodes(ilhkv_fem, inIlhkv)
            call rebuild_2d_nodes(hk_fem, inHkv)
            call rebuild_2d_elems(he_fem,inHev)
          else
            call rebuild_2d_nodes(ilhkv_fem)
            call rebuild_2d_nodes(hk_fem)
            call rebuild_2d_elems(he_fem)
          end if

          if(my_id .ne. 0) then
            allocate(inIlhkv(nkndi))
            allocate(inHkv(nkndi))
            allocate(inHev(neldi))
          end if

          hlv_dim = shympi_max(size(hlv_fem))

          if(size(hlv_fem) .lt. hlv_dim) then
            deallocate(hlv_fem)
            allocate(hlv_fem(hlv_dim))
            hlv_fem = 0.
            allocate(new_hlv_fem(hlv_dim))
          end if
          if(.not. allocated(new_hlv_fem)) then
            allocate(new_hlv_fem(size(hlv_fem)))
          end if


          call MPI_Bcast(inIlhkv(1),nkndi,MPI_INTEGER,0,MPI_COMM_WORLD ,ierr)

          call MPI_Bcast(inHkv(1),nkndi,MPI_REAL,0,MPI_COMM_WORLD ,ierr)

          call MPI_Bcast(inHev(1),neldi,MPI_REAL,0,MPI_COMM_WORLD ,ierr)

          call MPI_Allreduce(hlv_fem(1),new_hlv_fem(1),hlv_dim,MPI_REAL, &
                                MPI_MAX,MPI_COMM_WORLD,ierr)


          nkn_fem = nkndi
          nel_fem = neldi
          nlv_fem = hlv_dim

          deallocate(ilhkv_fem)
          allocate(ilhkv_fem(nkndi))
          deallocate(hk_fem)
          allocate(hk_fem(nkndi))
          ilhkv_fem=inIlhkv 
          hk_fem=inHkv 
          hlv_fem=new_hlv_fem

        return

        end subroutine

!******************************************************************

        subroutine rebuild_integ_2d(struct, size2, sizeArray, mode, newStruct)

        use mpi_common_struct

        implicit none

        integer size2, mode
        integer :: count_tag, i
        integer :: st(MPI_STATUS_SIZE), ierr
        integer :: sendbuffer
        integer,optional, dimension(size2) :: newStruct
        integer, dimension(size2) :: struct
        integer, dimension(n_threads) :: displsProc, sizeArray, recvbuffer, displs
        integer, allocatable, dimension(:) :: fullStruct

        if(my_id .eq. 0) then
          allocate(fullStruct(size2)) 
        end if

        sendbuffer = sizeArray(my_id+1)

        displs(1) = 0
        do i=2,n_threads
          displs(i) = displs(i-1) + sizeArray(i-1)
        end do

        call MPI_GATHERV(struct, sendbuffer, MPI_INTEGER, fullStruct, &
                sizeArray, displs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        if(my_id .eq. 0) then
          displsProc(1) = 0
          do i=2,n_threads
            displsProc(i) = displsProc(i-1) + sizeArray(i-1)
          end do
          if(mode .eq. 1) then
            call rebuild_integ_2d_int(fullStruct, displsProc, newStruct, allPartAssign, size2)
          else
            call rebuild_integ_2d_int(fullStruct, displsProc, newStruct, univocalNodesAssign, size2)
          end if
          deallocate(fullStruct)
        end if

      end subroutine

!******************************************************************

      subroutine rebuild_integ_2d_int(fullStruct, displsProc, newStruct, itemMap, size2)

        implicit none

        integer :: place, i, process, size2
        integer :: st(MPI_STATUS_SIZE), ierr
        integer, dimension(:) :: countProc(n_threads), displsProc(1)
        integer, dimension(:) :: fullStruct, newStruct
        integer itemMap(size2)

        ! countProc(1) = counter for items of the master process
        do i=1,n_threads
          countProc(i)=1
        end do

        do i=1, size2
          process = itemMap(i)
          place = displsProc(process+1) + countProc(process+1)
          newStruct(i) = fullStruct(place)
          countProc(process+1) = countProc(process+1) + 1
        end do
   
        return 

      end subroutine

!******************************************************************

      subroutine rebuild_real_2d(struct, size2, sizeArray, mode, newStruct)

        use mpi_common_struct

        implicit none

        integer size2, mode
        integer :: count_tag, i
        integer :: st(MPI_STATUS_SIZE), ierr
        integer :: sendbuffer
        real,optional, dimension(size2) :: newStruct
        real, dimension(size2) :: struct
        integer, dimension(n_threads) :: displsProc, sizeArray, recvbuffer, displs
        real, allocatable, dimension(:) :: fullStruct

        if(my_id .eq. 0) then
          allocate(fullStruct(size2)) 
        end if

        sendbuffer = sizeArray(my_id+1)

        displs(1) = 0
        do i=2,n_threads
          displs(i) = displs(i-1) + sizeArray(i-1)
        end do

        call MPI_GATHERV(struct, sendbuffer, MPI_REAL, fullStruct, &
                sizeArray, displs, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

        if(my_id .eq. 0) then
          displsProc(1) = 0
          do i=2,n_threads
            displsProc(i) = displsProc(i-1) + sizeArray(i-1)
          end do
          if(mode .eq. 1) then
            call rebuild_real_2d_int(fullStruct, displsProc, newStruct, allPartAssign, size2)
          else
            call rebuild_real_2d_int(fullStruct, displsProc, newStruct, univocalNodesAssign, size2)
          end if
          deallocate(fullStruct)
        end if

      end subroutine

!******************************************************************

      subroutine rebuild_real_2d_int(fullStruct, displsProc, newStruct, itemMap, size2)

        implicit none

        integer :: place, i, process, size2
        integer :: st(MPI_STATUS_SIZE), ierr
        integer, dimension(:) :: countProc(n_threads), displsProc(1)
        real, dimension(:) :: fullStruct, newStruct
        integer itemMap(size2)

        ! countProc(1) = counter for items of the master process
        do i=1,n_threads
          countProc(i)=1
        end do

        do i=1, size2
          process = itemMap(i)
          place = displsProc(process+1) + countProc(process+1)
          newStruct(i) = fullStruct(place)
          countProc(process+1) = countProc(process+1) + 1
        end do
   
        return 

      end subroutine

!******************************************************************

      subroutine rebuild_real_3d(struct, size1, size2, sizeArray, mode, newStruct)

        use mpi_common_struct

        implicit none

        integer size1, size2, mode
        integer :: count_tag, i
        integer :: st(MPI_STATUS_SIZE), ierr
        integer :: sendbuffer
        real, optional, dimension(size1,size2) :: newStruct
        real, dimension(size1,size2) :: struct
        integer, dimension(n_threads) :: displsProc, sizeArray, recvbuffer, displs
        real, allocatable, dimension(:,:) :: fullStruct

        if(my_id .eq. 0) then
          allocate(fullStruct(size1,size2)) 
          do i=1,n_threads
            recvbuffer(i) = sizeArray(i) * size1
          end do
        end if

        sendbuffer = size1 * sizeArray(my_id+1)

        displs(1) = 0
        do i=2,n_threads
          displs(i) = displs(i-1) + size1 * sizeArray(i-1)
        end do

        call MPI_GATHERV(struct, sendbuffer, MPI_REAL, fullStruct, &
                recvbuffer, displs, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

        if(my_id .eq. 0) then
          displsProc(1) = 0
          do i=2,n_threads
            displsProc(i) = displsProc(i-1) + sizeArray(i-1)
          end do
          if(mode .eq. 1) then
            call rebuild_real_3d_int(fullStruct,displsProc,newStruct,allPartAssign,size1,size2)
          else
            call rebuild_real_3d_int(fullStruct,displsProc,newStruct,univocalNodesAssign,size1,size2)
          end if
          deallocate(fullStruct)
        end if

      end subroutine

!******************************************************************

      subroutine rebuild_real_3d_int(fullStruct, displsProc, newStruct, itemMap, size1, size2)

        implicit none

        integer :: place, i, process, j, size1, size2
        integer :: st(MPI_STATUS_SIZE), ierr
        integer, dimension(:) :: countProc(n_threads), displsProc(1)
        real, dimension(:,:) :: fullStruct, newStruct
        integer itemMap(size2)
    
        ! countProc(1) = counter for items of the master process
        do i=1,n_threads
          countProc(i)=1
        end do

        do i=1, size2
          process = itemMap(i)
          place = displsProc(process+1) + countProc(process+1)
          do j=1,size1
            newStruct(j,i) = fullStruct(j,place)
          end do
          countProc(process+1) = countProc(process+1) + 1
        end do

      end subroutine

!******************************************************************

!==================================================================
        end module mod_mpi_io
!==================================================================



