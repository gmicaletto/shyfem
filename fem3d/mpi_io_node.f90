
!==================================================================
        module mod_mpi_io
!==================================================================

        use shympi
        use mpi

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
        	MODULE PROCEDURE & 
     			  rebuild_2d_nodes_i &
     			  ,rebuild_2d_nodes_r
        END INTERFACE

        INTERFACE rebuild_3d_elems3
        	MODULE PROCEDURE &
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
        implicit none

        integer array(nkn)
        integer, allocatable, optional, dimension(:) :: reb_array

        end subroutine rebuild_2d_nodes_i

!******************************************************************

        subroutine rebuild_2d_nodes_r(array, reb_array)

        use basin
        implicit none

        real array(nkn)
        real, allocatable, optional, dimension(:) :: reb_array

        end subroutine rebuild_2d_nodes_r

!******************************************************************

        subroutine rebuild_3d_elems3_r(array, reb_array)

        use basin
        implicit none

        real array(3,nel)
        real, allocatable, optional, dimension(:,:) :: reb_array

        end subroutine rebuild_3d_elems3_r

!******************************************************************

        subroutine rebuild_2d_elems_i(array,reb_array)

        use basin
        implicit none

        integer array(nel)
        integer, allocatable, optional, dimension(:) :: reb_array

        end subroutine rebuild_2d_elems_i

!******************************************************************

        subroutine rebuild_2d_elems_r(array,reb_array)

        use basin
        implicit none

        real array(nel)
        real, allocatable, optional, dimension(:) :: reb_array

        end subroutine rebuild_2d_elems_r

!******************************************************************

        subroutine rebuild_3d_nodes_r(array, reb_array)

        use basin
        use levels

        implicit none

        real array(nlvdi,nkn)
        real, allocatable, optional, dimension(:,:) :: reb_array

        end subroutine rebuild_3d_nodes_r


!******************************************************************

        subroutine rebuild_3d_elems_r(array, reb_array)

        use basin
        use levels

        implicit none

        real array(nlvdi,nel)
        real, allocatable, optional, dimension(:,:) :: reb_array

        end subroutine rebuild_3d_elems_r

!******************************************************************

        subroutine rebuild_structures

        use mod_hydro
        use levels
        use basin

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

        use levels
        use mod_depth
        use basin

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

        use levels
        use mod_depth
        use basin

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

        use levels
        use basin
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

        use levels
        use basin

        implicit none

        real fullStruct(nlvdi,nkndi)
        real myStruct(nlvdi,nkndi)

        end subroutine compute_my_tracer

!******************************************************************

        subroutine shympi_ts_init(ilhkv_fem,hk_fem,he_fem,hlv_fem, &
                                nkn_fem,nel_fem,nlv_fem)

        use levels, only: ilhkv,hlv,ilhv,nlvdi
        use mod_depth
        use basin, only: neldi,nkndi

        implicit none
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


          call MPI_Bcast(inIlhkv(1),nkndi,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

          call MPI_Bcast(inHkv(1),nkndi,MPI_REAL,0,MPI_COMM_WORLD ,ierr)

          call MPI_Bcast(inHev(1),neldi,MPI_REAL,0,MPI_COMM_WORLD ,ierr)

          call MPI_Allreduce(hlv_fem(1),new_hlv_fem(1),hlv_dim,MPI_REAL,&
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

!==================================================================
        end module mod_mpi_io
!==================================================================

