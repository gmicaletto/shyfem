
!==================================================================
        module mod_mpi_io
!==================================================================

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

        end subroutine rebuild_structures

!******************************************************************

        subroutine rebuild_ous_header

        use levels
        use mod_depth
        use basin

        implicit none
        
        end subroutine rebuild_ous_header

!******************************************************************

        subroutine rebuild_nos_header

        use levels
        use mod_depth
        use basin

        implicit none
        
        end subroutine rebuild_nos_header

!******************************************************************

        subroutine rebuild_scalar(what,array)

        use levels
        use basin
        !use mod_ts

        implicit none

        real array(nlvdi,nkn)
        
        character*(5) what

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

        return

        end subroutine

!==================================================================
        end module mod_mpi_io
!==================================================================

