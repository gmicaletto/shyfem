module mod_petsc

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscviewer.h>
#include <petsc/finclude/petscvec.h90>

   Mat                :: mat_solver
   Vec                :: rhs, sol
   KSP                :: ksp
   integer,allocatable,dimension(:),save :: rhs_glob

contains

   subroutine map_rhs(nkn,r)

     use shympi

     integer r(nkn)
     integer i,nkn

     if(bmpi) then        
       do i=1, nkn
          r(i) = mynodes%globalID(i)-1
       end do
     else           
       do i=1, nkn
          r(i) = i-1
       end do
     end if
    
     return             

   end subroutine map_rhs          

   subroutine petsc_init(nkn_glob,nkn_in)

      use shympi
      use basin

      implicit none

      integer,intent(in) :: nkn_glob, nkn_in
      PetscInt           :: ierr

      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

      call MatCreate( PETSC_COMM_WORLD, mat_solver, ierr) 
      call MatSetSizes(mat_solver,nkn_in,nkn_in,nkn_glob,nkn_glob,ierr )  
      call MatSetType(mat_solver, MATMPIAIJ, ierr)  
      call MatMPIAIJSetPreallocation(mat_solver,ngr+1,PETSC_NULL_INTEGER,ngr,PETSC_NULL_INTEGER,ierr)  
      
  
      call VecCreate( PETSC_COMM_WORLD, rhs, ierr)
      call VecSetSizes(rhs, nkn_in, nkn_glob,ierr)  
      call VecSetFromOptions(rhs,ierr)
      call VecDuplicate(rhs, sol, ierr)  

      call KSPCreate( PETSC_COMM_WORLD, ksp, ierr)
      call KSPSetOperators(ksp, mat_solver, mat_solver, ierr)
      call KSPSetFromOptions(ksp, ierr)

      allocate(rhs_glob(nkn_in))
      call map_rhs(nkn_in,rhs_glob)  

   end subroutine

   subroutine sol_petsc(nkn_in,n2nz,icoo,jcoo,ccoo,rvec,raux_pet)

      use shympi

      implicit none

      integer,intent(in)        :: n2nz, nkn_in
      integer                   :: i
      integer,intent(in)        :: icoo(n2nz),jcoo(n2nz)
      real*8, intent(in)        :: ccoo(n2nz)
      double precision          :: rvec(nkn_in)  
      double precision          :: raux_pet(nkn_in)

      double precision, pointer :: xx_v(:)
      PetscInt          :: ierr, its
      Vec               :: residual
  
      do i=1,n2nz 
        call MatSetValue(mat_solver,icoo(i)-1,jcoo(i)-1,ccoo(i),ADD_VALUES,ierr)
      end do  

      call MatAssemblyBegin(mat_solver,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(mat_solver,MAT_FINAL_ASSEMBLY,ierr)


      call MPI_Barrier(PETSC_COMM_WORLD,ierr)

      do i=1,nkn_in
        call VecSetValue(rhs,rhs_glob(i),rvec(i),INSERT_VALUES,ierr)
      end do    

      call MPI_Barrier(PETSC_COMM_WORLD,ierr)
        
      call VecAssemblyBegin(rhs, ierr)  
      call VecAssemblyEnd(rhs, ierr)  

      call KSPSolve(ksp, rhs, sol, ierr)
      call KSPGetIterationNumber(ksp, its, ierr)

      call VecGetArrayReadF90(sol,xx_v,ierr)

      do i=1,nkn_in
        raux_pet(i) = xx_v(i)     
      end do

      call VecRestoreArrayF90(sol,xx_v,ierr)

      call MatZeroEntries(mat_solver,ierr)

      if (b_use_mpi) rvec = raux_pet

   end subroutine sol_petsc        

   subroutine petsc_gather_2d_nodes_d(array,reb_array)

        use shympi
        use basin
        use mpi_common_struct

        implicit none

        integer sendbuffer,ierr,i
        integer, dimension(n_threads) :: recvbuffer, displs

        double precision array(nkn_inner)
        double precision,allocatable,dimension(:),optional :: reb_array(:)

        if ( .not. bmpi ) return

        if(my_id .eq. 0) then
          if(.not. allocated(reb_array)) then
            allocate(reb_array(nkndi))
          end if
        end if

        sendbuffer = numberNodes(my_id+1)

        displs(1) = 0
        do i=2,n_threads
          displs(i) = displs(i-1) + numberNodes(i-1)
        end do
 
        call MPI_GATHERV(array, sendbuffer, MPI_DOUBLE_PRECISION, &
                      reb_array, numberNodes, displs, &
                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        return

   end subroutine

   subroutine petsc_final

      implicit none
      PetscInt ierr
      
      call MatDestroy(mat_solver, ierr)  
      call VecDestroy(rhs, ierr)  
      call VecDestroy(sol, ierr)  
      call KSPDestroy(ksp, ierr)  
      call PetscFinalize(ierr)  

   end subroutine
!
!   subroutine mumps_init(nkn_glob, n2nz,icoo, jcoo)
!
!      implicit none
!      integer nkn_glob,n2nz
!      integer :: icoo(:),jcoo(:)
!
!      return
!
!   end subroutine
!
!   subroutine mumps_solve(iteration,acoo,rvec,rhs_map,raux)
!
!      implicit none
!      integer iteration
!      double precision :: acoo(:),rvec(:),raux(:)
!      integer :: rhs_map
!
!      return
!
!   end subroutine
!
!   subroutine mumps_finalize
!
!      implicit none
!
!      return
!
!   end subroutine
!
!   subroutine mumps_gather_2d_nodes_d(raux2d,total_raux2d)
!
!      implicit none
!      double precision :: raux2d(:)
!      double precision,optional :: total_raux2d(:)
!
!      return
!
!   end subroutine
!
end module mod_petsc
