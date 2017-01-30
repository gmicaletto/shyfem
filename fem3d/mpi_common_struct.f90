
!#######################################################################!
!*******************  start module mpi_common_struct  ******************!
!***********************************************************************!
!#######################################################################!
!------------------------------------------------------------------------
! mpi_common_struct - common structures are declared in this module
! allowing to eliminate the deprecate common block in order to use
! the more flexible allocatable structures 
!------------------------------------------------------------------------

module mpi_common_struct

   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: myieltv
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: numberElements
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: numberNodes
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: procNodes
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: numberLevels
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: allNodesAssign
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: univocalNodesAssign
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: scatterNodes

   integer,public,save :: maxelements,maxlevelsproc,totalnodes,maxNodes

end module mpi_common_struct

!------------------------------------------------------------------------
! end of mpi_common_struct module
!------------------------------------------------------------------------
