module mod_petsc


contains

   subroutine map_rhs(nkn,r)

     integer r(nkn)
     integer nkn

     return             

   end subroutine map_rhs          

   subroutine petsc_init(nkn_glob,nkn_in)

      use shympi

      implicit none

      integer,intent(in) :: nkn_glob, nkn_in

   end subroutine

   subroutine sol_petsc(nkn_in,n2nz,icoo,jcoo,ccoo,rvec,raux_pet)

      implicit none

      integer,intent(in)        :: n2nz, nkn_in
      integer                   :: i
      integer,intent(in)        :: icoo(n2nz),jcoo(n2nz)
      real*8, intent(in)        :: ccoo(n2nz)
      real*8, intent(in)        :: rvec(nkn_in)  
      double precision          :: raux_pet(nkn_in)

      return

   end subroutine sol_petsc        

   subroutine petsc_gather_2d_nodes_d(array,reb_array)

        use shympi

        implicit none

        double precision array(nkn_inner)
        double precision,optional :: reb_array(:)

        return

   end subroutine

   subroutine petsc_final

      implicit none
      
      return

   end subroutine

end module mod_petsc
