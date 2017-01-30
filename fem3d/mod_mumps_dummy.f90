module mod_mumps

   integer,allocatable,dimension(:),save :: loc_rhs_map,total_rhs_map

contains

  subroutine map_rhs(nkn,r)

    integer  r(nkn)
    integer  nkn

    return               

  end subroutine map_rhs          
                
  subroutine mumps_init(nkn_glob, n2nz,icoo, jcoo) 

    implicit none 
    integer,intent(in)          :: nkn_glob   ! matrix dimension    
    integer,intent(in)          :: n2nz       ! local non-zeros
    integer,intent(in)          :: icoo(*)    ! local array of global row indices
    integer,intent(in)          :: jcoo(*)    ! local array of global col indices

    return

  end subroutine mumps_init

  subroutine mumps_solve(iteration,acoo,rhs_map,rvec)

    implicit none 
    integer, intent (in)        :: iteration    
    real*8 ,intent(in)          :: acoo(*)    ! matrix entries
    integer,intent(in)          :: rhs_map(*) ! local <-> global mapping for rhs    
    real*8 ,optional :: rvec(*)    ! solution     

    return

  end subroutine mumps_solve

  subroutine mumps_gather_2d_nodes_d(array,reb_array)

    use shympi

    implicit none 
    double precision array(nkn_inner)
    double precision,allocatable,dimension(:),optional :: reb_array

    return

  end subroutine


  subroutine mumps_finalize

    implicit none    

    return

  end subroutine mumps_finalize            

end module
