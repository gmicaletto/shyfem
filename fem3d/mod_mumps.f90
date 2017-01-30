module mod_mumps

#include <dmumps_struc.h>
#include <mpif.h>
   type ( DMUMPS_STRUC) :: mumps_par      
   integer,allocatable,dimension(:),save :: loc_rhs_map,total_rhs_map

contains

  subroutine map_rhs(nkn,r)

    use shympi

    integer  r(nkn)
    integer  i,nkn

    if(bmpi) then        
      do i=1, nkn
         r(i) = mynodes%globalID(i)
      end do
    else           
      do i=1, nkn
         r(i) = i
      end do
    end if
    
    return               

  end subroutine map_rhs          
                
  subroutine mumps_init(nkn_glob, n2nz,icoo, jcoo) 

    use shympi
    use basin
 
    implicit none 
    integer                     :: IERR    
    integer,intent(in)          :: nkn_glob   ! matrix dimension    
    integer,intent(in)          :: n2nz       ! local non-zeros
    integer,intent(in)          :: icoo(*)    ! local array of global row indices
    integer,intent(in)          :: jcoo(*)    ! local array of global col indices

    call MPI_COMM_RANK(MPI_COMM_WORLD,mumps_par%MYID,IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,mumps_par%NPROCS,IERR)

    if ( mumps_par%MYID == 0 ) write(6,*) 'using MUMPS 5.0.1 '
    ! Define a communicator for the package.
    mumps_par%COMM = MPI_COMM_WORLD

    mumps_par%SYM =  0   ! matrix symmetry
    mumps_par%PAR =  1   ! host involved in factorization, analysis and solution
    mumps_par%JOB = -1   ! initialize package
        
    call DMUMPS(mumps_par)
    if (mumps_par%INFOG(1) .lt. 0) then
      write(6,'(A,A,I6,A,I9)') " ERROR RETURN: ", &
      &            "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
      &            "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2)
      ! DECIDE WHAT TO DO IN CASE OF ERROR
    end if

    mumps_par%ICNTL(1)  =  6
    mumps_par%ICNTL(2)  =  6
    mumps_par%ICNTL(3)  =  6
    mumps_par%ICNTL(4)  =  1   ! level of printing ( 1 = only error messages )
    mumps_par%ICNTL(5)  =  0   ! matrix input format (distributed, centalized...)
    mumps_par%ICNTL(7)  =  5   ! Metis Ordering
    mumps_par%ICNTL(11) =  2   ! (0 = no stats,1 = all statistics,2= main stats)
    mumps_par%ICNTL(13) =  0   ! use of scalapack if = 0
    mumps_par%ICNTL(14) =  50  ! % of increase of workspace
    mumps_par%ICNTL(18) =  3   ! matrix distribution strategy
    mumps_par%ICNTL(20) =  0   ! dense rhs
    mumps_par%ICNTL(28) =  0   ! automatic choice sequential or parallel
                               ! computation of the ordering
    mumps_par%ICNTL(29) =  2   ! ParMetis

    ! matrix size    
    mumps_par%N = nkn_glob      

    ! local non-zeros of matrix    
    mumps_par%NZ_loc = n2nz      

    ! allocate and set IRN_loc and JCN_loc  
    allocate( mumps_par%IRN_loc ( mumps_par%NZ_loc ) )
    allocate( mumps_par%JCN_loc ( mumps_par%NZ_loc ) )
    mumps_par%IRN_loc(1:mumps_par%NZ_loc) = icoo(1:mumps_par%NZ_loc)     
    mumps_par%JCN_loc(1:mumps_par%NZ_loc) = jcoo(1:mumps_par%NZ_loc)
 
    ! allocate A_loc
    allocate( mumps_par%A_loc   ( mumps_par%NZ_loc ) )    

    ! allocate RHS
    if ( mumps_par%MYID == 0 ) then 
      allocate( mumps_par%RHS ( mumps_par%N ) )
    end if    

    ! create mapping for RHS    
    if (.not. allocated(loc_rhs_map)) then
      allocate( loc_rhs_map(nkn_inner) )
    end if
    call map_rhs(nkn_inner,loc_rhs_map)
    if ( bmpi) then
      if ( shympi_is_master()) then
        call shympi_gather_2d_nodes_i(loc_rhs_map,total_rhs_map)
      else
        call shympi_gather_2d_nodes_i(loc_rhs_map)
      end if
    else
      allocate( total_rhs_map( nkn_glob ))
      total_rhs_map = loc_rhs_map
    end if

    return

  end subroutine mumps_init

  subroutine mumps_solve(iteration,acoo,rhs_map,rvec)

    use shympi    

    implicit none 
    integer                     :: IERR,i      
    integer, intent (in)        :: iteration    
    real*8 ,intent(in)          :: acoo(*)    ! matrix entries
    integer,intent(in)          :: rhs_map(*) ! local <-> global mapping for rhs    
    real*8,dimension(:) ,optional :: rvec(*)    ! solution     
    real*8                      :: set_time     
    real*8                      :: ana_time    
    real*8                      :: fac_time   
    real*8                      :: sol_time    
    real*8                      :: return_time  
    real*8                      :: tot_time    

    ! set matrix entries    
    mumps_par%A_loc(1:mumps_par%NZ_loc) = acoo(1:mumps_par%NZ_loc) 
    
    ! set RHS    
    if ( mumps_par%MYID == 0 ) then 
      do i=1,mumps_par%N
        mumps_par%RHS(rhs_map(i)) = rvec(i)
      end do  
    end if    

    ! analysis step 
    if ( iteration .lt. 2 ) then      
      mumps_par%JOB = 1    
      call DMUMPS(mumps_par)
      if (mumps_par%INFOG(1) .lt. 0) then
        write(6,'(A,A,I6,A,I9)') " ERROR RETURN: ", &
        &            "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
        &            "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2)
        call shympi_barrier    
        stop
      end if
    end if    

    ! factorization step    
    mumps_par%JOB = 2    
    call DMUMPS(mumps_par)
    if (mumps_par%INFOG(1) .lt. 0) then
      write(6,'(A,A,I6,A,I9)') " ERROR RETURN: ", &
      &            "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
      &            "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2)
        call shympi_barrier    
        stop
    end if
           
    ! solution step    
    mumps_par%JOB = 3    
    call DMUMPS(mumps_par)
    if (mumps_par%INFOG(1) .lt. 0) then
      write(6,'(A,A,I6,A,I9)') " ERROR RETURN: ", &
      &            "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
      &            "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2)
        call shympi_barrier    
        stop
    end if

    ! return solution    
    if ( mumps_par%MYID == 0 ) then
      do i = 1,mumps_par%N
        rvec(i) = mumps_par%RHS(i) 
      end do 
    end if            
    
    if ( mumps_par%MYID == -1 ) then
      write(6,*) 'mumps    set_time: ', set_time  
      write(6,*) 'mumps    ana_time: ', ana_time  
      write(6,*) 'mumps    fac_time: ', fac_time  
      write(6,*) 'mumps    sol_time: ', sol_time  
      write(6,*) 'mumps return_time: ', return_time  
      write(6,*) 'mumps  total_time: ', tot_time  
    end if                    

  end subroutine mumps_solve

        subroutine shympi_gather_2d_nodes_i(array, reb_array)

        use basin
        use mpi_common_struct
        use shympi

        implicit none

        integer sendbuffer,ierr,i
        integer, dimension(n_threads) :: recvbuffer, displs

        integer array(nkn_inner)
        integer,allocatable,optional,dimension(:) :: reb_array

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

        call MPI_GATHERV(array, sendbuffer, MPI_INTEGER, &
                      reb_array, numberNodes, displs, &
                      MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)


        end subroutine shympi_gather_2d_nodes_i

  subroutine mumps_gather_2d_nodes_d(array,reb_array)

     use shympi
     use basin
     use mpi_common_struct

     implicit none

     integer sendbuffer,ierr,i
     integer, dimension(n_threads) :: recvbuffer, displs

     double precision array(nkn_inner)
     double precision,allocatable,dimension(:),optional :: reb_array

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


  subroutine mumps_finalize

    use shympi

    implicit none    
    ! terminate package    
    mumps_par%JOB = -2
    if ( mumps_par%MYID == 0 ) write(6,*) 'finalize MUMPS 5.0.1 '
    call DMUMPS(mumps_par)
    if (mumps_par%INFOG(1) .lt. 0) then
      write(6,'(A,A,I6,A,I9)') " ERROR RETURN: ", &
      &            "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
      &            "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2)
      call shympi_barrier    
      stop
    end if

  end subroutine mumps_finalize            

end module
