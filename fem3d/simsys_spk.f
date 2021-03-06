!
! revision log :
!
! 12.01.2009	ggu	new file for system routines
! 31.03.2009	ggu	call renamed to spk_*
! 25.05.2015	ggu	some calls changed (pass array in)
! 09.12.2015	ggu	adapted to new pointers and 3d matrix
! 15.12.2015	ggu&deb	finsihed and validated
!
!******************************************************************

        subroutine system_initialize

! allocates data structure

	use mod_system
	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        write(6,*) '----------------------------------------'
        write(6,*) 'initializing matrix inversion routines'
        write(6,*) 'using Sparskit routines'
        write(6,*) '----------------------------------------'

        call mod_system_init(nkn,nel,ngr,mbw,nlv)

        end

!******************************************************************

	subroutine system_init

! on first call initializes pointers - sets arrays to zero
!
! must be called before every assembly of matrix

	use mod_system

	implicit none

	call spk_init_system

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine system_solve_z(n,z)

! solves system - z is used for initial guess

	use mod_system

	implicit none

	integer n
	real z(n)

	call spk_solve_system(.false.,n2max,n,z)
	!call spk_solve_system(n,z)

	end

!******************************************************************

	subroutine system_solve_3d(n,nlvdi,nlv,z)

! solves system - z is used for initial guess

	use mod_system

	implicit none

	integer n,nlvdi,nlv
	real z(nlvdi,n)

	integer i,k,l
	real p((nlv+2)*n)
	integer nn !DEB

	i = 0

	do k=1,n
	  i = i + 1
	  p(i) = 0.
	  do l=1,nlv
	    i = i + 1
	    p(i) = z(l,k)
	  end do
	  i = i + 1
	  p(i) = 0.
	end do

	nn = n*(nlv + 2) !DEB
	!call spk_solve_system(.true.,n3max,n,p)
	call spk_solve_system(.true.,n3max,nn,p) !DEB

	end

!******************************************************************

	subroutine system_assemble(ie,kn,mass,rhs)

! assembles element matrix into system matrix

	use mod_system

	implicit none

	integer ie
	integer kn(3)
	real mass(3,3)
	real rhs(3)

	integer i,j,kk

	integer loclp,loccoo
	external loclp,loccoo

        do i=1,3
          do j=1,3
            !kk=loccoo(kn(i),kn(j),n,m)		!COOGGU	
            !if(kk.gt.0) coo(kk) = coo(kk) + mass(i,j)
            kk=ijp_ie(i,j,ie)			!COOGGU
            if(kk.gt.0) c2coo(kk) = c2coo(kk) + mass(i,j)
          end do
          rvec2d(kn(i)) = rvec2d(kn(i)) + rhs(i)
        end do

	end

!******************************************************************

	subroutine system_assemble_3d(ie,l,nlv,kn,mass,rhs)

! assembles element matrix into system matrix

	use mod_system

	implicit none

	integer ie,l,nlv
	integer kn(3)
	real mass(-1:1,3,3)
	real rhs(3)

	integer i,j,kk

	integer loccoo3d
	external loccoo3d

        do i=1,3
          do j=1,3
	    kk = loccoo3d(i,j,kn,l,ie)
            if(kk.gt.0) then
	       c3coo(kk-1) = c3coo(kk-1) + mass(-1,i,j)
	       c3coo(kk) = c3coo(kk) + mass(0,i,j)
	       c3coo(kk+1) = c3coo(kk+1) + mass(+1,i,j)
	    end if
          end do
	  kk = (nlv+2)*(kn(i)-1) + l + 1
          rvec3d(kk) = rvec3d(kk) + rhs(i) !DEB
        end do
	
	end

!******************************************************************

        subroutine system_adjust_z(n,z)

! copies solution back to z

	use mod_system

        implicit none

	integer n
	real z(n)

        integer k

	z = real(rvec2d)

        end

!******************************************************************

        subroutine system_adjust_3d(n,nlvdi,nlv,z)

! copies solution back to z

	use mod_system

        implicit none

	integer n,nlvdi,nlv
	real z(nlvdi,n)

	integer i,k,l

	i = 0

	do k=1,n
	  i = i + 1
	  do l=1,nlv
	    i = i + 1
	    z(l,k) = real(rvec3d(i)) !DEB
	  end do
	  i = i + 1
	end do

        end

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine system_add_rhs(dt,n,array)

! adds right hand side to system array

	use mod_system

        implicit none

        real dt
	integer n
        real array(n)

        integer k

        do k=1,n
          rvec2d(k) = rvec2d(k) + dt * array(k)
        end do

        end

!******************************************************************

        subroutine system_adjust_matrix_3d

        implicit none

        call coo_adjust

        end

!******************************************************************

