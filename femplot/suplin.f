c
c $Id: suplin.f,v 1.5 2010-03-11 15:32:16 georg Exp $
c
c routines for section plot (scalars)
c
c revision log :
c
c 14.09.2009    ggu     routines written from scratch
c 09.10.2009    ggu     vertical plot nearly ready
c 14.10.2009    ggu     vertical plot finished (scalar and velocities)
c 26.03.2010    ggu     vertical plot for velocity finished 
c
c************************************************************************

	subroutine plot_sect(bvel,sv)

c plots section

	implicit none

        include 'param.h'

	logical bvel			!plot velocities
	real sv(nlvdim,nkndim)		!scalar to be plotted

	integer nldim
	parameter (nldim=200)

	real hev(1)
	common /hev/hev
	real xgv(1), ygv(1)
	common /xgv/xgv, /ygv/ygv
	real hlv(1)
	common /hlv/hlv
	integer ilhkv(1)
	common /ilhkv/ilhkv
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

c elems(1) is not used, etc..

	real val(0:2*nlvdim,nldim)	!scalar value along line
	real vel(3,0:2*nlvdim,nldim)	!projected velocities along line
	integer nodes(nldim)		!nodes along line
	integer elems(nldim)		!elements along line
	integer lelems(nldim)		!levels in element
	real helems(nldim)		!depth in elements
	real dxy(2,nldim)		!direction of projection line

	character*80 file
	character*80 string,line
	logical bhoriz
	integer it,i,k1,k2,ie,l,lbot,ltop,j
	integer ltot
	integer nr,nc,mode,ir
	real xmin,ymin,xmax,ymax
	real xrmin,yrmin,xrmax,yrmax
	real d,dd,dx,dy,htop,htot,hbot
	real x0,y0,x1,x2,y1,y2
	real vmin,vmax
	real vhmin,vhmax
	real rrl,rrd,x,y,xtick,ytick,ylast,rdist,rmax
	real xcm,ycm
	real fact
	integer ndec,nctick

	real x0s,y0s,x1s,y1s
	real xmid,hmid,umid,wmid
	real scale

	integer n
	save n,nodes,elems,helems,lelems,dxy

	integer ll,lvmax
	save ll,lvmax
	real rl,rd,hvmax
	save rl,rd,hvmax
	logical bgrid
	integer ivert
	save bgrid,ivert

	character*80 vtitle,xtitle,ytitle,ltitle,rtitle
	save vtitle,xtitle,ytitle,ltitle,rtitle
	real ascale,rscale,stip
	save ascale,rscale,stip

	logical inboxdim_noabs
	integer gettime
	integer ialfa
	real getpar
	real hlog,divdist,roundm

	integer icall
	save icall
	data icall / 0 /

c----------------------------------------------------------------
c initialize
c----------------------------------------------------------------

	if( icall .eq. 0 ) then
	  call getfnm('vsect',file)
	  call line_read_nodes(file,nldim,n,nodes)
	  call line_find_elements(n,nodes,hev,hlv,elems,helems,lelems)
	  call line_find_min_max(n,nodes,helems,lelems,xgv,ygv,rl,rd,ll)
	  call make_proj_dir(n,nodes,xgv,ygv,dxy)

	  bgrid = nint(getpar('ivgrid')) .ne. 0
	  ivert = nint(getpar('ivert'))
	  lvmax = nint(getpar('lvmax'))
	  hvmax = getpar('hvmax')

	  call set_max_dep_lay(nlv,hlv,rd,ll,hvmax,lvmax) !vertical range

	  call getfnm('vtitle',vtitle)
	  call getfnm('xtitle',xtitle)
	  call getfnm('ytitle',ytitle)
	  call getfnm('ltitle',ltitle)
	  call getfnm('rtitle',rtitle)

	  ascale = getpar('avscal')	!absolute scale
	  rscale = getpar('rvscal')	!relative scale
	  stip = getpar('svtip')	!arrow tip size
	end if

	icall = icall + 1

	it = gettime()

        call qstart
        call annotes(it,'vertical plot')
	call annote

c----------------------------------------------------------------
c prepare data
c----------------------------------------------------------------

	write(6,*) 'plotting section: ',it,n,rl,rd

	if( bvel ) then
	  call proj_velox(n,nodes,ilhkv,dxy,vel,val,vmin,vmax,vhmin,vhmax)
	else
	  call line_insert_scalars(n,nodes,ilhkv,nlvdim,sv,val,vmin,vmax)
	end if
	call colauto(vmin,vmax)
	write(6,*) 'min/max on line: ',vmin,vmax

c----------------------------------------------------------------
c set viewport
c----------------------------------------------------------------

	xrmin = 0.
	yrmin = 0.
	xrmax = 1.
	yrmax = 1.

	call qgetvp(xmin,ymin,xmax,ymax)
	write(6,*) 'plot_sect',xmin,ymin,xmax,ymax

	ymax = ymax / 2.			!empirical
	call qsetvp(xmin,ymin,xmax,ymax)
	call qworld(xrmin,yrmin,xrmax,yrmax)
	call pbox(xrmin,yrmin,xrmax,yrmax)	!plot outer box
	call bpixel

        bhoriz = inboxdim_noabs('sect',x0s,y0s,x1s,y1s)

        if( bhoriz ) then
	  xmin = xmin + 1.5
	  ymin = ymin + 1.5
	  xmax = xmax - 1.5
	  ymax = ymax - 1.5
	else
	  xmin = xmin + 1.5
	  ymin = ymin + 1.5
	  xmax = xmax - 2.5
	  ymax = ymax - 1.5
	end if

	call qsetvpnc(xmin,ymin,xmax,ymax)	!this is the plotting area

c-----------------------------------------------------------------
c set world coordinates
c-----------------------------------------------------------------

	xrmin = 0.
	xrmax = rl
	yrmin = -rd
	yrmax = 0.

	if( ivert .eq. 1 ) yrmin = -ll

	call qworld(xrmin,yrmin,xrmax,yrmax)

c--------------------------------------------------------------------
c plot scalar
c--------------------------------------------------------------------

	d = 0.
	do i=2,n
	  k1 = nodes(i-1)
	  k2 = nodes(i)
	  htot = helems(i)
	  ltot = lelems(i)
	  htot = min(htot,hvmax)
	  ltot = min(ltot,lvmax)
	  ie = elems(i)
	  dx = xgv(k2) - xgv(k1)
	  dy = ygv(k2) - ygv(k1)
	  dd = sqrt( dx**2 + dy**2 )
	  x1 = d
	  x2 = d + dd
	  y1 = yrmin
	  y2 = -htot
	  if( ivert .eq. 1 ) y2 = -ltot 
	  if( ivert .eq. 2 ) y2 = -hlog(htot,rd)
	  call qgray(0.5)
	  call qrfill(x1,y1,x2,y2)	!land (bottom)
	  htop = 0.
	  do l=1,ltot
	    hbot = hlv(l)
	    ltop = 2*l - 2
	    if( hbot .gt. htot ) hbot = htot
	    if( ivert .eq. 1 ) then
	      htop = l-1
	      hbot = l
	    else if( ivert .eq. 2 ) then
	      hbot = hlog(hbot,rd)
	    end if
	    call plot_rect(x1,-htop,x2,-hbot,val(ltop,i-1),val(ltop,i))
	    htop = hbot
	  end do
	  d = d + dd
	end do

c--------------------------------------------------------------------
c plot vector
c--------------------------------------------------------------------

	vhmax = max(abs(vhmin),abs(vhmax))
	scale = rscale*rl/(2.*vhmax*(n-1))
	if( ascale .gt. 0. ) scale = ascale * rscale / 2.
	write(6,*) 'arrow scale: ',vhmax,ascale,rscale,scale

	d = 0.
	do i=2,n
	  k1 = nodes(i-1)
	  k2 = nodes(i)
	  htot = helems(i)
	  ltot = lelems(i)
	  htot = min(htot,hvmax)
	  ltot = min(ltot,lvmax)
	  ie = elems(i)
	  dx = xgv(k2) - xgv(k1)
	  dy = ygv(k2) - ygv(k1)
	  dd = sqrt( dx**2 + dy**2 )
	  x1 = d
	  x2 = d + dd
	  y1 = yrmin
	  y2 = -htot
	  if( ivert .eq. 1 ) y2 = -ltot 
	  if( ivert .eq. 2 ) y2 = -hlog(htot,rd)
	  htop = 0.
	  do l=1,ltot
	    hbot = hlv(l)
	    ltop = 2*l - 2
	    if( hbot .gt. htot ) hbot = htot
	    if( ivert .eq. 1 ) then
	      htop = l-1
	      hbot = l
	    else if( ivert .eq. 2 ) then
	      hbot = hlog(hbot,rd)
	    end if
	    if( bvel .and. stip .ge. 0. ) then
	      hmid = 0.5*(htop+hbot)
	      xmid = 0.5*(x1+x2)
	      umid = 0.5*(vel(2,ltop+1,i-1)+vel(2,ltop+1,i))
	      wmid = 0.5*(vel(3,ltop+1,i-1)+vel(3,ltop+1,i))
	      call qgray(0.0)
	      call plot_arrow(xmid,-hmid,umid,wmid,scale,stip)
	    end if
	    if( bgrid ) then
	      call qgray(0.5)
	      call pbox(x1,-htop,x2,-hbot)
	    end if
	    htop = hbot
	  end do
	  d = d + dd
	end do

c--------------------------------------------------------------------
c labeling
c--------------------------------------------------------------------

        call qfont('Times-Roman')
        call qgray(0.0)
	call qtxts(9)
	call qcm(xcm,ycm)

	rrl = divdist(rl,5,0)
	nr = rl/rrl
	nc = -1		!do not write decimal point
	mode = -1	!left flushing
	y = yrmin
	ytick = 0.1 * ycm

	do i=0,nr			!x-axis
	  x = i * rrl
	  ir = ialfa(x,string,nc,mode)
	  call qline(x,y,x,y-ytick)
          call qtxtcr(0.,+2.5)
          call qtext(x,y,string(1:ir))
	end do

	ylast = -100.	!big negative number
	x = xrmin
	xtick = 0.1 * xcm
	rmax = -yrmin
	rrd = 0.
	rdist = max(1,nint(ll/7.))	!for ivert == 1

	do while( rrd .le. rmax )
	  y = rrd
	  if( ivert .eq. 2 ) y = hlog(y,rmax)
	  !write(6,*) rmax,rrd,y
          if( rrd .ne. 0.5 .and. y-ylast .gt. rmax/50. ) then
	    ir = ialfa(rrd,string,nc,mode)
	    call qline(x,-y,x-xtick,-y)
            call qtxtcr(+1.,0.)
            call qtext(x-2*xtick,-y,string(1:ir))
	    ylast = y
	  end if
	  if( ivert .eq. 1 ) then
	    rrd = rrd + rdist
	  else
	    rrd = roundm(rrd+0.5,1)
	  end if
	end do

c--------------------------------------------------------------------
c titles
c--------------------------------------------------------------------

	call qcomm('Titles')
        call qtxtcr(0.,0.)

	call qtxts(15)
        call qtext(rl/2.,-yrmin/15.,vtitle)	!title

	call qtxts(12)
        call qtext(rl/2.,yrmin*1.1,xtitle)	!x-axis

        call qtxtcr(+1.,0.)
        call qtxtcr(0.,0.)
        call qtext(0.,yrmin*1.07,ltitle)	!left point
        call qtxtcr(-1.,0.)
        call qtxtcr(0.,0.)
        call qtext(rl,yrmin*1.07,rtitle)	!right point

        call qtxtr(90.)
        call qtxtcr(0.,0.)
	if( ivert .eq. 1 .and. ytitle .eq. 'Depth [m]' ) then
		ytitle = 'Layers'
	end if
        call qtext(-rl/15.,yrmin*0.5,ytitle)	!y-axis

	call pbox(xrmin,yrmin,xrmax,yrmax)

c--------------------------------------------------------------------
c color bar
c--------------------------------------------------------------------

        if( bhoriz ) then
	  dx = xmax - xmin
	  dy = ymax - ymin
	  xmax = xmin + x1s * dx
	  ymax = ymin + y1s * dy
	  xmin = xmin + x0s * dx
	  ymin = ymin + y0s * dy
	else
	  dx = 0.2
	  dy = 1.
	  xmin = xmax + dx
	  ymin = ymin + dy
	  xmax = xmax + 2.5 - dx
	  ymax = ymax - dy
	end if

	xrmin = 0.
	yrmin = 0.
	xrmax = 1.
	yrmax = 1.

	call qsetvp(xmin,ymin,xmax,ymax)	!this is the plotting area
	call qworld(xrmin,yrmin,xrmax,yrmax)

	call qcomm('Color Bar')
	call pbox(x0,y0,x1,y1)

        fact = getpar('faccol')
        ndec = nint(getpar('ndccol'))
        nctick = nint(getpar('nctick'))
        call getfnm('legcol',line)

	call qgray(0.0)
	call colbar(line,bhoriz,nctick,ndec,fact,xrmin,yrmin,xrmax,yrmax)

c--------------------------------------------------------------------
c end of routine
c--------------------------------------------------------------------

        call qend

	end

c************************************************************************
c************************************************************************
c************************************************************************

	function hlog(x,rd)

	implicit none

	real hlog
	real x,rd

	hlog = rd * log(x+1.)/log(rd+1.)

	end

c************************************************************************

	subroutine plot_rect(x1,y1,x2,y2,v1,v2)

	implicit none

	real x1,y1,x2,y2
	real v1(3),v2(3)

	include 'color.h'

	logical bfirst
	integer icsave
	real xm,ym,vm
	real x(3),y(3),f(3)

	bfirst = .false.	!what is this? -> old way: no middle point

	xm = (x1+x2)/2.
	ym = (y1+y2)/2.
	vm = (v1(2)+v2(2))/2.

	call set_auto_color_table

	if( bfirst ) then

	call setxyf(x1,x1,x2,y1,ym,y1,v1(1),v1(2),v2(1),x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)
	call setxyf(x1,x2,x2,ym,ym,y1,v1(2),v2(2),v2(1),x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)

	call setxyf(x1,x1,x2,ym,y2,ym,v1(2),v1(3),v2(2),x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)
	call setxyf(x1,x2,x2,y2,y2,ym,v1(3),v2(3),v2(2),x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)

	else

	!first plot upper triangles, than lower ones

	call setxyf(x1,x1,xm,y1,ym,ym,v1(1),v1(2),vm,x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)
	call setxyf(x1,xm,x2,y1,ym,y1,v1(1),vm,v2(1),x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)
	call setxyf(xm,x2,x2,ym,ym,y1,vm,v2(2),v2(1),x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)

	call setxyf(x1,x1,xm,ym,y2,ym,v1(2),v1(3),vm,x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)
	call setxyf(x1,xm,x2,y2,ym,y2,v1(3),vm,v2(3),x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)
	call setxyf(xm,x2,x2,ym,y2,ym,vm,v2(3),v2(2),x,y,f)
	call plcol(x,y,f,ciso,fiso,isoanz+1,fnull)

	end if

	call reset_auto_color_table

	end

c************************************************************************

	subroutine setxyf(x1,x2,x3,y1,y2,y3,f1,f2,f3,x,y,f)
	real x(3),y(3),f(3)
	x(1)=x1
	x(2)=x2
	x(3)=x3
	y(1)=y1
	y(2)=y2
	y(3)=y3
	f(1)=f1
	f(2)=f2
	f(3)=f3
	end

c************************************************************************

	subroutine plot_arrow(x,y,u,w,scale,stip)

	implicit none

	real x,y,u,w,scale,stip

	real uu,ww,s

	uu = x + u*scale
	ww = y + w*scale
	s = stip

	call fcmpfeil(x,y,uu,ww,s)

	end

c************************************************************************

	subroutine make_proj_dir(n,nodes,xgv,ygv,dxy)

c computes projection line for nodes

	implicit none

	integer n
	integer nodes(n)
	real xgv(1), ygv(1)
	real dxy(2,n)			!direction of line (for projection)

	integer i,i1,i2,k1,k2
	real dx,dy,dd

	do i=1,n
	  i1 = max(1,i-1)
	  i2 = min(n,i+1)
	  k1 = nodes(i1)
	  k2 = nodes(i2)
	  dx = xgv(k2) - xgv(k1)
	  dy = ygv(k2) - ygv(k1)
	  dd = sqrt( dx*dx + dy*dy )
	  dxy(1,i) = dx/dd
	  dxy(2,i) = dy/dd
	end do

	end

c************************************************************************

	subroutine proj_velox(n,nodes,ilhkv,dxy,vel,val
     +				,vmin,vmax,vhmin,vhmax)

	implicit none

	include 'param.h'

	integer n
	integer nodes(n)
	integer ilhkv(1)		!number of layers in node
	real dxy(2,n)			!direction of line (for projection)
	real vel(3,0:2*nlvdim,n)	!projected velocities along line
	real val(0:2*nlvdim,n)
	real vmin,vmax
	real vhmin,vhmax

        real uprv(nlvdim,nkndim)
        common /uprv/uprv
        real vprv(nlvdim,nkndim)
        common /vprv/vprv
        real wprv(nlvdim,nkndim)
        common /wprv/wprv

	integer i,k,l,lmax,llayer
	real ut,un,w

	do i=1,n
	  k = nodes(i)
	  lmax = ilhkv(k)
	  do l=1,lmax
	    llayer = 2*l-1
	    ut =  dxy(1,i)*uprv(l,k) + dxy(2,i)*vprv(l,k)	!tangent
	    un = -dxy(2,i)*uprv(l,k) + dxy(1,i)*vprv(l,k)	!normal
	    w  =  wprv(l,k)
	    val(llayer,i) = un
	    vel(1,llayer,i) = un
	    vel(2,llayer,i) = ut
	    vel(3,llayer,i) = w
	  end do
	  call insert_between_layers(3,2*lmax,vel(1,0,i))
	  call insert_between_layers(1,2*lmax,val(0,i))
	end do

	i=5
	k = nodes(i)
	lmax = ilhkv(k)
	write(66,*) (vel(1,l,i),l=0,2*lmax)
	write(66,*) (vel(2,l,i),l=0,2*lmax)
	write(66,*) (vel(3,l,i),l=0,2*lmax)
	write(66,*) '---- proj_velox ---------'

	vmin = val(0,1)
	vmax = val(0,1)
	vhmin = vel(2,0,1)
	vhmax = vel(2,0,1)
	do i=1,n
	  k = nodes(i)
	  lmax = ilhkv(k)
	  do l=0,2*lmax
	    vmin = min(vmin,val(l,i))
	    vmax = max(vmax,val(l,i))
	    vhmin = min(vhmin,vel(2,l,i))
	    vhmax = max(vhmax,vel(2,l,i))
	  end do
	end do

	end

c************************************************************************

	subroutine insert_between_layers(nd,nl,vals)

	implicit none

	integer nd,nl
	real vals(nd,0:2*nl)

	integer i,l

	do i=1,nd
	  vals(i,0) = vals(i,1)
	  do l=2,2*(nl-1),2
	    vals(i,l) = 0.5 * (vals(i,l-1)+vals(i,l+1))
	  end do
	  vals(i,2*nl) = vals(i,2*nl-1)
	end do

	end

c************************************************************************

	subroutine line_insert_scalars(n,nodes,ilhkv,nlvdim,sv
     +					,values,vmin,vmax)

c inserts scalar values into matrix section

	implicit none

	integer n
	integer nodes(n)
	integer ilhkv(1)		!number of layers in node
	integer nlvdim
	real sv(nlvdim,1)
	real values(0:2*nlvdim,n)
	real vmin,vmax

	integer i,k,lmax,l,llayer
	integer ipext

	vmin = 1.e+30
	vmax = -1.e+30

	do i=1,n
	  k = nodes(i)
	  lmax = ilhkv(k)
	  values(0,i) = sv(1,k)
	  values(1,i) = sv(1,k)
	  vmin = min(vmin,sv(1,k))
	  vmax = max(vmax,sv(1,k))
	  do l=2,lmax
	    llayer = 2*l-1
	    values(llayer,i) = sv(l,k)
	    values(llayer-1,i) = 0.5 * ( sv(l,k) + sv(l-1,k) )
	    vmin = min(vmin,sv(l,k))
	    vmax = max(vmax,sv(l,k))
	  end do
	  values(2*lmax,i) = sv(lmax,k)
	  !write(78,*) k,ipext(k),lmax
	  !write(78,*) (sv(l,k),l=1,lmax)
	  !write(79,*) k,ipext(k),lmax,2*lmax
	  !write(79,*) (values(l,i),l=0,2*lmax)
	end do

	end

c************************************************************************

	subroutine line_find_min_max(n,nodes,helems,lelems
     +					,xgv,ygv,rl,rd,ll)

c finds length and max depth of line

	implicit none

	integer n
	integer nodes(n)
	real helems(n)		!depth in chosen elements
	integer lelems(n)
	real xgv(1), ygv(1)
	real rl,rd		!length and depth (return)
	integer ll

	integer i,k1,k2
	real dx,dy

	rl = 0.
	rd = 0.
	ll = 0

	do i=2,n
	  k1 = nodes(i-1)
	  k2 = nodes(i)
	  dx = xgv(k2) - xgv(k1)
	  dy = ygv(k2) - ygv(k1)
	  rl = rl + sqrt( dx**2 + dy**2 )
	  rd = max(rd,helems(i))
	  ll = max(ll,lelems(i))
	end do

	end

c************************************************************************

	subroutine line_find_elements(n,nodes,hev,hlv
     +					,elems,helems,lelems)

c finds elements along line given by nodes
c
c deepest element is chosen

	implicit none

	integer n
	integer nodes(n)
	real hev(1)		!depth in elements
	real hlv(1)		!layer structure
	integer elems(n)	!element number of chosen elements
	real helems(n)		!depth in chosen elements
	integer lelems(n)	!layers in element number of 

	integer i,k1,k2,ie1,ie2,l
	real h
	integer ipext,ieext

	do i=2,n
	  k1 = nodes(i-1)
	  k2 = nodes(i)
	  call find_elems_to_segment(k1,k2,ie1,ie2)
	  elems(i) = ie1
	  if( ie2 .gt. 0 ) then
	    if( hev(ie2) .gt. hev(ie1) ) elems(i) = ie2
	  end if
	  helems(i) = hev(elems(i))
	  !write(6,*) 'line_find_elements: ',k1,k2,ie1,ie2,hev(ie1)
	  !write(6,*) '                    ',ipext(k1),ipext(k2)
	  !write(6,*) '                    ',ieext(ie1),ieext(ie2)
	end do

	do i=2,n
	  h = helems(i)
	  l = 1
	  do while( hlv(l) .lt. h ) 
	    l = l + 1
	  end do
	  lelems(i) = l
	end do
	  
	end

c************************************************************************

	subroutine line_read_nodes(file,ndim,n,nodes)

c reads line given by nodes

	implicit none

	character*(*) file
	integer ndim
	integer n
	integer nodes(ndim)

	integer iunit,nn
	integer ifileo,ipint

	n = 0

	iunit = ifileo(0,file,'form','old')
	if( iunit .le. 0 ) goto 99 

    1	continue
	  read(iunit,*,end=2) nn
	  if( nn .le. 0 ) goto 2
	  n = n + 1
	  if( n .gt. ndim ) stop 'error stop line_read_nodes: ndim'
	  nodes(n) = ipint(nn)
	goto 1
    2	continue

	write(6,*) 'finished reading line section: ',n

	close(iunit)

	return
   99	continue
	write(6,*) 'Cannot open section file: '
	write(6,*) file
	stop 'error stop line_read_nodes: no such file'
	end

c************************************************************************

	subroutine set_max_dep_lay(nlv,hlv,rd,ll,hvmax,lvmax)

c sets hvmax and lvmax

	implicit none

	integer nlv
	real hlv(1)
	real rd
	integer ll
	real hvmax
	integer lvmax

	integer l

        if( hvmax .gt. 0. .and. lvmax .gt. 0 ) then
          write(6,*) 'hvmax,lvmax: ',hvmax,lvmax
          write(6,*) 'only one of hvmax - lvmax can be given'
          stop 'error stop set_max_dep_lay: hvmax/lvmax'
        end if

        if( hvmax .gt. 0. ) then
          do l=1,nlv
            if( hlv(l) .gt. hvmax ) goto 7
          end do
    7     continue
          lvmax = l             !could be one higher than max layer
          ll = lvmax
          rd = hvmax
        else if( lvmax .gt. 0 ) then
          if( lvmax .gt. nlv .and. nlv .gt. 1 ) then
            hvmax = hlv(nlv) + hlv(nlv) - hlv(nlv-1)
          else
            hvmax = hlv(lvmax)
          end if
          ll = lvmax
          rd = hvmax
        else
          hvmax = rd
          lvmax = ll
        end if

	end

c************************************************************************

