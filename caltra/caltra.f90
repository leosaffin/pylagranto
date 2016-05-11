module caltra
  use inter, only: get_index4, int_index4

  implicit none

  contains
  !-----------------------------------------------------------------------------
  subroutine main(xx0, yy0, pp0, leftflag,                                      &
                  ts, nsubs, imethod, numit, jflag, wfactor, fbflag,            &
                  spt0, spt1, p3t0, p3t1,                                       &
                  uut0, uut1, vvt0, vvt1, wwt0,wwt1,                            &
                  xmin, ymin, dx, dy, per, hem,                                 &
                  nx, ny, nz, ntra)

  !f2py real, intent(inout) :: xx0(ntra), yy0(ntra), pp0(ntra)
  !f2py integer, intent(inout) :: leftflag(ntra)
  !f2py real, intent(in) :: ts
  !f2py integer, intent(in) :: nsubs
  !f2py integer, intent(in) :: imethod
  !f2py integer, intent(in) :: numit
  !f2py integer, intent(in) :: jflag
  !f2py real, intent(in) :: wfactor
  !f2py integer, intent(in) :: fbflag
  !f2py real, intent(in) :: spt0(nx*ny), spt1(nx*ny)
  !f2py real, intent(in) :: p3t0(nx*ny*nz), p3t1(nx*ny*nz)
  !f2py real, intent(in) :: uut0(nx*ny*nz), uut1(nx*ny*nz)
  !f2py real, intent(in) :: vvt0(nx*ny*nz), vvt1(nx*ny*nz)
  !f2py real, intent(in) :: wwt0(nx*ny*nz), wwt1(nx*ny*nz)
  !f2py real, intent(in) :: xmin, ymin
  !f2py real, intent(in) :: dx ,dy
  !f2py real, intent(in) :: per
  !f2py integer, intent(in) :: hem
  !f2py integer, intent(hide) :: nx, ny, nz, ntra

  !-Input variables-------------------------------------------------------------
  ! Positions along trajectories
  real, intent(inout) :: xx0(ntra), yy0(ntra), pp0(ntra)

  ! Flag if trajectory leaves domain
  integer, intent(inout) :: leftflag(ntra)

  ! Computational time step
  real, intent(in) :: ts

  ! Number of sub steps between files
  integer, intent(in) :: nsubs

  ! Integration method (1=Euler, 2=Runge-Kutta)
  integer, intent(in) :: imethod

  ! Number of iterations of the Euler method
  integer, intent(in) :: numit

  ! Flag to whether trajectories jump on hitting the ground
  integer, intent(in) :: jflag

  ! Vertical velocity scaling factor
  real, intent(in) :: wfactor

  ! Forward or reverse trajectories (+1 or -1)
  integer, intent(in) :: fbflag

  ! Surface pressure at each timestep
  real, intent(in) :: spt0(nx*ny), spt1(nx*ny)

  ! Pressure
  real, intent(in) :: p3t0(nx*ny*nz), p3t1(nx*ny*nz)

  ! Wind fields
  real, intent(in) :: uut0(nx*ny*nz), uut1(nx*ny*nz)
  real, intent(in) :: vvt0(nx*ny*nz), vvt1(nx*ny*nz)
  real, intent(in) :: wwt0(nx*ny*nz), wwt1(nx*ny*nz)

  ! Domain boundaries
  real, intent(in) :: xmin, ymin

  ! Grid spacing (in degrees)
  real, intent(in) :: dx ,dy

  ! Periodic data flag
  real, intent(in) :: per

  ! Hemispheric data flag
  integer, intent(in) :: hem

  ! Array dimensions
  integer, intent(in) :: nx, ny, nz, ntra

  !-Local variables-------------------------------------------------------------
  ! Missing data value
  real, parameter :: mdv = -1000

  ! Loop counters
  integer :: iloop, i

  ! Relative time positioning in loop
  real :: reltpos0, reltpos1

  ! Trajectory position update
  real :: xx1, yy1, pp1

  !-----------------------------------------------------------------------------
  !Split the interval <timeinc> into computational time steps <ts>
  do iloop=1, nsubs

    ! Calculate relative time position in the interval timeinc
    ! 0=beginning, 1=end
    reltpos0 = (real(iloop)-1.) * ts
    reltpos1 =  real(iloop) * ts

    ! Timestep for all trajectories
    do i=1,ntra
      ! Check if trajectory has already left the data domain
      if (leftflag(i).ne.1) then
        ! Iterative Euler timestep (x0,y0,p0 -> x1,y1,p1)
        if (imethod.eq.1) then
          call euler(xx1, yy1, pp1, leftflag(i), xx0(i), yy0(i), pp0(i),        &
                     reltpos0, reltpos1, ts, numit, jflag, mdv, wfactor,        &
                     fbflag, spt0, spt1, p3t0, p3t1, uut0, uut1, vvt0, vvt1,    &
                     wwt0,wwt1, xmin, ymin, dx, dy, per, hem, nx, ny, nz)

        ! Runge-Kutta timestep (x0,y0,p0 -> x1,y1,p1)
        else if (imethod.eq.2) then
          call runge(xx1, yy1, pp1, leftflag(i), xx0(i), yy0(i), pp0(i),        &
                     reltpos0, reltpos1, ts, numit, jflag, mdv, wfactor,        &
                     fbflag, spt0, spt1, p3t0, p3t1, uut0, uut1, vvt0, vvt1,    &
                     wwt0, wwt1, xmin, ymin, dx, dy, per, hem, nx, ny, nz)
        endif

        ! Update trajectory position, or increase number of trajectories leaving
        ! domain
        if (leftflag(i).ne.1) then
          xx0(i)=xx1
          yy0(i)=yy1
          pp0(i)=pp1
        endif

      ! Trajectory has already left data domain (mark as <mdv>)
      else
        xx0(i)=mdv
        yy0(i)=mdv
        pp0(i)=mdv
      endif
    end do
  end do

  end subroutine main


!*******************************************************************
!* Time step : either Euler or Runge-Kutta                         *
!*******************************************************************

!Time-step from (x0,y0,p0) to (x1,y1,p1)

!(x0,y0,p0) input	coordinates (long,lat,p) for starting point
!(x1,y1,p1) output	coordinates (long,lat,p) for end point
!deltat	 input	timestep in seconds
!numit	 input	number of iterations
!jump	 input  flag (=1 trajectories don't enter the ground)
!left	 output	flag (=1 if trajectory leaves data domain)

!-------------------------------------------------------------------
!Iterative Euler time step
!-------------------------------------------------------------------

      subroutine euler(x1,y1,p1,left,x0,y0,p0,reltpos0,reltpos1,                &
                       deltat,numit,jump,mdv,wfactor,fbflag,                    &
                       spt0,spt1,p3d0,p3d1,uut0,uut1,vvt0,vvt1,wwt0,wwt1,       &
                        xmin,ymin,dx,dy,per,hem,nx,ny,nz)

      implicit none

!Declaration of subroutine parameters
      integer :: nx,ny,nz
      real :: x1,y1,p1
      integer :: left
      real :: x0,y0,p0
      real :: reltpos0,reltpos1
      real :: deltat
      integer :: numit
      integer :: jump
      real :: wfactor
      integer :: fbflag
      real :: spt0(nx*ny)   ,spt1(nx*ny)
      real :: uut0(nx*ny*nz),uut1(nx*ny*nz)
      real :: vvt0(nx*ny*nz),vvt1(nx*ny*nz)
      real :: wwt0(nx*ny*nz),wwt1(nx*ny*nz)
      real :: p3d0(nx*ny*nz),p3d1(nx*ny*nz)
      real :: xmin,ymin,dx,dy
      real :: per
      integer :: hem
      real :: mdv

!Numerical and physical constants
      real :: deltay
      parameter(deltay=1.112E5)  ! Distance in m between 2 lat circles
      real :: pi
      parameter(pi=3.1415927)    ! Pi

!Auxiliary variables
      real :: xmax,ymax
      real :: xind,yind,pind
      real :: u0,v0,w0,u1,v1,w1,u,v,w,sp
      integer :: icount
      character :: ch

!Reset the flag for domain-leaving
      left=0

!Set the east-north boundary of the domain
      xmax = xmin+real(nx-1)*dx
      ymax = ymin+real(ny-1)*dy

!Interpolate wind fields to starting position (x0,y0,p0)
      call get_index4 (xind,yind,pind,x0,y0,p0,reltpos0,&
                      p3d0,p3d1,spt0,spt1,3,&
                      nx,ny,nz,xmin,ymin,dx,dy,mdv)
      u0 = int_index4(uut0,uut1,nx,ny,nz,xind,yind,pind,reltpos0,mdv)
      v0 = int_index4(vvt0,vvt1,nx,ny,nz,xind,yind,pind,reltpos0,mdv)
      w0 = int_index4(wwt0,wwt1,nx,ny,nz,xind,yind,pind,reltpos0,mdv)

!Force the near-surface wind to zero
      if (pind.lt.1.) w0=w0*pind

!For first iteration take ending position equal to starting position
      x1=x0
      y1=y0
      p1=p0

!Iterative calculation of new position
      do icount=1,numit

!   Calculate new winds for advection
         call get_index4 (xind,yind,pind,x1,y1,p1,reltpos1,&
                         p3d0,p3d1,spt0,spt1,3,&
                         nx,ny,nz,xmin,ymin,dx,dy,mdv)
         u1 = int_index4(uut0,uut1,nx,ny,nz,xind,yind,pind,reltpos1,mdv)
         v1 = int_index4(vvt0,vvt1,nx,ny,nz,xind,yind,pind,reltpos1,mdv)
         w1 = int_index4(wwt0,wwt1,nx,ny,nz,xind,yind,pind,reltpos1,mdv)

!   Force the near-surface wind to zero
         if (pind.lt.1.) w1=w1*pind
 
!   Get the new velocity in between
         u=(u0+u1)/2.
         v=(v0+v1)/2.
         w=(w0+w1)/2.
         
!   Calculate new positions
         x1 = x0 + fbflag*u*deltat/(deltay*cos(y0*pi/180.))
         y1 = y0 + fbflag*v*deltat/deltay
         p1 = p0 + fbflag*wfactor*w*deltat/100.

!  Handle pole problems (crossing and near pole trajectory)
        if ((hem.eq.1).and.(y1.gt.90.)) then
          y1=180.-y1
          x1=x1+per/2.
        endif
        if ((hem.eq.1).and.(y1.lt.-90.)) then
          y1=-180.-y1
          x1=x1+per/2.
        endif
        if (y1.gt.89.99) then
           y1=89.99
        endif

!  Handle crossings of the dateline
        if ((hem.eq.1).and.(x1.gt.xmin+per-dx)) then
           x1=xmin+amod(x1-xmin,per)
        endif
        if ((hem.eq.1).and.(x1.lt.xmin)) then
           x1=xmin+per+amod(x1-xmin,per)
        endif

!  Interpolate surface pressure to actual position
        call get_index4 (xind,yind,pind,x1,y1,1050.,reltpos1,&
                        p3d0,p3d1,spt0,spt1,3,&
                       nx,ny,nz,xmin,ymin,dx,dy,mdv)
        sp = int_index4 (spt0,spt1,nx,ny,1,xind,yind,1.,reltpos1,mdv)

!  Handle trajectories which cross the lower boundary (jump flag)
        if ((jump.eq.1).and.(p1.gt.sp)) p1=sp-10.
 
!  Check if trajectory leaves data domain
        if ( ( (hem.eq.0).and.(x1.lt.xmin)    ).or.&
            ( (hem.eq.0).and.(x1.gt.xmax-dx) ).or.&
              (y1.lt.ymin).or.(y1.gt.ymax).or.(p1.gt.sp) )&
       then
          left=1
          goto 100
        endif

      enddo

!Exit point for subroutine
 100  continue

      return

      end

!-------------------------------------------------------------------
!Runge-Kutta (4th order) time-step
!-------------------------------------------------------------------

      subroutine runge(x1,y1,p1,left,x0,y0,p0,reltpos0,reltpos1,&
                      deltat,numit,jump,mdv,wfactor,fbflag,&
                      spt0,spt1,p3d0,p3d1,uut0,uut1,vvt0,vvt1,wwt0,wwt1,&
                      xmin,ymin,dx,dy,per,hem,nx,ny,nz)

      implicit none

!Declaration of subroutine parameters
      integer :: nx,ny,nz
      real :: x1,y1,p1
      integer :: left
      real :: x0,y0,p0
      real :: reltpos0,reltpos1
      real :: deltat
      integer :: numit
      integer :: jump
      real :: wfactor
      integer :: fbflag
      real :: spt0(nx*ny)   ,spt1(nx*ny)
      real :: uut0(nx*ny*nz),uut1(nx*ny*nz)
      real :: vvt0(nx*ny*nz),vvt1(nx*ny*nz)
      real :: wwt0(nx*ny*nz),wwt1(nx*ny*nz)
      real :: p3d0(nx*ny*nz),p3d1(nx*ny*nz)
      real :: xmin,ymin,dx,dy
      real :: per
      integer :: hem
      real :: mdv

!Numerical and physical constants
      real :: deltay
      parameter(deltay=1.112E5)  ! Distance in m between 2 lat circles
      real :: pi
      parameter(pi=3.1415927)    ! Pi

!Auxiliary variables
      real :: xmax,ymax
      real :: xind,yind,pind
      real :: u0,v0,w0,u1,v1,w1,u,v,w,sp
      integer :: icount,n
      real :: xs,ys,ps,xk(4),yk(4),pk(4)
      real :: reltpos

!Reset the flag for domain-leaving
      left=0

!Set the esat-north bounray of the domain
      xmax = xmin+real(nx-1)*dx
      ymax = ymin+real(ny-1)*dy

!Apply the Runge Kutta scheme
      do n=1,4
 
!  Get intermediate position and relative time
        if (n.eq.1) then
          xs=0.
          ys=0.
          ps=0.
          reltpos=reltpos0
        else if (n.eq.4) then
          xs=xk(3)
          ys=yk(3)
          ps=pk(3)
          reltpos=reltpos1
        else
          xs=xk(n-1)/2.
          ys=yk(n-1)/2.
          ps=pk(n-1)/2.
          reltpos=(reltpos0+reltpos1)/2.
        endif
        
!  Calculate new winds for advection
        call get_index4 (xind,yind,pind,x0+xs,y0+ys,p0+ps,reltpos,&
                        p3d0,p3d1,spt0,spt1,3,&
                        nx,ny,nz,xmin,ymin,dx,dy,mdv)
        u = int_index4 (uut0,uut1,nx,ny,nz,xind,yind,pind,reltpos,mdv)
        v = int_index4 (vvt0,vvt1,nx,ny,nz,xind,yind,pind,reltpos,mdv)
        w = int_index4 (wwt0,wwt1,nx,ny,nz,xind,yind,pind,reltpos,mdv)
         
!  Force the near-surface wind to zero
        if (pind.lt.1.) w1=w1*pind
 
!  Update position and keep them
        xk(n)=fbflag*u*deltat/(deltay*cos(y0*pi/180.))
        yk(n)=fbflag*v*deltat/deltay
        pk(n)=fbflag*w*deltat*wfactor/100.

      enddo
 
!Calculate new positions
      x1=x0+(1./6.)*(xk(1)+2.*xk(2)+2.*xk(3)+xk(4))
      y1=y0+(1./6.)*(yk(1)+2.*yk(2)+2.*yk(3)+yk(4))
      p1=p0+(1./6.)*(pk(1)+2.*pk(2)+2.*pk(3)+pk(4))

!Handle pole problems (crossing and near pole trajectory)
      if ((hem.eq.1).and.(y1.gt.90.)) then
         y1=180.-y1
         x1=x1+per/2.
      endif
      if ((hem.eq.1).and.(y1.lt.-90.)) then
         y1=-180.-y1
         x1=x1+per/2.
      endif
      if (y1.gt.89.99) then
         y1=89.99
      endif
      
!Handle crossings of the dateline
      if ((hem.eq.1).and.(x1.gt.xmin+per-dx)) then
         x1=xmin+amod(x1-xmin,per)
      endif
      if ((hem.eq.1).and.(x1.lt.xmin)) then
         x1=xmin+per+amod(x1-xmin,per)
      endif
      
!Interpolate surface pressure to actual position
      call get_index4 (xind,yind,pind,x1,y1,1050.,reltpos1,&
                      p3d0,p3d1,spt0,spt1,3,&
                      nx,ny,nz,xmin,ymin,dx,dy,mdv)
      sp = int_index4 (spt0,spt1,nx,ny,1,xind,yind,1.,reltpos,mdv)

!Handle trajectories which cross the lower boundary (jump flag)
      if ((jump.eq.1).and.(p1.gt.sp)) p1=sp-10.
      
!Check if trajectory leaves data domain
      if ( ( (hem.eq.0).and.(x1.lt.xmin)    ).or.&
          ( (hem.eq.0).and.(x1.gt.xmax-dx) ).or.&
            (y1.lt.ymin).or.(y1.gt.ymax).or.(p1.gt.sp) )&
     then
         left=1
         goto 100
      endif
      
!Exit point fdor subroutine
 100  continue

      return
      end

end module caltra
