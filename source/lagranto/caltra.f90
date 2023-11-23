module caltra
  use inter, only: get_index4, int_index4

  implicit none

  !Numerical and physical constants
  real, parameter :: deltay=1.112E5  ! Distance in m between 2 lat circles
  real, parameter :: pi=3.1415927
  real, parameter :: mdv = -1000     ! Missing data value

  contains
  !-----------------------------------------------------------------------------
  subroutine main(xx0, yy0, pp0, leftflag,                                      &
                  ts, nsubs, imethod, numit, jflag, fbflag,                     &
                  spt0, spt1, p3t0, p3t1,                                       &
                  uut0, uut1, vvt0, vvt1, wwt0,wwt1,                            &
                  xmin, ymin, dx, dy, per, hem,                                 &
                  nx, ny, nz, ntra,                                             &
                  x0, y0, p0, left)

  !f2py real, intent(in) :: xx0(ntra), yy0(ntra), pp0(ntra)
  !f2py integer, intent(in) :: leftflag(ntra)
  !f2py real, intent(in) :: ts
  !f2py integer, intent(in) :: nsubs
  !f2py integer, intent(in) :: imethod
  !f2py integer, intent(in) :: numit
  !f2py integer, intent(in) :: jflag
  !f2py integer, intent(in) :: fbflag
  !f2py real, intent(in) :: spt0(nx,ny), spt1(nx,ny)
  !f2py real, intent(in) :: p3t0(nx,ny,nz), p3t1(nx,ny,nz)
  !f2py real, intent(in) :: uut0(nx,ny,nz), uut1(nx,ny,nz)
  !f2py real, intent(in) :: vvt0(nx,ny,nz), vvt1(nx,ny,nz)
  !f2py real, intent(in) :: wwt0(nx,ny,nz), wwt1(nx,ny,nz)
  !f2py real, intent(in) :: xmin, ymin
  !f2py real, intent(in) :: dx ,dy
  !f2py real, intent(in) :: per
  !f2py integer, intent(in) :: hem
  !f2py integer, intent(hide) :: nx, ny, nz, ntra
  !f2py real, intent(out) :: x0(ntra), y0(ntra), p0(ntra)
  !f2py integer, intent(out) :: left(ntra)

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

  !-Output variables------------------------------------------------------------
  ! Positions along trajectories
  real, intent(out) :: x0(ntra), y0(ntra), p0(ntra)

  ! Flag if trajectory leaves domain
  integer, intent(out) :: left(ntra)

  !-Local variables-------------------------------------------------------------
  ! Loop counters
  integer :: iloop, i

  ! Input wind relative time positioning in loop
  real :: reltpos0, reltpos1

  ! Trajectory position update
  real :: xx1, yy1, pp1

  !-----------------------------------------------------------------------------
  ! Copy input to output
  x0 = xx0
  y0 = yy0
  p0 = pp0
  left = leftflag
  !Split the interval between files (ts) into (nsubs) computational time steps
  do iloop=1, nsubs

    ! Calculate relative time position in the interval timeinc
    ! 0=beginning, 1=end
    reltpos0 =  real(iloop - 1) / real(nsubs)
    reltpos1 =  real(iloop)     / real(nsubs)

    ! Timestep for all trajectories
    do i=1,ntra
      ! Check if trajectory has already left the data domain
      if (left(i).ne.1) then
        ! Iterative Euler timestep (x0,y0,p0 -> x1,y1,p1)
        if (imethod.eq.1) then
          call euler(xx1, yy1, pp1, left(i), x0(i), y0(i), p0(i),               &
                     reltpos0, reltpos1, ts, numit, jflag,                      &
                     fbflag, spt0, spt1, p3t0, p3t1, uut0, uut1,                &
                     vvt0, vvt1, wwt0, wwt1, xmin, ymin, dx, dy, per, hem,      &
                     nx, ny, nz)

        ! Runge-Kutta timestep (x0,y0,p0 -> x1,y1,p1)
        else if (imethod.eq.2) then
          call runge(xx1, yy1, pp1, left(i), x0(i), y0(i), p0(i),               &
                     reltpos0, reltpos1, ts, numit, jflag,                      &
                     fbflag, spt0, spt1, p3t0, p3t1, uut0, uut1,                &
                     vvt0, vvt1, wwt0, wwt1, xmin, ymin, dx, dy, per, hem,      &
                     nx, ny, nz)
        endif

        ! Update trajectory position, or increase number of trajectories leaving
        ! domain
        if (left(i).ne.1) then
          x0(i)=xx1
          y0(i)=yy1
          p0(i)=pp1
        endif

      ! Trajectory has already left data domain (mark as <mdv>)
      else
        x0(i)=mdv
        y0(i)=mdv
        p0(i)=mdv
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
                       deltat,numit,jump,fbflag,                                &
                       spt0,spt1,p3d0,p3d1,uut0,uut1,vvt0,vvt1,wwt0,wwt1,       &
                        xmin,ymin,dx,dy,per,hem,nx,ny,nz)

      ! Declaration of subroutine parameters
      integer :: nx,ny,nz
      real :: x1,y1,p1
      integer :: left
      real :: x0,y0,p0
      real :: reltpos0,reltpos1
      real :: deltat
      integer :: numit
      integer :: jump
      integer :: fbflag
      real :: spt0(nx*ny)   ,spt1(nx*ny)
      real :: uut0(nx*ny*nz),uut1(nx*ny*nz)
      real :: vvt0(nx*ny*nz),vvt1(nx*ny*nz)
      real :: wwt0(nx*ny*nz),wwt1(nx*ny*nz)
      real :: p3d0(nx*ny*nz),p3d1(nx*ny*nz)
      real :: xmin,ymin,dx,dy
      real :: per
      integer :: hem

      ! Auxiliary variables
      real :: xind,yind,pind
      real :: u0,v0,w0,u1,v1,w1,u,v,w
      integer :: icount
      character :: ch

      ! Reset the flag for domain-leaving
      left=0

      ! Interpolate wind fields to starting position (x0,y0,p0)
      call get_index4 (xind,yind,pind,x0,y0,p0,reltpos0,&
                      p3d0,p3d1,spt0,spt1,1,&
                      nx,ny,nz,xmin,ymin,dx,dy,mdv)
      u0 = int_index4(uut0,uut1,nx,ny,nz,xind,yind,pind,reltpos0,mdv)
      v0 = int_index4(vvt0,vvt1,nx,ny,nz,xind,yind,pind,reltpos0,mdv)
      w0 = int_index4(wwt0,wwt1,nx,ny,nz,xind,yind,pind,reltpos0,mdv)

      ! Force the near-surface wind to zero
      if (pind.lt.1.) w0=w0*pind

      ! For first iteration take ending position equal to starting position
      x1=x0
      y1=y0
      p1=p0


      ! Iterative calculation of new position
      do icount=1,numit
          ! Calculate new winds for advection
         call get_index4 (xind,yind,pind,x1,y1,p1,reltpos1,&
                         p3d0,p3d1,spt0,spt1,1,&
                         nx,ny,nz,xmin,ymin,dx,dy,mdv)
         u1 = int_index4(uut0,uut1,nx,ny,nz,xind,yind,pind,reltpos1,mdv)
         v1 = int_index4(vvt0,vvt1,nx,ny,nz,xind,yind,pind,reltpos1,mdv)
         w1 = int_index4(wwt0,wwt1,nx,ny,nz,xind,yind,pind,reltpos1,mdv)

         ! Force the near-surface wind to zero
         if (pind.lt.1.) w1=w1*pind

         ! Get the new velocity in between
         u=(u0+u1)/2.
         v=(v0+v1)/2.
         w=(w0+w1)/2.

         ! Calculate new positions
         x1 = x0 + fbflag*u*deltat/(deltay*cos(y0*pi/180.))
         y1 = y0 + fbflag*v*deltat/deltay
         p1 = p0 + fbflag*w*deltat

         ! Check if trajectory leaves data domain
         call check_boundaries(x1,y1,p1,left,reltpos1,jump,&
                 spt0,spt1,p3d0,p3d1,&
                 xmin,ymin,dx,dy,per,hem,nx,ny,nz)
         if (left==1) then
             goto 100
         end if
      end do

!Exit point for subroutine
 100  continue

      return

      end subroutine euler

!-------------------------------------------------------------------
!Runge-Kutta (4th order) time-step
!-------------------------------------------------------------------
      subroutine runge(x1,y1,p1,left,x0,y0,p0,reltpos0,reltpos1,&
                      deltat,numit,jump,fbflag,&
                      spt0,spt1,p3d0,p3d1,uut0,uut1,vvt0,vvt1,wwt0,wwt1,&
                      xmin,ymin,dx,dy,per,hem,nx,ny,nz)

      ! Declaration of subroutine parameters
      integer :: nx,ny,nz
      real :: x1,y1,p1
      integer :: left
      real :: x0,y0,p0
      real :: reltpos0,reltpos1
      real :: deltat
      integer :: numit
      integer :: jump
      integer :: fbflag
      real :: spt0(nx*ny)   ,spt1(nx*ny)
      real :: uut0(nx*ny*nz),uut1(nx*ny*nz)
      real :: vvt0(nx*ny*nz),vvt1(nx*ny*nz)
      real :: wwt0(nx*ny*nz),wwt1(nx*ny*nz)
      real :: p3d0(nx*ny*nz),p3d1(nx*ny*nz)
      real :: xmin,ymin,dx,dy
      real :: per
      integer :: hem

      ! Auxiliary variables
      real :: xind,yind,pind
      real :: u0,v0,w0,u1,v1,w1,u,v,w
      integer :: icount,n
      real :: xs,ys,ps,xk(4),yk(4),pk(4)
      real :: reltpos

      ! Reset the flag for domain-leaving
      left=0

      ! Apply the Runge Kutta scheme
      do n=1,4
          ! Get intermediate position and relative time
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
                        p3d0,p3d1,spt0,spt1,1,&
                        nx,ny,nz,xmin,ymin,dx,dy,mdv)
        u = int_index4 (uut0,uut1,nx,ny,nz,xind,yind,pind,reltpos,mdv)
        v = int_index4 (vvt0,vvt1,nx,ny,nz,xind,yind,pind,reltpos,mdv)
        w = int_index4 (wwt0,wwt1,nx,ny,nz,xind,yind,pind,reltpos,mdv)


        !  Force the near-surface wind to zero
        if (pind.lt.1.) w1=w1*pind


        !  Update position and keep them
        xk(n)=fbflag*u*deltat/(deltay*cos(y0*pi/180.))
        yk(n)=fbflag*v*deltat/deltay
        pk(n)=fbflag*w*deltat

      end do

      ! Calculate new positions
      x1=x0+(1./6.)*(xk(1)+2.*xk(2)+2.*xk(3)+xk(4))
      y1=y0+(1./6.)*(yk(1)+2.*yk(2)+2.*yk(3)+yk(4))
      p1=p0+(1./6.)*(pk(1)+2.*pk(2)+2.*pk(3)+pk(4))


      ! Check if trajectory leaves data domain
      call check_boundaries(x1,y1,p1,left,reltpos1,jump,&
              spt0,spt1,p3d0,p3d1,&
              xmin,ymin,dx,dy,per,hem,nx,ny,nz)
      if (left==1) then
          goto 100
      end if

      !Exit point for subroutine
      100  continue

      return
      end subroutine runge

    subroutine check_boundaries(x1,y1,p1,left,reltpos,jump,&
                      spt0,spt1,p3d0,p3d1,&
                      xmin,ymin,dx,dy,per,hem,nx,ny,nz)

        !f2py real, intent(in) :: x1 y1, p1
        !f2py real, intent(in) :: reltpos
        !f2py integer, intent(in) :: jump
        !f2py real, intent(in) :: spt0(nx, ny), spt1(nx, ny)
        !f2py  real, intent(in) :: p3d0(nx, ny, nz), p3d1(nx, ny, nz)
        !f2py real, intent(in) :: xmin,ymin,dx,dy
        !f2py real, intent(in) :: per
        !f2py integer, intent(in) :: hem
        !f2py integer, intent(hide) :: nx,ny,nz
        !f2py integer, intent(inout) :: left

        ! Declaration of subroutine parameters
        real, intent(inout) :: x1,y1,p1
        integer, intent(inout) :: left
        real, intent(in) :: reltpos
        integer, intent(in) :: jump
        real, intent(in) :: spt0(nx*ny)   ,spt1(nx*ny)
        real, intent(in) :: p3d0(nx*ny*nz),p3d1(nx*ny*nz)
        real, intent(in) :: xmin,ymin,dx,dy
        real, intent(in) :: per
        integer, intent(in) :: hem
        integer, intent(in) :: nx,ny,nz

        ! Auxiliary variables
        real :: xmax,ymax
        real :: xind,yind,pind
        real :: sp, pb

        ! Set the east-north boundary of the domain
        xmax = xmin+real(nx-1)*dx
        ymax = ymin+real(ny-1)*dy

        ! Handle pole problems (crossing and near pole trajectory)
        if ((hem==1).and.(y1>90.)) then
            y1=180.-y1
            x1=x1+per/2.
        end if

        if ((hem==1).and.(y1<-90.)) then
            y1=-180.-y1
            x1=x1+per/2.
        end if

        ! Not sure about this check. Would be better to just say it leaves the
        ! domain for non-hemispheric data
        !if (y1>89.99) then
        !    y1=89.99
        !end if

        ! Handle crossings of the dateline
        if ((hem==1).and.(x1>xmin+per-dx)) then
            x1=xmin+amod(x1-xmin,per)
        end if

        if ((hem==1).and.(x1<xmin)) then
            x1=xmin+per+amod(x1-xmin,per)
        end if

        ! Interpolate surface pressure and lowest pressure to actual position
        ! Vertical position (p1) unimportant as we are interpolating a 2d
        ! field but using get_index4 and int_index4 to interpolate in time
        call get_index4 (&
                xind,yind,pind,x1,y1,0.,reltpos,&
                p3d0,p3d1,spt0,spt1,1,&
                nx,ny,nz,xmin,ymin,dx,dy,mdv)
        sp = int_index4 (spt0,spt1,nx,ny,1, xind,yind,1.,reltpos,mdv)
        pb = int_index4 (p3d0,p3d1,nx,ny,nz,xind,yind,1.,reltpos,mdv)

        ! Check if trajectory leaves data domain
        ! Check if vertical position is off the lower boundary in a way that works for
        ! decreasing (e.g. pressure) and increasing (e.g. height) coordinates
        if (p1<sp .and. sp<pb .or. p1>sp .and. sp>pd) then
            if (jump==1) then
                p1=pb
            else
                left=1
            end if
        else if (x1<xmin .or. x1>xmax-dx .or. y1<ymin .or. y1>ymax) then
            left=1
        end if
    end subroutine check_boundaries
end module caltra
