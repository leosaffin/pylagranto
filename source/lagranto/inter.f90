module inter

  implicit none

  ! Numerical epsilon
  real, parameter :: eps=1.e-8

  contains

! *************************************************************
! * This package provides a general interpolaton routine      *
! *************************************************************

!The main interface routines are:
!    get_index3,4 : get the grid indices for interpolation
!    int_index3,4 : interpolate to the grid position

!-------------------------------------------------------------
!Get index in grid space for interpolation
!-------------------------------------------------------------

      subroutine get_index4 (rid,rjd,rkd,xpo,ypo,ppo,rtp,&
                            vert0,vert1,surf0,surf1,mode,&
                            nx,ny,nz,lonw,lats,dlon,dlat,misdat)

!Purpose:
!   This subroutine determines the indices (rid,rjd,rkd) in grid
!   space for a point in physical space (xpo,ypo,ppo). The
!   horizontal grid is specified by the south-west point (lats,lonw)
!   and the grid spacing (dlat,dlon). The vertical grid is given
!   by <vert(n1,n2,n3)>. The lower boundary (typicall surface
!   pressure) is given by <surf(n1,n2)>.
!Arguments:
!   rid,rjd,rkd  real  output   grid location to be interpolated to
!   xpo,ypo,ppo  real  input    physical coordinates
!   rtp          real  input    relative time position (0=beginning, 1=end)
!   n1,n2,n3     int   input    grid dimensions in x-, y- and p-direction
!   lats,lonw    real  input    south and west boundary of grid space
!   vert         real  input    vertical coordinate grid
!   surf         real  input    lower boundary (surface pressure)
!   mode         int   input    direction of vertical axis (p=1,th=-1)
!                                   1: linear, 1 -> nz (th)
!                                   2: linear, nz -> 1 (pv)
!                                   3: binary (p)

!f2py real, intent(out) :: rid,rjd,rkd
!f2py real, intent(in) :: xpo, ypo, ppo
!f2py real, intent(in) :: rtp
!f2py real, intent(in) :: vert0(nx*ny*nz), vert1(nx*ny*nz)
!f2py real, intent(in) :: surf0(nx*ny), surf1(nx*ny)
!f2py integer, intent(in) :: mode
!f2py integer, intent(in) :: nx, ny, nz
!f2py real, intent(in) :: lonw, lats, dlon, dlat
!f2py real, intent(in) :: misdat


!Declartion of function parameters
      integer   nx,ny,nz
      real      xpo,ypo,ppo,rtp
      real      vert0(nx*ny*nz),vert1(nx*ny*nz)
      real      surf0(nx*ny)   ,surf1(nx*ny)
      real      rid,rjd,rkd
      real      dlat,dlon,lats,lonw
      real      misdat
      integer   mode

!Auxiliary variables
      real      rid0,rjd0,rkd0,rid1,rjd1,rkd1

!Get the inidices
      if (abs(rtp).lt.eps) then
         call  get_index3 (rid,rjd,rkd,xpo,ypo,ppo,mode,&
                          vert0,surf0,nx,ny,nz,lonw,lats,dlon,dlat)
      elseif (abs(rtp-1.).lt.eps) then
         call  get_index3 (rid,rjd,rkd,xpo,ypo,ppo,mode,&
                          vert1,surf1,nx,ny,nz,lonw,lats,dlon,dlat)
      else
         call  get_index3 (rid0,rjd0,rkd0,xpo,ypo,ppo,mode,&
                         vert0,surf0,nx,ny,nz,lonw,lats,dlon,dlat)
         call  get_index3 (rid1,rjd1,rkd1,xpo,ypo,ppo,mode,&
                          vert1,surf1,nx,ny,nz,lonw,lats,dlon,dlat)
         rid = int_time (rid0,rid1,rtp,misdat)
         rjd = int_time (rjd0,rjd1,rtp,misdat)
         rkd = int_time (rkd0,rkd1,rtp,misdat)

      endif

      end subroutine get_index4

!-------------------------------------------------------------
!Interpolate to an arbitrary position in grid space and time
!-------------------------------------------------------------

      real function int_index4 (ar0,ar1,n1,n2,n3,rid,rjd,rkd,rtp,misdat)

!Purpose:
!   This subroutine interpolates a 3d-array to an arbitrary
!   location within the grid including a linear time-interpolation.
!Arguments:
!   rid,rjd,rkd  real  output   grid location to be interpolated to
!   xpo,ypo,ppo  real  input    physical coordinates
!   n1,n2,n3     int   input    grid dimensions in x-, y- and p-direction
!   lats,lonw    real  input    south and west boundary of grid space
!   vert         real  input    vertical coordinate grid
!   surf         real  input    lower boundary (surface pressure)

!Declartion of function parameters
      integer   n1,n2,n3
      real      ar0(n1*n2*n3),ar1(n1*n2*n3)
      real      rid,rjd,rkd
      real      rtp
      real      misdat

!Auxiliary variables
      real      val0,val1,val

!Do the 3d-interpolation
      if (abs(rtp).lt.eps) then
         val = int_index3 (ar0,n1,n2,n3,rid,rjd,rkd,misdat)
      elseif (abs(rtp-1.).lt.eps) then
         val = int_index3 (ar1,n1,n2,n3,rid,rjd,rkd,misdat)
      else
         val0 = int_index3 (ar0,n1,n2,n3,rid,rjd,rkd,misdat)
         val1 = int_index3 (ar1,n1,n2,n3,rid,rjd,rkd,misdat)
         val  = int_time (val0,val1,rtp,misdat)
      endif

!Return value
      int_index4 = val

      return
      end function int_index4


!-------------------------------------------------------------
!Interpolate to an arbitrary position in grid space
!-------------------------------------------------------------

      real function int_index3 (ar,n1,n2,n3,rid,rjd,rkd,misdat)

!Purpose:
!   This subroutine interpolates a 3d-array to an arbitrary
!   location within the grid. The interpolation includes the
!   testing of the missing data flag 'misdat'. If one dimension
!   is 1, a 2d-interpolation is performed; if two dimensions
!   are 1, it is a 1d-interpolation; if all three dimensions are
!   1, no interpolation is performed and the input value is
!   returned.
!Arguments:
!   ar        real  input   input data array
!   n1,n2,n3  int   input   dimensions of ar
!   ri,rj,rk  real  input   grid location to be interpolated to
!   misdat    real  input   missing data flag (on if misdat<>0)

!Declartion of function parameters
      integer   n1,n2,n3
      real      ar(n1*n2*n3)
      real      rid,rjd,rkd
      real      misdat

!Local variables
      integer   i,j,k,ip1,jp1,kp1
      real      frac0i,frac0j,frac0k,frac1i,frac1j,frac1k
      real      ri,rj,rk
      real      val000,val001,val010,val011,val100,val101,val110,val111
      real      frc000,frc001,frc010,frc011,frc100,frc101,frc110,frc111
      real      frc
      real      mdv
      real      val

!Elementary test for dimensions
      if ( (n1.lt.1).or.(n2.lt.1).or.(n3.lt.1) ) then
         print*,'Invalid grid dimensions ',n1,n2,n3
         stop
      endif

!Activate or inactive the missing data check (quick and dirty)
      if (misdat.ne.0.) then
         mdv = misdat
      else
         mdv = 257.22725394015
      endif

!Bring the indices into the grid space
      ri = amax1(1.,amin1(float(n1),rid))
      rj = amax1(1.,amin1(float(n2),rjd))
      rk = amax1(1.,amin1(float(n3),rkd))

!Get the index of the west-south-bottom corner of the box
      i   = min0(int(ri),n1-1)
      ip1 = i+1
      j   = min0(int(rj),n2-1)
      jp1 = j+1
      k   = min0(int(rk),n3-1)
      kp1 = k+1

!Special handling for 2d arrays
      if (n3.eq.1) then
         k=1
         kp1=1
      endif

!Get location relative to grid box
      if ( i.ne.ip1 ) then
         frac0i = ri-float(i)
         frac1i = 1.-frac0i
      else
         frac0i = 0.
         frac1i = 1.
      endif
      if ( j.ne.jp1 ) then
         frac0j = rj-float(j)
         frac1j = 1.-frac0j
      else
         frac0j = 0.
         frac1j = 1.
      endif
      if ( k.ne.kp1 ) then
         frac0k = rk-float(k)
         frac1k = 1.-frac0k
      else
         frac0k = 0.
         frac1k = 1.
      endif

!On a grid point - take the grid point value
      if ( ( abs(frac0i).lt.eps ).and.&
          ( abs(frac0j).lt.eps ).and.&
          ( abs(frac0k).lt.eps ) ) then
         
         val = ar( i + n1*(j -1) + n1*n2*(k -1) )
         goto 100
         
      endif

!Init the fractions
      frc000 = frac1i * frac1j * frac1k
      frc001 = frac0i * frac1j * frac1k
      frc010 = frac1i * frac0j * frac1k
      frc011 = frac0i * frac0j * frac1k
      frc100 = frac1i * frac1j * frac0k
      frc101 = frac0i * frac1j * frac0k
      frc110 = frac1i * frac0j * frac0k
      frc111 = frac0i * frac0j * frac0k

!Init the values
      val000 = ar( i   + n1*(j  -1) + n1*n2*(k  -1) )
      val001 = ar( ip1 + n1*(j  -1) + n1*n2*(k  -1) )
      val010 = ar( i   + n1*(jp1-1) + n1*n2*(k  -1) )
      val011 = ar( ip1 + n1*(jp1-1) + n1*n2*(k  -1) )
      val100 = ar( i   + n1*(j  -1) + n1*n2*(kp1-1) )
      val101 = ar( ip1 + n1*(j  -1) + n1*n2*(kp1-1) )
      val110 = ar( i   + n1*(jp1-1) + n1*n2*(kp1-1) )
      val111 = ar( ip1 + n1*(jp1-1) + n1*n2*(kp1-1) )

!Handle missing data
      if ( abs(val000-mdv).lt.eps ) frc000 = 0.
      if ( abs(val001-mdv).lt.eps ) frc001 = 0.
      if ( abs(val010-mdv).lt.eps ) frc010 = 0.
      if ( abs(val011-mdv).lt.eps ) frc011 = 0.
      if ( abs(val100-mdv).lt.eps ) frc100 = 0.
      if ( abs(val101-mdv).lt.eps ) frc101 = 0.
      if ( abs(val110-mdv).lt.eps ) frc110 = 0.
      if ( abs(val111-mdv).lt.eps ) frc111 = 0.

!Build the final value
      frc = frc000 + frc001 + frc010 + frc011 + &
           frc100 + frc101 + frc110 + frc111
      if ( frc.gt.0. ) then
         val = 1./frc * ( frc000 * val000 + frc001 * val001 +&
                         frc010 * val010 + frc011 * val011 +&
                         frc100 * val100 + frc101 * val101 +&
                         frc110 * val110 + frc111 * val111 )
      else
         val = misdat
      endif

!Return the value
 100  continue

      int_index3 = val

      end function int_index3


!-------------------------------------------------------------
!Time interpolation
!-------------------------------------------------------------

      real function int_time (val0,val1,reltpos,misdat)

!Purpose:
!   This subroutine interpolates linearly in time between two
!   values.
!Arguments:
!   val0      real  input   value at time 0
!   val1      real  input   value at time 1
!   reltpos   real  input   relative time (between 0 and 1)
!   misdat    real  input   missing data flag (on if misdat<>0)

!Declaration of parameters
      real      val0
      real      val1
      real      reltpos
      real      misdat

!Local variables
      real      val
      real      mdv

!Activate or inactive the missing data check (quick and dirty)
      if (misdat.ne.0.) then
         mdv = misdat
      else
         mdv = 257.22725394015
      endif

!Do the linear interpolation
      if ( abs(reltpos).lt.eps ) then
         val = val0
      elseif ( abs(reltpos-1.).lt.eps ) then
         val = val1
      elseif ( (abs(val0-mdv).gt.eps).and.&
               (abs(val1-mdv).gt.eps) ) then
         val = (1.-reltpos)*val0+reltpos*val1
      else
         val = mdv
      endif

!Return value
      int_time = val

      end function int_time


!-------------------------------------------------------------
!Get the position of a physical point in grid space
!-------------------------------------------------------------

      subroutine get_index3 (rid,rjd,rkd,xpo,ypo,ppo,mode,&
                            vert,surf,nx,ny,nz,lonw,lats,dlon,dlat)

!Purpose:
!   This subroutine determines the indices (rid,rjd,rkd) in grid
!   space for a point in physical space (xpo,ypo,ppo). The
!   horizontal grid is specified by the south-west point (lats,lonw)
!   and the grid spacing (dlat,dlon). The vertical grid is given
!   by <vert(n1,n2,n3)>. The lower boundary (typicall surface
!   pressure) is given by <surf(n1,n2)>.
!Arguments:
!   rid,rjd,rkd  real  output   grid location to be interpolated to
!   xpo,ypo,ppo  real  input    physical coordinates
!   n1,n2,n3     int   input    grid dimensions in x-, y- and p-direction
!   lats,lonw    real  input    south and west boundary of grid space
!   vert         real  input    vertical coordinate grid
!   surf         real  input    lower boundary (surface pressure)
!   mode         int   input    direction of vertical axis
!                                   1: linear, 1 -> nz (th)
!                                   2: linear, nz -> 1 (pv)
!                                   3: binary (p)
 
!f2py real, intent(out) :: rid,rjd,rkd
!f2py real, intent(in) :: xpo, ypo, ppo
!f2py real, intent(in) :: vert(nx,ny,nz)
!f2py real, intent(in) :: surf(nx,ny)
!f2py integer, intent(in) :: mode
!f2py integer, intent(hide) :: nx, ny, nz
!f2py real, intent(in) :: lonw, lats, dlon, dlat

!Declartion of function parameters
      integer   nx,ny,nz
      real      vert(nx*ny*nz)
      real      surf(nx*ny)
      real      rid,rjd,rkd
      real      xpo,ypo,ppo
      real      dlat,dlon,lats,lonw
      integer   mode

      !Local variables
      integer   i,j,k
      real      ppo0,ppo1,ppom,psur
      integer   i0,im,i1

      !Get the horizontal grid indices
      rid=(xpo-lonw)/dlon+1.
      rjd=(ypo-lats)/dlat+1.

      !Two-dimensional interpolation on horizontal plane: return 1.
      if (nz == 1) then
         rkd = 1.
      else
        !Get the pressure at the lowest level and at the surface
        ppo0 = int_index3(vert,nx,ny,nz,rid,rjd,real(1),0.)
        psur = int_index3(surf,nx,ny, 1,rid,rjd,real(1),0.)

        !The point is between the surface and the lowest level: return
        if ( between(ppo, psur, ppo0) ) then
          rkd  = (psur-ppo)/(psur-ppo0)
        else

          !Full-level search (TH): linear ascending scanning through all levels
          if ( mode == 1 ) then
            rkd=0
            i = 1
            do while (i<nz .and. rkd == 0)
              ppo1 = int_index3(vert,nx,ny,nz,rid,rjd,real(i+1),0.)
              if ( between(ppo, ppo0, ppo1) ) then
                rkd=real(i)+(ppo0-ppo)/(ppo0-ppo1)
              endif
              ppo0 = ppo1
              i = i+1
            end do

          !Full-level search (PV): linear descending scanning through all levels
          elseif ( mode == 2 ) then
            ppo1 = int_index3(vert,nx,ny,nz,rid,rjd,real(nz),0.)
            rkd=0
            i = nz - 1
            do while (i >= 1 .and. rkd == 0)
              ppo0 = int_index3(vert,nx,ny,nz,rid,rjd,real(i),0.)
              if ( between(ppo, ppo0, ppo1) ) then
                rkd=real(i)+(ppo0-ppo)/(ppo0-ppo1)
              endif
              ppo1 = ppo0
              i = i-1
            end do

          end if
        end if
      end if
      end subroutine get_index3

    logical function between(value, bound_1, bound_2)
        real, intent(in) :: value, bound_1, bound_2

        between = ((bound_1 <= value .AND. value <= bound_2) .OR. &
                   (bound_2 <= value .AND. value <= bound_1))

    end function between
end module inter
