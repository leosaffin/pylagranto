module trace
  use inter, only: get_index3, int_index3

  implicit none

  contains
  !-----------------------------------------------------------------------------
  subroutine interp_to(tracer, xt, yt, zt, leftflag, z, surf,                   &
                       xmin, ymin, dx, dy, nx, ny, nz, ntra, values)


  !f2py integer, intent(in) :: nx, ny, nz
  !f2py integer, intent(in) :: ntra
  !f2py real, intent(in) :: tracer(nx*ny*nz)
  !f2py real, intent(in) :: xt(ntra), yt(ntra), zt(ntra)
  !f2py integer, intent(in) :: leftflag(ntra)
  !f2py real, intent(in) :: z(nx*ny*nz)
  !f2py real, intent(in) :: surf(nx*ny)
  !f2py real, intent(in) :: xmin, ymin
  !f2py real, intent(in) :: dx, dy
  !f2py real, intent(out) :: values(ntra)

  !-Input variables-------------------------------------------------------------
  ! Grid dimensions
  integer, intent(in) :: nx, ny, nz

  ! Number of trajectories
  integer, intent(in) :: ntra

  ! Variable to be traced
  real, intent(in) :: tracer(nx*ny*nz)

  ! Trajectory positions
  real, intent(in) :: xt(ntra), yt(ntra), zt(ntra)

  ! Trajectory in domain flag
  integer, intent(in) :: leftflag(ntra)

  ! Vertical coordinate
  real, intent(in) :: z(nx*ny*nz)

  ! Surface boundary
  real, intent(in) :: surf(nx*ny)

  ! Horizontal coordinate boundaries
  real, intent(in) :: xmin, ymin

  ! Horizontal grid spacing
  real, intent(in) :: dx, dy

  !-Output variables------------------------------------------------------------
  real, intent(out) :: values(ntra)

  !-Local variables-------------------------------------------------------------
  ! Interpolation indices
  real :: i,j,k

  ! loop counter
  integer :: n
  !-----------------------------------------------------------------------------

  do n=1,ntra
    call get_index3 (i, j, k, xt(n), yt(n), zt(n), 1, z, surf,                           &
                     nx, ny, nz, xmin, ymin, dx, dy)

    values(n) = int_index3 (tracer, nx, ny, nz, i, j, k, -1000.)
  end do

  end subroutine interp_to

end module trace
