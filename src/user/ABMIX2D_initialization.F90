module ABMIX2D_initialization
! This module sets the topography and initial thicknesses for the
! Abyssal mixing 2D simulations
!
! Ryan Holmes, June 2016
  
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_io, only : close_file, fieldtype, file_exists
use MOM_io, only : open_file, read_data, read_axis_data, SINGLE_FILE
use MOM_io, only : write_field, slasher, vardesc
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

! Public functions
public ABMIX2D_initialize_topography
public ABMIX2D_initialize_thickness

character(len=40) :: mod = "ABMIX2D_initialization" !< This module's name.

contains

!> Initialize topography with a shelf and slope in a 2D domain
subroutine ABMIX2D_initialize_topography ( D, G, param_file, max_depth )
  ! Arguments
  type(dyn_horgrid_type),             intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out) :: D !< Ocean bottom depth in m
  type(param_file_type),              intent(in)  :: param_file !< Parameter file structure
  real,                               intent(in)  :: max_depth  !< Maximum depth of model in m
  ! Local variables
  integer :: i, j
  real    :: y, yb, a, b, c, basin_width, idepth

  call get_param(param_file, mod, "ABMIX2D_SLOPE_DEPTH", idepth, &
                 'Depth of slope, as fraction of full depth, at y=0 in ABMIX2D configuration.', &
                 units='nondim',default=0.1)
  call get_param(param_file, mod, "ABMIX2D_BASIN_WIDTH", basin_width, &
                 'Width of deep ocean basin, as fraction of domain, in ABMIX2D configuration.', &
                 units='nondim',default=0.4)
  call get_param(param_file, mod, "ABMIX2D_SLOPE_CURV", c, &
                 'Quadratic coefficient of slope (sets curvature), in ABMIX2D configuration.', &
                 units='nondim',default=0.0)

  !! The following code sets a simple slope topography, which can be
  !! curved if abs(c)>0
  !!
  !! a = depth at y=0
  !! b = slope of bathymetry, chosen such that bathymetry intersects
  !! bottom at y = yb = 1 - basin_width.
  !! c = normalized curvature. c = 1 implies that the bathymetry slope
  !! is zero at y = 0
  
  yb = 1.0 - basin_width

  a = idepth * max_depth
  if (yb .gt. 0.0) then
     c = c * (max_depth - a) / yb / yb
     b = (max_depth - a - c * yb * yb) / yb 
  else
     b = 0.0
     c = 0.0
  endif
  
  do i=G%isc,G%iec
    do j=G%jsc,G%jec

      ! Compute normalized meridional coordinate
      y = ( G%geoLatT(i,j) - G%south_lat ) / G%len_lat;

      if ( y .lt. yb ) then
        D(i,j) = a + b * y + c * y * y
      else
        D(i,j) = max_depth
      end if

    enddo
  enddo
end subroutine ABMIX2D_initialize_topography

!> Initialize thicknesses according to coordinate mode
subroutine ABMIX2D_initialize_thickness ( h, G, GV, param_file )
  type(ocean_grid_type),                     intent(in)  :: G !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)  :: GV !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: h !< Layer thicknesses
  type(param_file_type),                     intent(in)  :: param_file !< Parameter file structure
  ! Local variables
  real :: e0(SZK_(G))     ! The resting interface heights, in m, usually !
                          ! negative because it is positive upward.      !
  real :: eta1D(SZK_(G)+1)! Interface height relative to the sea surface !
                          ! positive upward, in m.                       !
  integer :: i, j, k, is, ie, js, je, nz

  real    :: botlay_thickness, excess

  call get_param(param_file, mod, "ABMIX2D_INITIAL_BOTTH", botlay_thickness, &
                 'Thickness of bottom layer above max depth, in ABMIX2D configuration.', &
                 units='m',default=0.0)

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call MOM_mesg("MOM_initialization.F90, ABMIX2D_initialize_thickness: setting thickness")

  !! The following code sets a simple constant stratification

  ! WARNING: this routine specifies the interface heights so that the last layer
  !          is vanished, even at maximum depth. In order to have a uniform
  !          layer distribution, use this line of code within the loop:
  !          e0(k) = -G%max_depth * real(k-1) / real(nz)
  !          To obtain a thickness distribution where the last layer is
  !          vanished and the other thicknesses uniformly distributed, use:
  !          e0(k) = -G%max_depth * real(k-1) / real(nz-1)
  do k=1,nz
    e0(k) = -G%max_depth * real(k-1) / real(nz)
  enddo
  
  ! Thicken bottom layer:
  excess = botlay_thickness - (e0(nz)+G%max_depth)
  if (excess .gt. 0.0) then
     e0(nz) = e0(nz) + excess
     do k=1,nz-1
        e0(k) = e0(k) * (1.0 - excess / G%max_depth)
     enddo
  endif

  do j=js,je ; do i=is,ie
     eta1D(nz+1) = -1.0*G%bathyT(i,j)
     do k=nz,1,-1
        eta1D(k) = e0(k)
        if (eta1D(k) < (eta1D(k+1) + GV%Angstrom_z)) then
           eta1D(k) = eta1D(k+1) + GV%Angstrom_z
           h(i,j,k) = GV%Angstrom_z
        else
           h(i,j,k) = eta1D(k) - eta1D(k+1)
        endif
     enddo
  end do; end do

end subroutine ABMIX2D_initialize_thickness

end module ABMIX2D_initialization
