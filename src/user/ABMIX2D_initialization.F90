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
  real    :: y, l1
  real    :: idepth, basin_width

  call get_param(param_file, mod, "ABMIX2D_SLOPE_DEPTH", idepth, &
                 'Depth of slope, as fraction of full depth, at y=0 in ABMIX2D configuration.', &
                 units='nondim',default=0.1)
  call get_param(param_file, mod, "ABMIX2D_BASIN_WIDTH", basin_width, &
                 'Width of deep ocean basin, as fraction of domain, in ABMIX2D configuration.', &
                 units='nondim',default=0.4)

  !! The following code sets a simple slope topography

  ! location where downslope reaches maximum depth
  l1 = 1.0 - basin_width

  do i=G%isc,G%iec
    do j=G%jsc,G%jec

      ! Compute normalized meridional coordinate
      y = ( G%geoLatT(i,j) - G%south_lat ) / G%len_lat;

      if ( y .lt. l1 ) then
        D(i,j) = idepth * max_depth + (1.0-idepth) * max_depth * &
                 ( y / l1)
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
