! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009-2012 Chris Brady <C.S.Brady@warwick.ac.uk>
! Copyright (C) 2011-2012 Martin Ramsay <M.G.Ramsay@warwick.ac.uk>
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE calc_df

  USE boundary

  IMPLICIT NONE

CONTAINS

  SUBROUTINE calc_boundary(data_array, species)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN), OPTIONAL :: species

    CALL processor_summation_bcs(data_array, ng, species=species)

  END SUBROUTINE calc_boundary



  SUBROUTINE calc_mass_density(data_array, current_species)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_m
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, idx
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    data_array = 0.0_num
    part_m = 0.0_num
    fac = 0.0_num

    idx = 1.0_num / dx / dy

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    DO ispecies = spec_start, spec_end
      IF (spec_sum .AND. &
          io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_m  = io_list(ispecies)%mass
      fac = io_list(ispecies)%weight
      wdata = part_m * fac

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_m  = current%mass
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
#endif
        wdata = part_m * fac
#else
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
        wdata = part_m * fac
#endif
#endif

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
        ENDDO
        ENDDO

        current => current%next
      ENDDO
      CALL calc_boundary(data_array, ispecies)
    ENDDO

    CALL calc_boundary(data_array)

    data_array = data_array * idx
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

  END SUBROUTINE calc_mass_density



  SUBROUTINE calc_ekbar(data_array, current_species)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_ux, part_uy, part_uz, part_mc, part_u2
    ! The weight of a particle
    REAL(num) :: part_w
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, gamma_rel, gamma_rel_m1
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: wt
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    ALLOCATE(wt(1-ng:nx+ng,1-ng:ny+ng))
    data_array = 0.0_num
    wt = 0.0_num
    part_mc  = 1.0_num
    part_w = 1.0_num

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    DO ispecies = spec_start, spec_end
      IF (spec_sum .AND. &
          io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_mc = c * io_list(ispecies)%mass
      part_w = io_list(ispecies)%weight
      fac = part_mc * part_w * c

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        fac = part_mc * part_w * c
#else
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
        fac = part_mc * part_w * c
#endif
#endif

        IF (io_list(ispecies)%species_type /= c_species_id_photon) THEN
          part_ux = current%part_p(1) / part_mc
          part_uy = current%part_p(2) / part_mc
          part_uz = current%part_p(3) / part_mc

          part_u2 = part_ux**2 + part_uy**2 + part_uz**2
          gamma_rel = SQRT(part_u2 + 1.0_num)
          gamma_rel_m1 = part_u2 / (gamma_rel + 1.0_num)

          wdata = gamma_rel_m1 * fac
        ELSE
#ifdef PHOTONS
          wdata = current%particle_energy * part_w
#else
          wdata = 0.0_num
#endif
        ENDIF

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
          wt(cell_x+ix, cell_y+iy) = &
              wt(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * part_w
        ENDDO
        ENDDO

        current => current%next
      ENDDO
      CALL calc_boundary(data_array, ispecies)
      CALL calc_boundary(wt, ispecies)
    ENDDO

    CALL calc_boundary(data_array)
    CALL calc_boundary(wt)

    data_array = data_array / MAX(wt, c_tiny)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

    DEALLOCATE(wt)

  END SUBROUTINE calc_ekbar



  SUBROUTINE calc_ekflux(data_array, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species, direction
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_ux, part_uy, part_uz, part_mc, part_u2
    ! The weight of a particle
    REAL(num) :: part_w
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, gamma_rel, gamma_rel_m1, part_flux, xfac, yfac, zfac
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: wt
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    ALLOCATE(wt(1-ng:nx+ng,1-ng:ny+ng))
    data_array = 0.0_num
    wt = 0.0_num
    part_mc  = 1.0_num
    part_w = 1.0_num

    xfac = c * dy
    yfac = c * dx
    zfac = c * dx * dy

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    DO ispecies = spec_start, spec_end
      IF (spec_sum .AND. &
          io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_mc = c * io_list(ispecies)%mass
      part_w = io_list(ispecies)%weight
      fac = part_mc * part_w * c

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        fac = part_mc * part_w * c
#else
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
        fac = part_mc * part_w * c
#endif
#endif

        IF (io_list(ispecies)%species_type /= c_species_id_photon) THEN
          part_ux = current%part_p(1) / part_mc
          part_uy = current%part_p(2) / part_mc
          part_uz = current%part_p(3) / part_mc

          part_u2 = part_ux**2 + part_uy**2 + part_uz**2
          gamma_rel = SQRT(part_u2 + 1.0_num)
          gamma_rel_m1 = part_u2 / (gamma_rel + 1.0_num)

          wdata = gamma_rel_m1 * fac
        ELSE
#ifdef PHOTONS
          fac = c / current%particle_energy
          part_ux = current%part_p(1) * fac
          part_uy = current%part_p(2) * fac
          part_uz = current%part_p(3) * fac
          gamma_rel = 1.0_num
          wdata = current%particle_energy * part_w
#else
          gamma_rel = 1.0_num
          wdata = 0.0_num
#endif
        ENDIF

        SELECT CASE(direction)
        CASE(-c_dir_x)
          ! negative flux in x
          part_flux = xfac * part_ux / gamma_rel
          wdata = -wdata * MIN(part_flux, 0.0_num)
        CASE( c_dir_x)
          ! positive flux in x
          part_flux = xfac * part_ux / gamma_rel
          wdata =  wdata * MAX(part_flux, 0.0_num)
        CASE(-c_dir_y)
          ! negative flux in y
          part_flux = yfac * part_uy / gamma_rel
          wdata = -wdata * MIN(part_flux, 0.0_num)
        CASE( c_dir_y)
          ! positive flux in y
          part_flux = yfac * part_uy / gamma_rel
          wdata =  wdata * MAX(part_flux, 0.0_num)
        CASE(-c_dir_z)
          ! negative flux in z
          part_flux = zfac * part_uz / gamma_rel
          wdata = -wdata * MIN(part_flux, 0.0_num)
        CASE( c_dir_z)
          ! positive flux in z
          part_flux = zfac * part_uz / gamma_rel
          wdata =  wdata * MAX(part_flux, 0.0_num)
        END SELECT

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
          wt(cell_x+ix, cell_y+iy) = &
              wt(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * part_w
        ENDDO
        ENDDO

        current => current%next
      ENDDO
      CALL calc_boundary(data_array, ispecies)
      CALL calc_boundary(wt, ispecies)
    ENDDO

    CALL calc_boundary(data_array)
    CALL calc_boundary(wt)

    data_array = data_array / MAX(wt, c_tiny)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

    DEALLOCATE(wt)

  END SUBROUTINE calc_ekflux



  SUBROUTINE calc_poynt_flux(data_array, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: direction
    INTEGER :: ix, iy
    REAL(num) :: ex_cc, ey_cc, ez_cc, bx_cc, by_cc, bz_cc

    SELECT CASE(direction)
    CASE(c_dir_x)
      DO iy = 1, ny
      DO ix = 1, nx
        ey_cc = 0.5_num  * (ey(ix  , iy-1) + ey(ix, iy  ))
        ez_cc = ez(ix, iy)
        by_cc = 0.5_num  * (by(ix-1, iy  ) + by(ix, iy  ))
        bz_cc = 0.25_num * (bz(ix-1, iy-1) + bz(ix, iy-1) &
                         +  bz(ix-1, iy  ) + bz(ix, iy  ))
        data_array(ix,iy) = (ey_cc * bz_cc - ez_cc * by_cc) / mu0
      ENDDO
      ENDDO
    CASE(c_dir_y)
      DO iy = 1, ny
      DO ix = 1, nx
        ex_cc = 0.5_num  * (ex(ix-1, iy  ) + ex(ix, iy  ))
        ez_cc = ez(ix, iy)
        bx_cc = 0.5_num  * (bx(ix  , iy-1) + bx(ix, iy  ))
        bz_cc = 0.25_num * (bz(ix-1, iy-1) + bz(ix, iy-1) &
                         +  bz(ix-1, iy  ) + bz(ix, iy  ))
        data_array(ix,iy) = (ez_cc * bx_cc - ex_cc * bz_cc) / mu0
      ENDDO
      ENDDO
    CASE(c_dir_z)
      DO iy = 1, ny
      DO ix = 1, nx
        ex_cc = 0.5_num  * (ex(ix-1, iy  ) + ex(ix, iy))
        ey_cc = 0.5_num  * (ey(ix  , iy-1) + ey(ix, iy))
        bx_cc = 0.5_num  * (bx(ix  , iy-1) + bx(ix, iy))
        by_cc = 0.5_num  * (by(ix-1, iy  ) + by(ix, iy))
        data_array(ix,iy) = (ex_cc * by_cc - ey_cc * bx_cc) / mu0
      ENDDO
      ENDDO
    END SELECT

  END SUBROUTINE calc_poynt_flux



  SUBROUTINE calc_charge_density(data_array, current_species)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_q
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, idx
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    data_array = 0.0_num
    part_q = 0.0_num
    fac = 0.0_num

    idx = 1.0_num / dx / dy

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    DO ispecies = spec_start, spec_end
      IF (spec_sum .AND. &
          io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_q  = io_list(ispecies)%charge
      fac = io_list(ispecies)%weight
      wdata = part_q * fac

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_q  = current%charge
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
#endif
        wdata = part_q * fac
#else
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
        wdata = part_q * fac
#endif
#endif

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
        ENDDO
        ENDDO

        current => current%next
      ENDDO
      CALL calc_boundary(data_array, ispecies)
    ENDDO

    CALL calc_boundary(data_array)

    data_array = data_array * idx
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

  END SUBROUTINE calc_charge_density



  SUBROUTINE calc_number_density(data_array, current_species)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: idx
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    data_array = 0.0_num

    idx = 1.0_num / dx / dy

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    DO ispecies = spec_start, spec_end
      IF (spec_sum .AND. &
          io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      wdata = io_list(ispecies)%weight

      DO WHILE (ASSOCIATED(current))
#ifndef PER_SPECIES_WEIGHT
        wdata = current%weight
#endif

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
        ENDDO
        ENDDO

        current => current%next
      ENDDO
      CALL calc_boundary(data_array, ispecies)
    ENDDO

    CALL calc_boundary(data_array)

    data_array = data_array * idx
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

  END SUBROUTINE calc_number_density



  SUBROUTINE calc_ppc(data_array, current_species)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    INTEGER :: ispecies, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
    REAL(num) :: cell_x_r, cell_y_r
    INTEGER :: cell_x, cell_y

    data_array = 0.0_num

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    DO ispecies = spec_start, spec_end
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head

      DO WHILE (ASSOCIATED(current))
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = (current%part_pos(1) - x_grid_min_local) / dx
        cell_y_r = (current%part_pos(2) - y_grid_min_local) / dy
#else
        cell_x_r = (current%part_pos(1) - x_grid_min_local) / dx + 0.5_num
        cell_y_r = (current%part_pos(2) - y_grid_min_local) / dy + 0.5_num
#endif
        cell_x = FLOOR(cell_x_r) + 1
        cell_y = FLOOR(cell_y_r) + 1

        data_array(cell_x,cell_y) = data_array(cell_x,cell_y) + 1.0_num

        current => current%next
      ENDDO
    ENDDO

  END SUBROUTINE calc_ppc



  SUBROUTINE calc_average_weight(data_array, current_species)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: part_count
    INTEGER :: ispecies, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
    REAL(num) :: cell_x_r, cell_y_r
    INTEGER :: cell_x, cell_y

    data_array = 0.0_num

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    ALLOCATE(part_count(1-ng:nx+ng,1-ng:ny+ng))
    part_count = 0.0_num

    DO ispecies = spec_start, spec_end
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      wdata = io_list(ispecies)%weight

      DO WHILE (ASSOCIATED(current))
#ifndef PER_SPECIES_WEIGHT
        wdata = current%weight
#endif

#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = (current%part_pos(1) - x_grid_min_local) / dx
        cell_y_r = (current%part_pos(2) - y_grid_min_local) / dy
#else
        cell_x_r = (current%part_pos(1) - x_grid_min_local) / dx + 0.5_num
        cell_y_r = (current%part_pos(2) - y_grid_min_local) / dy + 0.5_num
#endif
        cell_x = FLOOR(cell_x_r) + 1
        cell_y = FLOOR(cell_y_r) + 1

        data_array(cell_x,cell_y) = data_array(cell_x,cell_y) + wdata
        part_count(cell_x,cell_y) = part_count(cell_x,cell_y) + 1.0_num

        current => current%next
      ENDDO
    ENDDO

    data_array = data_array / MAX(part_count, c_tiny)

    DEALLOCATE(part_count)

  END SUBROUTINE calc_average_weight



  SUBROUTINE calc_temperature(sigma, current_species)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: sigma
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_pmx, part_pmy, part_pmz, sqrt_part_m
    ! The weight of a particle
    REAL(num) :: part_w
    REAL(num) :: gf
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: part_count, meanx, meany, meanz
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    ALLOCATE(meanx(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(meany(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(meanz(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(part_count(1-ng:nx+ng,1-ng:ny+ng))
    meanx = 0.0_num
    meany = 0.0_num
    meanz = 0.0_num
    part_count = 0.0_num
    sigma = 0.0_num

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    DO ispecies = spec_start, spec_end
      IF (spec_sum .AND. &
          io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      sqrt_part_m  = SQRT(io_list(ispecies)%mass)
      part_w = io_list(ispecies)%weight

      DO WHILE(ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
        sqrt_part_m  = SQRT(current%mass)
#endif
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        ! Copy the particle properties out for speed
        part_pmx = current%part_p(1) / sqrt_part_m
        part_pmy = current%part_p(2) / sqrt_part_m
        part_pmz = current%part_p(3) / sqrt_part_m

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          gf = gx(ix) * gy(iy) * part_w
          meanx(cell_x+ix, cell_y+iy) = &
              meanx(cell_x+ix, cell_y+iy) + gf * part_pmx
          meany(cell_x+ix, cell_y+iy) = &
              meany(cell_x+ix, cell_y+iy) + gf * part_pmy
          meanz(cell_x+ix, cell_y+iy) = &
              meanz(cell_x+ix, cell_y+iy) + gf * part_pmz
          part_count(cell_x+ix, cell_y+iy) = &
              part_count(cell_x+ix, cell_y+iy) + gf
        ENDDO
        ENDDO
        current => current%next
      ENDDO
      CALL calc_boundary(meanx, ispecies)
      CALL calc_boundary(meany, ispecies)
      CALL calc_boundary(meanz, ispecies)
      CALL calc_boundary(part_count, ispecies)
    ENDDO

    CALL calc_boundary(meanx)
    CALL calc_boundary(meany)
    CALL calc_boundary(meanz)
    CALL calc_boundary(part_count)

    part_count = MAX(part_count, 1.e-6_num)

    meanx = meanx / part_count
    meany = meany / part_count
    meanz = meanz / part_count

    part_count = 0.0_num
    DO ispecies = spec_start, spec_end
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      sqrt_part_m  = SQRT(io_list(ispecies)%mass)

      DO WHILE(ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
        sqrt_part_m  = SQRT(current%mass)
#endif
        ! Copy the particle properties out for speed
        part_pmx = current%part_p(1) / sqrt_part_m
        part_pmy = current%part_p(2) / sqrt_part_m
        part_pmz = current%part_p(3) / sqrt_part_m

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          gf = gx(ix) * gy(iy)
          sigma(cell_x+ix, cell_y+iy) = sigma(cell_x+ix, cell_y+iy) + gf &
              * ((part_pmx - meanx(cell_x+ix, cell_y+iy))**2 &
              + (part_pmy - meany(cell_x+ix, cell_y+iy))**2 &
              + (part_pmz - meanz(cell_x+ix, cell_y+iy))**2)
          part_count(cell_x+ix, cell_y+iy) = &
              part_count(cell_x+ix, cell_y+iy) + gf
        ENDDO
        ENDDO
        current => current%next
      ENDDO
      CALL calc_boundary(sigma, ispecies)
      CALL calc_boundary(part_count, ispecies)
    ENDDO

    CALL calc_boundary(sigma)
    CALL calc_boundary(part_count)

    ! 3/2 kT = <p^2>/(2m)
    sigma = sigma / MAX(part_count, 1.e-6_num) / kb / 3.0_num

    DEALLOCATE(part_count, meanx, meany, meanz)

  END SUBROUTINE calc_temperature



  SUBROUTINE calc_on_grid_with_evaluator(data_array, current_species, evaluator)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    INTERFACE
      FUNCTION evaluator(a_particle, species_eval)
        USE shared_data
        TYPE(particle), POINTER :: a_particle
        INTEGER, INTENT(IN) :: species_eval
        REAL(num) :: evaluator
      END FUNCTION evaluator
    END INTERFACE

    data_array = 0.0_num

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    DO ispecies = spec_start, spec_end
      IF (spec_sum .AND. &
          io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head

      DO WHILE (ASSOCIATED(current))
#include "particle_to_grid.inc"

        wdata = evaluator(current, ispecies)

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
        ENDDO
        ENDDO

        current => current%next
      ENDDO
      CALL calc_boundary(data_array, ispecies)
    ENDDO

    CALL calc_boundary(data_array)

    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

  END SUBROUTINE calc_on_grid_with_evaluator



  SUBROUTINE calc_per_species_current(data_array, current_species, direction)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species, direction
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_q, part_mc
    REAL(num) :: part_px, part_py, part_pz
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, idx, root
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    data_array = 0.0_num
    part_q = 0.0_num
    fac = 0.0_num

    idx = 1.0_num / dx / dy

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species <= 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    DO ispecies = spec_start, spec_end
      IF (spec_sum .AND. &
          io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_mc = c * io_list(ispecies)%mass
      part_q  = io_list(ispecies)%charge
      fac = io_list(ispecies)%weight
      wdata = part_q * fac

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
        part_q  = current%charge
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
#endif
        wdata = part_q * fac
#else
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
        wdata = part_q * fac
#endif
#endif

        ! Copy the particle properties out for speed
        part_px = current%part_p(1)
        part_py = current%part_p(2)
        part_pz = current%part_p(3)
        root = 1.0_num / SQRT(part_mc**2 + part_px**2 + part_py**2 + part_pz**2)
        SELECT CASE (direction)
          CASE(c_dir_x)
            wdata = wdata * part_px * root
          CASE(c_dir_y)
            wdata = wdata * part_py * root
          CASE(c_dir_z)
            wdata = wdata * part_pz * root
        END SELECT

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
        ENDDO
        ENDDO

        current => current%next
      ENDDO
      CALL calc_boundary(data_array, ispecies)
    ENDDO

    CALL calc_boundary(data_array)

    idx = c * idx
    data_array = data_array * idx
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

  END SUBROUTINE calc_per_species_current



  SUBROUTINE calc_per_species_jx(data_array, current_species)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    CALL calc_per_species_current(data_array, current_species, c_dir_x)

  END SUBROUTINE calc_per_species_jx



  SUBROUTINE calc_per_species_jy(data_array, current_species)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    CALL calc_per_species_current(data_array, current_species, c_dir_y)

  END SUBROUTINE calc_per_species_jy



  SUBROUTINE calc_per_species_jz(data_array, current_species)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    CALL calc_per_species_current(data_array, current_species, c_dir_z)

  END SUBROUTINE calc_per_species_jz



  SUBROUTINE calc_total_energy_sum

    REAL(num) :: particle_energy, field_energy
    REAL(num) :: part_ux, part_uy, part_uz, part_u2
    REAL(num) :: part_mc, part_w, fac, gamma_rel, gamma_rel_m1
    REAL(num) :: sum_out(2), sum_in(2)
    REAL(num), PARAMETER :: c2 = c**2
    INTEGER :: ispecies, i, j
    TYPE(particle), POINTER :: current

    particle_energy = 0.0_num

    ! Sum over all particles to calculate total kinetic energy
    DO ispecies = 1, n_species
    !It is strange to skip the photon calculation.
      IF (io_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (species_list(ispecies)%tracer) CYCLE
#endif
      current => species_list(ispecies)%attached_list%head
      part_mc = c * species_list(ispecies)%mass
      part_w = species_list(ispecies)%weight
      fac = part_mc * part_w * c

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        fac = part_mc * part_w * c
#else
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
        fac = part_mc * part_w * c
#endif
#endif

        IF (species_list(ispecies)%species_type /= c_species_id_photon) THEN
          part_ux = current%part_p(1) / part_mc
          part_uy = current%part_p(2) / part_mc
          part_uz = current%part_p(3) / part_mc

          part_u2 = part_ux**2 + part_uy**2 + part_uz**2
          gamma_rel = SQRT(part_u2 + 1.0_num)
          gamma_rel_m1 = part_u2 / (gamma_rel + 1.0_num)

          particle_energy = particle_energy + gamma_rel_m1 * fac
#ifdef PHOTONS
        ELSE
          particle_energy = particle_energy + current%particle_energy * part_w
#endif
        ENDIF

        current => current%next
      ENDDO
    ENDDO

    ! EM field energy
    field_energy = 0.0_num
    DO j = 1, ny
    DO i = 1, nx
      field_energy = field_energy + ex(i,j)**2 + ey(i,j)**2 &
          + ez(i,j)**2 + c2 * (bx(i,j)**2 + by(i,j)**2 + bz(i,j)**2)
    ENDDO
    ENDDO
    field_energy = 0.5_num * epsilon0 * field_energy * dx * dy

    sum_out(1) = particle_energy
    sum_out(2) = field_energy
    CALL MPI_REDUCE(sum_out, sum_in, 2, mpireal, MPI_SUM, 0, comm, errcode)
    total_particle_energy = sum_in(1)
    total_field_energy = sum_in(2)

  END SUBROUTINE calc_total_energy_sum



  SUBROUTINE calc_initial_current

    REAL(num), DIMENSION(:,:), ALLOCATABLE :: jx, jy, jz
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_q, part_mc
    REAL(num) :: part_px, part_py, part_pz
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: part_jx, part_jy, part_jz
    REAL(num) :: fac, idx, root
    REAL(num) :: sum_in(3), sum_out(3)
    INTEGER :: ispecies, ix, iy
    TYPE(particle), POINTER :: current
#include "particle_head.inc"

    ALLOCATE(jx(1-jng:nx+jng, 1-jng:ny+jng))
    ALLOCATE(jy(1-jng:nx+jng, 1-jng:ny+jng))
    ALLOCATE(jz(1-jng:nx+jng, 1-jng:ny+jng))

    jx = 0.0_num
    jy = 0.0_num
    jz = 0.0_num

    idx = 1.0_num / dx / dy

    DO ispecies = 1, n_species
      IF (species_list(ispecies)%species_type == c_species_id_photon) CYCLE
#ifndef NO_TRACER_PARTICLES
      IF (species_list(ispecies)%tracer) CYCLE
#endif
      current => species_list(ispecies)%attached_list%head
      part_mc = c * species_list(ispecies)%mass
      part_q  = species_list(ispecies)%charge
      fac = species_list(ispecies)%weight
      wdata = part_q * fac

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
        part_q  = current%charge
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
#endif
        wdata = part_q * fac
#else
#ifndef PER_SPECIES_WEIGHT
        fac = current%weight
        wdata = part_q * fac
#endif
#endif

        ! Copy the particle properties out for speed
        part_px = current%part_p(1)
        part_py = current%part_p(2)
        part_pz = current%part_p(3)
        root = 1.0_num / SQRT(part_mc**2 + part_px**2 + part_py**2 + part_pz**2)

        part_jx = wdata * part_px * root
        part_jy = wdata * part_py * root
        part_jz = wdata * part_pz * root

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          jx(cell_x+ix, cell_y+iy) = &
              jx(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * part_jx
          jy(cell_x+ix, cell_y+iy) = &
              jy(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * part_jy
          jz(cell_x+ix, cell_y+iy) = &
              jz(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * part_jz
        ENDDO
        ENDDO

        current => current%next
      ENDDO
    ENDDO

    sum_out(1) = SUM(jx)
    sum_out(2) = SUM(jy)
    sum_out(3) = SUM(jz)
    DEALLOCATE(jx, jy, jz)

    CALL MPI_ALLREDUCE(sum_out, sum_in, 3, mpireal, MPI_SUM, comm, errcode)

    fac = c * idx / nx_global / ny_global

    initial_jx = sum_in(1) * fac
    initial_jy = sum_in(2) * fac
    initial_jz = sum_in(3) * fac

  END SUBROUTINE calc_initial_current

END MODULE calc_df
