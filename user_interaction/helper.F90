! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009-2010 Chris Brady <C.S.Brady@warwick.ac.uk>
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

MODULE helper

  USE balance
  USE boundary
  USE strings
  USE partlist
  USE calc_df
  USE simple_io
  USE deltaf_loader

  IMPLICIT NONE

CONTAINS

  SUBROUTINE set_thermal_bcs

    INTEGER :: ispecies
    TYPE(particle_species), POINTER :: species

    DO ispecies = 1, n_species
      species => species_list(ispecies)

      ! Set temperature at boundary for thermal bcs.

      IF (species%bc_particle(c_bd_x_min) == c_bc_thermal) THEN
        species_list(ispecies)%ext_temp_x_min(:,:) = &
            species_list(ispecies)%initial_conditions%temp(1,:,:)
      ENDIF
      IF (species%bc_particle(c_bd_x_max) == c_bc_thermal) THEN
        species_list(ispecies)%ext_temp_x_max(:,:) = &
            species_list(ispecies)%initial_conditions%temp(nx,:,:)
      ENDIF
      IF (species%bc_particle(c_bd_y_min) == c_bc_thermal) THEN
        species_list(ispecies)%ext_temp_y_min(:,:) = &
            species_list(ispecies)%initial_conditions%temp(:,1,:)
      ENDIF
      IF (species%bc_particle(c_bd_y_max) == c_bc_thermal) THEN
        species_list(ispecies)%ext_temp_y_max(:,:) = &
            species_list(ispecies)%initial_conditions%temp(:,ny,:)
      ENDIF
    ENDDO

  END SUBROUTINE set_thermal_bcs



  SUBROUTINE auto_load

    INTEGER :: ispecies
    TYPE(particle_species), POINTER :: species

    CALL set_thermal_bcs

    DO ispecies = 1, n_species
      species => species_list(ispecies)

#ifndef PER_SPECIES_WEIGHT
      CALL setup_particle_density(&
          species_list(ispecies)%initial_conditions%density, species, &
          species_list(ispecies)%initial_conditions%density_min, &
          species_list(ispecies)%initial_conditions%density_max)
#else
      CALL non_uniform_load_particles(&
          species_list(ispecies)%initial_conditions%density, species, &
          species_list(ispecies)%initial_conditions%density_min, &
          species_list(ispecies)%initial_conditions%density_max)
#endif
!   !Xiey  only drift in z-direction.
    IF (species_list(ispecies)%species_type == c_species_id_photon) THEN
      CALL setup_particle_temperature(&
          species_list(ispecies)%initial_conditions%temp(:,:,1), c_dir_x, &
          species, species_list(ispecies)%initial_conditions%drift(:,:,1))
      CALL setup_particle_temperature(&
          species_list(ispecies)%initial_conditions%temp(:,:,2), c_dir_y, &
          species, species_list(ispecies)%initial_conditions%drift(:,:,2))
      CALL setup_particle_temperature(&
          species_list(ispecies)%initial_conditions%temp(:,:,3), c_dir_z, &
          species, species_list(ispecies)%initial_conditions%drift(:,:,3))
    ELSE
      CALL setup_particle_temperature(&
          species_list(ispecies)%initial_conditions%temp(:,:,3), c_dir_z, &
          species, species_list(ispecies)%initial_conditions%drift(:,:,3))
    ENDIF

    ENDDO

    IF (rank == 0) THEN
      DO ispecies = 1, n_species
        species => species_list(ispecies)
        IF (species%count < 0) THEN
          WRITE(*,*) 'No particles specified for species ', &
              '"' // TRIM(species%name) // '"'
#ifndef NO_IO
          WRITE(stat_unit,*) 'No particles specified for species ', &
              '"' // TRIM(species%name) // '"'
#endif
          species%count = 0
        ENDIF
      ENDDO
    ENDIF

    CALL deltaf_load

  END SUBROUTINE auto_load



  SUBROUTINE allocate_ic

    INTEGER :: ispecies

    DO ispecies = 1, n_species
      ALLOCATE(species_list(ispecies)%initial_conditions&
          %density(1-ng:nx+ng,1-ng:ny+ng))
      ALLOCATE(species_list(ispecies)%initial_conditions&
          %temp(1-ng:nx+ng,1-ng:ny+ng,1:3))
      ALLOCATE(species_list(ispecies)%initial_conditions&
          %drift(1-ng:nx+ng,1-ng:ny+ng,1:3))

      species_list(ispecies)%initial_conditions%density = 1.0_num
      species_list(ispecies)%initial_conditions%temp = 0.0_num
      species_list(ispecies)%initial_conditions%drift = 0.0_num
      species_list(ispecies)%initial_conditions%density_min = EPSILON(1.0_num)
      species_list(ispecies)%initial_conditions%density_max = HUGE(1.0_num)
      species_list(ispecies)%initial_conditions%density_back = 0.0_num
      species_list(ispecies)%initial_conditions%temp_back = 0.0_num
      species_list(ispecies)%initial_conditions%drift_back = 0.0_num
    ENDDO

  END SUBROUTINE allocate_ic



  SUBROUTINE deallocate_ic

    INTEGER :: ispecies

    DO ispecies = 1, n_species
      DEALLOCATE(species_list(ispecies)%initial_conditions%density)
      DEALLOCATE(species_list(ispecies)%initial_conditions%temp)
      DEALLOCATE(species_list(ispecies)%initial_conditions%drift)
    ENDDO

  END SUBROUTINE deallocate_ic



#ifdef PER_SPECIES_WEIGHT
  SUBROUTINE non_uniform_load_particles(density, species, density_min, &
      density_max)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(INOUT) :: density
    TYPE(particle_species), POINTER :: species
    REAL(num), INTENT(INOUT) :: density_min, density_max
    INTEGER(i8) :: num_valid_cells_local, num_valid_cells_global
    INTEGER(i8) :: npart_per_cell
    REAL(num) :: density_total, density_total_global, density_average
    REAL(num) :: npart_per_cell_average
    INTEGER(i8) :: npart_this_proc_new, ipart, npart_this_species
    INTEGER :: ix, iy
    CHARACTER(LEN=15) :: string
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current, next

    partlist => species%attached_list

    num_valid_cells_local = 0
    density_total = 0.0_num

    DO iy = 1-ng, ny+ng
    DO ix = 1-ng, nx+ng
      IF (density(ix,iy) > density_max) density(ix,iy) = density_max
    ENDDO ! ix
    ENDDO ! iy

    DO iy = 1, ny
    DO ix = 1, nx
      IF (density(ix,iy) >= density_min) THEN
        num_valid_cells_local = num_valid_cells_local + 1
        density_total = density_total + density(ix,iy)
      ENDIF
    ENDDO ! ix
    ENDDO ! iy

    CALL MPI_ALLREDUCE(num_valid_cells_local, num_valid_cells_global, 1, &
        MPI_INTEGER8, MPI_SUM, comm, errcode)

    IF (species%npart_per_cell >= 0) THEN
      npart_per_cell_average = FLOOR(species%npart_per_cell, num)
    ELSE
      npart_per_cell_average = REAL(species%count, num) &
          / REAL(num_valid_cells_global, num)
    ENDIF

    IF (npart_per_cell_average <= 0) RETURN

    CALL MPI_ALLREDUCE(density_total, density_total_global, 1, mpireal, &
        MPI_SUM, comm, errcode)
    density_average = density_total_global / REAL(num_valid_cells_global, num)

    npart_this_proc_new = 0
    DO iy = 1, ny
    DO ix = 1, nx
      npart_per_cell = NINT(density(ix, iy) / density_average &
          * npart_per_cell_average)
      npart_this_proc_new = npart_this_proc_new + npart_per_cell
    ENDDO ! ix
    ENDDO ! iy

    CALL destroy_partlist(partlist)
    CALL create_allocated_partlist(partlist, npart_this_proc_new)

    ! Randomly place npart_per_cell particles into each valid cell
    current => partlist%head
    DO iy = 1, ny
    DO ix = 1, nx
      npart_per_cell = NINT(density(ix, iy) / density_average &
          * npart_per_cell_average)

      ipart = 0
      DO WHILE(ASSOCIATED(current) .AND. ipart < npart_per_cell)
#ifdef PER_PARTICLE_CHARGE_MASS
        ! Even if particles have per particle charge and mass, assume
        ! that initially they all have the same charge and mass (user
        ! can easily over_ride)
        current%charge = species%charge
        current%mass = species%mass
#endif
#ifdef DELTAF_METHOD
        ! Store the number of particles per cell to allow calculation
        ! of phase space volume later
        current%pvol = npart_per_cell
#endif
        current%part_pos(1) = x(ix) + (random() - 0.5_num) * dx
        current%part_pos(2) = y(iy) + (random() - 0.5_num) * dy

        ipart = ipart + 1
        current => current%next
      ENDDO
    ENDDO ! ix
    ENDDO ! iy

    ! Remove any unplaced particles from the list. This should never be
    ! called if the above routines worked correctly.
    DO WHILE(ASSOCIATED(current))
      next => current%next
      CALL remove_particle_from_partlist(partlist, current)
      DEALLOCATE(current)
      current => next
    ENDDO

    CALL MPI_ALLREDUCE(partlist%count, npart_this_species, 1, MPI_INTEGER8, &
        MPI_SUM, comm, errcode)

    species%count = npart_this_species
    species%weight = density_total_global * dx * dy / npart_this_species

    IF (rank == 0) THEN
      CALL integer_as_string(npart_this_species, string)
      WRITE(*,*) 'Loaded ', TRIM(ADJUSTL(string)), &
          ' particles of species ', '"' // TRIM(species%name) // '"'
#ifndef NO_IO
      WRITE(stat_unit,*) 'Loaded ', TRIM(ADJUSTL(string)), &
          ' particles of species ', '"' // TRIM(species%name) // '"'
#endif
    ENDIF

    CALL particle_bcs

  END SUBROUTINE non_uniform_load_particles
#endif



  ! This subroutine automatically loads a uniform density of pseudoparticles
  SUBROUTINE load_particles(species, load_list)

    TYPE(particle_species), POINTER :: species
    LOGICAL, DIMENSION(1-ng:,1-ng:), INTENT(IN) :: load_list
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: valid_cell_list
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current, next
    INTEGER(i8) :: ipart, npart_per_cell, num_int, num_total, idx
    INTEGER(i8) :: num_valid_cells_local, num_valid_cells_global
    INTEGER(i8) :: npart_this_species, num_new_particles, npart_left
    INTEGER(i8), ALLOCATABLE :: num_valid_cells_all(:), num_idx(:)
    REAL(num) :: valid_cell_frac, num_real, f0, f1
    REAL(num), ALLOCATABLE :: num_frac(:)
    INTEGER(i8) :: cell_x
    INTEGER(i8) :: cell_y
    INTEGER(i8) :: i, ipos
    INTEGER :: ix, iy, nx_e
    INTEGER :: ix_min, ix_max, iy_min, iy_max
    CHARACTER(LEN=15) :: string
    LOGICAL :: sweep

    npart_this_species = species%count
    IF (npart_this_species <= 0) RETURN

    ix_min = 1
    ix_max = nx

    iy_min = 1
    iy_max = ny

    IF (species%fill_ghosts) THEN
      IF (x_min_boundary) ix_min = ix_min - png
      IF (x_max_boundary) ix_max = ix_max + png

      IF (y_min_boundary) iy_min = iy_min - png
      IF (y_max_boundary) iy_max = iy_max + png
    ENDIF

    nx_e = ix_max - ix_min + 1

    num_valid_cells_local = 0
    DO iy = iy_min, iy_max
    DO ix = ix_min, ix_max
      IF (load_list(ix, iy)) num_valid_cells_local = num_valid_cells_local + 1
    ENDDO ! ix
    ENDDO ! iy

    IF (species%npart_per_cell >= 0) THEN
      npart_per_cell = FLOOR(species%npart_per_cell, KIND=i8)
      num_new_particles = &
          FLOOR(species%npart_per_cell * num_valid_cells_local, KIND=i8)
    ELSE
      ALLOCATE(num_valid_cells_all(nproc), num_idx(nproc), num_frac(nproc))

      ! Calculate global number of particles per cell
      CALL MPI_ALLGATHER(num_valid_cells_local, 1, MPI_INTEGER8, &
          num_valid_cells_all, 1, MPI_INTEGER8, comm, errcode)

      num_valid_cells_global = 0
      DO i = 1,nproc
        num_valid_cells_global = num_valid_cells_global + num_valid_cells_all(i)
      ENDDO

      IF (num_valid_cells_global == 0) THEN
        IF (rank == 0) THEN
          WRITE(*,*) '*** ERROR ***'
          WRITE(*,*) 'Intial condition settings mean that there are no cells ' &
              // 'where particles may'
          WRITE(*,*) 'validly be placed for species "' // TRIM(species%name) &
              // '". ', 'Code will now terminate.'
          CALL abort_code(c_err_bad_setup)
        ENDIF
      ENDIF

      valid_cell_frac = REAL(num_valid_cells_local, num) &
          / REAL(num_valid_cells_global, num)
      num_real = npart_this_species * valid_cell_frac
      num_new_particles = FLOOR(num_real, KIND=i8)

      ! Work out which processors get the remaining fractional numbers
      ! of particles

      ! Get a list of the fractional part on each processor, along with
      ! the total
      num_total = 0
      DO i = 1,nproc
        valid_cell_frac = REAL(num_valid_cells_all(i), num) &
            / REAL(num_valid_cells_global, num)
        num_real = npart_this_species * valid_cell_frac
        num_int = FLOOR(num_real, KIND=i8)
        num_frac(i) = num_real - num_int
        num_idx (i) = i - 1
        num_total = num_total + num_int
      ENDDO
      num_total = npart_this_species - num_total

      IF (num_total > 0) THEN
        ! Sort the list of fractions into decreasing order using bubble sort
        sweep = .TRUE.
        DO WHILE(sweep)
          sweep = .FALSE.
          f0 = num_frac(1)
          DO i = 2,nproc
            f1 = num_frac(i)
            IF (f1 > f0) THEN
              num_frac(i-1) = f1
              num_frac(i) = f0
              f1 = f0
              idx = num_idx(i-1)
              num_idx(i-1) = num_idx(i)
              num_idx(i) = idx
              sweep = .TRUE.
            ENDIF
            f0 = f1
          ENDDO
        ENDDO

        ! Accumulate fractional particles until they have all been accounted
        ! for. If any of them have been assigned to the current processor,
        ! add them and exit the loop.

        DO i = 1,nproc
          IF (num_idx(i) == rank) THEN
            num_new_particles = num_new_particles + 1
            EXIT
          ENDIF
          num_total = num_total - 1
          IF (num_total <= 0) EXIT
        ENDDO
      ENDIF

      DEALLOCATE(num_valid_cells_all, num_idx, num_frac)

      species%npart_per_cell = &
          REAL(npart_this_species,num) / num_valid_cells_global
      npart_per_cell = FLOOR(species%npart_per_cell, KIND=i8)
    ENDIF

    partlist => species%attached_list

    CALL destroy_partlist(partlist)
    CALL create_allocated_partlist(partlist, num_new_particles)

    ! Randomly place npart_per_cell particles into each valid cell
    npart_left = num_new_particles
    current => partlist%head
    IF (npart_per_cell > 0) THEN

      DO iy = iy_min, iy_max
      DO ix = ix_min, ix_max
        IF (.NOT. load_list(ix, iy)) CYCLE

        ipart = 0
        DO WHILE(ASSOCIATED(current) .AND. ipart < npart_per_cell)
#ifdef PER_PARTICLE_CHARGE_MASS
          ! Even if particles have per particle charge and mass, assume
          ! that initially they all have the same charge and mass (user
          ! can easily over_ride)
          current%charge = species%charge
          current%mass = species%mass
#endif
#ifdef DELTAF_METHOD
          ! Store the number of particles per cell to allow calculation
          ! of phase space volume later
          current%pvol = npart_per_cell
#endif
          current%part_pos(1) = x(ix) + (random() - 0.5_num) * dx
          current%part_pos(2) = y(iy) + (random() - 0.5_num) * dy

          ipart = ipart + 1
          current => current%next

          ! One particle sucessfully placed
          npart_left = npart_left - 1
        ENDDO
      ENDDO ! ix
      ENDDO ! iy

    ENDIF

    ! When num_new_particles does not equal
    ! npart_per_cell * num_valid_cells_local there will be particles left
    ! over that didn't get placed.
    ! The following loop randomly place remaining particles into valid cells.
    IF (npart_left > 0) THEN
      ALLOCATE(valid_cell_list(num_valid_cells_local))

      ipos = 0
      DO iy = iy_min, iy_max
      DO ix = ix_min, ix_max
        IF (load_list(ix,iy)) THEN
          ipos = ipos + 1
          valid_cell_list(ipos) = ix - ix_min + nx_e * (iy - iy_min)
        ENDIF
      ENDDO ! ix
      ENDDO ! iy

      DO i = 1, npart_left
        ipos = INT(random() * (num_valid_cells_local - 1)) + 1
        ipos = valid_cell_list(ipos)

        cell_y = ipos / nx_e
        ipos = ipos - nx_e * cell_y
        cell_y = cell_y + iy_min

        cell_x = ipos + ix_min

#ifdef PER_PARTICLE_CHARGE_MASS
        ! Even if particles have per particle charge and mass, assume
        ! that initially they all have the same charge and mass (user
        ! can easily over_ride)
        current%charge = species%charge
        current%mass = species%mass
#endif
#ifdef DELTAF_METHOD
        ! Store the number of particles per cell to allow calculation
        ! of phase space volume later
        current%pvol = npart_per_cell
#endif
        current%part_pos(1) = x(cell_x) + (random() - 0.5_num) * dx
        current%part_pos(2) = y(cell_y) + (random() - 0.5_num) * dy

        current => current%next
      ENDDO

      DEALLOCATE(valid_cell_list)
    ENDIF

    ! Remove any unplaced particles from the list. This should never be
    ! called if the above routines worked correctly.
    DO WHILE(ASSOCIATED(current))
      next => current%next
      CALL remove_particle_from_partlist(partlist, current)
      DEALLOCATE(current)
      current => next
    ENDDO

    CALL MPI_ALLREDUCE(partlist%count, npart_this_species, 1, MPI_INTEGER8, &
        MPI_SUM, comm, errcode)

    species%count = npart_this_species

    IF (rank == 0) THEN
      CALL integer_as_string(npart_this_species, string)
      WRITE(*,*) 'Loaded ', TRIM(ADJUSTL(string)), &
          ' particles of species ', '"' // TRIM(species%name) // '"'
#ifndef NO_IO
      WRITE(stat_unit,*) 'Loaded ', TRIM(ADJUSTL(string)), &
          ' particles of species ', '"' // TRIM(species%name) // '"'
#endif
    ENDIF

    CALL particle_bcs

  END SUBROUTINE load_particles



#ifndef PER_SPECIES_WEIGHT
  SUBROUTINE setup_particle_density(density_in, species, density_min, &
      density_max)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(IN) :: density_in
    TYPE(particle_species), POINTER :: species
    REAL(num), INTENT(IN) :: density_min, density_max
    TYPE(particle), POINTER :: current, next
    INTEGER(i8) :: ipart
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: npart_in_cell
#ifdef PARTICLE_SHAPE_TOPHAT
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: rpart_in_cell
#endif
    REAL(num) :: wdata, x0, x1, y0, y1
    TYPE(particle_list), POINTER :: partlist
    INTEGER :: ix, iy, i, j, isubx, isuby
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: density
    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: density_map
#include "particle_head.inc"

    ALLOCATE(density(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(density_map(1-ng:nx+ng,1-ng:ny+ng))
    density = density_in
    density_map = .FALSE.

    CALL field_bc(density, ng)

    DO iy = 1-ng, ny+ng
    DO ix = 1-ng, nx+ng
      IF (density(ix,iy) > density_max) density(ix,iy) = density_max
      IF (density(ix,iy) >= density_min) THEN
        density_map(ix,iy) = .TRUE.
      ELSE
        density(ix,iy) = 0.0_num
      ENDIF
    ENDDO ! ix
    ENDDO ! iy

    ! Uniformly load particles in space
    CALL load_particles(species, density_map)

    ALLOCATE(npart_in_cell(1-ng:nx+ng,1-ng:ny+ng))
    npart_in_cell = 0

    partlist => species%attached_list
    ! If using per particle weighing then use the weight function to match the
    ! uniform pseudoparticle density to the real particle density
    current => partlist%head
    ipart = 0
    DO WHILE(ipart < partlist%count)
      IF (.NOT. ASSOCIATED(current)) PRINT *, 'Bad Particle'

#include "particle_to_grid.inc"

      ! Calculate density at the particle position
      wdata = 0.0_num
      DO isuby = sf_min, sf_max
        i = cell_x
        j = cell_y + isuby
#ifdef PARTICLE_SHAPE_TOPHAT
        IF (.NOT. density_map(i,j)) j = cell_y + 1 - isuby
#else
        IF (.NOT. density_map(i,j)) THEN
          j = cell_y + isuby / 2
#ifdef PARTICLE_SHAPE_BSPLINE3
          IF (.NOT. density_map(i,j)) j = cell_y - isuby / 2
#endif
        ENDIF
#endif
        DO isubx = sf_min, sf_max
          i = cell_x + isubx
#ifdef PARTICLE_SHAPE_TOPHAT
          IF (.NOT. density_map(i,j)) i = cell_x + 1 - isubx
#else
          IF (.NOT. density_map(i,j)) THEN
            i = cell_x + isubx / 2
#ifdef PARTICLE_SHAPE_BSPLINE3
            IF (.NOT. density_map(i,j)) i = cell_x - isubx / 2
#endif
          ENDIF
#endif
          wdata = wdata + gx(isubx) * gy(isuby) * density(i,j)
        ENDDO ! isubx
      ENDDO ! isuby

      current%weight = wdata
      npart_in_cell(cell_x,cell_y) = npart_in_cell(cell_x,cell_y) + 1

      current => current%next
      ipart = ipart + 1
    ENDDO
    DEALLOCATE(density_map)
    DEALLOCATE(density)

    wdata = dx * dy

#ifdef PARTICLE_SHAPE_TOPHAT
    ! For the TOPHAT shape function, particles can be located on a
    ! neighbouring process
    ALLOCATE(rpart_in_cell(1-ng:nx+ng,1-ng:ny+ng))

    rpart_in_cell = npart_in_cell
    CALL processor_summation_bcs(rpart_in_cell, ng)
    npart_in_cell = INT(rpart_in_cell)

    DEALLOCATE(rpart_in_cell)
#endif

    partlist => species%attached_list
    ! Second loop renormalises particle weights
    current => partlist%head
    ipart = 0
    DO WHILE(ipart < partlist%count)
#ifdef PARTICLE_SHAPE_TOPHAT
      cell_x = FLOOR((current%part_pos(1) - x_grid_min_local) / dx) + 1
      cell_y = FLOOR((current%part_pos(2) - y_grid_min_local) / dy) + 1
#else
      cell_x = FLOOR((current%part_pos(1) - x_grid_min_local) / dx + 1.5_num)
      cell_y = FLOOR((current%part_pos(2) - y_grid_min_local) / dy + 1.5_num)
#endif

      current%weight = current%weight * wdata / npart_in_cell(cell_x,cell_y)

      current => current%next
      ipart = ipart + 1
    ENDDO

    DEALLOCATE(npart_in_cell)

    ! If you are filling ghost cells to meet an injector
    ! Then you have overfilled by half a cell but need those particles
    ! To calculate weights correctly. Now delete those particles that
    ! Overlap with the injection region
    IF (species%fill_ghosts) THEN
      x1 = 0.5_num * dx * png
      x0 = x_min - x1
      x1 = x_max + x1

      y1 = 0.5_num * dy * png
      y0 = y_min - y1
      y1 = y_max + y1

      current => partlist%head
      DO WHILE(ASSOCIATED(current))
        next => current%next
        IF (current%part_pos(1) < x0 .OR. current%part_pos(1) >= x1 &
            .OR. current%part_pos(2) < y0 .OR. current%part_pos(2) >= y1) THEN
          CALL remove_particle_from_partlist(partlist, current)
          DEALLOCATE(current)
        ENDIF
        current => next
      ENDDO
    ENDIF

  END SUBROUTINE setup_particle_density
#endif



  FUNCTION sample_dist_function(axis, dist_fn)

    REAL(num), DIMENSION(:), INTENT(IN) :: axis, dist_fn
    REAL(num), DIMENSION(:), ALLOCATABLE :: cdf
    REAL(num) :: position, d_cdf
    INTEGER :: n_points, ipoint, start, endpoint, current
    REAL(num) :: sample_dist_function

    n_points = SIZE(dist_fn)
    ALLOCATE(cdf(n_points))

    cdf(1) = dist_fn(1)
    DO ipoint = 2, n_points
      cdf(ipoint) = cdf(ipoint-1) + dist_fn(ipoint)
    ENDDO

    cdf = cdf / cdf(n_points)

    position = random()
    sample_dist_function = 0.0_num

    start = 1
    endpoint = n_points
    current = (start + endpoint) / 2

    DO current = 1, n_points-1
      IF (cdf(current) <= position .AND. cdf(current+1) >= position) THEN
        d_cdf = cdf(current+1) - cdf(current)
        sample_dist_function = (axis(current) * (position - cdf(current)) &
            + axis(current+1) * (cdf(current+1) - position)) / d_cdf
        EXIT
      ENDIF
    ENDDO

    DEALLOCATE(cdf)

  END FUNCTION sample_dist_function



  SUBROUTINE custom_particle_load

    LOGICAL :: file_inconsistencies
    INTEGER :: current_loader_num
    INTEGER :: part_count, read_count
    CHARACTER(LEN=string_length) :: stra
    REAL(num), DIMENSION(:), POINTER :: xbuf, ybuf
    REAL(num), DIMENSION(:), POINTER :: pxbuf, pybuf, pzbuf
#if !defined(PER_SPECIES_WEIGHT) || defined (PHOTONS)
    REAL(num), DIMENSION(:), POINTER :: wbuf
#endif
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
    INTEGER(KIND=i4), DIMENSION(:), POINTER :: idbuf4
    INTEGER(KIND=i8), DIMENSION(:), POINTER :: idbuf8
    INTEGER :: i, id_offset
    INTEGER, DIMENSION(:), POINTER :: part_counts
#endif
    TYPE(particle_species), POINTER :: species
    TYPE(custom_particle_loader), POINTER :: curr_loader
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: new_particle

    ! Only needed if there are custom loaders to act on
    IF (n_custom_loaders < 1) RETURN

    ! For every custom loader
    DO current_loader_num = 1, n_custom_loaders
      file_inconsistencies = .FALSE.

      curr_loader => custom_loaders_list(current_loader_num)

      ! Grab associated particle lists
      species => species_list(curr_loader%species_id)
      partlist => species%attached_list

      ! Just to be sure
      CALL destroy_partlist(partlist)
      CALL create_empty_partlist(partlist)

      ! MPI read files
      part_count = load_1d_real_array(curr_loader%x_data, xbuf, &
          curr_loader%x_data_offset, errcode)

      read_count = load_1d_real_array(curr_loader%y_data, ybuf, &
          curr_loader%y_data_offset, errcode)
      IF (part_count /= read_count) file_inconsistencies = .TRUE.

#if !defined(PER_SPECIES_WEIGHT) || defined (PHOTONS)
      read_count = load_1d_real_array(curr_loader%w_data, wbuf, &
          curr_loader%w_data_offset, errcode)
      IF (part_count /= read_count) file_inconsistencies = .TRUE.
#endif
      IF (curr_loader%px_data_given) THEN
        read_count = load_1d_real_array(curr_loader%px_data, pxbuf, &
            curr_loader%px_data_offset, errcode)
        IF (part_count /= read_count) file_inconsistencies = .TRUE.
      ENDIF

      IF (curr_loader%py_data_given) THEN
        read_count = load_1d_real_array(curr_loader%py_data, pybuf, &
            curr_loader%py_data_offset, errcode)
        IF (part_count /= read_count) file_inconsistencies = .TRUE.
      ENDIF

      IF (curr_loader%pz_data_given) THEN
        read_count = load_1d_real_array(curr_loader%pz_data, pzbuf, &
            curr_loader%pz_data_offset, errcode)
        IF (part_count /= read_count) file_inconsistencies = .TRUE.
      ENDIF

#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
      IF (curr_loader%id_data_given) THEN
        IF (curr_loader%id_data_4byte) THEN
          read_count = load_1d_integer4_array(curr_loader%id_data, idbuf4, &
              curr_loader%id_data_offset, errcode)
        ELSE
          read_count = load_1d_integer8_array(curr_loader%id_data, idbuf8, &
              curr_loader%id_data_offset, errcode)
        ENDIF
        IF (part_count /= read_count) file_inconsistencies = .TRUE.
      ENDIF
#endif

      CALL MPI_ALLREDUCE(MPI_IN_PLACE, file_inconsistencies, 1, MPI_LOGICAL, &
          MPI_LOR, comm, errcode)

      IF (file_inconsistencies) THEN
        IF (rank == 0) THEN
          WRITE(*,*) '*** ERROR ***'
          WRITE(*,*) 'Error while loading particles_from_file for species ', &
              TRIM(species%name)
        ENDIF
        CALL abort_code(c_err_bad_setup)
      ENDIF

! This is needed to get the IDs assigned properly
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
      IF (.NOT.curr_loader%id_data_given) THEN
        ALLOCATE(part_counts(0:nproc-1))
        CALL MPI_ALLGATHER(part_count, 1, MPI_INTEGER4, part_counts, 1, &
            MPI_INTEGER4, comm, errcode)
        id_offset = 0
        DO i = 0, rank
          id_offset = id_offset + part_counts(i)
        ENDDO
      ENDIF
#endif

      DO read_count = 1, part_count
        CALL create_particle(new_particle)
        CALL add_particle_to_partlist(partlist, new_particle)

        ! Insert data to particle
        new_particle%part_pos(1) = xbuf(read_count)
        new_particle%part_pos(2) = ybuf(read_count)
#if !defined(PER_SPECIES_WEIGHT) || defined (PHOTONS)
        new_particle%weight = wbuf(read_count)
#endif
        IF (curr_loader%px_data_given) THEN
          new_particle%part_p(1) = pxbuf(read_count)
        ENDIF
        IF (curr_loader%py_data_given) THEN
          new_particle%part_p(2) = pybuf(read_count)
        ENDIF
        IF (curr_loader%pz_data_given) THEN
          new_particle%part_p(3) = pzbuf(read_count)
        ENDIF
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
        IF (curr_loader%id_data_given) THEN
#if defined(PARTICLE_ID4)
          IF (curr_loader%id_data_4byte) THEN
            new_particle%id = idbuf4(read_count)
          ELSE
            new_particle%id = INT(idbuf8(read_count), i4)
          ENDIF
#else
          IF (curr_loader%id_data_4byte) THEN
            new_particle%id = INT(idbuf4(read_count), i8)
          ELSE
            new_particle%id = idbuf8(read_count)
          ENDIF
#endif
        ELSE
#if defined(PARTICLE_ID4)
          new_particle%id = INT(id_offset + read_count, i4)
#else
          new_particle%id = INT(id_offset + read_count, i8)
#endif
        ENDIF
#endif
        ! Just being careful
        NULLIFY(new_particle)
      ENDDO

      ! Need to keep totals accurate
      CALL MPI_ALLREDUCE(partlist%count, species%count, 1, MPI_INTEGER8, &
          MPI_SUM, comm, errcode)

      IF (rank == 0) THEN
        CALL integer_as_string(species%count, stra)
        WRITE(*,*) 'Inserted ', TRIM(stra), &
            ' custom particles of species "', TRIM(species%name), '"'
#ifndef NO_IO
        WRITE(stat_unit,*) 'Inserted ', TRIM(stra), &
            ' custom particles of species "', TRIM(species%name), '"'
#endif
      ENDIF
    ENDDO

    DEALLOCATE(custom_loaders_list)

    CALL distribute_particles

  END SUBROUTINE custom_particle_load

END MODULE helper
