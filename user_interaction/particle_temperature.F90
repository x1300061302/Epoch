! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009      Chris Brady <C.S.Brady@warwick.ac.uk>
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

MODULE particle_temperature

  USE shared_data
  USE random_generator

  IMPLICIT NONE

CONTAINS

  ! Subroutine to initialise a thermal particle distribution
  ! Assumes linear interpolation of temperature between cells
  SUBROUTINE setup_particle_temperature(temperature, direction, part_species, &
      drift)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(IN) :: temperature
    INTEGER, INTENT(IN) :: direction
    TYPE(particle_species), POINTER :: part_species
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(IN) :: drift
    TYPE(particle_list), POINTER :: partlist
    REAL(num) :: mass, temp_local, drift_local,px,py,pz
    TYPE(particle), POINTER :: current
    INTEGER(i8) :: ipart
    INTEGER :: ix, iy
#include "particle_head.inc"

    partlist => part_species%attached_list
    current => partlist%head
    ipart = 0
    DO WHILE(ipart < partlist%count)
#ifdef PER_PARTICLE_CHARGE_MASS
      mass = current%mass
#else
      mass = part_species%mass
#endif

      ! Assume that temperature is cell centred
#include "particle_to_grid.inc"

      temp_local = 0.0_num
      drift_local = 0.0_num
      DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          temp_local = temp_local &
              + gx(ix) * gy(iy) * temperature(cell_x+ix, cell_y+iy)
          drift_local = drift_local &
              + gx(ix) * gy(iy) * drift(cell_x+ix, cell_y+iy)
        ENDDO
      ENDDO

      !Not Use Juttner Distribution
    IF (.NOT. use_juttner) THEN
        IF (direction == c_dir_x) current%part_p(1) = &
            momentum_from_temperature(mass, temp_local, drift_local)

        IF (direction == c_dir_y) current%part_p(2) = &
            momentum_from_temperature(mass, temp_local, drift_local)

        IF (direction == c_dir_z) current%part_p(3) = &
            momentum_from_temperature(mass, temp_local, drift_local)
       ELSE
      !Use Juttner distribution 
!      IF (direction == c_dir_x) THEN 
!          call momentum_from_temperature(mass, temp_local, drift_local)
!          !call rel_momentum_from_temperature_x(mass, temp_local, drift_local,px,py,pz)
!          current%part_p(1) = px
!          current%part_p(2) = py
!          current%part_p(3) = pz
!      END IF
!
!      IF (direction == c_dir_y) THEN
!          call momentum_from_temperature(mass, temp_local, drift_local)
!          !call rel_momentum_from_temperature_y(mass, temp_local, drift_local, px, py, pz)
!          current%part_p(1) = px
!          current%part_p(2) = py
!          current%part_p(3) = pz
!      END IF
!!     Xiey Only drift in z-direction
      !test write(*,*) temp_local
      IF (direction == c_dir_z) THEN 
          !call momentum_from_temperature(mass, temp_local, drift_local)
          call rel_momentum_from_temperature_z(mass, temp_local, drift_local,px,py,pz)
          current%part_p(1) = px
          current%part_p(2) = py
          current%part_p(3) = pz
      END IF ! direction == c_dir_z
  ENDIF!use_juttner
      current => current%next
      ipart = ipart + 1
    ENDDO

  END SUBROUTINE setup_particle_temperature



  FUNCTION momentum_from_temperature(mass, temperature, drift)

    REAL(num), INTENT(IN) :: mass, temperature, drift
    REAL(num) :: momentum_from_temperature
    DOUBLE PRECISION :: stdev, mu

    stdev = DBLE(SQRT(temperature * kb * mass))
    mu = DBLE(drift)
    momentum_from_temperature = random_box_muller(stdev, mu)

  END FUNCTION momentum_from_temperature



  ! Function for generating momenta of thermal particles in a particular
  ! direction, e.g. the +x direction.
  ! These satisfy a Rayleigh distribution, formed by combining two
  ! normally-distributed (~N(0,sigma)) random variables as follows:
  ! Z = SQRT(X**2 + Y**2)
  FUNCTION flux_momentum_from_temperature(mass, temperature, drift)

    REAL(num), INTENT(IN) :: mass, temperature, drift
    REAL(num) :: flux_momentum_from_temperature
    REAL(num) :: mom1, mom2

    mom1 = momentum_from_temperature(mass, temperature, 0.0_num)
    mom2 = momentum_from_temperature(mass, temperature, 0.0_num)

    flux_momentum_from_temperature = SQRT(mom1**2 + mom2**2) + drift

  END FUNCTION flux_momentum_from_temperature

  ! Add By Xu.Z and W P. Yao 
    SUBROUTINE rel_momentum_from_temperature_z(mass, temperature, driftz, px, py, pz)
        REAL(num), INTENT(IN)  :: mass, temperature, driftz
        REAL(num), INTENT(OUT) :: px, py, pz
       
        REAL(num) :: mu,tp, ur, eta, gammaz, betaz, ga
        REAL(num) :: x1, x2, x3, x4, x5, x6, x7
        mu =(mass * c * c)/kb/temperature     ! normalized temperature
        tp = 1.0/mu
        gammaz = sqrt(1+(driftz/mass/c)**2.0_num)
        betaz =  sqrt(1-1/gammaz**2.0_num)
        ! use Gamma distribution as comparision function for the 
        ! rejection method
        IF (tp < 1)THEN
            px = momentum_from_temperature(mass,temperature,0.0_num)
            py = momentum_from_temperature(mass,temperature,0.0_num)
            pz = momentum_from_temperature(mass,temperature,driftz)
            RETURN
        ENDIF
        DO 
            x1 = random()
            x2 = random()
            x3 = random()
            x4 = random()
            IF(abs(x1*x2*x3*x4) <= c_tiny) cycle
            ur = -tp * log(x1*x2*x3)
            eta = ur - tp * log(x4)
            IF( eta*eta - ur *ur > 1) EXIT
        END DO
        x5 = random()
        x6 = random()
        x7 = random()
        px = 2.0_num * ur * sqrt(x5*(1-x5))*cos(2*pi*x6)
        py = 2.0_num * ur * sqrt(x5*(1-x5))*sin(2*pi*x6)
        ! Xiey
        pz = ur * (2*x5-1)
        ga = sqrt(1+ur*ur)
        if(-betaz*pz/ga > x7) pz = -pz
        pz = gammaz * (pz + betaz*sqrt(1+ur*ur))
        px = px * mass * c
        py = py * mass * c
        pz = pz * mass * c
    END SUBROUTINE rel_momentum_from_temperature_z

END MODULE particle_temperature
