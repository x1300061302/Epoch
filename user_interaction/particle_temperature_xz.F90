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

        REAL(num), DIMENSION(-2:,-2:), INTENT(IN) :: temperature
        INTEGER, INTENT(IN) :: direction
        TYPE(particle_species), POINTER :: part_species
        REAL(num), DIMENSION(-2:,-2:), INTENT(IN) :: drift
        TYPE(particle_list), POINTER :: partlist
        REAL(num) :: mass, temp_local, drift_local
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

        IF (direction == c_dir_x) current%part_p(1) = &
            momentum_from_temperature(mass, temp_local, drift_local)

        IF (direction == c_dir_y) current%part_p(2) = &
            momentum_from_temperature(mass, temp_local, drift_local)

        IF (direction == c_dir_z) current%part_p(3) = &
            momentum_from_temperature(mass, temp_local, drift_local)

        current => current%next
        ipart = ipart + 1
        ENDDO

    END SUBROUTINE setup_particle_temperature



    FUNCTION momentum_from_temperature(mass, temperature, drift)

        REAL(num), INTENT(IN) :: mass, temperature, drift
        REAL(num) :: momentum_from_temperature

        REAL(num) :: stdev
        REAL(num) :: rand1, rand2, w
        REAL(num), SAVE :: val
        LOGICAL, SAVE :: cached = .FALSE.

        ! This is a basic polar Box-Muller transform
        ! It generates gaussian distributed random numbers
        ! The standard deviation (stdev) is related to temperature

        stdev = SQRT(temperature * kb * mass)

        IF (cached) THEN
            cached = .FALSE.
            momentum_from_temperature = val * stdev + drift
        ELSE
            cached = .TRUE.

            DO
            rand1 = random()
            rand2 = random()

            rand1 = 2.0_num * rand1 - 1.0_num
            rand2 = 2.0_num * rand2 - 1.0_num

            w = rand1**2 + rand2**2

            IF (w > c_tiny .AND. w < 1.0_num) EXIT
            ENDDO

            w = SQRT((-2.0_num * LOG(w)) / w)

            momentum_from_temperature = rand1 * w * stdev + drift
            val = rand2 * w
        ENDIF

    END FUNCTION momentum_from_temperature

    SUBROUTINE rel_momentum_from_temperature_x(mass, temperature, driftx, px, py, pz)
        REAL(num), INTENT(IN)  :: mass, temperature, driftx
        REAL(num), INTENT(OUT) :: px, py, pz
       
        REAL(num) :: mu,tp, ur, eta, gammax, betax, ga
        REAL(num) :: x1, x2, x3, x4, x5, x6, x7
        mu = mass * c*c /(kb*temperature)     ! normalized temperature
        tp = 1.0/mu
        gammax = sqrt(1+(driftx/mass/c)**2.0_num)
        betax =  sqrt(1-1/gammax**2.0_num)
        ! use Gamma distribution as comparision function for the 
        ! rejection method ( more efficient than our comparision function
        ! (1/(1+(x-x0)^2/a0^2)
        
        !Avoiding ineffient sampling when temperature is non-relativistic.
        !Add by HH 2017.09.30
        IF (tp < 0.1)THEN
            px = momentum_from_temperature(mass,temperature,driftx)
            py = momentum_from_temperature(mass,temperature,0.0_num)
            pz = momentum_from_temperature(mass,temperature,0.0_num)
            RETURN
        ENDIF
        DO 
            x1 = random()
            x2 = random()
            x3 = random()
            x4 = random()
            IF((x1*x3*x3) <= c_tiny .OR. x4 <= c_tiny) cycle
            ur = -tp * log(x1*x3*x3)
            eta = ur - tp * log(x4)
            IF( eta*eta - ur *ur > 1) EXIT
        END DO
        x5 = random()
        x6 = random()
        x7 = random()
        px = 2.0_num * ur * (2*x5-1)
        py = 2.0_num * ur * sqrt(x5*(1-x5))*cos(2*pi*x6)
        pz = 2.0_num * ur * sqrt(x5*(1-x5))*sin(2*pi*x6)
        ! flipping method 
        ! see Zenitani, Phys. Plasma, 22, 042116 (2015)
        ga = sqrt(1+ur*ur)
        if(-betax*px/ga > x7) px = -px
        px = gammax * (px + betax*sqrt(1+ur*ur))
        px = px * mass * c
        py = py * mass * c
        pz = pz * mass * c
    END SUBROUTINE rel_momentum_from_temperature_x

    SUBROUTINE rel_momentum_from_temperature_y(mass, temperature, drifty, px, py, pz)
        REAL(num), INTENT(IN)  :: mass, temperature, drifty
        REAL(num), INTENT(OUT) :: px, py, pz
       
        REAL(num) :: mu,tp, ur, eta, gammay, betay, ga
        REAL(num) :: x1, x2, x3, x4, x5, x6, x7
        mu = mass * c /(kb*temperature)     ! normalized temperature
        tp = 1.0/mu
        gammay = sqrt(1+(drifty/mass/c)**2.0_num)
        betay =  sqrt(1-1/gammay**2.0_num)
        ! use Gamma distribution as comparision function for the 
        ! rejection method 
        DO 
            x1 = random()
            x2 = random()
            x3 = random()
            x4 = random()
            IF(abs(x1*x3*x3*x4) <= c_tiny) cycle
            ur = -tp * log(x1*x3*x3)
            eta = ur - tp * log(x4)
            IF( eta*eta - ur *ur > 1) EXIT
        END DO
        x5 = random()
        x6 = random()
        x7 = random()
        px = 2.0_num * ur * sqrt(x5*(1-x5))*sin(2*pi*x6)
        py = 2.0_num * ur * (2*x5-1)
        pz = 2.0_num * ur * sqrt(x5*(1-x5))*cos(2*pi*x6)
        ga = sqrt(1+ur*ur)
        if(-betay*py/ga > x7) py = -py
        py = gammay * (py + betay*sqrt(1+ur*ur))
        px = px * mass * c
        py = py * mass * c
        pz = pz * mass * c
    END SUBROUTINE rel_momentum_from_temperature_y
    
    SUBROUTINE rel_momentum_from_temperature_z(mass, temperature, driftz, px, py, pz)
        REAL(num), INTENT(IN)  :: mass, temperature, driftz
        REAL(num), INTENT(OUT) :: px, py, pz
       
        REAL(num) :: mu,tp, ur, eta, gammaz, betaz, ga
        REAL(num) :: x1, x2, x3, x4, x5, x6, x7
        mu = mass * c /(kb*temperature)     ! normalized temperature
        tp = 1.0/mu
        gammaz = sqrt(1+(driftz/mass/c)**2.0_num)
        betaz =  sqrt(1-1/gammaz**2.0_num)
        ! use Gamma distribution as comparision function for the 
        ! rejection method
        DO 
            x1 = random()
            x2 = random()
            x3 = random()
            x4 = random()
            IF(abs(x1*x3*x3*x4) <= c_tiny) cycle
            ur = -tp * log(x1*x3*x3)
            eta = ur - tp * log(x4)
            IF( eta*eta - ur *ur > 1) EXIT
        END DO
        x5 = random()
        x6 = random()
        x7 = random()
        px = 2.0_num * ur * sqrt(x5*(1-x5))*cos(2*pi*x6)
        py = 2.0_num * ur * sqrt(x5*(1-x5))*sin(2*pi*x6)
        pz = 2.0_num * ur * (2*x5-1)
        ga = sqrt(1+ur*ur)
        if(-betaz*pz/ga > x7) pz = -pz
        pz = gammaz * (pz + betaz*sqrt(1+ur*ur))
        px = px * mass * c
        py = py * mass * c
        pz = pz * mass * c
    END SUBROUTINE rel_momentum_from_temperature_z
END MODULE particle_temperature
