PROGRAM bessel

IMPLICIT NONE 
REAL(kind = 16), PARAMETER :: pi = 3.141592653589791638462643383279503
REAL(kind = 16) :: n
REAL(kind = 16) :: z
REAL(kind = 16) :: values
INTEGER :: iz 
values = 0.0
Do iz = 1,2
    z = REAL(iz)*1
    CALL bessel_kn(n,z,values)
WRITE(*,*) values
END DO 

END PROGRAM bessel

SUBROUTINE bessel_kn(n,z,values)
    !Calculate bessel function values = In(z)
    REAL(kind = 16) , INTENT(IN) :: z,n
    REAL(kind = 16), PARAMETER :: pi = 3.141592653589793238462643383279503
    INTEGER :: step,it
    REAL(kind = 16) :: dt,t,f1,f2,f3
    REAL(kind = 16), INTENT(OUT):: values

    step = 1000
    dt = 1.0
    t = 0.0
    values = 0.0

    !Simpson integral formula:
    !I(f) = S(h) = h/6*sigma^n-1_i=0(f(xi)+4*f(xi+1/2*h)+f(xi+h))
    Do it = 0 , step - 1 
       t = t +REAL(it*dt)
       f1 = exp(-z*cosh(t))*cosh(2*t)
       f2 = exp(-z*cosh((t+dt/2)))*cosh(2*(t+dt/2));
       f3 = exp(-z*cosh((t+dt)))*cosh(2*(t+dt));
       values = values + (f1 + 4*f2 + f3)*dt/6.0
       !write(*,*) it,theta, f1 ,f2 ,f3 ,values
    END DO
END SUBROUTINE bessel_kn
