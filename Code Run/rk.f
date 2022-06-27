c	Runge Kutta for set of first order differential equations
c
 	PROGRAM oscillator
 	IMPLICIT none
c
c 		declarations
c N:number of equations, nsteps:number of steps, tstep:length of steps
c 		y(1): initial position, y(2):initial velocity
	REAL*8 t, tstep, y(5)
	INTEGER	N, i, j, nsteps
	N=2
        nsteps=300
        tstep=0.1
        y(1)=1.0
        y(2)=0.0
c
c          	open file
        OPEN(6, FILE='rungef.dat')
        WRITE (6,*) 0, y(1)
c
c 		do loop n steps of Runga-Kutta algorithm
	DO 60 j = 1, nsteps
	t=j*tstep
   	call rk4(t, y, tstep, N)
        WRITE (6,*) t, y(1)
 60	CONTINUE
c
	CLOSE(6)
	STOP
        END
c------------------------end of main program------------------------
c
c 		fourth-order Runge-Kutta subroutine 
	SUBROUTINE rk4(t, y, tstep, N)
	IMPLICIT none
c
c		declarations
	REAL*8 DERIV, h, t, tstep, y(5) 
        REAL*8 k1(5), k2(5),k3(5), k4(5), temp1(5), temp2(5), temp3(5)
     	INTEGER i, N
     	h=tstep/2.0
c
	DO 10 i = 1,N
		k1(i) = tstep * DERIV(t, y, i)
		temp1(i) = y(i) + 0.5*k1(i)
 10	CONTINUE
c
	DO 20 i = 1,N
		k2(i) = tstep * DERIV(t+h, temp1, i)
		temp2(i) = y(i) + 0.5*k2(i)
 20	CONTINUE
c
	DO 30 i = 1,N
		k3(i) = tstep * DERIV(t+h, temp2, i)
		temp3(i) = y(i) + k3(i)
 30	CONTINUE
c
	DO 40 i = 1,N
		k4(i) = tstep * DERIV(t+tstep, temp3, i)
		y(i) = y(i) + (k1(i) + (2.*(k2(i) + k3(i))) + k4(i))/6.0
 40	CONTINUE
c 
	RETURN
	END
c
c		function which returns the derivatives
	FUNCTION DERIV(t, temp, i)
	IMPLICIT none
c
c		declarations
	REAL*8 DERIV, t, temp(5)
	INTEGER i
c
	IF (i .EQ. 1) DERIV=temp(2)
	IF (i .EQ. 2) DERIV=-temp(1)
	RETURN
	END
