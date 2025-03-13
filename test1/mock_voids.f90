module common
    implicit none
    !cosmology
    real, parameter :: omegaM = 0.31
    real, parameter :: rhoB = 3.9*1e10 !background density in Msun/Mpc^3
    real, parameter :: L0 = 147.5 !box side in Mpc
    real, parameter :: zeta = 0.00
    real, parameter :: total_mass = rhoB*L0**3
    integer, parameter :: npart = 1000000 ! number of particles to throw into the box
    real, parameter :: masspart = total_mass / npart
    !voids
    real, parameter :: deltae = -0.9 !mean density contrast inside void
    real, parameter :: delta_shell = 0.5 !mean density contrast inside shell
    real, parameter :: fshell = 0.1 !fraction of Re to consider the shell
    real, parameter :: pi = acos(-1.0)
    real, parameter :: max_rad = 50.
    real, parameter :: min_rad = 12.
    real, parameter :: v0 = 100. !km/s
    integer, parameter :: nvoid = 50 !number of voids
    integer, parameter :: seed = 14
    integer, parameter :: ellipsoids = 1
end module


!/////////////////////////////////
program montecarlo_voids
!/////////////////////////////////
use common
implicit none
real :: num_rand, theta_rand, phi_rand, bass_rand
real :: rx,ry,rz,mod,X,rxprime,ryprime,rzprime
real :: xc(nvoid), yc(nvoid), zc(nvoid), Re(nvoid), a(nvoid), b(nvoid), c(nvoid)
real :: uax(nvoid), uay(nvoid), uaz(nvoid), ubx(nvoid), uby(nvoid), ubz(nvoid), ucx(nvoid), ucy(nvoid), ucz(nvoid)
integer :: npvoids(nvoid), npshells(nvoid)
integer :: npvoidstot, nfield, npshellstot
integer :: i, j, low
integer :: overlap,outside
real :: xpart(npart), ypart(npart), zpart(npart), vx(npart), vy(npart), vz(npart), mass(npart)
integer :: idvoid(npart), oripa(npart)
real :: buffer1 = 0.


!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MOCK DESIGNING PART
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Here we design "nvoid" spherical voids that
! cannot overlap with each other. We use the mass profile
! provided by Colberg (2005)
!!!!!!!!!!!!!!!!!!!!!

write(*,*) '//////////////////////////////'
write(*,*) '/// Mock voids for testing ///'
write(*,*) '//////////////////////////////'
write(*,*)
write(*,*)

!file to save mock voids
open(unit=1, file='mocks', form = 'formatted')
write(1,*) nvoid

!random generate nvoids with center between -L0/2 and L0/2
!with different axis a,b,c and centers 

write(*,*)
write(*,*) 'Placing voids in the box'
call srand(seed)
num_rand = 0.
do i=1,nvoid
    !generate random radius and center
    call exp_rand(2, num_rand)
    Re(i) = max_rad*num_rand
    do while(Re(i) < min_rad)
        call exp_rand(2, num_rand)
        Re(i) = max_rad*num_rand
    enddo

    xc(i) = L0*(rand() - 0.5)
    yc(i) = L0*(rand() - 0.5)
    zc(i) = L0*(rand() - 0.5)
    
    !ensure no overlap!!!!!!!!!!!!!!!!
    overlap = 0
    do j=1,i
        if (i==j) cycle
        if (sqrt((xc(i)-xc(j))**2 + (yc(i)-yc(j))**2 + (zc(i)-zc(j))**2) < (Re(i) + Re(j))*(1.+fshell) ) then
            overlap = 1
            exit
        endif
    enddo

    do while (overlap==1)

        call exp_rand(2, num_rand)
        Re(i) = max_rad*num_rand
        do while(Re(i) < min_rad)
            call exp_rand(2, num_rand)
            Re(i) = max_rad*num_rand
        enddo
        
        xc(i) = L0*(rand() - 0.5)
        yc(i) = L0*(rand() - 0.5)
        zc(i) = L0*(rand() - 0.5)

        overlap = 0
        do j=1,i
            if (i==j) cycle
            if (sqrt((xc(i)-xc(j))**2 + (yc(i)-yc(j))**2 + (zc(i)-zc(j))**2) < (Re(i) + Re(j))*(1.+fshell) ) then
                overlap = 1
                exit
            endif
        enddo
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !ensure inside box
    outside = 0
    if (xc(i) + Re(i)*(1.+fshell) > L0/2 .or. xc(i) - Re(i)*(1.+fshell) < -L0/2) outside = 1
    if (yc(i) + Re(i)*(1.+fshell) > L0/2 .or. yc(i) - Re(i)*(1.+fshell) < -L0/2) outside = 1
    if (zc(i) + Re(i)*(1.+fshell) > L0/2 .or. zc(i) - Re(i)*(1.+fshell) < -L0/2) outside = 1

    do while (outside == 1)

        !generate random radius and center
        call exp_rand(2, num_rand)
        Re(i) = max_rad*num_rand
        do while(Re(i) < min_rad)
            call exp_rand(2, num_rand)        
            if (abs(xc(i)) > L0/2) outside = 1
            if (abs(yc(i)) > L0/2) outside = 1
            if (abs(zc(i)) > L0/2) outside = 1
            Re(i) = max_rad*num_rand
        enddo
    
        xc(i) = L0*(rand() - 0.5)
        yc(i) = L0*(rand() - 0.5)
        zc(i) = L0*(rand() - 0.5)
        
        !ensure no overlap!!!!!!!!!!!!!!!!
        overlap = 0
        do j=1,i
            if (i==j) cycle
            if (sqrt((xc(i)-xc(j))**2 + (yc(i)-yc(j))**2 + (zc(i)-zc(j))**2) < (Re(i) + Re(j))*(1.+fshell)) then
                overlap = 1
                exit
            endif
        enddo
    
        do while (overlap==1)
            call exp_rand(2, num_rand)
            Re(i) = max_rad*num_rand
            do while(Re(i) < min_rad)
                call exp_rand(2, num_rand)
                Re(i) = max_rad*num_rand
            enddo
            xc(i) = L0*(rand() - 0.5)
            yc(i) = L0*(rand() - 0.5)
            zc(i) = L0*(rand() - 0.5)
    
            overlap = 0
            do j=1,i
                if (i==j) cycle
                if (sqrt((xc(i)-xc(j))**2 + (yc(i)-yc(j))**2 + (zc(i)-zc(j))**2) < (Re(i) + Re(j))*(1.+fshell)) then
                    overlap = 1
                    exit
                endif
            enddo
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        outside = 0

        if (xc(i) + Re(i)*(1.+fshell) > L0/2 .or. xc(i) - Re(i)*(1.+fshell) < -L0/2) outside = 1
        if (yc(i) + Re(i)*(1.+fshell) > L0/2 .or. yc(i) - Re(i)*(1.+fshell) < -L0/2) outside = 1
        if (zc(i) + Re(i)*(1.+fshell) > L0/2 .or. zc(i) - Re(i)*(1.+fshell) < -L0/2) outside = 1
    enddo

enddo
write(*,*)'done..'

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! METROPOLIS MOTECARLO PART
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now we plance "npart" particles in the box with some
! mass using the montecarlo method
! we can also add noise to test the void finder accuracy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(*,*)
write(*,*) 'Placing particles inside voids'

!divide between particles inside voids and outside
npvoidstot = 0
nfield = 0
idvoid = 0

!Distribute particles according to void sizes
call part_distribution(Re, npvoidstot, npvoids, npshellstot, npshells)

!place particles in voids according to the standard profile
!given in Colberg 2005
!velocity according to ..
low = 0
do i=1,nvoid
    do j=low+1,low+npvoids(i)

        !random Colberg position
        ! call dens_rand(num_rand)
        !random Ricciardelli position 
        call dens_rand_ricciardelli(num_rand)
        num_rand = Re(i)*num_rand !radius
        phi_rand = 2*pi*rand()
        theta_rand = pi*rand()
        bass_rand = rand()
        do while(abs(sin(theta_rand)) < bass_rand)
            theta_rand = pi*rand()
            bass_rand = rand()
        enddo
        xpart(j) = xc(i) + num_rand*sin(theta_rand)*cos(phi_rand)
        ypart(j) = yc(i) + num_rand*sin(theta_rand)*sin(phi_rand)
        zpart(j) = zc(i) + num_rand*cos(theta_rand)

        mass(j) = masspart
        oripa(j) = j
        idvoid(j) = i

        !velocity
        !unitary radial vector
        rx = xpart(j) - xc(i)
        ry = ypart(j) - yc(i)
        rz = zpart(j) - zc(i)
        mod = (rx**2 + ry**2 + rz**2)**0.5
        rx = rx/mod
        ry = ry/mod
        rz = rz/mod

        X=mod/Re(i)
        vx(j) = 10*12./15.*v0*X*(2 - 3./4.*X )*rx
        vy(j) = 10*12./15.*v0*X*(2 - 3./4.*X )*ry
        vz(j) = 10*12./15.*v0*X*(2 - 3./4.*X )*rz

        ! !randomize as we get close to the border
        ! vx(j) = vx(j) + 0.5*v0*(rand()-0.5)*(mod/Re(i))
        ! vy(j) = vy(j) + 0.5*v0*(rand()-0.5)*(mod/Re(i))
        ! vz(j) = vz(j) + 0.5*v0*(rand()-0.5)*(mod/Re(i))

    enddo
    low = low + npvoids(i)
enddo

write(*,*) 'done..'
write(*,*)
write(*,*) 'Placing particles in shells'
!place particles in shells, random uniform
low = npvoidstot
do i=1,nvoid
    do j=low+1,low+npshells(i)

        !random radial position
        num_rand = rand()
        num_rand = Re(i)*fshell*num_rand + Re(i) ! rand radius between Re and Re*(1+fshell)
        phi_rand = 2*pi*rand()
        theta_rand = pi*rand()
        bass_rand = rand()
        do while(abs(sin(theta_rand)) < bass_rand)
            theta_rand = pi*rand()
            bass_rand = rand()
        enddo
        xpart(j) = xc(i) + num_rand*sin(theta_rand)*cos(phi_rand)
        ypart(j) = yc(i) + num_rand*sin(theta_rand)*sin(phi_rand)
        zpart(j) = zc(i) + num_rand*cos(theta_rand)

        mass(j) = masspart
        oripa(j) = j
        idvoid(j) = i

        !velocity
        vx(i) = v0*(rand()-0.5)
        vy(i) = v0*(rand()-0.5)
        vz(i) = v0*(rand()-0.5)

        ! !randomize as we get close to the border
        ! vx(j) = vx(j) + 0.5*v0*(rand()-0.5)*(mod/Re(i))
        ! vy(j) = vy(j) + 0.5*v0*(rand()-0.5)*(mod/Re(i))
        ! vz(j) = vz(j) + 0.5*v0*(rand()-0.5)*(mod/Re(i))

    enddo
    low = low + npshells(i)
enddo

!now rescale positions to ellipsoid with a=Re, obtain b and c randomly
if (ellipsoids .eq. 1) then
    !restart the seed, otherwise, different the number of different calls to rand()
    !will give different results
    call srand(seed)
    !voids
    low = 0
    do i=1,nvoid
        a(i) = Re(i)
        b(i) = 0.25*a(i)*rand() + 0.75*a(i)
        c(i) = (b(i) - 0.6*a(i))*rand() + 0.6*a(i)
        write(1,*) i,xc(i),yc(i),zc(i),a(i),b(i),c(i)
        !random a unitary vector, random sign
        uax(i) = rand()*(2*rand()-1.)
        uay(i) = rand()*(2*rand()-1.)
        uaz(i) = rand()*(2*rand()-1.)
        !random b unitary vector, perpendicular to a
        ubx(i) = rand()*(2*rand()-1.)
        uby(i) = rand()*(2*rand()-1.)
        ubz(i) = -(uax(i)*ubx(i) + uay(i)*uby(i))/uaz(i)
        !random c unitary vector, perpendicular to a and b
        ucx(i) = rand()*(2*rand()-1.)
        ucz(i) = (uby(i)*uax(i)*ucx(i)/uay(i) - ubx(i)*ucx(i))/(ubz(i) - uby(i)*uaz(i)/uay(i))
        ucy(i) = -(uax(i)*ucx(i) + uaz(i)*ucz(i))/uay(i) 
        !normalize
        mod = (uax(i)**2 + uay(i)**2 + uaz(i)**2)**0.5
        uax(i) = uax(i)/mod
        uay(i) = uay(i)/mod
        uaz(i) = uaz(i)/mod
        mod = (ubx(i)**2 + uby(i)**2 + ubz(i)**2)**0.5
        ubx(i) = ubx(i)/mod
        uby(i) = uby(i)/mod
        ubz(i) = ubz(i)/mod
        mod = (ucx(i)**2 + ucy(i)**2 + ucz(i)**2)**0.5
        ucx(i) = ucx(i)/mod
        ucy(i) = ucy(i)/mod
        ucz(i) = ucz(i)/mod
        !rescale to ellipsoid
        do j=low+1,low+npvoids(i)
            rx = xpart(j) - xc(i)
            ry = ypart(j) - yc(i)
            rz = zpart(j) - zc(i)
            rxprime = (uax(i)*rx + uay(i)*ry + uaz(i)*rz)*(a(i)/Re(i))
            ryprime = (ubx(i)*rx + uby(i)*ry + ubz(i)*rz)*(b(i)/Re(i))
            rzprime = (ucx(i)*rx + ucy(i)*ry + ucz(i)*rz)*(c(i)/Re(i))
            xpart(j) = xc(i) + uax(i)*rxprime + ubx(i)*ryprime + ucx(i)*rzprime
            ypart(j) = yc(i) + uay(i)*rxprime + uby(i)*ryprime + ucy(i)*rzprime
            zpart(j) = zc(i) + uaz(i)*rxprime + ubz(i)*ryprime + ucz(i)*rzprime
        enddo
    low = low + npvoids(i)
    enddo
    !shells
    low = npvoidstot
    do i=1,nvoid
        !rescale to ellipsoid
        do j=low+1,low+npshells(i)
            rx = xpart(j) - xc(i)
            ry = ypart(j) - yc(i)
            rz = zpart(j) - zc(i)
            rxprime = (uax(i)*rx + uay(i)*ry + uaz(i)*rz)*(a(i)/Re(i))
            ryprime = (ubx(i)*rx + uby(i)*ry + ubz(i)*rz)*(b(i)/Re(i))
            rzprime = (ucx(i)*rx + ucy(i)*ry + ucz(i)*rz)*(c(i)/Re(i))
            xpart(j) = xc(i) + uax(i)*rxprime + ubx(i)*ryprime + ucx(i)*rzprime
            ypart(j) = yc(i) + uay(i)*rxprime + uby(i)*ryprime + ucy(i)*rzprime
            zpart(j) = zc(i) + uaz(i)*rxprime + ubz(i)*ryprime + ucz(i)*rzprime
        enddo
    low = low + npshells(i)
    enddo


    !all field particles, completely random position out of voids
    !and random velocity field
    write(*,*)
    write(*,*) 'Placing particles outside voids'
    do i=npvoidstot+npshellstot+1,npart
        xpart(i) = L0*(rand() - 0.5)
        ypart(i) = L0*(rand() - 0.5)
        zpart(i) = L0*(rand() - 0.5)
        overlap = 0
        !check particle is not inside ellipsoids
        do j=1,nvoid
            if ( ( (xpart(i)-xc(j))*uax(j) + (ypart(i)-yc(j))*uay(j) + (zpart(i)-zc(j))*uaz(j) )**2/a(j)**2 + &
                 ( (xpart(i)-xc(j))*ubx(j) + (ypart(i)-yc(j))*uby(j) + (zpart(i)-zc(j))*ubz(j) )**2/b(j)**2 + &
                 ( (xpart(i)-xc(j))*ucx(j) + (ypart(i)-yc(j))*ucy(j) + (zpart(i)-zc(j))*ucz(j) )**2/c(j)**2 < 1 ) then
                overlap = 1
                exit
            endif
        enddo

        do while(overlap==1)
            xpart(i) = L0*(rand() - 0.5)
            ypart(i) = L0*(rand() - 0.5)
            zpart(i) = L0*(rand() - 0.5)
            overlap = 0
            do j=1,nvoid
                if ( ( (xpart(i)-xc(j))*uax(j) + (ypart(i)-yc(j))*uay(j) + (zpart(i)-zc(j))*uaz(j) )**2/a(j)**2 + &
                     ( (xpart(i)-xc(j))*ubx(j) + (ypart(i)-yc(j))*uby(j) + (zpart(i)-zc(j))*ubz(j) )**2/b(j)**2 + &
                     ( (xpart(i)-xc(j))*ucx(j) + (ypart(i)-yc(j))*ucy(j) + (zpart(i)-zc(j))*ucz(j) )**2/c(j)**2 < 1 ) then
                    overlap = 1
                    exit
                endif
            enddo
        enddo

        mass(i) = masspart
        oripa(i) = i
        vx(i) = v0*(rand()-0.5)
        vy(i) = v0*(rand()-0.5)
        vz(i) = v0*(rand()-0.5)
    enddo
    write(*,*) 'done..'

else
    do i=1,nvoid
        write(1,*) i,xc(i),yc(i),zc(i),Re(i),Re(i),Re(i)
    enddo

    !all field particles, completely random position out of voids
    !and random velocity field
    write(*,*)
    write(*,*) 'Placing particles outside voids'
    do i=npvoidstot+npshellstot+1,npart
        xpart(i) = L0*(rand() - 0.5)
        ypart(i) = L0*(rand() - 0.5)
        zpart(i) = L0*(rand() - 0.5)
        overlap = 0
        do j=1,nvoid
            if (sqrt((xpart(i)-xc(j))**2 + (ypart(i)-yc(j))**2 + (zpart(i)-zc(j))**2) < Re(j)*(1.+fshell)) then
                overlap = 1
                exit
            endif
        enddo

        do while(overlap==1)
            xpart(i) = L0*(rand() - 0.5)
            ypart(i) = L0*(rand() - 0.5)
            zpart(i) = L0*(rand() - 0.5)
            overlap = 0
            do j=1,nvoid
                if (sqrt((xpart(i)-xc(j))**2 + (ypart(i)-yc(j))**2 + (zpart(i)-zc(j))**2) < Re(j)*(1.+fshell)) then
                    overlap = 1
                    exit
                endif
            enddo
        enddo

        mass(i) = masspart
        oripa(i) = i
        vx(i) = v0*(rand()-0.5)
        vy(i) = v0*(rand()-0.5)
        vz(i) = v0*(rand()-0.5)
    enddo
    write(*,*) 'done..'

endif

close(1)



!file to save mock voids
write(*,*)
write(*,*) 'writing output'
open(unit=2, file='bin_file_part00001', form='unformatted')
write(2) int(npart, kind=8), zeta
write(2) (xpart(i), i=1,npart)
write(2) (ypart(i), i=1,npart)
write(2) (zpart(i), i=1,npart)
write(2) (vx(i), i=1,npart)
write(2) (vy(i), i=1,npart)
write(2) (vz(i), i=1,npart)
write(2) (mass(i), i=1,npart)
! write(2) (oripa(i), i=1,npart)
! write(2) (idvoid(i), i=1,npart)
close(2)

write(*,*)
write(*,*) '//////////////////////////////'
write(*,*) '/// Done ///'
write(*,*) '//////////////////////////////'

!/////////////////////////////////
end program
!/////////////////////////////////


!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine exp_rand(a, num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! produce a random number according to P(x) = exp(-a*x)
! num in [0,1]
implicit none
!in
integer :: a
!out
real :: num, num2, num3

!const = a / (1-exp(real(-a)))
num = rand()
num3 = exp(-a*num) !between 0 and 1
num2 = rand()     !between 0 and 1
do while (num2 > num3)
    num = rand()
    num3 = exp(-a*num)
    num2 = rand()
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dens_rand(num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! produce a random number according to Colberg density profile
! REJECTION METHOD
! num in [0,1]
implicit none
!in
real :: alpha = 1.85
!out
real :: num, num2, num3

!const = alpha*(1 - 1/exp(1.))
num = rand()
num3 = num**2*(1+(alpha/3)*num**alpha)*exp(num**alpha - 1)/(1+alpha/3) !between 0 and 1
num2 = rand()             !between 0 and 1
do while (num2 > num3)
    num = rand()
    num3 = num**2*(1+(alpha/3)*num**alpha)*exp(num**alpha - 1)/(1+alpha/3)
    num2 = rand()
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dens_rand_ITF(x)
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! produce a random number according to Colberg density profile
! INVERSE TRANSFORM METHOD
! num in [0,1]
implicit none
!in
real :: alpha = 1.85
!out
real :: x
!locals
real :: u,CDF,CDF_der
integer :: nit = 40, it

!u such that x=CDF-1(u)
u = rand()

!x init
x = 0.5
CDF = x**3*exp(x**alpha - 1)
CDF_der = x**2*(3 + alpha*x**alpha)*exp(x**alpha - 1)
!Apply newton raphson to find x
do it=1,nit
    x = x - (CDF - u) / CDF_der
    CDF = x**3*exp(x**alpha - 1)
    CDF_der = x**2*(3 + alpha*x**alpha)*exp(x**alpha - 1)
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dens_rand_ricciardelli(num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! produce a random number according to Ricciardelli et al 2014 density profile
! REJECTION METHOD
! num in [0,1]
implicit none
!in
real :: alpha2 = 0.08
real :: beta2 = 1.29
!out
real :: num, num2, num3

!const = alpha*(1 - 1/exp(1.))
num = rand()
num3 = num**(2+alpha2)*exp(num**beta2 - 1)*(1+(beta2/(alpha2+3.))*num**beta2)/(1+beta2/(alpha2+3.)) !between 0 and 1
num2 = rand()             !between 0 and 1
do while (num2 > num3)
    num = rand()
    num3 = num**(2+alpha2)*exp(num**beta2 - 1)*(1+(beta2/(alpha2+3.))*num**beta2)/(1+beta2/(alpha2+3.))
    num2 = rand()
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine part_distribution(Re, npvoidstot, npvoids, npshellstot, npshells)
!!!!!!!!!!!!!!!!!!!!!!!!!!!
use common
implicit none
!in 
real :: Re(nvoid)
real :: rhoe, rhoshell
!out
integer :: npvoidstot, npshellstot
integer :: npvoids(nvoid), npshells(nvoid)
!locals
integer :: i

!mean delta inside the void is -0.8
rhoe = rhoB * (deltae + 1.)

!how many particles per void: depends on Re
npvoidstot = 0
do i=1,nvoid
    npvoids(i) = int( npart*(4./3.)*pi*Re(i)**3*rhoe/total_mass ) 
    npvoidstot = npvoidstot + npvoids(i)
enddo

!mean delta inside the shell is 100
rhoshell = rhoB * (delta_shell + 1.)

!how many particles per shell: depends on Re
npshellstot = 0
do i=1,nvoid
    npshells(i) = int( npart*(4./3.)*pi*Re(i)**3*( (1+fshell)**3 - 1.)*rhoshell/total_mass )
    npshellstot = npshellstot + npshells(i)
enddo

write(*,*) 'Particles in voids:', npvoidstot
write(*,*) 'Particles in shells:', npshellstot

!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!
