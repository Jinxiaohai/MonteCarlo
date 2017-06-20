module definitions
  implicit none
  integer :: i, j, k, m, n, p ! dummy integers
  integer, allocatable :: A(:,:,:) ! matrix containing spins
  integer :: nx, ny, nz ! number of rows and cols of A
  DOUBLE PRECISION :: ntotal ! ntotal=DBLE(nx*ny*nz). The total number of lattices
  DOUBLE PRECISION :: temp, beta ! temperature, inverse temperature
  integer :: npass ! number of passes for MC algorithm
  integer :: ipass ! the current pass number
  integer :: nequil ! number of equilibration steps
  integer :: nskip ! sampling every nskip steps, reduce correlation error
  integer :: trial_spin ! values of changed spin
  DOUBLE PRECISION, external :: ran1
  DOUBLE PRECISION :: high_temp ! starting temp for scan
  DOUBLE PRECISION :: low_temp ! final temp for scan
  DOUBLE PRECISION :: temp_interval ! interval between scan points
  integer :: nscans ! number of scans (each at diff T)
  integer :: iscan ! current scan number
  DOUBLE PRECISION :: deltaU ! change in energy between 2 configs
  DOUBLE PRECISION :: log_eta ! log of random number to compare to
  DOUBLE PRECISION :: magnetization ! magnetization of all spins in lattice
  DOUBLE PRECISION :: magnetization_ave ! cumulative average magnetization
  DOUBLE PRECISION :: magnetization_old ! magnetization before or after a spin flip
  DOUBLE PRECISION :: magnetization2_ave ! cumulative average of mag. squared
  DOUBLE PRECISION :: energy ! energy of all spins in lattice
  DOUBLE PRECISION :: energy_old ! energy before or after a spin-flip
  DOUBLE PRECISION :: energy_ave ! cumulative average of energy
  DOUBLE PRECISION :: energy2_ave ! cumulative average of energy squared
  DOUBLE PRECISION :: J1 ! nearest exchange constant
  DOUBLE PRECISION :: J2 ! next nearest exchange constant
  integer :: output_count ! # times things have been added to averages
  save  i, j, k, m, n, p, A, nx, ny, nz, temp, beta, npass,    &
       ipass, nequil, nskip, trial_spin, high_temp, low_temp, &
       temp_interval, nscans, iscan, deltaU, log_eta,         &
       magnetization, magnetization_ave, magnetization_old,   &
       magnetization2_ave, energy, energy_old, energy_ave,    &
       energy2_ave, J1, J2, output_count
end module definitions

program ising ! 3D Monte Carlo Simulation of Ising Model
  use definitions
  implicit none

  print*, "The critical temperature is approximately 2.3"

  call read_input()
  call initial_spin()
  call temperature_scan()

  write(32,*)"finial spin configuration"
  do k = 1, nz+1
     write(32,*), "k=", k
     do i = 1, nx + 2
        write(32,*), A(i,:,k)
     enddo
  enddo

  close(32)
  close(33)
  close(34)

  print*, "Program ising.f90 complete!"

end program ising

double precision function ran1(idum)
  !_______RANDOM NUMBER GENERATING FUNCTION______!
  implicit none
  double precision :: r(97)
  integer, intent(IN) :: idum
  save
  integer, parameter :: M1=259200,IA1=7141,IC1=54773
  DOUBLE PRECISION, parameter :: RM1=1.0d0/M1
  integer, parameter :: M2=134456,IA2=8121,IC2=28411
  DOUBLE PRECISION, parameter :: RM2=1.0d0/M2
  integer, parameter :: M3=243000,IA3=4561,IC3=51349
  integer :: IX1, IX2, IX3, jjj
  integer :: iff=0
  if (idum < 0 .or. iff == 0) then
     iff = 1
     IX1 = mod(IC1-idum,M1)
     IX1 = mod(IA1*IX1+IC1,M1)
     IX2 = mod(IX1,M2)
     IX1 = mod(IA1*IX1+IC1,M1)
     IX3 = mod(IX1,M3)
     do jjj = 1,97
        IX1 = mod(IA1*IX1+IC1,M1)
        IX2 = mod(IA2*IX2+IC2,M2)
        r(jjj) = (dfloat(IX1)+dfloat(IX2)*RM2)*RM1
     end do
  end if
  IX1 = mod(IA1*IX1+IC1,M1)
  IX2 = mod(IA2*IX2+IC2,M2)
  IX3 = mod(IA3*IX3+IC3,M3)
  jjj = 1+(97*IX3)/M3
  if (jjj > 97 .or. jjj < 1) PAUSE
  ran1 = r(jjj)
  r(jjj) = (dfloat(IX1)+dfloat(IX2)*RM2)*RM1
end function ran1

subroutine read_input()
  use definitions
  implicit none
  open(unit=11,file='ising.in',status='old',action='read')
  read(11,*)
  read(11,*)nx
  read(11,*)
  read(11,*)ny
  read(11,*)
  read(11,*)nz
  read(11,*)
  read(11,*)npass
  read(11,*)
  read(11,*)nequil
  read(11,*)
  read(11,*)nskip
  read(11,*)
  read(11,*)high_temp
  read(11,*)
  read(11,*)low_temp
  read(11,*)
  read(11,*)temp_interval
  read(11,*)
  read(11,*)J1
  read(11,*)
  read(11,*)J2
  close(11)

  J1=11.6045*J1
  J2=11.6045*J2 !from mev to K
  ntotal = DBLE(nx*ny*nz)

  open(unit=32,file='spin-array',status='replace',action='write')
  open(unit=34,file='energy',status='replace',action='write')
  open(unit=33,file='magnetization',status='replace',action='write')

  write(33,*) "temp magnetization/lattice ave_magnetization^2/lattice susceptibility/lattice"
  write(34,*) "nx,ny,nz", nx, ny, nz
  write(34,*) "total MC steps,  warm up steps", npass, nequil
  write(34,*) "temp energy/lattice ave_energy^2/lattice C_v/lattice"
end subroutine read_input

subroutine initial_spin()
  use definitions
  implicit none

  allocate(A(nx+2,ny+2,nz+2))
  !initial spin configuration
  A(1,1,1) = 1
  do i = 1, nx+1
     A(i+1,1,1) = A(i,1,1)
  enddo
  do j = 1, ny+1
     A(:,j+1,1) = A(:,j,1)
  enddo
  do k = 1, nz+1
     A(:,:,k+1) = A(:,:,k)
  enddo

  write(32,*)"initial spin configuration"
  do k = 1, nz+1
     write(32,*), "k=", k
     do i = 1, nx + 2
        write(32,*), A(i,:,k)
     enddo
  enddo
end subroutine initial_spin

subroutine temperature_scan()
  use definitions
  implicit none
  nscans = int((high_temp - low_temp)/temp_interval) + 1
  scan_loop: do iscan = 1, nscans
     temp = high_temp - temp_interval*DBLE(iscan-1)
     print*, "Running program for T =", temp

     beta = 1.0/temp
     output_count = 0
     energy_ave = 0.0
     energy2_ave = 0.0
     magnetization_ave = 0.0
     magnetization2_ave = 0.0

     call MC_steps()

     write(33, *) temp, abs(magnetization_ave/(DBLE(output_count))), &
          magnetization2_ave/(DBLE(output_count))*ntotal,        &
          beta*(magnetization2_ave/(DBLE(output_count)) -        &
          (magnetization_ave/(DBLE(output_count)))**2)*ntotal
     write(34, *) temp, energy_ave/(DBLE(output_count)),             & 
          energy2_ave/(DBLE(output_count))*ntotal,               &
          (beta**2)*(energy2_ave/(DBLE(output_count)) -           &
          (energy_ave/(DBLE(output_count)))**2)*ntotal
     print*, "output_count", output_count
  enddo scan_loop
end subroutine temperature_scan

subroutine MC_steps
  use definitions
  implicit none
  MC_passes: do ipass = 0, npass
     if (ipass > nequil) then
        output_count = output_count + 1
        if (output_count .EQ. 1) then
           magnetization = DBLE(sum(A(2:nx+1,2:ny+1,2:nz+1)))
           magnetization_old = magnetization
        else
           if (-beta*deltaU > log_eta) then
              magnetization = magnetization_old + DBLE(trial_spin*2)
              magnetization_old  = magnetization
           endif
        endif
        if (MOD(ipass,nskip) == 0) then
           magnetization_ave = magnetization_ave + magnetization/ntotal
           magnetization2_ave = magnetization2_ave + (magnetization/ntotal)**2
        endif

        if (output_count .EQ. 1) then
           energy = 0.0
           call current_energy()
           energy_old = energy
        else
           if (deltaU < 0) then
              energy = energy_old + deltaU
              energy_old = energy
           elseif (-beta*deltaU > log_eta) then
              energy = energy_old + deltaU
              energy_old = energy
           endif
        endif
        if (MOD(ipass,nskip) == 0) then
           energy_ave = energy_ave + energy/ntotal
           energy2_ave = energy2_ave + (energy/ntotal)**2
        endif

        if (MOD(ipass,nskip) /= 0) then
           output_count = output_count - 1
        endif
     endif

     m = nint((nx-1)*ran1(5) + 2) ! choose a random row, from 2 to (nx +1)
     n = nint((ny-1)*ran1(5) + 2) ! choose a random column, from 2 to (ny +1)
     p = nint((nz-1)*ran1(5) + 2) ! choose a random column, from 2 to (nz +1)
     trial_spin = -A(m,n,p) ! trial spin value

     call deltaU_flip

     if (deltaU < 0) then
        call spin_flip
     else
        log_eta = dlog(ran1(5) + 1.0d-4) ! random number 0-1 (+ tiny offset)
        if (-beta*deltaU > log_eta) then
           call spin_flip
        endif
     endif

  enddo MC_passes
end subroutine MC_steps

subroutine current_energy()
  use definitions
  implicit none
  do i = 2, nx + 1
     do j = 2, ny + 1
        do k = 2, nz + 1
           ! six nearest neighboring
           energy = energy - J1*DBLE(A(i,j,k)*&
                ( A(i-1,j,k)+A(i+1,j,k)+A(i,j-1,k)+A(i,j+1,k)+A(i,j,k-1)+A(i,j,k+1) ))

           ! tweleve  next nearest neighboring
           energy = energy - J2*DBLE(A(i,j,k)*&
                ( A(i-1,j-1,k)+A(i-1,j+1,k)+A(i+1,j-1,k)+A(i+1,j+1,k) + &
                A(i-1,j,k-1)+A(i-1,j,k+1)+A(i+1,j,k-1)+A(i+1,j,k+1) + &
                A(i,j-1,k-1)+A(i,j-1,k+1)+A(i,j+1,k-1)+A(i,j+1,k+1) ))
        enddo
     enddo
  enddo
  energy = energy/2  !double counting
end subroutine current_energy

subroutine deltaU_flip
  use definitions
  implicit none

  ! six nearest neighboring
  deltaU = -J1*DBLE(trial_spin*2*&
       ( A(m-1,n,p)+A(m+1,n,p)+A(m,n-1,p)+A(m,n+1,p)+A(m,n,p-1)+A(m,n,p+1) ))

  ! tweleve  next nearest neighboring
  deltaU = deltaU - J2*DBLE(trial_spin*2*&
       ( A(m-1,n-1,p)+A(m-1,n+1,p)+A(m+1,n-1,p)+A(m+1,n+1,p) + &
       A(m-1,n,p-1)+A(m-1,n,p+1)+A(m+1,n,p-1)+A(m+1,n,p+1) +   &
       A(m,n-1,p-1)+A(m,n-1,p+1)+A(m,n+1,p-1)+A(m,n+1,p+1)   ))
end subroutine deltaU_flip

subroutine spin_flip
  use definitions
  implicit none
  A(m,n,p) = trial_spin
  ! six faces
  if (m == 2) A(nx+2,n,p) = trial_spin
  if (n == 2) A(m,ny+2,p) = trial_spin
  if (p == 2) A(m,n,nz+2) = trial_spin
  if (m == nx+1) A(1,n,p) = trial_spin
  if (n == ny+1) A(m,1,p) = trial_spin
  if (p == nz+1) A(m,n,1) = trial_spin

  ! tweleve edges
  if (m == 2 .AND. n == 2) A(nx+2,ny+2,p) = trial_spin
  if (m == 2 .AND. p == 2) A(nx+2,n,nz+2) = trial_spin
  if (n == 2 .AND. p == 2) A(m,ny+2,nz+2) = trial_spin
  if (m == 2 .AND. n == ny+1) A(nx+2,1,p) = trial_spin
  if (m == nx+1 .AND. n == 2) A(1,ny+2,p) = trial_spin
  if (m == 2 .AND. p == nz+1) A(nx+2,n,1) = trial_spin
  if (m == nx+1 .AND. p == 2) A(1,n,nz+2) = trial_spin
  if (n == 2 .AND. p == nz+1) A(m,ny+2,1) = trial_spin
  if (n == ny+1 .AND. p == 2) A(m,1,nz+2) = trial_spin
  if (m == nx+1 .AND. n == ny+1) A(1,1,p) = trial_spin
  if (m == nx+1 .AND. p == nz+1) A(1,n,1) = trial_spin
  if (n == ny+1 .AND. p == nz+1) A(m,1,1) = trial_spin

  ! eight corners
  if (m == 2 .AND. n == 2 .AND. p == 2) A(nx+2,ny+2,nz+2) = trial_spin
  if (m == 2 .AND. n == 2 .AND. p == nz+1) A(nx+2,ny+2,1) = trial_spin
  if (m == 2 .AND. n == ny+1 .AND. p == 2) A(nx+2,1,nz+2) = trial_spin
  if (m == nx+1 .AND. n == 2 .AND. p == 2) A(1,ny+2,nz+2) = trial_spin
  if (m == nx+1 .AND. n == ny+1 .AND. p == 2) A(1,1,nz+2) = trial_spin
  if (m == nx+1 .AND. n == 2 .AND. p == nz+1) A(1,ny+2,1) = trial_spin
  if (m == 2 .AND. n == ny+1 .AND. p == nz+1) A(nx+2,1,1) = trial_spin
  if (m == nx+1 .AND. n == ny+1 .AND. p == nz+1) A(1,1,1) = trial_spin
end subroutine spin_flip

