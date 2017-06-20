program ising ! 2D Monte Carlo Simulation of Ising Model

  implicit none

  integer :: i, j, m, n, m2, n2 ! dummy integers
  integer, allocatable :: A(:,:) ! matrix containing spins
  integer :: nrows, ncols ! number of rows and cols of A
  DOUBLE PRECISION :: temp, beta ! temperature, inverse temperature
  integer :: ConfigType ! starting configuration type
  integer :: npass ! number of passes for MC algorithm
  integer :: ipass ! the current pass number
  integer :: nequil ! number of equilibration steps
  integer :: trial_spin ! values of changed spin
  DOUBLE PRECISION :: high_temp ! starting temp for scan
  DOUBLE PRECISION :: low_temp ! final temp for scan
  DOUBLE PRECISION :: temp_interval ! interval between scan points
  integer :: nscans ! number of scans (each at diff T)
  integer :: iscan ! current scan number
  logical :: MovieOn ! set to .true. to make movie of 1 temp
  DOUBLE PRECISION :: deltaU ! change in energy between 2 configs
  DOUBLE PRECISION :: energy_old ! energy before or after a spin-flip
  DOUBLE PRECISION :: log_eta ! log of random number to compare to
  DOUBLE PRECISION :: magnetization ! magnetization of all spins in lattice
  DOUBLE PRECISION :: magnetization_old ! magnetization before or after a spin flip
  DOUBLE PRECISION :: magnetization_ave ! cumulative average magnetization
  DOUBLE PRECISION :: magnetization2_ave ! cumulative average of mag. squared
  DOUBLE PRECISION :: energy ! energy of all spins in lattice
  DOUBLE PRECISION :: energy_ave ! cumulative average of energy
  DOUBLE PRECISION :: energy2_ave ! cumulative average of energy squared
  integer :: output_count ! # times things have been added to averages

  print*, "The critical temperature is approximately 2.3"

  open(unit=11,file='ising.in')
  read(11,*)
  read(11,*) nrows
  read(11,*)
  read(11,*) ncols
  read(11,*)
  read(11,*) npass
  read(11,*)
  read(11,*) nequil
  read(11,*)
  read(11,*) high_temp
  read(11,*)
  read(11,*) low_temp
  read(11,*)
  read(11,*) temp_interval
  close(11)
  allocate(A(nrows+2,ncols+2))

  open(unit=32,file='spin-array',status='replace',action='write')
  write(32,*) "nrows", nrows
  write(32,*) "ncols", ncols
  nscans = int((high_temp - low_temp)/temp_interval) + 1
  write(32,*) "nscans", nscans

  open(unit=33,file='magnetization',status='replace',action='write')
  write(33,*) "temp ave_magnetization ave_magnetization^2 susceptibility"
  open(unit=34,file='energy',status='replace',action='write')
  write(34,*) "temp ave_energy ave_energy^2 C_v"

  scan_loop: do iscan = 1, nscans
     temp = high_temp - temp_interval*(iscan-1)
     print*, "Running program for T =", temp

     beta = 1.0/temp
     output_count = 0
     energy_ave = 0.0
     energy2_ave = 0.0
     magnetization_ave = 0.0
     magnetization2_ave = 0.0

     A(1,1) = 1
     ! chushihua
     do i = 1, nrows+1
        A(i+1,1) = -A(i,1)
     enddo
     do j = 1, ncols+1
        A(:,j+1) = -A(:,j)
     enddo

     write(32,*)"initial spin configuration"
     do i = 1, nrows + 2
        write(32,*), A(i,:)
     enddo

     MC_passes: do ipass = 0, npass
        if (ipass > nequil) then
           output_count = output_count + 1
           if (output_count .EQ. 1) then
              magnetization = sum(A(2:nrows+1,2:ncols+1))
              magnetization_old = magnetization
           else
              if (-beta*deltaU > log_eta) then
                 magnetization = magnetization_old + trial_spin*2
                 magnetization_old  = magnetization
              endif
           endif
           magnetization_ave = magnetization_ave + magnetization/(ncols*nrows*1.0)
           magnetization2_ave = magnetization2_ave + (magnetization/(ncols*nrows*1.0))**2

           if (output_count .EQ. 1) then
              energy = 0.0
              do i = 2, nrows + 1
                 do j = 2, ncols + 1
                    energy = energy - A(i,j)*(A(i-1,j)+A(i+1,j)+A(i,j-1)+A(i,j+1))
                 enddo
              enddo
              energy = energy/2
              energy_old = energy
           else
              if (-beta*deltaU > log_eta) then
                 energy = energy_old + deltaU
                 energy_old = energy
              endif
           endif
           energy_ave = energy_ave + energy/(ncols*nrows)
           energy2_ave = energy2_ave + (energy/(ncols*nrows))**2
        endif

        m = nint((nrows-1)*ran1(5) + 2) ! choose a random row, from 2 to (nrows +1)
        n = nint((ncols-1)*ran1(5) + 2) ! choose a random column, from 2 to (ncols +1)
        trial_spin = -A(m,n) ! trial spin value

        deltaU = -trial_spin*(A(m-1,n)+A(m+1,n)+A(m,n-1)+A(m,n+1))*2
        if (deltaU < 0) then
           A(m,n) = trial_spin
           if (m == 2) A(nrows+2,n) = trial_spin
           if (m == nrows+1) A(1,n) = trial_spin
           if (n == 2) A(m,ncols+2) = trial_spin
           if (n == ncols+1) A(m,1) = trial_spin
        else
           log_eta = dlog(ran1(5) + 1.0d-4) ! random number 0-1 (+ tiny offset)
           if (-beta*deltaU > log_eta) then
              A(m,n) = trial_spin
              if (m == 2) A(nrows+2,n) = trial_spin
              if (m == nrows+1) A(1,n) = trial_spin
              if (n == 2) A(m,ncols+2) = trial_spin
              if (n == ncols+1) A(m,1) = trial_spin
           endif
        endif

     enddo MC_passes

     write(32,*)"iscan",  iscan
     do i = 1, nrows + 2
        write(32,*), A(i,:)
     enddo

     write(33, *) temp, abs(magnetization_ave/output_count), &
          magnetization2_ave/output_count, &
          beta*(magnetization2_ave/output_count - &
          (magnetization_ave/output_count)**2)
     write(34, *) temp, energy_ave/output_count, & 
          energy2_ave/output_count, &
          (beta**2)*(energy2_ave/output_count - &
          (energy_ave/output_count)**2)

  enddo scan_loop

  close(32)
  close(33)
  close(34)

  print*, "Program ising.f90 complete!"

contains
  !_______RANDOM NUMBER GENERATING FUNCTION______!

  double precision function ran1(idum)
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

end program ising


