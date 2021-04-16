!==============================================================================
! PROGRAM trRIXS is the main driving program. 
!==============================================================================
program trRIXS
use mod_parameters
use mod_hilbert_space
use mod_hamiltonian
use mod_random
use mod_lanczos
use mod_chebyshev
implicit none
integer(KIND=8), parameter :: numtime = 1200
integer(KIND=8), parameter :: ntskip = 10
double precision, parameter :: time0 = 0.0d0
double precision, parameter :: start_time = 0.0d0
double precision, parameter :: end_time = 12000.0d0
double precision, parameter :: kb = 8.617e-5
logical compute_RIXS
logical orthogonalize
logical, parameter :: verbose = .false.
character filename1*200
character filename2*200
character screen_output_file*200
character filesuffix*100
integer istate
integer, parameter :: screen = 6      !output units
integer, parameter :: nfwrk1 = 70
integer, parameter :: nfwrk2 = 71
integer(KIND=8), parameter :: numstates = 3
integer(KIND=8), parameter :: numomega = 500
integer(KIND=8) spin, site, counter
integer(KIND=8) n_lanczos_states
integer(KIND=8) nt, i, j, ch, nw
integer(KIND=8) maxHUDlen, maxHUm1Dlen, maxHUDm1len
integer(KIND=8) HUDlen, HUm1Dlen, HUDm1len
integer(KIND=8), allocatable, dimension(:) :: H0i, H0j
integer(KIND=8), allocatable, dimension(:) :: HUDi, HUDj
integer(KIND=8), allocatable, dimension(:) :: HUm1Di, HUm1Dj
integer(KIND=8), allocatable, dimension(:) :: HUDm1i, HUDm1j
double precision, dimension(1:numstates) :: bose
double precision dnorm2, displacement, heaviside, domega, integral, I0
double precision, parameter :: omegab = 9.0d0
double precision, parameter :: omegaa = -0.5d0
double precision, dimension(1:numstates) :: SiSip1, Mz, Stot
double precision, dimension(0:numomega,1:numstates) :: RIXS
double precision, dimension(0:numomega) :: RIXS0
double precision, allocatable, dimension(:) :: H0
double precision, allocatable, dimension(:,:) :: States   
double precision, allocatable, dimension(:) :: Ei
double precision, allocatable, dimension(:) :: Ei_t
double precision, dimension(1:2) :: VecPotential
double precision DME_xin, DME_yin
double precision DME_xout, DME_yout
double precision win, omega, Stmp, SiSip1tmp, Mtmp
double precision dtime, time, temp, Z
double precision Eg
double precision pi
double precision rseed, fill
double precision, dimension(1:Ncu) :: FillCu
double precision, dimension(1:Nox) :: FillOx
double precision delta, arg
double precision, dimension(1:nfmax) :: EE
double complex cfac1, cfac2, ctmp
double complex, parameter :: eye = (0.0d0,1.0d0)
double complex, allocatable, dimension(:) :: HUD, HUm1D, HUDm1
double complex, allocatable, dimension(:) :: Psi_t, Psi_tp   
double complex, dimension(1:nfmax) :: final
double complex, allocatable, dimension(:) :: psi1, psi2,  cvec, psi0
double complex, allocatable, dimension(:) :: fvecs, fvecs2, cgs, gs


rseed = 0.42542323d0
pi = 2.0d0*asin(1.0d0)
win = 0.50d0
domega = (omegab-omegaa)/dfloat(numomega)

if(nox.lt.10)then
 if(Amp.lt.10.0d0)then
  if(win.lt.0.0d0)then
   write(unit=filesuffix,fmt=10) Ncu,Nox,amp, abs(win), gz, dimag(igam), nquanta
  else
   write(unit=filesuffix,fmt=11) Ncu,Nox,amp, win, gz, dimag(igam), nquanta
  endif
 else
  if(win.lt.0.0d0)then
   write(unit=filesuffix,fmt=12) Ncu,Nox,amp, abs(win), gz, dimag(igam), nquanta
  else
   write(unit=filesuffix,fmt=13) Ncu,Nox,amp, win, gz, dimag(igam), nquanta
  endif
 endif
else
 if(Amp.lt.10.0d0)then
  if(win.lt.0.0d0)then
   write(unit=filesuffix,fmt=14) Ncu,Nox,amp, abs(win), gz, dimag(igam), nquanta
  else
   write(unit=filesuffix,fmt=15) Ncu,Nox,amp, win, gz, dimag(igam), nquanta
  endif
 else
  if(win.lt.0.0d0)then
   write(unit=filesuffix,fmt=16) Ncu,Nox,amp, abs(win), gz, dimag(igam), nquanta
  else
   write(unit=filesuffix,fmt=17) Ncu,Nox,amp, win, gz, dimag(igam), nquanta
  endif
 endif
endif

!Notation is I incoming, outgoing. 
! sigma is in the y-direction.
10 format('_Cu',i1,'O',i1,'_Ips_Spump_Amp',f5.3,'_winm',f5.3,'_geq',f5.3,'_gam',f5.3,'_',i3.3,'.dat')
11 format('_Cu',i1,'O',i1,'_Ips_Spump_Amp',f5.3,'_win',f5.3,'_geq',f5.3,'_gam',f5.3,'_',i3.3,'.dat')
12 format('_Cu',i1,'O',i1,'_Ips_Spump_Amp',f5.2,'_winm',f5.3,'_geq',f5.3,'_gam',f5.3,'_',i3.3,'.dat')
13 format('_Cu',i1,'O',i1,'_Ips_Spump_Amp',f5.2,'_win',f5.3,'_geq',f5.3,'_gam',f5.3,'_',i3.3,'.dat')
14 format('_Cu',i1,'O',i2,'_Ips_Spump_Amp',f5.3,'_winm',f5.3,'_geq',f5.3,'_gam',f5.3,'_',i3.3,'.dat')
15 format('_Cu',i1,'O',i2,'_Ips_Spump_Amp',f5.3,'_win',f5.3,'_geq',f5.3,'_gam',f5.3,'_',i3.3,'.dat')
16 format('_Cu',i1,'O',i2,'_Ips_Spump_Amp',f5.2,'_winm',f5.3,'_geq',f5.3,'_gam',f5.3,'_',i3.3,'.dat')
17 format('_Cu',i1,'O',i2,'_Ips_Spump_Amp',f5.2,'_win',f5.3,'_geq',f5.3,'_gam',f5.3,'_',i3.3,'.dat')

!define the DME
DME_xin = cos(45.0d0*pi/180.0d0)
DME_yin = 0.0d0
DME_yout = cos(56.0d0*pi/180.0d0)
DME_xout = 0.0d0

write(unit=screen,fmt=*) 'trRIXS for a Cu3O8 Cluster.'
write(unit=screen,fmt=*) 'Random seed = ', rseed
write(unit=screen,fmt=*) 'Outputting file to suffix: ', trim(filesuffix)

call init_ranLM(rseed)

! Step 1: Determine the size of the relevant Hilbert Space
! First determine the size of the electronic subspace
NelUD   = size_hilbert_space(Nup,Ndn,N)
NelUm1D = size_hilbert_space(Nup-1,Ndn,N)
NelUDm1 = size_hilbert_space(Nup,Ndn-1,N)

! Determine the total hilbert space size
hs_size_UD = NelUD*(nquanta+1)
hs_size_Um1D = NelUm1D*(nquanta+1)
hs_size_UDm1 = NelUDm1*(nquanta+1)

! Output the size of the space to the user. 

write(unit=screen,fmt=*) ' '
write(unit=screen,fmt=*) 'Number up carriers: ', nup
write(unit=screen,fmt=*) 'Number dn carriers: ', ndn
write(unit=screen,fmt=*) 'Number of phonon quanta: ', nquanta
write(unit=screen,fmt=*) ' '
write(unit=screen,fmt=*) 'Electronic HS size for (Nup,Ndn) problem = ', NelUD
write(unit=screen,fmt=*) 'Electronic HS size for (Nup-1,Ndn) problem = ', NelUm1D
write(unit=screen,fmt=*) 'Electronic HS size for (Nup,Ndn-1) problem = ', NelUDm1
write(unit=screen,fmt=*) ' '
write(unit=screen,fmt=*) 'Total HS size for (Nup,Ndn) problem = ', hs_size_UD
write(unit=screen,fmt=*) 'Total HS size for (Nup-1,Ndn) problem = ', hs_size_Um1D
write(unit=screen,fmt=*) 'Total HS size for (Nup,Ndn-1) problem = ', hs_size_UDm1

! Step 2: Construct the basis for the problems
write(unit=screen,fmt=*) ' '
write(unit=screen,fmt=*) 'Allocating space for basis arrays.'
call init_basis()
write(unit=screen,fmt=*) 'constructing basis.'
call construct_basis(nup,ndn,N,NelUD,basis_UD)
call construct_basis(nup-1,ndn,N,NelUm1D,basis_Um1D)
call construct_basis(nup,ndn-1,N,NelUDm1,basis_UDm1)

!Output the basis states
if(verbose)then
 write(unit=filename1,fmt="('basis_',i1,'u',i1,'d.dat')") nup, ndn
 open(file=FILENAME1,unit=60,action='write')
 do i = 1,NelUD
  write(unit=60,fmt=*) i, basis_UD(i)
 enddo
 close(unit=60)

 write(unit=filename1,fmt="('basis_',i1,'u',i1,'d.dat')") nup-1, ndn
 open(file=FILENAME1,unit=60,action='write')
 do i = 1,NelUm1D
  write(unit=60,fmt=*) i, basis_Um1D(i)
 enddo
 close(unit=60)

 write(unit=filename1,fmt="('basis_',i1,'u',i1,'d.dat')") nup, ndn-1
 open(file=FILENAME1,unit=60,action='write')
 do i = 1,NelUDm1
  write(unit=60,fmt=*) i, basis_UDm1(i)
 enddo
 close(unit=60)
endif

! Step 3: Setup all the onsite energies and hopping elements for time 
call init_cluster()


!First alloccate space for the Hamiltonian in the (nup,ndn) problem.
HUDlen = 0
maxHUDLen = 30*NelUD
allocate(HUDi(1:maxHUDLen))
allocate(HUDj(1:maxHUDLen))
allocate(HUD(1:maxHUDLen))

!Now construct the Hamiltonian for the (Nup, Ndn) problem. 
ch = 0                  !Generate the Hamiltonian with no core hole
time = start_time       !assume time t = start_time 
displacement = X0 + XAmp*sin(omegamode*(time-time0))*exp(-(time-time0)/tau)*heaviside(time-time0)
print*, ' '
print*, 'Constructing H for (nup,ndn) problem.'
call construct_hamiltonian(HUDi,HUDj,HUD,HUDlen,maxHUDlen,basis_UD,NelUD,ch,displacement)

!output the matricies for testing purpose
if(verbose)then
 open(file='Hij_el_2u1d.txt',unit=60,action='write')
 do i = 1,HUDlen
  write(unit=60,fmt="(i10,' ',i10,' ',E15.8,' ',E15.8)") & 
                    HUDi(i), HUDj(i), dreal(HUD(i)), dimag(HUD(i))
 enddo
 close(unit=60)
endif

! Now compute the ground and first few exited states
! Allocate space to store the states and their energies
allocate(states(1:hs_size_UD,1:numstates))
allocate(Ei(1:numstates))
!convert the complex valued Hamiltonian to real
allocate(H0(1:HUDlen))
allocate(H0i(1:HUDlen))
allocate(H0j(1:HUDlen))
do i = 1,HUDlen
 H0(i) = dreal(HUD(i))
 H0i(i) = HUDi(i)
 H0j(i) = HUDj(i)
enddo
n_lanczos_states = 200
orthogonalize = .true.

call get_ground_state(states,Ei,numstates,basis_UD,hs_size_UD,NelUD,& 
                 H0i,H0j,H0,HUDlen,maxHUDlen,n_lanczos_states,orthogonalize) 

!Output the energies of the eigenstates to be sure.
print*, 'Energies of the lowest lying states: '
do i = 1,numstates
 print*, i, Ei(i)
enddo

dtime = (end_time-start_time)/dfloat(numtime)

!Define the vector potential as a function of time
filename1 = 'At'//trim(filesuffix)
open(file=filename1,unit=102,action='write')
do nt = 1,numtime
 time = start_time + dtime*dfloat(nt-1)
 !Assume it points in the x-direction
 call get_vec_potential(VecPotential,time)
 write(102,*) time*hbar, VecPotential(1), VecPotential(2)
enddo
close(unit=102)

allocate(psi_t(1:hs_size_UD))
allocate(psi_tp(1:hs_size_UD))
allocate(Ei_t(1:numtime))
allocate(Psi0(1:hs_size_UD))

!Start from the ground state
Psi_t(:) = states(:,1) 

filename1 = "projections"//trim(filesuffix)
open(unit=nfwrk1,action='write',file=filename1,status='replace')
close(unit=nfwrk1)
filename2 = 'rixs'//trim(filesuffix)
open(file=filename2,unit=nfwrk1,action='write',status='replace')
close(unit=nfwrk1)

!Declare arrays to hold the initial state of the RIXS calculation and 
!set the Lanczos parameters
allocate(gs(1:hs_size_UD))

counter = 0
do nt = 1,numtime-1,ntskip
 time = start_time + dtime*dfloat(nt-1)

 displacement = X0 + XAmp*sin(omegamode*(time-time0))*exp(-(time-time0)/tau)*heaviside(time-time0)

 !Compute the ground state
 n_lanczos_states = 200
 orthogonalize = .true.
 ch = 0 
 call construct_hamiltonian(HUDi,HUDj,HUD,HUDlen,maxHUDlen,basis_UD,NelUD,ch,displacement)
 do i = 1,HUDlen
  H0(i) = dreal(HUD(i))
  H0i(i) = HUDi(i)
  H0j(i) = HUDj(i)
 enddo
 call get_ground_state(states,Ei,numstates,basis_UD,hs_size_UD,NelUD,&
                       H0i,H0j,H0,HUDlen,maxHUDlen,n_lanczos_states,orthogonalize)
 
 RIXS = 0.0d0      !Init the RIXS intensity to zero.
 do istate = 1,numstates
  gs(:) = states(:,istate)  
  call get_filling_and_spin(gs,hs_size_UD,basis_UD,NelUD,Stmp,SiSip1tmp,Mtmp,Fillcu,&
                             Fillox,Fill)
  Stot(istate) = Stmp
  SiSip1(istate) = SiSip1tmp
  Mz(istate) = Mtmp

  orthogonalize = .false.
  !now compute the RIXS spectra for all of the points. 
  write(unit=screen,fmt=*) 'Calculating RIXS spectra at time ', time
  Psi0(:) = dcmplx(0.0d0,0.0d0)
  Eg = Ei(istate)
  !Loop over all possible core hole locations
  do ch = 1,nox
   !Destroy a spin up hole in the valence band
   spin = spinup
   HUm1Dlen = 0
   maxHUm1Dlen = 30*NelUm1D
   allocate(HUm1Di(1:maxHUm1DLen))
   allocate(HUm1Dj(1:maxHUm1DLen))
   allocate(HUm1D(1:maxHUm1DLen))
   allocate( cgs(1:hs_size_Um1D))
   allocate(psi1(1:hs_size_Um1D))
   allocate(psi2(1:hs_size_UD))
   allocate(cvec(1:hs_size_Um1D))
   allocate(fvecs(1:hs_size_Um1D))

   !Build the intermediate state Hamiltonian at time t
   call construct_hamiltonian(HUm1Di,HUm1Dj,HUm1D,HUm1Dlen,maxHUm1Dlen,& 
                              basis_Um1D,NelUm1D,ch,displacement)
   !first apply the destruction operator to the gs and store it in cgs
   psi1 = dcmplx(0.0d0,0.0d0)
   !Oxygen x
   site = 3*ncu+0*nox+ch
   call remove_particle(gs,hs_size_UD,basis_UD,NelUD,& 
                        cgs,hs_size_Um1D,basis_Um1d,NelUm1D,site,spin,N)
   psi1 = dcmplx(DME_yin,0.0d0)*cgs
 
   !Oxygen y
   site = 3*ncu+1*nox+ch
   call remove_particle(gs,hs_size_UD,basis_UD,NelUD,& 
                        cgs,hs_size_Um1D,basis_Um1d,NelUm1D,site,spin,N)
   psi1 = psi1 + dcmplx(DME_xin,0.0d0)*cgs
 
 
   !here |psi_1> = [D_x + D_y] |gs>  
   call get_frequ_vec_lanczos(psi1,hs_size_Um1D,NelUm1D,basis_Um1D, & 
                              HUm1Di,HUm1Dj,HUm1D,HUm1Dlen,maxHUm1Dlen,&
                              win,Eg,fvecs,n_lanczos_states,orthogonalize)
   !now add a paticle
   !|cvec> = |i><i| [D_x + D_y] |gs>
   cvec = fvecs
 
   site = 3*ncu+0*nox+ch
   call add_particle(cvec,hs_size_Um1D,basis_Um1D,NelUm1D,&  
                     psi2,hs_size_UD,  basis_UD,NelUD,site,spin,N)
   psi0 = psi0 + dcmplx(DME_yout,0.0d0)*psi2
 
   site = 3*ncu+1*nox+ch
   call add_particle(cvec,hs_size_Um1D,basis_Um1D,NelUm1D,&  
                     psi2,hs_size_UD,  basis_UD,NelUD,site,spin,N)
   psi0 = psi0 + dcmplx(DME_xout,0.0d0)*psi2
 
   deallocate(HUm1Di)
   deallocate(HUm1Dj)
   deallocate(HUm1D)
   deallocate(cgs)
   deallocate(psi1)
   deallocate(psi2)
   deallocate(cvec)
   deallocate(fvecs)
 
   !Destroy a spin up hole in the valence band
   spin = spindn
   HUDm1len = 0
   maxHUDm1len = 30*NelUDm1
   allocate(HUDm1i(1:maxHUDm1Len))
   allocate(HUDm1j(1:maxHUDm1Len))
   allocate(HUDm1(1:maxHUDm1Len))
   allocate( cgs(1:hs_size_UDm1))
   allocate(psi1(1:hs_size_UDm1))
   allocate(psi2(1:hs_size_UD))
   allocate(cvec(1:hs_size_UDm1))
   allocate(fvecs(1:hs_size_UDm1))
 
   !Build the intermediate state Hamiltonian at time t
   call construct_hamiltonian(HUDm1i,HUDm1j,HUDm1,HUDm1len,maxHUDm1len,& 
                              basis_UDm1,NelUDm1,ch,displacement)
   !first apply the destruction operator to the gs and store it in cgs
   psi1 = dcmplx(0.0d0,0.0d0)
   !Oxygen x
   site = 3*ncu+0*nox+ch
   call remove_particle(gs,hs_size_UD,basis_UD,NelUD,& 
                        cgs,hs_size_UDm1,basis_UDm1,NelUDm1,site,spin,N)
   psi1 = dcmplx(DME_yin,0.0d0)*cgs
 
   !Oxygen y
   site = 3*ncu+1*nox+ch
   call remove_particle(gs,hs_size_UD,basis_UD,NelUD,& 
                        cgs,hs_size_UDm1,basis_UDm1,NelUDm1,site,spin,N)
   psi1 = psi1 + dcmplx(DME_xin,0.0d0)*cgs
 
 
   !here |psi_1> = [D_x + D_y] |gs>  
   call get_frequ_vec_lanczos(psi1,hs_size_UDm1,NelUDm1,basis_UDm1, & 
                              HUDm1i,HUDm1j,HUDm1,HUDm1len,maxHUDm1len,&
                              win,Eg,fvecs,n_lanczos_states,orthogonalize)
   !now add a paticle
   !|cvec> = |i><i| [D_x + D_y] |gs>
   cvec = fvecs
 
   site = 3*ncu+0*nox+ch
   call add_particle(cvec,hs_size_UDm1,basis_UDm1,NelUDm1,&  
                     psi2,hs_size_UD,  basis_UD,NelUD,site,spin,N)
   psi0 = psi0 + dcmplx(DME_yout,0.0d0)*psi2
 
   site = 3*ncu+1*nox+ch
   call add_particle(cvec,hs_size_UDm1,basis_UDm1,NelUDm1,&  
                     psi2,hs_size_UD,  basis_UD,NelUD,site,spin,N)
   psi0 = psi0 + dcmplx(DME_xout,0.0d0)*psi2
 
   deallocate(HUDm1i)
   deallocate(HUDm1j)
   deallocate(HUDm1)
   deallocate(cgs)
   deallocate(psi1)
   deallocate(psi2)
   deallocate(cvec)
   deallocate(fvecs)
  enddo
 
  !now we need to do the final projection on to the final states
  ch = 0
  call construct_hamiltonian(HUDi,HUDj,HUD,HUDlen,maxHUDlen,basis_UD,NelUD,ch,displacement)
  call project_on_final(psi0,hs_size_UD,basis_UD,NelUD,EE,Final,& 
                        n_lanczos_states,HUDi,HUDj,HUD,HUDlen,maxHUDlen,& 
                        j,orthogonalize)
  !and build the RIXS spectra
  open(unit=nfwrk2,file=filename2,access='append')
  integral = 0.0d0
  do nw = 0,numomega
   omega = dfloat(nw)*(omegab-omegaa)/dfloat(numomega) + omegaa
   do i = 1,j
    ctmp = final(i)
    arg = omega + Eg - EE(i)
    delta = (sig*0.5d0/pi)/(arg*arg+0.25d0*sig*sig)
    RIXS(nw,istate) = RIXS(nw,istate) + delta*dreal(ctmp*conjg(ctmp))
   enddo
  enddo
 enddo !ito state = 1,3

 !output the data after constructing the total RIXS
 RIXS0 = 0.0d0
 temp = 21.3013d0 + 181.2d0*(1.0d0-exp(-time*0.001d0/1.949d0))
 bose = (2.0d0*Stot+1.0d0)*exp(-(Ei-Ei(1))/(temp*kb))
 bose = bose/sum(bose)
 open(unit=nfwrk2,file=filename2,access='append')
 integral = 0.0d0
 do nw = 0,numomega
  omega = dfloat(nw)*domega + omegaa  
  RIXS0(nw) = RIXS0(nw) + RIXS(nw,1)*bose(1) & 
                        + RIXS(nw,2)*bose(2) & 
                        + RIXS(nw,3)*bose(3)
  write(nfwrk2,*) time*hbar, omega, RIXS0(nw)
  if(omega.gt.3.0d0.and.omega.le.4.1d0)then
   integral = integral + RIXS0(nw)*domega
  endif
 enddo
 SiSip1tmp = sum(SiSip1*bose)
 write(nfwrk2,*) ' '
 close(unit=nfwrk2)
 if(nt.eq.1) I0 = integral
enddo !do nt = 0,numtime


stop
end program trRIXS

double precision function heaviside(x)
implicit none
double precision x
if(x.gt.0.0d0)then
  heaviside = 1.0d0
elseif(x.lt.0.0d0)then
  heaviside = 0.0d0
else
  heaviside = 0.5d0
endif
return
end function heaviside
