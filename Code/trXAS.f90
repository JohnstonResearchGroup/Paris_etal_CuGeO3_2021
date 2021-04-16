!==============================================================================
! PROGRAM trXAS is the main driving program. 
!==============================================================================
program trXAS
use mod_parameters
use mod_hilbert_space
use mod_hamiltonian
use mod_random
use mod_lanczos
use mod_chebyshev
implicit none
logical compute_xas
integer(KIND=8), parameter :: numtime = 3000
integer(KIND=8), parameter :: ntskip = 30
double precision, parameter :: start_time = -500.0d0
double precision, parameter :: end_time = 2500.0d0
logical orthogonalize
logical, parameter :: verbose = .false.
character filename1*200
character filename2*200
character filename3*200
character filesuffix*100
integer, parameter :: screen = 6      !output units
integer, parameter :: nfwrk1 = 70
integer, parameter :: nfwrk2 = 71
integer seed
integer(KIND=8), parameter :: numstates = 6
integer(KIND=8), parameter :: numomega = 1000
integer(KIND=8) spin, site, counter
integer(KIND=8) n_lanczos_states
integer(KIND=8) nt, i, j, ch, nw
integer(KIND=8) maxHUDlen, maxHUm1Dlen, maxHUDm1len
integer(KIND=8) HUDlen, HUm1Dlen, HUDm1len
integer(KIND=8), allocatable, dimension(:) :: H0i, H0j
integer(KIND=8), allocatable, dimension(:) :: HUDi, HUDj
integer(KIND=8), allocatable, dimension(:) :: HUm1Di, HUm1Dj
integer(KIND=8), allocatable, dimension(:) :: HUDm1i, HUDm1j
double precision dnorm2
double precision, parameter :: omegab = 10.0d0
double precision, parameter :: omegaa = -10.0d0
double precision, dimension(0:numomega) :: XAS 
double precision, allocatable, dimension(:) :: H0
double precision, allocatable, dimension(:,:) :: States   
double precision, allocatable, dimension(:) :: Ei
double precision, allocatable, dimension(:) :: Ei_t
double precision, dimension(1:2) :: VecPotential
double precision DME_xin, DME_yin, Mz, Stot, SiSip1
double precision, dimension(1:Ncu) :: FillCu
double precision, dimension(1:Nox) :: FillOx
double precision omega, Fill
double precision dtime, time, displacement, heaviside
double precision Eg
double precision pi
double precision rseed
double precision, dimension(1:nfmax) :: EE
double complex cfac1, cfac2, ctmp
double complex, parameter :: eye = (0.0d0,1.0d0)
double complex, allocatable, dimension(:) :: HUD, HUm1D, HUDm1
double complex, allocatable, dimension(:) :: Psi_t   
double complex, allocatable, dimension(:) :: Psi_tp   
double complex, dimension(1:nfmax) :: final
double complex, allocatable, dimension(:) :: psi1, psi0
double complex, allocatable, dimension(:) :: cgs, gs

seed = 54884
call srand(seed)
pi = 2.0d0*asin(1.0d0)

if(Nox.lt.10)then
 if(Amp.lt.10.0d0)then
  write(unit=filesuffix,fmt=11) Ncu, Nox, Amp, gz, dimag(igam), nquanta
 else
  write(unit=filesuffix,fmt=12) Ncu, Nox, Amp, gz, dimag(igam), nquanta
 endif
11 format('_Cu',i1,'O',i1,'_Spump_Amp',f5.3,'_geq',f5.3,'_gam',f5.3,'_',i3.3,'.dat')
12 format('_Cu',i1,'O',i1,'_Spump_Amp',f5.2,'_geq',f5.3,'_gam',f5.3,'_',i3.3,'.dat')
else
 if(Amp.lt.10.0d0)then
  write(unit=filesuffix,fmt=13) Ncu, Nox, Amp, gz, dimag(igam), nquanta
 else
  write(unit=filesuffix,fmt=14) Ncu, Nox, Amp, gz, dimag(igam), nquanta
 endif
13 format('_Cu',i1,'O',i2,'_Spump_Amp',f5.3,'_geq',f5.3,'_gam',f5.3,'_',i3.3,'.dat')
14 format('_Cu',i1,'O',i2,'_Spump_Amp',f5.2,'_geq',f5.3,'_gam',f5.3,'_',i3.3,'.dat')
endif

!define the DME
DME_xin = 0.0d0
DME_yin = cos(56.0d0*pi/180.0d0)

write(unit=screen,fmt=*) 'trXAS for a Cu3O8 Cluster.'
write(unit=screen,fmt=*) 'Random seed = ', seed
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
maxHUDLen = 40*NelUD
allocate(HUDi(1:maxHUDLen))
allocate(HUDj(1:maxHUDLen))
allocate(HUD(1:maxHUDLen))

!Now construct the Hamiltonian for the (Nup, Ndn) problem. 
ch = 0                  !Generate the Hamiltonian with no core hole
time = start_time       !assume time t = start_time 
displacement = X0 + XAmp*sin(omegamode*time)*exp(-time/tau)*heaviside(time)

print*, ' '
print*, 'Constructing H for (nup,ndn) problem.'
call construct_hamiltonian(HUDi,HUDj,HUD,HUDlen,maxHUDlen,basis_UD,NelUD,ch,displacement)

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
n_lanczos_states = 350
orthogonalize = .true.

call get_ground_state(states,Ei,numstates,basis_UD,hs_size_UD,NelUD,& 
                 H0i,H0j,H0,HUDlen,maxHUDlen,n_lanczos_states,orthogonalize) 

!Output the energies of the eigenstates to be sure.
print*, 'Energies of the lowest lying states: '
allocate(psi_t(1:hs_size_UD))    !Psi at time t

do i = 1,numstates
 Psi_t(:) = states(:,i) 
 call get_filling_and_spin(Psi_t,hs_size_UD,basis_UD,NelUD,Stot,SiSip1,Mz,Fillcu,&
                           Fillox,Fill) 
 print*, i, Ei(i),Stot,Mz
enddo

print*, 'Energy splitting:'
print*, 'E(2) - E(1) = ', Ei(2) - Ei(1)
print*, 'E(3) - E(2) = ', Ei(3) - Ei(1)

dtime = (end_time-start_time)/dfloat(numtime)
allocate(Psi0(1:hs_size_UD))

filename1 = "projections"//trim(filesuffix)
open(unit=nfwrk1,action='write',file=filename1,status='replace')
close(unit=nfwrk1)
filename2 = 'xas'//trim(filesuffix)
open(file=filename2,unit=nfwrk1,action='write',status='replace')
close(unit=nfwrk1)
filename3 = "occupations"//trim(filesuffix)
open(unit=nfwrk1,action='write',file=filename3,status='replace')
close(unit=nfwrk1)

!Set the lanczos parameters for the XAS calculations
!orthogonalize = .false.
!n_lanczos_states = 400
!Allocate an array to store the initial state for the XAS calculation
allocate(gs(1:hs_size_UD))

! Now do the time propagation 
counter = 0
open(file='X_of_t.dat',unit=999,action='write')
do nt = 1,numtime-1,ntskip
 time = start_time + dtime*dfloat(nt-1)
 ch = 0
 displacement = X0 + XAmp*sin(omegamode*time)*exp(-time/tau)*heaviside(time)

 call construct_hamiltonian(HUDi,HUDj,HUD,HUDlen,maxHUDlen,basis_UD,NelUD,ch,displacement)
 do i = 1,HUDlen
  H0(i) = dreal(HUD(i))
  H0i(i) = HUDi(i)
  H0j(i) = HUDj(i)
 enddo
 n_lanczos_states = 350
 orthogonalize = .true.

 !compute the ground state
 call get_ground_state(states,Ei,numstates,basis_UD,hs_size_UD,NelUD,& 
                       H0i,H0j,H0,HUDlen,maxHUDlen,n_lanczos_states,orthogonalize) 
 Psi_t(:) = states(:,1) 
 call get_filling_and_spin(Psi_t,hs_size_UD,basis_UD,NelUD,Stot,SiSip1,Mz,Fillcu,&
                           Fillox,Fill) 
 write(  6,*) time, displacement, Ei(1), stot, Mz, SiSip1
 write(999,*) time, displacement, Ei(1), stot, Mz, SiSip1
enddo

!
! !Chebyshev propagate the wavefunction forward in time.
! call cheby_propagate(psi_t(:), psi_tp(:), &
!                      HUDi, HUDj, HUD, HUDlen, maxHUDlen, basis_UD, NelUD )
!
! ! This output is actually one time slice behind
! dnorm2 = dreal( dot_product(Psi_t(:),Psi_t(:)))
!
! !Get the energy of the state  
! ! Psi0 = H * Psi_t(:,nt)
! call multiply_cvector_by_cH(Psi_t(:),Psi0,hs_size_UD,HUDi,HUDj,HUD,HUDlen,&
!                             maxHUDlen,basis_UD,NelUD)
! Ei_t(nt) = dreal(dot_product(Psi_t(:),Psi0))
!
! call get_filling_and_spin(Psi_t(:),hs_size_UD,basis_UD,NelUD,Stot,SiSip1,Mz,Fillcu,&
!                           Fillox,Fill) 
!
! !Output the time evolution of the energy (output time in fs) to the output file 
! !and the screen
! open(unit=nfwrk1,file=filename1,access='append')
! open(unit=nfwrk2,file=filename3,access='append')
! write(nfwrk1,*) time*hbar, Ei_t(nt), dreal(dot_product(states(:,1),psi_t(:))), &
!                                      dreal(dot_product(states(:,2),psi_t(:))), &
!                                      dreal(dot_product(states(:,3),psi_t(:))), &
!                                      dnorm2, &
!                                      dreal(psi_t(1)), dimag(psi_t(1))
! write(nfwrk2,*) time*hbar, Stot, Mz, SiSip1, (FillCu(i),i=1,Ncu), & 
!                                              FillOx(1),FillOx(2), Fill
! write(unit=screen,fmt=20) time,dnorm2,Ei_t(nt),Stot,SiSip1,Mz,Fill,sum(FillCu),sum(FillOx)
!20 format('t=',f6.1,' norm=',f8.6,' E(t)=',f8.5,' Stot=', f8.5,' Si*Si+1=',f8.5,&
!          ' Mz=', f8.5,' Fill=', f8.5,' NCu=',f8.5,' Nox=',f8.5)
! close(unit=nfwrk1)
! close(unit=nfwrk2)
! 
! if(counter.eq.ntskip)then
!  counter = 0
!  compute_xas = .true.
! else
!  counter = counter + 1
!  compute_xas = .false.
! endif
!
! if(compute_xas)then
!  print*, 'Computing XAS for time ', nt
!  XAS = 0.0d0      !Init the XAS intensity to zero.
! 
!  Psi0(:) = dcmplx(0.0d0,0.0d0)
!  !Set the current Psi(t) to be equal to the "ground state" of the XAS 
!  !calculation
!  gs(:) = psi_t(:) 
!  Eg = Ei_t(nt)
! 
!  !Loop over all possible core hole locations
!  do ch = 1,nox 
!   !Destroy a spin up hole in the valence band
!   spin = spinup
!   HUm1Dlen = 0
!   maxHUm1Dlen = 30*NelUm1D
!   allocate(HUm1Di(1:maxHUm1DLen))
!   allocate(HUm1Dj(1:maxHUm1DLen))
!   allocate(HUm1D(1:maxHUm1DLen))
!   allocate( cgs(1:hs_size_Um1D))
!   allocate(psi1(1:hs_size_Um1D))
! 
!   !Build the intermediate state Hamiltonian at time t
!   call construct_hamiltonian(HUm1Di,HUm1Dj,HUm1D,HUm1Dlen,maxHUm1Dlen,& 
!                              basis_Um1D,NelUm1D,ch,time)
!  !first apply the destruction operator to the gs and store it in cgs
!   psi1 = dcmplx(0.0d0,0.0d0)
!  !Oxygen y
!   site = 3*ncu+0*nox+ch
!   call remove_particle(gs,hs_size_UD,basis_UD,NelUD,& 
!                        cgs,hs_size_Um1D,basis_Um1d,NelUm1D,site,spin,N)
!   psi1 = dcmplx(DME_yin,0.0d0)*cgs
!
!   !Oxygen x
!   site = 3*ncu+1*nox+ch
!   call remove_particle(gs,hs_size_UD,basis_UD,NelUD,& 
!                        cgs,hs_size_Um1D,basis_Um1d,NelUm1D,site,spin,N)
!   psi1 = psi1 + dcmplx(DME_xin,0.0d0)*cgs
!  
!   call project_on_final(psi1,hs_size_Um1D,basis_Um1D,NelUm1D,EE,Final,& 
!                        n_lanczos_states,HUm1Di,HUm1Dj,HUm1D,HUm1Dlen,maxHUm1Dlen,& 
!                        j,orthogonalize)
! 
!   do i = 1,j
!    do nw = 0,numomega
!     omega = dfloat(nw)*(omegab-omegaa)/dfloat(numomega) + omegaa
!     XAS(nw) = XAS(nw) + dreal(final(i)*conjg(final(i)))* & 
!                  aimag(1.0d0/(EE(i)-Eg-omega-igam))
!    enddo
!   enddo
!
!   deallocate(HUm1Di)
!   deallocate(HUm1Dj)
!   deallocate(HUm1D)
!   deallocate(cgs)
!   deallocate(psi1)
!  
!   !Destroy a spin dn hole in the valence band
!   spin = spindn
!   HUDm1len = 0
!   maxHUDm1len = 30*NelUDm1
!   allocate(HUDm1i(1:maxHUDm1Len))
!   allocate(HUDm1j(1:maxHUDm1Len))
!   allocate(HUDm1(1:maxHUDm1Len))
!   allocate( cgs(1:hs_size_UDm1))
!   allocate(psi1(1:hs_size_UDm1))
! 
!   !Build the intermediate state Hamiltonian at time t
!   call construct_hamiltonian(HUDm1i,HUDm1j,HUDm1,HUDm1len,maxHUDm1len,& 
!                              basis_UDm1,NelUDm1,ch,time)
!   !first apply the destruction operator to the gs and store it in cgs
!   psi1 = dcmplx(0.0d0,0.0d0)
!   !Oxygen y
!   site = 3*ncu+0*nox+ch
!   call remove_particle(gs,hs_size_UD,basis_UD,NelUD,& 
!                        cgs,hs_size_UDm1,basis_Udm1,NelUDm1,site,spin,N)
!   psi1 = dcmplx(DME_yin,0.0d0)*cgs
! 
!   !Oxygen x
!   site = 3*ncu+1*nox+ch
!   call remove_particle(gs,hs_size_UD,basis_UD,NelUD,& 
!                        cgs,hs_size_UDm1,basis_Udm1,NelUDm1,site,spin,N)
!   psi1 = psi1 + dcmplx(DME_xin,0.0d0)*cgs
!  
!   call project_on_final(psi1,hs_size_UDm1,basis_UDm1,NelUDm1,EE,Final,& 
!                        n_lanczos_states,HUDm1i,HUDm1j,HUDm1,HUDm1len,maxHUDm1len,& 
!                        j,orthogonalize)
! 
!   do i = 1,j
!    do nw = 0,numomega
!     omega = dfloat(nw)*(omegab-omegaa)/dfloat(numomega) + omegaa
!     XAS(nw) = XAS(nw) + dreal(final(i)*conjg(final(i)))* & 
!                  aimag(1.0d0/(EE(i)-Eg-omega-igam))
!    enddo
!   enddo
! 
!   deallocate(HUDm1i)
!   deallocate(HUDm1j)
!   deallocate(HUDm1)
!   deallocate(cgs)
!   deallocate(psi1)
!  enddo !loop over core-hole positions
! 
!  !and output the XAS spectra at this time to a file
!  open(unit=nfwrk2,file=filename2,access='append')
!  do nw = 0,numomega
!   omega = dfloat(nw)*(omegab-omegaa)/dfloat(numomega) + omegaa
!   write(nfwrk2,*) time*hbar, omega, XAS(nw)
!  enddo
!  write(nfwrk2,*) ' '
!  close(unit=nfwrk2)
! endif !If(compute_xas)then
!enddo !do nt = 0,numtime
close(unit=999)

stop
end program trXAS

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

