!==============================================================================
! module mod_hamiltonian defines the routines needed to construct the 
! hamiltonian and multiple by it, which is needed by the Lanczos routines. 
!==============================================================================
module mod_hamiltonian
 use mod_parameters, only: N
 double precision, dimension(1:N) :: esite
 double precision, dimension(1:N,1:N) :: hopping
 double precision, dimension(1:N,1:N,1:2) :: Rij 
 double precision, dimension(1:3,1:3) :: Uprime
 double precision, dimension(1:3,1:3) :: Jd 
contains
 !=============================================================================
 subroutine init_cluster()
 use mod_parameters
 implicit none
 integer(KIND=8) i, j, itl, ibl, itr, ibr, orb

 Uprime = 0.0d0
 jd = 0.0d0
 hopping = 0.0d0
 include 'hopping.f90'
 include 'coulomb.f90'
 include 'rij.f90'

 !check that the hopping table is hermitian
 do i = 1,N
  do j = 1,N
   if(hopping(i,j).ne.hopping(j,i))then
    print*, 'Error: Hopping table is not hermitian.'
    print*, 'Element: ', i, j
    stop
   endif
   if((rij(i,j,1).ne.-rij(j,i,1)).or. & 
      (rij(i,j,1).ne.-rij(j,i,1)))then
    print*, 'Error: Displacement table is not correct.'
    print*, 'Element: ', i, j
    stop
   endif
  enddo
 enddo

 do i = 1,ncu
  esite(i) = edxy
  esite(i+ncu) = edx2_y2
  esite(i+2*ncu) = edz2_r2
 enddo
 do i = 3*Ncu+1,3*Ncu+nox
  esite(i) = epy
  esite(i+nox) = epx
 enddo

 esite(1+0*ncu) = esite(1+0*ncu) + deltaE_Cu !Cu
 esite(1+1*ncu) = esite(1+1*ncu) + deltaE_Cu !Cu
 esite(1+2*ncu) = esite(1+2*ncu) + deltaE_cu !Cu

 esite(  ncu) = esite(  ncu) + deltaE_cu
 esite(2*ncu) = esite(2*ncu) + deltaE_cu
 esite(3*ncu) = esite(3*ncu) + deltaE_cu

 esite(3*Ncu    +1) = esite(3*Ncu + 1) + deltaE_O
 esite(3*Ncu+nox+1) = esite(3*Ncu + nox + 1) + deltaE_O
 esite(3*Ncu    +nox/2) = esite(3*Ncu    +nox/2) + deltaE_O
 esite(3*Ncu+nox+nox/2) = esite(3*Ncu+nox+nox/2) + deltaE_O
 esite(3*Ncu    +nox/2+1) = esite(3*Ncu    +nox/2+1) + deltaE_O
 esite(3*Ncu+nox+nox/2+1) = esite(3*Ncu+nox+nox/2+1) + deltaE_O
 esite(3*Ncu+  Nox) = esite(3*Ncu+  nox) + deltaE_O
 esite(3*Ncu+2*nox) = esite(3*Ncu+2*nox) + deltaE_O

 return
 end subroutine init_cluster

 subroutine get_vec_potential(A,time)
 use mod_parameters, only: amp, pulsewidth, pump_freq, eps_pumpAx, eps_pumpAy
 implicit none
 double precision sigma, time, pi
 double precision, dimension(1:2) :: A
 pi = 2.0d0*asin(1.0d0)
 sigma = pulsewidth/(2.355d0)   !pulse width is the FWHM so divide by 2.355 to get sigma
 A(1) = eps_pumpAx*Amp*exp(-time*time/(2.0d0*sigma*sigma))/sqrt(2*pi*sigma*sigma)
 A(2) = eps_pumpAy*Amp*exp(-time*time/(2.0d0*sigma*sigma))/sqrt(2*pi*sigma*sigma)
 !Note time is actually t/hbar and pump_freq is in units of eV
 A = A*cos(time*pump_freq)
 return
 end subroutine get_vec_potential
 !=============================================================================
 ! subroutine construct_hamiltonian
 !=============================================================================
 subroutine construct_hamiltonian(Hi,Hj,Hel,Hlen,maxHLen,basis,hs_size,ch,disp)
 use mod_parameters
 use mod_hilbert_space, only: return_idx_for_state
 implicit none
 integer(KIND=8) spin, i, j, newstate, ii, jj
 integer(KIND=8) spinp, orb1, orb2, idx1, idx2
 integer(KIND=8) nn, state
 integer(KIND=8) c1, c2
 integer(KIND=8) Hlen, maxHLen, hs_size, ch
 integer(KIND=8), dimension(1:hs_size) :: basis
 integer(KIND=8), dimension(1:maxHLen) :: Hi, Hj
 double precision disp, jsp
 double precision, dimension(1:N) :: E
 double precision rtmp, rtmp2, fsgn
 double precision, dimension(1:2) :: A  !The complex part due to the gauge field
 double complex, parameter :: eye = (0.0d0,1.0d0)
 double complex ctmp
 double complex, dimension(1:maxHlen) :: Hel

 !call get_vec_potential(A,time)

 !Init everything to zero 
 Hlen = 0
 do i = 1,maxHlen
  Hi(i) = 0
  Hj(i) = 0
  Hel(i) = dcmplx(0.0d0,0.0d0)
 enddo
 
 !set the site energies
 E = esite

 !if ch is not zero then we need to add the core-hole potential.
 if(ch.ne.0)then
  E(3*ncu + 0*nox + ch) = E(3*ncu + 0*nox + ch) + Uq
  E(3*ncu + 1*nox + ch) = E(3*ncu + 1*nox + ch) + Uq
 endif

 do orb1 = 0,2
  E(1+orb1*ncu) = E(1+orb1*ncu) + decu*disp
  E(2+orb1*ncu) = E(2+orb1*ncu) + decu*disp
  E(3+orb1*ncu) = E(3+orb1*ncu) + decu*disp
 enddo

 !loop over basis states and begin constructing the hamiltonian
 do nn = 1,hs_size
  state = basis(nn)
  rtmp = 0.0d0

  !First apply the site energies
  do j = 1,N
   do spin = 0,1
    if(btest(state,j-1+spin*N))then
     rtmp = rtmp + E(j)
    endif
   enddo
  enddo

  !On-site Cu interactions
  do j = 1,3*Ncu
   if(btest(state,j-1).and.btest(state,j+N-1).and.Ud.ne.0.0d0)then
    rtmp = rtmp + Ud
   endif
  enddo

  !oxygen on-site interactions
  do j = 3*Ncu+1,N
   if(btest(state,j-1).and.btest(state,j+N-1).and.Up.ne.0.0d0)then
    rtmp = rtmp + Up
   endif
  enddo

  !apply the U' term for copper
  do i = 1,Ncu
   do spin = 0,1
    do spinp = 0,1
     do orb1 = 0,2
      do orb2 = 0,2
       if(orb1.ne.orb2)then
        idx1 = i + ncu*orb1 + spin*N
        idx2 = i + ncu*orb2 + spinp*N
        if(btest(state,idx1-1).and.btest(state,idx2-1))then
         rtmp = rtmp + Uprime(orb1+1,orb2+1)/2.0d0
        endif
       endif
      enddo
     enddo
    enddo
   enddo
  enddo

  !U' for the oxygen terms
  do i = 3*Ncu+1,3*Ncu+nox   !sum over sites
   do spin = 0,1
    do spinp = 0,1
     do orb1 = 0,1
      do orb2 = 0,1
       if(orb1.ne.orb2)then
        idx1 = i + nox*orb1 + spin*N
        idx2 = i + nox*orb2 + spinp*N
        if(btest(state,idx1-1).and.btest(state,idx2-1))then
         rtmp = rtmp + (Up-2.0d0*Jp)/2.0d0
        endif
       endif
      enddo
     enddo
    enddo
   enddo
  enddo

  !copper-copper interactions
  do i = 1,Ncu-1
   do spin = 0,1
    do spinp = 0,1
     do orb1 = 0,2
      do orb2 = 0,2
       idx1 = i + orb1*Ncu + spin*N
       idx2 = i + 1 + orb2*Ncu + spinp*N
       if(btest(state,idx1-1).and.btest(state,idx2-1))then
        rtmp = rtmp + Vdd
       endif
      enddo
     enddo
    enddo
   enddo
  enddo

  !finally, sum over the Upd terms
  rtmp2 = 0.0d0
  do i = 1,Ncu
   do spin = 0,1
    do spinp = 0,1
     do orb1 = 0,2 !copper orbitals
      do orb2 = 0,1 !oxygen orbitals
       idx1 = i + Ncu*orb1 + spin*N
       !oxygen in the bottom left
       idx2 = i + 3*Ncu + orb2*nox + spinp*N
       if(btest(state,idx1-1).and.btest(state,idx2-1))then
        rtmp2 = rtmp2 + Upd
       endif

       !oxygen in the top left
       idx2 = (i+Ncu+1) + 3*Ncu + orb2*nox + spinp*N
       if(btest(state,idx1-1).and.btest(state,idx2-1))then
        rtmp2 = rtmp2 + Upd
       endif

       !oxygen bottom right
       idx2 = i+1 + 3*Ncu + orb2*nox + spinp*N
       if(btest(state,idx1-1).and.btest(state,idx2-1))then
        rtmp2 = rtmp2 + Upd
       endif

       !oxygen in the top right
       idx2 = (i + Ncu + 2) + 3*Ncu + orb2*nox + spinp*N
       if(btest(state,idx1-1).and.btest(state,idx2-1))then
        rtmp2 = rtmp2 + Upd
       endif

      enddo
     enddo
    enddo
   enddo
  enddo
  ctmp = dcmplx(rtmp + rtmp2,0.0d0)
  !and insert the diagonal element to the hamiltonian array.
  if(rtmp.ne.0.0d0)then
   call insert_element(Hi,Hj,Hel,nn,nn,ctmp,maxHLen,Hlen)
  endif
 
  !That is the end of the diagonal terms.  
  !Now we have to do all the off-diagonal parts.
  do i = 1,N  !loop over starting sites
   do j = 1,N  !loop over end sites
    if(hopping(i,j).ne.0.0d0)then
     do spin = 0,1
      call apply_hop(state,basis,hs_size,N,spin,i,j,hopping(i,j),ii,jj,rtmp)
      if(ii.gt.0.and.jj.gt.0)then  !if ii and jj are > 0 this is a valid hop
       ctmp = dcmplx(rtmp,0.0d0)
       call insert_element(Hi,Hj,Hel,ii,jj,ctmp,maxHLen,Hlen)
      endif
     enddo
    endif
   enddo ! do j = 1,N
  enddo !do i = 1,N

  !=============================================================
  ! Apply the Hunds (J) coupling terms
  !=============================================================
  !Copper atoms first
  do i = 1,Ncu
   do spin = 0,1
    do spinp = 0,1
     do orb1 = 0,2
      do orb2 = 0,2
       if(orb1.ne.orb2)then
        c1 = i + orb1*Ncu
        c2 = i + orb2*Ncu
        call apply_hunds(c1,c2,spin,spinp,state,newstate,N,fsgn)
        if(newstate.gt.0)then
         call return_idx_for_state(state,ii,basis,hs_size)
         call return_idx_for_state(newstate,jj,basis,hs_size)
         ctmp = dcmplx(fsgn*jd(orb1+1,orb2+1)*0.5d0,0.0d0)
         call insert_element(Hi,Hj,Hel,ii,jj,ctmp,maxHLen,Hlen)
        endif
       endif
      enddo
     enddo
    enddo
   enddo
  enddo
 
  !coupling for the oxgyen
  do i = 3*ncu+1,3*ncu+nox
   do spin = 0,1
    do spinp = 0,1
     do orb1 = 0,1
      do orb2 = 0,1
       if(orb1.ne.orb2)then
        c1 = i+orb1*Nox
        c2 = i+orb2*nox
        call apply_hunds(c1,c2,spin,spinp,state,newstate,N,fsgn)
        if(newstate.gt.0)then
         call return_idx_for_state(state,ii,basis,hs_size)
         call return_idx_for_state(newstate,jj,basis,hs_size)
         ctmp = dcmplx(fsgn*jp*0.5d0,0.0d0)
         call insert_element(Hi,Hj,Hel,ii,jj,ctmp,maxHLen,Hlen)
        endif
       endif
      enddo
     enddo
    enddo
   enddo
  enddo

  !=============================================================
  ! Apply the Hunds (J') coupling terms
  !=============================================================
  !Copper
  do i = 1,Ncu
   do spin = 0,1
    if(spin.eq.1)then
     spinp = 0
    else
     spinp = 1
    endif

    do orb1 = 0,2
     do orb2 = 0,2
      if(orb1.ne.orb2)then
       c1 = i+ncu*orb1
       c2 = i+ncu*orb2
       call apply_hunds_2(c1, c2,spin,spinp,state,newstate,N,fsgn)
       if(newstate.gt.0)then
        call return_idx_for_state(state,ii,basis,hs_size)
        call return_idx_for_state(newstate,jj,basis,hs_size)
        ctmp = dcmplx(fsgn*jd(orb1+1,orb2+1)*0.5d0,0.0d0)
        call insert_element(Hi,Hj,Hel,ii,jj,ctmp,maxHLen,Hlen)
       endif
      endif
     enddo
    enddo
   enddo
  enddo

  !oxygen 
  do i = 3*Ncu+1,3*Ncu+nox
   do spin = 0,1
    if(spin.eq.1)then
     spinp = 0
    else
     spinp = 1
    endif
    do orb1 = 0,1
     do orb2 = 0,1
      if(orb1.ne.orb2)then
       c1 = i+nox*orb1
       c2 = i+nox*orb2
       call apply_hunds_2(c1, c2,spin,spinp,state,newstate,N,fsgn)
       if(newstate.gt.0)then
        call return_idx_for_state(state,ii,basis,hs_size)
        call return_idx_for_state(newstate,jj,basis,hs_size)
        ctmp = dcmplx(fsgn*jp*0.5d0,0.0d0)
        call insert_element(Hi,Hj,Hel,ii,jj,ctmp,maxHLen,Hlen)
       endif
      endif
     enddo
    enddo
   enddo
  enddo

  !now the J_{pd} interactions by looping over copper sites
  do i = 1,Ncu
   !sum over sigma and sigmap
   do spin = 0,1
    do spinp = 0,1
     do orb1 = 0,2  !copper orbitals
      do orb2 = 0,1 !oxygen orbitals
 
       idx1 = i + Ncu*orb1
       !oxygen in the bottom left
       idx2 = i + 3*Ncu + orb2*nox

       call apply_hunds(idx1,idx2,spin,spinp,state,newstate,N,fsgn)
       if(newstate.gt.0)then
        call return_idx_for_state(state,ii,basis,hs_size)
        call return_idx_for_state(newstate,jj,basis,hs_size)
        ctmp = dcmplx(fsgn*jpd,0.0d0)
        call insert_element(Hi,Hj,Hel,ii,jj,ctmp,maxHLen,Hlen)
       endif

       !oxygen in the top left
       idx2 = (i+Ncu+1) + 3*Ncu + orb2*nox
       call apply_hunds(idx1,idx2,spin,spinp,state,newstate,N,fsgn)
       if(newstate.gt.0)then
        call return_idx_for_state(state,ii,basis,hs_size)
        call return_idx_for_state(newstate,jj,basis,hs_size)
        ctmp = dcmplx(fsgn*jpd,0.0d0)
        call insert_element(Hi,Hj,Hel,ii,jj,ctmp,maxHLen,Hlen)
       endif

       !oxygen in the bottom right
       idx2 = i+1 + 3*Ncu + orb2*nox
       call apply_hunds(idx1,idx2,spin,spinp,state,newstate,N,fsgn)
       if(newstate.gt.0)then
        call return_idx_for_state(state,ii,basis,hs_size)
        call return_idx_for_state(newstate,jj,basis,hs_size)
        ctmp = dcmplx(fsgn*jpd,0.0d0)
        call insert_element(Hi,Hj,Hel,ii,jj,ctmp,maxHLen,Hlen)
       endif

       !oxygen in the top right
       idx2 = (i + Ncu + 2) + 3*Ncu + orb2*nox
       call apply_hunds(idx1,idx2,spin,spinp,state,newstate,N,fsgn)
       if(newstate.gt.0)then
        call return_idx_for_state(state,ii,basis,hs_size)
        call return_idx_for_state(newstate,jj,basis,hs_size)
        ctmp = dcmplx(fsgn*jpd,0.0d0)
        call insert_element(Hi,Hj,Hel,ii,jj,ctmp,maxHLen,Hlen)
       endif
      enddo
     enddo
    enddo
   enddo
  enddo

  !finally we apply the second Hund's term for Jpd
  do i = 1,Ncu
   !sum over sigma and sigmap
   do spin = 0,1
    !set the value of spinp
    if(spin.eq.0)then
     spinp = 1
    elseif(spin.eq.1)then
     spinp = 0
    endif

    do orb1 = 0,2  !copper orbitals
     do orb2 = 0,1 !oxygen orbitals
      idx1 = i + Ncu*orb1
      !oxygen in the bottom left
      idx2 = i + 3*Ncu + orb2*nox
      call apply_hunds_2(idx1,idx2,spin,spinp,state,newstate,N,fsgn)
      if(newstate.gt.0)then
       call return_idx_for_state(state,ii,basis,hs_size)
       call return_idx_for_state(newstate,jj,basis,hs_size)
       ctmp = dcmplx(fsgn*jpd*0.5d0,0.0d0)
       call insert_element(Hi,Hj,Hel,ii,jj,ctmp,maxHLen,Hlen)
      endif

      call apply_hunds_2(idx2,idx1,spin,spinp,state,newstate,N,fsgn)
      if(newstate.gt.0)then
       call return_idx_for_state(state,ii,basis,hs_size)
       call return_idx_for_state(newstate,jj,basis,hs_size)
       ctmp = dcmplx(fsgn*jpd*0.5d0,0.0d0)
       call insert_element(Hi,Hj,Hel,ii,jj,ctmp,maxHLen,Hlen)
      endif

      !oxygen in the top left
      idx2 = (i+Ncu+1) + 3*Ncu + orb2*nox
      call apply_hunds_2(idx1,idx2,spin,spinp,state,newstate,N,fsgn)
      if(newstate.gt.0)then
       call return_idx_for_state(state,ii,basis,hs_size)
       call return_idx_for_state(newstate,jj,basis,hs_size)
       ctmp = dcmplx(fsgn*jpd*0.5d0,0.0d0)
       call insert_element(Hi,Hj,Hel,ii,jj,ctmp,maxHLen,Hlen)
      endif

      call apply_hunds_2(idx2,idx1,spin,spinp,state,newstate,N,fsgn)
      if(newstate.gt.0)then
       call return_idx_for_state(state,ii,basis,hs_size)
       call return_idx_for_state(newstate,jj,basis,hs_size)
       ctmp = dcmplx(fsgn*jpd*0.5d0,0.0d0)
       call insert_element(Hi,Hj,Hel,ii,jj,ctmp,maxHLen,Hlen)
      endif
      !oxygen in the bottom right
      idx2 = i+1 + 3*Ncu + orb2*nox
      call apply_hunds_2(idx1,idx2,spin,spinp,state,newstate,N,fsgn)
      if(newstate.gt.0)then
       call return_idx_for_state(state,ii,basis,hs_size)
       call return_idx_for_state(newstate,jj,basis,hs_size)
       ctmp = dcmplx(fsgn*jpd*0.5d0,0.0d0)
       call insert_element(Hi,Hj,Hel,ii,jj,ctmp,maxHLen,Hlen)
      endif

      call apply_hunds_2(idx2,idx1,spin,spinp,state,newstate,N,fsgn)
      if(newstate.gt.0)then
       call return_idx_for_state(state,ii,basis,hs_size)
       call return_idx_for_state(newstate,jj,basis,hs_size)
       ctmp = dcmplx(fsgn*jpd*0.5d0,0.0d0)
       call insert_element(Hi,Hj,Hel,ii,jj,ctmp,maxHLen,Hlen)
      endif

      !oxygen in the top right
      idx2 = (i + Ncu + 2) + 3*Ncu + orb2*nox
      call apply_hunds_2(idx1,idx2,spin,spinp,state,newstate,N,fsgn)
      if(newstate.gt.0)then
       call return_idx_for_state(state,ii,basis,hs_size)
       call return_idx_for_state(newstate,jj,basis,hs_size)
       ctmp = dcmplx(fsgn*jpd*0.5d0,0.0d0)
       call insert_element(Hi,Hj,Hel,ii,jj,ctmp,maxHLen,Hlen)
      endif

      call apply_hunds_2(idx2,idx1,spin,spinp,state,newstate,N,fsgn)
      if(newstate.gt.0)then
       call return_idx_for_state(state,ii,basis,hs_size)
       call return_idx_for_state(newstate,jj,basis,hs_size)
       ctmp = dcmplx(fsgn*jpd*0.5d0,0.0d0)
       call insert_element(Hi,Hj,Hel,ii,jj,ctmp,maxHLen,Hlen)
      endif
     enddo
    enddo
   enddo
  enddo
  
  !now spin - peierls term
  jsp =Jdimer*disp
  do i = 2,Ncu
   jsp = -jsp
   !sum over sigma and sigmap
   do spin = 0,1
    do spinp = 0,1
     do orb1 = 0,0  !copper orbitals
      do orb2 = 0,0 !oxygen orbitals 
       idx1 = i + Ncu*orb1
       !oxygen in the bottom left
       idx2 = i - 1 + Ncu*orb2
       call apply_hunds(idx1,idx2,spin,spinp,state,newstate,N,fsgn)
       if(newstate.gt.0)then
        call return_idx_for_state(state,ii,basis,hs_size)
        call return_idx_for_state(newstate,jj,basis,hs_size)
        ctmp = dcmplx(fsgn*jsp,0.0d0)
        call insert_element(Hi,Hj,Hel,ii,jj,ctmp,maxHLen,Hlen)
       endif
      enddo
     enddo
    enddo
   enddo
  enddo

  !finally, apply the second spin-peierls
  jsp =Jdimer*disp
  do i = 2,Ncu
   jsp =  -jsp
   !sum over sigma and sigmap
   do spin = 0,1
    !set the value of spinp
    if(spin.eq.0)then
     spinp = 1
    elseif(spin.eq.1)then
     spinp = 0
    endif

    do orb1 = 0,0  !copper orbitals
     do orb2 = 0,0 !oxygen orbitals
      idx1 = i + Ncu*orb1
      !oxygen in the bottom left
      idx2 = i - 1 + Ncu*orb2
      call apply_hunds_2(idx1,idx2,spin,spinp,state,newstate,N,fsgn)
      if(newstate.gt.0)then
       call return_idx_for_state(state,ii,basis,hs_size)
       call return_idx_for_state(newstate,jj,basis,hs_size)
       ctmp = dcmplx(fsgn*jsp*0.5d0,0.0d0)
       call insert_element(Hi,Hj,Hel,ii,jj,ctmp,maxHLen,Hlen)
      endif
     enddo
    enddo
   enddo
  enddo
 enddo !nn = 1,hs_size
 return
 end subroutine construct_hamiltonian
 !==============================================================================
 ! Subroutine insert_element
 !==============================================================================
 subroutine insert_element(Hi,Hj,H,i,j,value,siz,point)
 implicit none
 logical found, loop
 integer(KIND=8) i, j, siz, point, ii
 double complex, dimension(1:siz) :: H
 integer(KIND=8), dimension(1:siz) :: Hi, Hj
 double complex value
 !
 !first we need to check if this term already has an entry.  We only need to
 !search the last few terms as all "i"th elements are grouped together.
 found = .false.
 ii = point      !set the counter to the current array
 if(ii.gt.0)then
  loop =.true.
  do while(loop)
   if(Hi(ii).eq.i.and.Hj(ii).eq.j)then
    found = .true.
    loop = .false.
   else
    ii = ii - 1
    if(ii.eq.0)then
     loop = .false.
    else
     if(Hi(ii).ne.i) loop = .false.
    endif
   endif
  enddo
 endif

 !if found is still false then we increment "point" and add the new matrix
 !element to the list
 if(.not.found)then
  point = point + 1
  if(point.ge.siz)then
   print*, 'Array is too small for storage of H, please increase the size.'
   print*, 'point = ', point, 'siz = ', siz
   stop
  else
   Hi(point) = i
   Hj(point) = j
   H(point) = value
  endif
 else !if found is true then increment the matrix element we found in the list.
  H(ii) = H(ii) + value
 endif
 return
 end subroutine insert_element
 
 !==============================================================================
 ! Computes the Fermion Sign for an operator
 !==============================================================================
 double precision function sgn(creat,annil,spin,N,state)
 implicit none
 integer(KIND=8) creat, annil, spin, idx, N, state, itmp, j
 itmp = state
 sgn = 1.0d0
 idx = N*spin
 do j = N+idx,annil+1+idx,-1
  if(btest(state,j-1)) sgn = -sgn
 enddo
 state = ibclr(state,annil+idx-1)
 do j = N+idx,creat+1+idx,-1
  if(btest(state,j-1)) sgn = -sgn
 enddo
 state = itmp
 return
 end function sgn

 !==============================================================================
 subroutine apply_hop(state,psi,hs_size,N,spin,annil,creat,t,ii,jj,rtmp)
 use mod_hilbert_space, only: return_idx_for_state
 implicit none
 integer(KIND=8) newstate
 integer(KIND=8) state, hs_size, annil, creat, ii, jj, spin, N
 integer(KIND=8), dimension(1:hs_size) :: psi
 double precision t, rtmp

 rtmp = 0.0d0
 call return_idx_for_state(state,ii,psi,hs_size)
 if(btest(state,annil+spin*N-1).and..not.btest(state,creat+spin*N-1) &
                                    .and.t.ne.0.0d0)then
  newstate = ibclr(state,annil+spin*N-1)
  newstate = ibset(newstate,creat+spin*N-1)
  call return_idx_for_state(newstate,jj,psi,hs_size)
  rtmp = t*(sgn(creat,annil,spin,N,state))
 else
  jj = -1
 endif
 return
 end subroutine apply_hop
 !==============================================================================
 ! function double precision sgn_hunds
 ! This routine calculates the sign for the hunds coupling term (J)
 ! H = (J/2) \sum_{sigma,sigmap,gamma.ne.gammap}
 !            d^\dagger_{i,gamma,sigma}d^\dagger_{i,gammap,sigmap} * 
 !            d_{i,gamma,sigmap}d_{i,gammap,sigma}
 !==============================================================================
 subroutine apply_hunds(gam,gamp,spin,spinp,state,newstate,N,fsgn)
 implicit none
 integer(KIND=8) tmpstate, gam, gamp, spin, spinp,state,newstate,N
 integer(KIND=8) k
 double precision fsgn
 !default the sign to 1 and new state to zero
 fsgn = 1.0d0
 newstate = 0
 tmpstate = state

 !okay, first destroy the electron at site (gamp,spin) 
 if(btest(tmpstate,gamp+spin*N-1))then
  do k = 2*N,gamp+spin*N+1,-1
   if(btest(tmpstate,k-1)) fsgn = -fsgn
  enddo
  tmpstate = ibclr(tmpstate,gamp+spin*N-1)
 else
  newstate = 0
  return
 endif

 !Next, destroy the electron at site (gam,spinp) 
 if(btest(tmpstate,gam+spinp*N-1))then
  do k = 2*N,gam+spinp*N+1,-1
   if(btest(tmpstate,k-1)) fsgn = -fsgn
  enddo
  tmpstate = ibclr(tmpstate,gam+spinp*N-1)
 else
  newstate = 0
  return
 endif
 !now create the electron at site (gamp,spinp)
 if(.not.btest(tmpstate,gamp+spinp*N-1))then
  do k = 2*N,gamp+spinp*N+1,-1
   if(btest(tmpstate,k-1)) fsgn = -fsgn
  enddo
  tmpstate = ibset(tmpstate,gamp+spinp*N-1)
 else
  newstate = 0
  return
 endif

 !now create the electron at site (gam,spin)
 if(.not.btest(tmpstate,gam+spin*N-1))then
  do k = 2*N,gam+spin*N+1,-1
   if(btest(tmpstate,k-1)) fsgn = -fsgn
  enddo
  tmpstate = ibset(tmpstate,gam+spin*N-1)
 else
  newstate = 0
  return
 endif
 newstate = tmpstate
 return
 end subroutine apply_hunds

 !==============================================================================
 subroutine apply_hunds_2(gam,gamp,spin,spinp,state,newstate,N,fsgn)
 implicit none
 integer(KIND=8) tmpstate, gam, gamp, spin, spinp,state,newstate,N
 integer(KIND=8) k
 double precision fsgn
 !default the sign to 1 and new state to zero
 fsgn = 1.0d0
 newstate = 0
 tmpstate = state
 !okay, first destroy the electron at site (gamp,spin) 
 if(btest(tmpstate,gamp+spin*N-1))then
  do k = 2*N,gamp+spin*N+1,-1
   if(btest(tmpstate,k-1)) fsgn = -fsgn
  enddo
  tmpstate = ibclr(tmpstate,gamp+spin*N-1)
 else
  newstate = 0
  return
 endif

 !Next, destroy the electron at site (gamp,spinp) 
 if(btest(tmpstate,gamp+spinp*N-1))then
  do k = 2*N,gamp+spinp*N+1,-1
   if(btest(tmpstate,k-1)) fsgn = -fsgn
  enddo
  tmpstate = ibclr(tmpstate,gamp+spinp*N-1)
 else
  newstate = 0
  return
 endif

 !now create the electron at site (gamp,spinp)
 if(.not.btest(tmpstate,gam+spinp*N-1))then
  do k = 2*N,gam+spinp*N+1,-1
   if(btest(tmpstate,k-1)) fsgn = -fsgn
  enddo
  tmpstate = ibset(tmpstate,gam+spinp*N-1)
 else
  newstate = 0
  return
 endif

 !now create the electron at site (gam,spin)
 if(.not.btest(tmpstate,gam+spin*N-1))then
  do k = 2*N,gam+spin*N+1,-1
   if(btest(tmpstate,k-1)) fsgn = -fsgn
  enddo
  tmpstate = ibset(tmpstate,gam+spin*N-1)
 else
  newstate = 0
  return
 endif
 newstate = tmpstate
 return
 end subroutine apply_hunds_2
 !============================================================================= 
 subroutine get_filling_and_spin(Psi,Nhs,basis,Nbasis,Stot,SiSip1,Mz,Fillcu,Fillox,Fill)
 use mod_hilbert_space, only: ph_sector, return_idx_for_state
 use mod_parameters
 logical iszero
 integer(KIND=8) i, j, nn, mm, k, idx, orb, site, spin, state, newstate
 integer(KIND=8) Nhs, Nbasis
 integer(KIND=8), dimension(1:Nbasis) :: basis
 double complex, dimension(1:Nhs) :: Psi, Psip, Psip2
 double precision Stot, Mz, Fill, Niup, Nidn, Njip, Njdn, SiSip1
 double precision fsgn
 double precision, dimension(1:Ncu) :: FillCu
 double precision, dimension(1:nox) :: FillOx
 Stot = 0.0d0
 Mz = 0.0d0
 FillCu = 0.0d0
 FillOx = 0.0d0
 Fill = 0.0d0
 SiSip1 = 0.0d0
 Psip = 0.0d0
 Psip2= 0.0d0
 !loop over the electronic space
 do nn = 1,Nbasis
  do k = 0,nquanta
   idx = nn + ph_sector(k)*Nbasis
   state = Basis(nn)
   !measure the Cu occupancy  
   do i = 1,Ncu
    do orb = 0,2
     site = i + ncu*orb
     do spin = 0,1
      if(btest(state,site+spin*N-1))then
       FillCu(orb+1) = FillCu(orb+1) + conjg(Psi(idx))*Psi(idx)
      endif
     enddo
    enddo
   enddo

   !Now measure the Oxygen Occupancy
   do i = 1,Nox 
    do orb = 0,1
     site = i + 3*Ncu + nox*orb 
     do spin = 0,1
      if(btest(state,site+spin*N-1))then
       FillOx(orb+1) = FillOx(orb+1) + conjg(Psi(idx))*Psi(idx)
      endif
     enddo
    enddo
   enddo
   
   !Now measure the total Occupancy
   do i = 1,N 
    site = i  
    do spin = 0,1
     if(btest(state,site+spin*N-1))then
      Fill = Fill + conjg(Psi(idx))*Psi(idx)
     endif
    enddo
   enddo

   !Now measure the total Spin and Mz quantum numbers 
   do i = 1,N 
    Niup = 0.0d0
    Nidn = 0.0d0
    if(btest(state,i+N-1)) Niup = 1.0d0
    if(btest(state,i  -1)) Nidn = 1.0d0
    Mz = Mz + 0.5d0*(Niup-Nidn)*conjg(Psi(idx))*Psi(idx)

    do j = 1,N
     !First, compute the z-components
     Njup = 0.0d0
     Njdn = 0.0d0
     if(btest(state,j+N-1)) Njup = 1.0d0
     if(btest(state,j  -1)) Njdn = 1.0d0
     psip(idx) = psip(idx) + 0.25d0*psi(idx)*(Niup-Nidn)*(Njup-Njdn)

     !now apply the Siy*Sjy terms
     newstate = state
     iszero = .false.
     fsgn = 1.0d0
     call apply_op(newstate,fsgn,j,spindn,annhil,iszero)
     call apply_op(newstate,fsgn,j,spinup,create,iszero)
     call apply_op(newstate,fsgn,i,spinup,annhil,iszero)
     call apply_op(newstate,fsgn,i,spindn,create,iszero)
     if(.not.iszero)then
      call return_idx_for_state(newstate,mm,basis,nbasis)
      psip(mm+ph_sector(k)*nbasis) = psip(mm+ph_sector(k)*nbasis) + 0.5d0*fsgn*psi(idx)
     endif
       
     !now apply the Six*Sjx terms
     newstate = state
     iszero = .false.
     fsgn = 1.0d0
     call apply_op(newstate,fsgn,j,spinup,annhil,iszero)
     call apply_op(newstate,fsgn,j,spindn,create,iszero)
     call apply_op(newstate,fsgn,i,spindn,annhil,iszero)
     call apply_op(newstate,fsgn,i,spinup,create,iszero)
     if(.not.iszero)then
      call return_idx_for_state(newstate,mm,basis,nbasis)
      psip(mm+ph_sector(k)*nbasis) = psip(mm+ph_sector(k)*nbasis) + 0.5d0*fsgn*psi(idx)
     endif

    enddo !j = 1,N
   enddo !i = 1,N
   
   !Finally, measure the S_i*S_{i+1} correlations
   do i = 1,Ncu-1 
    j = i + 1
    Niup = 0.0d0
    Nidn = 0.0d0
    Njup = 0.0d0
    Njdn = 0.0d0
    if(btest(state,i+N-1)) Niup = 1.0d0
    if(btest(state,i  -1)) Nidn = 1.0d0
    if(btest(state,j+N-1)) Njup = 1.0d0
    if(btest(state,j  -1)) Njdn = 1.0d0
    psip2(idx) = psip2(idx) + 0.25d0*psi(idx)*(Niup-Nidn)*(Njup-Njdn)

    !now apply the Siy*Sjy terms
    newstate = state
    iszero = .false.
    fsgn = 1.0d0
    call apply_op(newstate,fsgn,j,spindn,annhil,iszero)
    call apply_op(newstate,fsgn,j,spinup,create,iszero)
    call apply_op(newstate,fsgn,i,spinup,annhil,iszero)
    call apply_op(newstate,fsgn,i,spindn,create,iszero)
    if(.not.iszero)then
     call return_idx_for_state(newstate,mm,basis,nbasis)
     psip2(mm+ph_sector(k)*nbasis) = psip2(mm+ph_sector(k)*nbasis) + 0.5d0*fsgn*psi(idx)
    endif
       
    !now apply the Six*Sjx terms
    newstate = state
    iszero = .false.
    fsgn = 1.0d0
    call apply_op(newstate,fsgn,j,spinup,annhil,iszero)
    call apply_op(newstate,fsgn,j,spindn,create,iszero)
    call apply_op(newstate,fsgn,i,spindn,annhil,iszero)
    call apply_op(newstate,fsgn,i,spinup,create,iszero)
    if(.not.iszero)then
     call return_idx_for_state(newstate,mm,basis,nbasis)
     psip2(mm+ph_sector(k)*nbasis) = psip2(mm+ph_sector(k)*nbasis) + 0.5d0*fsgn*psi(idx)
    endif
   enddo !i = 1,Ncu-1
  enddo
 enddo

 SiSip1 = dot_product(psi,psip2)
 Stot = dot_product(psi,psip)
 Stot = 0.5d0*(sqrt(1+4.0d0*Stot)-1.0d0)

 return
 end subroutine get_filling_and_spin
 !============================================================================= 
 subroutine apply_op(state,fsgn,site,spin,operation,iszero)
 use mod_parameters, only: create, annhil, N
 implicit none
 logical iszero
 double precision fsgn
 integer(KIND=8) k
 integer(KIND=8) state, site, spin,operation
 
 !If the iszero flag is set, then we do nothing
 if(.not.iszero)then
  if(operation.eq.create)then
   !now create the electron at site (gamp,spinp)
   if(.not.btest(state,site+spin*N-1))then
    do k = 2*N,site+spin*N+1,-1
     if(btest(state,k-1)) fsgn = -fsgn
    enddo
    state = ibset(state,site+spin*N-1)
    iszero = .false. 
   else
    iszero=.true.
    return
   endif
  elseif(operation.eq.annhil)then
   if(btest(state,site+spin*N-1))then
    do k = 2*N,site+spin*N+1,-1
     if(btest(state,k-1)) fsgn = -fsgn
    enddo
    state = ibclr(state,site+spin*N-1)
    iszero = .false. 
   else
    iszero=.true.
    return
   endif
  else
   print*, 'Invalid operation in subroutine apply_op.'
   stop
  endif
 else
  return
 endif

 return
 end subroutine apply_op
 !============================================================================= 
 subroutine multiply_dvector_by_dH(vec,newvec,Nhs,Hi,Hj,Hval,Hlen,maxHLen,& 
                                                                psi,hs_size) 
 use mod_hilbert_space, only: ph_sector
 use mod_parameters
 implicit none
 integer omp_get_thread_num, omp_get_num_threads, tid
 integer(KIND=8) idx1, idx2, k,i
 integer(KIND=8) nn, i1, i2, i3, i4, j1, j2, j3, j4
 integer(KIND=8) ii, jj, Nel, Hlen, maxHLen, Nhs, hs_size
 integer(kind=8), dimension(1:hs_size) :: psi
 integer(KIND=8), dimension(1:maxHLen) :: Hi, Hj
 double precision, dimension(1:Nhs) :: vec, newvec
 double precision, dimension(1:maxHLen) :: Hval
 double precision rtmp

 !init the vector to zero
 newvec = 0.0d0
 !loop over the phonon subspace
 !$OMP parallel do default(none) & 
 !$OMP private(idx1,idx2,ii,jj,rtmp,i1,j1,k) &
 !$OMP shared(newvec,vec,Hi,Hj,Hlen,Hval,ph_sector,hs_size,psi)
 do i1 = 0,nquanta
  idx1 = ph_sector(i1)
  if(idx1.ne.-1)then
   !do the electronic part first
   do nn = 1,Hlen
    ii = Hi(nn) + idx1*Hs_size
    jj = Hj(nn) + idx1*Hs_Size
    rtmp = Hval(nn)
    if(ii.eq.jj)then
     rtmp = rtmp + dfloat(i1)*wph
    endif
    !$OMP ATOMIC 
    newvec(ii) = newvec(ii) + rtmp*vec(jj)
   enddo

   !the electron-phonon coupling terms
   j1 = i1 + 1    !subtract one quanta from oxygen site 1
   idx2 = ph_sector(j1)
   if(idx2.ne.-1)then
    rtmp = 0.0d0
    do nn = 1,hs_size
     ii = nn + idx1*hs_size
     jj = nn + idx2*hs_size
     k = 0
     do i = 1,9
      if(btest(psi(nn),i  -1)) k = k + 1
      if(btest(psi(nn),i+N-1)) k = k + 1
     enddo
     rtmp = gz*dfloat(k)*sqrt(dfloat(i1+1))
     !$OMP ATOMIC 
     newvec(ii) = newvec(ii) + rtmp*vec(jj)
    enddo
   endif
   j1 = i1 - 1    !subtract one quanta from oxygen site 1
   idx2 = ph_sector(j1)
   if(idx2.ne.-1)then
    rtmp = 0.0d0
    do nn = 1,hs_size
     ii = nn + idx1*hs_size
     jj = nn + idx2*hs_size
     k = 0
     do i = 1,9
      if(btest(psi(nn),i  -1)) k = k + 1
      if(btest(psi(nn),i+N-1)) k = k + 1
     enddo
     rtmp = gz*dfloat(k)*sqrt(dfloat(i1))
     !$OMP ATOMIC 
     newvec(ii) = newvec(ii) + rtmp*vec(jj)
    enddo
   endif
  endif
 enddo !i1 = 0,nquanta
 !$OMP END PARALLEL DO
 return
 end subroutine multiply_dvector_by_dH
 !============================================================================= 
 subroutine multiply_cvector_by_cH(vec,newvec,Nhs,Hi,Hj,Hval,Hlen,maxHLen,& 
                                                                psi,hs_size) 
 use mod_hilbert_space, only: ph_sector
 use mod_parameters
 implicit none
 integer(KIND=8) idx1, idx2, k,i
 integer(KIND=8) nn, i1, i2, i3, i4, j1, j2, j3, j4
 integer(KIND=8) ii, jj, Nel, Hlen, maxHLen, Nhs, hs_size
 integer(kind=8), dimension(1:hs_size) :: psi
 integer(KIND=8), dimension(1:maxHLen) :: Hi, Hj
 double complex, dimension(1:Nhs) :: vec, newvec
 double complex, dimension(1:maxHLen) :: Hval
 double complex ctmp

 !init the vector to zero
 newvec = dcmplx(0.0d0,0.0d0)
 !loop over the phonon subspace
 do i1 = 0,nquanta
  idx1 = ph_sector(i1)
  if(idx1.ne.-1)then
   !do the electronic part first
   do nn = 1,Hlen
    ii = Hi(nn) + idx1*Hs_size
    jj = Hj(nn) + idx1*Hs_Size
    ctmp = Hval(nn)
    if(ii.eq.jj)then
     ctmp = ctmp + dcmplx(dfloat(i1)*wph,0.0d0)
    endif
    newvec(ii) = newvec(ii) + ctmp*vec(jj)
   enddo

   !the electron-phonon coupling terms
   j1 = i1 + 1    !subtract one quanta from oxygen site 1
   idx2 = ph_sector(j1)
   if(idx2.ne.-1)then
    ctmp = dcmplx(0.0d0,0.0d0)
    do nn = 1,hs_size
     ii = nn + idx1*hs_size
     jj = nn + idx2*hs_size
     k = 0
     do i = 1,9
      if(btest(psi(nn),i  -1)) k = k + 1
      if(btest(psi(nn),i+N-1)) k = k + 1
     enddo
     ctmp = dcmplx(gz*dfloat(k)*sqrt(dfloat(i1+1)),0.0d0)
     newvec(ii) = newvec(ii) + ctmp*vec(jj)
    enddo
   endif
   j1 = i1 - 1    !subtract one quanta from oxygen site 1
   idx2 = ph_sector(j1)
   if(idx2.ne.-1)then
    ctmp = dcmplx(0.0d0, 0.0d0)
    do nn = 1,hs_size
     ii = nn + idx1*hs_size
     jj = nn + idx2*hs_size
     k = 0
     do i = 1,9
      if(btest(psi(nn),i  -1)) k = k + 1
      if(btest(psi(nn),i+N-1)) k = k + 1
     enddo
     ctmp = dcmplx(gz*dfloat(k)*sqrt(dfloat(i1)),0.0d0)
     newvec(ii) = newvec(ii) + ctmp*vec(jj)
    enddo
   endif
  endif
 enddo !i1 = 0,nquanta
 return
 end subroutine multiply_cvector_by_cH

 !=============================================================================
 subroutine remove_particle(vecin,Nin,basisin,hs_size_in,&
                            vecout,Nout,basisout,hs_size_out,site,spin,nsites)
 use mod_parameters, only: nquanta
 use mod_hilbert_space, only: return_idx_for_state
 implicit none
 integer(KIND=8) state, newstate
 integer(KIND=8) nn,i, j, k, Nin,hs_size_in, spin, nsites
 integer(KIND=8) Nout, hs_size_out, site
 double complex, dimension(1:Nin) :: vecin
 double complex, dimension(1:Nout) :: vecout
 integer(KIND=8), dimension(1:hs_size_in) :: basisin
 integer(KIND=8), dimension(1:hs_size_out) :: basisout
 double precision sgn

 vecout(:) = 0.0d0
 do i = 1,hs_size_in
  state = basisin(i)
  if(btest(state,site+spin*Nsites-1))then
   sgn = 1.0d0
   do k = 2*Nsites,site+spin*Nsites+1,-1
    if(btest(state,k-1)) sgn = -sgn
   enddo

   newstate = ibclr(state,site+spin*nsites-1)
   call return_idx_for_state(newstate,k,basisout,hs_size_out)

   do nn = 0,Nquanta
    vecout(k+nn*hs_size_out)=vecout(k+nn*hs_size_out)+sgn*vecin(i+nn*hs_size_in)
   enddo
  endif
 enddo
 return
 end subroutine remove_particle
  !=======================================================================================
 subroutine add_particle(vecin,Nin,basisin,hs_size_in, &
                         vecout,Nout,basisout,hs_size_out,site,spin,nsites)
 use mod_parameters, only: nquanta
 use mod_hilbert_space, only: return_idx_for_state
 implicit none
 integer(KIND=8) i,k,j,nn, state, newstate
 integer(KIND=8) Nin, Nout, hs_size_in, hs_size_out, site, spin, nsites
 integer(KIND=8), dimension(1:hs_size_in) :: basisin
 integer(KIND=8), dimension(1:hs_size_out):: basisout
 double complex, dimension(1:Nin) :: vecin
 double complex, dimension(1:Nout) :: vecout
 double complex sgn

 vecout = dcmplx(0.0d0,0.0d0)
 do i = 1,hs_size_in
  state = basisin(i)
  if(.not.btest(state,site+spin*Nsites-1))then
   sgn = dcmplx(1.0d0,0.0d0)
   do k = 2*nsites,site+spin*nsites+1,-1
    if(btest(state,k-1)) sgn = -sgn
   enddo
   newstate = ibset(state,site+spin*nsites-1)
   call return_idx_for_state(newstate,k,BasisOut,hs_size_out)
   do nn = 0,Nquanta
    vecout(k+nn*hs_size_out)=vecout(k+nn*hs_size_out)+sgn*vecin(i+nn*hs_size_in)
   enddo
  endif
 enddo
 return
 end subroutine add_particle
 
end module mod_hamiltonian
