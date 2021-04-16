module mod_hilbert_space
 integer(KIND=8) hs_size_UD         !Total size of the (Nup,Ndn) hilbert space
 integer(KIND=8) hs_size_Um1D       !Total size of the (Nup-1,Ndn) space 
 integer(KIND=8) hs_size_UDm1       !Total size of the (Nup,Ndn-1) space
 integer(KIND=8) NelUD              !Size of the (Nup,Ndn) electronic subspace 
 integer(KIND=8) NelUm1D            !Size of the (Nup-1,Ndn) electronic subspace
 integer(KIND=8) NelUDm1            !Size of the (Nup,Ndn-1) electronic subspace
 integer(KIND=8), allocatable, dimension(:) :: ph_sector
 integer(KIND=8), allocatable, dimension(:) :: basis_UD
 integer(KIND=8), allocatable, dimension(:) :: basis_Um1D
 integer(KIND=8), allocatable, dimension(:) :: basis_UDm1
contains
 !=============================================================================
 ! Subroutine init_basis allocates space for the basis arrays
 !=============================================================================
 subroutine init_basis()
 use mod_parameters, only: nquanta
 implicit none
 integer(KIND=8) i
 allocate(ph_sector(-1:nquanta+1))
 allocate(basis_UD(1:NelUD))
 allocate(basis_Um1D(1:NelUm1D))
 allocate(basis_UDm1(1:NelUDm1))
 ph_sector = -1
 do i = 0,nquanta
  ph_sector(i) = i
 enddo
 return
 end subroutine init_basis

 !=============================================================================
 ! Subroutine construct_basis generates the basis as an ordered array.
 !=============================================================================
 subroutine construct_basis(nup, ndn, N, hs_size, psi)
 use mod_quicksort
 implicit none
 integer(KIND=8) idx, state
 integer(KIND=8) nup, ndn, N, i,j, hs_size, hsup, hsdn
 integer(KIND=8), dimension(1:hs_size) :: psi
 integer(KIND=8), allocatable, dimension(:) :: pup, pdn
 integer(KIND=8), dimension(1:nup) :: vup
 integer(KIND=8), dimension(1:ndn) :: vdn
 i = 0
 hsup = size_hilbert_space(nup,i,N)
 hsdn = size_hilbert_space(i,ndn,N)
 allocate(pup(1:hsup))
 allocate(pdn(1:hsdn))

  i = 1
 idx = 0
 call combinations(pup,vup,i,N,i,nup,idx)
 !check idx is correct
 if(idx.ne.hsup)then
  print*, 'Incorrect number of spin up states.'
  print*, 'Expected: ', hsup, ' but generated ', idx,'.'
  stop
 endif
 idx = 0
 call combinations(pdn,vdn,i,N,i,ndn,idx)
 if(idx.ne.hsdn)then
  print*, 'Incorrect number of spin dn states.'
  print*, 'Expected: ', hsdn, ' but generated ', idx,'.'
  stop
 endif

 !now generate the basis set
 idx = 0
 do i = 1,hsup
  do j = 1,hsdn
   state = pup(i)*(2**N) + pdn(j)
   idx = idx + 1
   psi(idx) = state
  enddo
 enddo

 !check the total number of states
 if(idx.ne.hs_size)then
  print*, 'Incorrect total number of states.'
  print*, 'Expected: ', hs_size, ' but generated ', idx,'.'
  stop
 endif

 call quicksort(psi)
 deallocate(pup)
 deallocate(pdn)
 return
 end subroutine construct_basis


 !=============================================================================
 recursive subroutine combinations(array,v,start,n,k,maxk,idx)
 implicit none
 integer(KIND=8) i, start, n, k, maxk, state, idx
 integer(KIND=8), dimension(1:maxk) :: V
 integer(KIND=8), intent(in out), dimension(:) :: array
 if(k.gt.maxk)then
  state = 0
  idx = idx + 1
  do i = 1,maxk
   state = ibset(state,v(i)-1)
  enddo
  array(idx) = state
  return
 endif

 do i = start, n
  v(k) = i;
  call combinations(array,v,i+1,n,k+1,maxk,idx);
 enddo

 return
 end subroutine combinations
 !==============================================================================
 ! This routine returns the state index for a given state.  The search is
 ! performed using a binary search algoritm. 
 !==============================================================================
 subroutine return_idx_for_state(state,idx,psi,hs_size)
 implicit none
 integer(KIND=8) state, hs_size
 integer(KIND=8) idx, mid, low, hi, maxiter, iter
 integer(KIND=8), dimension(1:hs_size) :: psi
 logical loop
 maxiter = dlog(dfloat(hs_size))/dlog(dfloat(2))+10
 low = 0
 hi = hs_size+1
 mid = (hi+low)/2
 iter = 0
 loop = .true.
 do while(loop)
  iter = iter + 1
  if(state.eq.psi(mid))then
   loop = .false.
  elseif(state.lt.psi(mid))then
   hi=mid
   mid=(hi+low)/2
  else
   low = mid
   mid=(hi+low)/2
  endif
  if(iter.ge.maxiter)then
   print*, 'State not found in hilbert space after', maxiter, ' lookups.'
   print*, 'Looking for state ', state
   stop
  endif
 enddo
 idx = mid
 return
 end subroutine return_idx_for_state
 !=============================================================================
 ! function size_hilbert_space 
 ! returns the size of the hilbert space for nup and ndn fermions.
 !=============================================================================
 integer(KIND=8) function size_hilbert_space(nup,ndn,nsites)
 implicit none
 integer(KIND=8) nup, ndn, nsites
 size_hilbert_space = nchoosek(nsites,nup)*nchoosek(nsites,ndn)
 return
 end function size_hilbert_space
 !=============================================================================
 ! integer function nchoosek 
 ! integer implemenation of a combination.
 !=============================================================================
 integer(KIND=8) function nchoosek(n,k)
 implicit none
 integer(KIND=8) i, n, k, itmp1, itmp2, itmp3
 itmp1 = 1
 do i = n,n-k+1,-1
  itmp1 = itmp1*i
 enddo
 call factorial(k,itmp2)
 nchoosek = itmp1/itmp2
 return
 end function nchoosek
 !=============================================================================
 ! RECURSIVE SUBROUTINE FACTORIAL(n,fac)
 ! Recursive implementation of n factorial where n is an integer.
 !=============================================================================
 recursive subroutine factorial(n,fac)
 implicit none
 integer(KIND=8), intent(in) :: n
 integer(KIND=8), intent(out) :: fac
 integer(KIND=8) itmp
 if(n.ge.1)then
  call factorial(n-1,itmp)
  fac = n*itmp
 else
  fac = 1
 endif
 return
 end subroutine factorial
 
end module mod_hilbert_space 
