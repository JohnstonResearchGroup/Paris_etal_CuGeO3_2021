module mod_chebyshev
  integer, parameter :: cheby_M = 15
  double precision, parameter :: cheby_alpha = 0.01
  double complex, dimension(0:cheby_M-1) :: ck
  double precision :: cheby_a, cheby_b, cheby_Emin, cheby_Emax, cheby_dt
  double complex :: ipow
  integer(kind=8) :: Hsize

  double complex, allocatable, target, dimension(:) :: tempstate1, tempstate2, tempstate3

contains
  !==========================
  subroutine cheby_setup(Emin, Emax, dtime, Hsize_in)
    implicit none
    double precision :: Emin, Emax
    double precision :: dtime
    double complex, parameter :: eye = (0.0d0, 1.0d0)
    integer :: i
    integer(kind=8) :: Hsize_in

    Hsize = Hsize_in

    cheby_Emin = Emin
    cheby_Emax = Emax
    cheby_dt = dtime

    cheby_a = 0.5d0 * (1.0d0 + cheby_alpha) * (cheby_Emax - cheby_Emin)
    cheby_b = 0.5d0 * (cheby_Emax + cheby_Emin)
    !print *,cheby_Emin, cheby_Emax,cheby_a, cheby_b
    !print *,dtime,cheby_a*dtime

    print *,
    print *,'Chebyshev coefficients'
    print *,'Ensure these are sufficiently small'

    ! Calculate chebyshev expansion parameters
    ipow = 1.0d0
    do i = 0,cheby_M-1
        ck(i) = bessel_jn(i, cheby_a*dtime) * ipow
        ipow = -ipow * eye
        print*, 'ck(',i,') =',ck(i)
    end do


    print*, "Allocating temp states"
    allocate(tempstate1(0:Hsize-1))
    allocate(tempstate2(0:Hsize-1))
    allocate(tempstate3(0:Hsize-1))

  end subroutine

  !==========================

  subroutine cheby_cleanup()
      implicit none

      deallocate(tempstate1)
      deallocate(tempstate2)
      deallocate(tempstate3)

  end subroutine

  !==========================

  subroutine cheby_propagate(groundstate, newstate,        &
                             Hi, Hj, Hz, Hlen, maxHlen,    &
                             basis, nel)
    use mod_hamiltonian, only: multiply_cvector_by_cH

    implicit none
    double complex, parameter :: eye = (0.0d0, 1.0d0)

    ! Hamiltonian stuff for pass-through
    integer(kind=8) :: maxHlen, Hlen, nel
    integer(kind=8), dimension(1:maxHlen) :: Hi, Hj
    double complex, dimension(1:maxHlen)  :: Hz
    integer(kind=8), dimension(1:nel)   :: basis

    ! Input and output states
    double complex, target, dimension(0:Hsize-1) :: groundstate, newstate
    
    ! Pointers that we'll actually use for addressing
    double complex, pointer :: v0(:), v1(:)
    double complex, pointer :: vkmin1(:), vk(:), vkplus1(:)
    double complex, pointer :: tmp(:)

    ! Local variables
    double precision :: scalar, shift
    integer(kind=8) :: ii
    integer         :: k

    scalar = 1.0d0/cheby_a
    shift  = cheby_b

    ! Start by copying the ground state over since in this version of
    ! the algorithm we cannot overwrite it
    tempstate3 = groundstate

    v0(0:) => tempstate3
    v1(0:) => tempstate1



    ! 0th term
    newstate = ck(0) * v0

    ! Calculate v1
    call multiply_cvector_by_cH(v0, v1, Hsize,             &
                                Hi, Hj, Hz, Hlen, maxHlen, &
                                basis, nel)

    ! Scale by a and shift
    v1 = scalar * (v1 - shift * v0)
                                
    ! Add 1st term
    ! Factor of 2 here is from the chebyshev expansion
    newstate = newstate + 2.0d0 * ck(1) * v1


    ! Set my pointers
    vkmin1  => v0
    vk      => v1
    vkplus1 => tempstate2


    ! Recurrence relation: v(k+1) = 2 * Htilde * v(k) - v(k-1)
    do k=2,cheby_M-1
        
        ! Calculate product
        call multiply_cvector_by_cH(vk, vkplus1, Hsize,        &
                                    Hi, Hj, Hz, Hlen, maxHlen, &
                                    basis, nel)

        ! Scale and shift
        ! Factor of 2 here is from recurrence
        vkplus1 = 2.0d0 * scalar * (vkplus1 - shift * vk) 

        ! Subtract vkmin1
        vkplus1 = vkplus1 - vkmin1

        ! Shift vectors around
        tmp     => vkmin1
        vkmin1  => vk
        vk      => vkplus1
        vkplus1 => tmp

        ! Add kth term to cumulant with 2.0 prefactor from Cheb. expansion
        newstate = newstate + 2.0d0 * ck(k) * vk

    end do

    newstate = newstate * exp(-eye * cheby_b * cheby_dt)

  end subroutine

end module mod_chebyshev












