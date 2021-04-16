module mod_lanczos
 integer, parameter :: northo = 10
 integer, parameter :: ndiag =  10
 integer, parameter :: nimax = 400
 integer, parameter :: nfmax = 400
 integer, parameter :: nconvergef = 400
 integer, parameter :: nconvergei = 400
 double precision, parameter :: orth_limit = 1e-12
contains
 !=============================================================================
 subroutine project_on_final(vecin,N,basis,hs_size,Ef,Overlap,nfinal,Hi,hj,hval,&
                             Hlen,maxHlen,reached,orthogonalize)
 use mod_random
 use mod_hamiltonian, only: multiply_cvector_by_cH
 implicit none
 logical orthogonalize
 logical loop
 integer(KIND=8) N, hs_size, nfinal, reached, info, i, j, k, cnt
 integer(KIND=8) maxHlen, Hlen, nlanczos
 integer(KIND=8), dimension(1:maxHlen) :: Hi
 integer(KIND=8), dimension(1:maxHlen) :: Hj
 integer(KIND=8), dimension(1:hs_size) :: basis
 double precision, dimension(1:nfmax) :: Ef
 double precision, dimension(1:nfmax) :: Enew,Eold
 double complex, dimension(1:N) :: vecin, C0, C1, C2
 double complex, dimension(1:nfmax) :: overlap
 double complex, dimension(1:maxHlen) :: Hval
 double precision norm, rtmp
 double complex ctmp
 double complex, parameter :: czero = (0.0d0,0.0d0)
 double precision, allocatable, dimension(:) :: alpha, beta, tmp1d
 double precision, allocatable, dimension(:) :: alpha_temp, beta_temp
 complex(kind=8), dimension(1:N) :: tmpvec
 double precision, allocatable, dimension(:) :: work
 complex(kind=8), allocatable, dimension(:,:) :: l_vecs_c, tmp2d
 double precision, allocatable, dimension(:,:) :: tmpmat

 if(nfinal.gt.nfmax)then
  print*, 'Error, nfinal larger than nmax.'
  stop
 else
  nlanczos = min(nfinal,nfmax)
  allocate(alpha(1:nlanczos))
  allocate(beta(1:nlanczos))
 endif
 Enew = 0.0d0
 Eold = 1e4

 c0 = vecin
 norm = sqrt(dot_product(c0,c0))
 c0 = c0/norm

 if(ORTHOGONALIZE)then
  allocate(l_vecs_c(1:N,1:nlanczos))
  l_vecs_c(:,1) = c0(:)
  cnt = 0
 endif

 call multiply_cvector_by_cH(c0,c1,N,Hi,Hj,Hval,Hlen,maxHlen,basis,hs_size)

 alpha(1) = dot_product(c0,c1)
 beta(1) = czero
 c1 = c1 - alpha(1)*c0
 c1 = c1/sqrt(dot_product(c1,c1))

 i = 1
 loop = .true.
 do while(loop)
  i = i + 1

  call multiply_cvector_by_cH(c1,c2,N,Hi,Hj,Hval,Hlen,maxHlen,basis,hs_size)
  alpha(i) = dot_product(c1,c2)
  beta(i) = dot_product(c0,c2)
  c2 = c2 - alpha(i)*c1 - beta(i)*c0
  c0 = c1

  if(ORTHOGONALIZE)then
   cnt = cnt + 1
   l_vecs_c(:,i) = c1
   rtmp = conjg(dot_product(c1,VecIn))
   if(cnt.ge.northo.or.rtmp.gt.orth_limit)then
    do j = 1,i
     ctmp = dot_product(l_vecs_c(:,j),C2)
     c2 = c2 - ctmp*l_vecs_c(:,j)
     enddo
   endif
  endif
  c1 = c2/sqrt(dot_product(c2,c2))

  !Now check the eigenvalues for convergence
  if(i.eq.NLanczos)then
   !diagonalize the problem
   !first allocate space for the problem
   allocate(work(1:2*Nlanczos-2))
   allocate(tmpmat(1:nlanczos,1:nlanczos))
   allocate(beta_temp(1:nlanczos))
   allocate(alpha_temp(1:nlanczos))
   alpha_temp = alpha
   beta_temp = beta
   !lapack call
   call dstev('N',nlanczos,alpha_temp,beta_temp(2:),tmpmat,nlanczos,work,info)
   Eold(:) = Enew(:)
   Enew(1:nlanczos) = alpha_temp(1:nlanczos)   !store new eigenvaules for the ground states
   deallocate(work)       !get rid of memory we need to reuse
   deallocate(tmpmat)
   deallocate(beta_temp)
   deallocate(alpha_temp)

   if(maxval(abs(enew(1:nconvergef)-eold(1:nconvergef))).le.1e-8)then
    loop = .false.
   else
    k = nlanczos
    nlanczos = min(nlanczos+ndiag,nfmax)
    !reallocate memory
    !redimension alpha
    allocate(tmp1d(1:k))
    tmp1d(:) = alpha(:)
    deallocate(alpha)
    allocate(alpha(1:nlanczos))
    alpha(1:k) = tmp1d(:);
    alpha(k+1:) = 0.0d0

    !redimension beta
    tmp1d(1:k) = beta(:)
    deallocate(beta)
    allocate(beta(1:nlanczos))
    beta(1:k) = tmp1d;
    beta(k+1:) = 0.0d0
    deallocate(tmp1d)

    !redimension l_vecs_c
    if ( orthogonalize ) then
     allocate(tmp2d(1:N,1:k))
     do j = 1,k
      tmp2d(:,j) = l_vecs_c(:,j)
     enddo
     deallocate(l_vecs_c)
     allocate(l_vecs_c(1:N,1:nlanczos))
     l_vecs_c = dcmplx(0.0d0, 0.0d0)
     do j = 1,k
      l_vecs_c(:,j) = tmp2d(:,j)
     enddo
     deallocate(tmp2d)
    endif

   endif !if(maxval(abs(enew-eold)).le.1e-8)then                         
  endif

  if(i.eq.nfmax)then
   loop = .false.
  endif
 enddo

 !diagonalize
 allocate(work(1:2*nlanczos-2))
 allocate(tmpmat(1:nlanczos,1:nlanczos))
 call dstev('V',nlanczos,alpha,beta(2:),tmpmat,nlanczos,work,info)
 Ef(1:nlanczos) = alpha(:)

! do i = 1,20 !nlanczos
!  write(6,*) alpha(i)
! enddo

 overlap = dcmplx(0.0d0,0.0d0)
 do i = 1,nlanczos
  overlap(i) = tmpmat(1,i)*norm
 enddo

 deallocate(work)
 deallocate(tmpmat)

 if(ORTHOGONALIZE)then
  deallocate(l_vecs_c)
 endif
 !THIS IS J
 reached = nlanczos

 return
 end subroutine project_on_final
 !=============================================================================
 subroutine get_ground_state(states,Eg,num,basis,N,hs_size,Hi,Hj,Hel,Hlen,& 
                             maxHlen,nstates,orthogonalize)
 use mod_random
 use mod_hamiltonian, only: multiply_dvector_by_dH
 implicit none
 logical orthogonalize
 integer(KIND=8) num, N, hs_size, Hlen, maxHLen, nstates
 integer(KIND=8) i, j, cnt, info
 integer(KIND=8), dimension(1:hs_size) :: basis
 integer(KIND=8), dimension(1:maxHLen) :: Hi, Hj
 double precision, dimension(1:maxHLen) :: Hel
 double precision rtmp
 double precision, dimension(1:num) ::Eg
 double precision, dimension(1:N) :: c0,c1,ctmp,W,C2,ci
 double precision, dimension(1:N,1:num) :: states 
 double precision, dimension(1:nstates) :: alpha, beta
 double precision, dimension(1:nstates) :: alpha_temp, beta_temp
 double precision, allocatable, dimension(:,:) :: l_vecs_c, vec
 double precision, allocatable, dimension(:) :: work

 print*, 'Fetching the ground state and first ', num-1, ' exited states.'
 !set the vector c1 to a random vector with norm 1 and set the vector c2 to zero
 do i = 1,N
  c1(i) = rand() - 0.5d0
  c0(i) = 0.0d0
 enddo
 c1 = c1/sqrt(dot_product(c1,c1))
 ci = c1

 if(ORTHOGONALIZE)then
  allocate(l_vecs_c(1:N,1:nstates))
  l_vecs_c(:,1) = c1(:)
  cnt = 0
 endif

 !beta(1) is zero by design
 beta(1) = 0.0d0

 do i = 1,nstates
  !here c1 plays the role of V_i, c0 is V_{i-1}, and c2 is V_{i+1}
  call multiply_dvector_by_dH(c1,ctmp,N,Hi,Hj,Hel,Hlen,maxHlen,basis,hs_size)

  W = ctmp - beta(i)*c0
  alpha(i) = dot_product(c1,W)
  W = W - alpha(i)*c1
  if(i.ne.nstates)then
   beta(i+1) = sqrt(dot_product(W,W))
   C2 = W/beta(i+1)
  endif

  if(ORTHOGONALIZE)then
   l_vecs_c(:,i) = c1(:)
   cnt = cnt + 1
   rtmp = abs(dot_product(c2,l_vecs_c(:,1)))
   !here we need to reorthogonalize 
   !if(rtmp.ge.orth_limit)then
   if(cnt.ge.northo.or.rtmp.gt.orth_limit)then
    do j = 1,i-1
     rtmp = dot_product(l_vecs_c(:,j),c2)
     c2 = c2 - rtmp*l_vecs_c(:,j)
    enddo
   endif
  endif

  c0 = c1
  c1 = c2/sqrt(dot_product(c2,c2))
 enddo

  !now we diagonalize the problem using a lapack call.
 allocate(work(1:2*nstates-2))
 allocate(vec(1:nstates,1:nstates))
 call dstev('V',nstates,alpha,beta(2:),vec,nstates,work,info)
 do i = 1,num
  Eg(i) = alpha(i)
 enddo

 !construct the ground state
 c1 = ci
 c0 = 0.d0
 states = 0.0d0

 do i = 1,nstates
  if( orthogonalize ) then
   do j = 1,num
    states(:,j) = states(:,j) + vec(i,j)*l_vecs_c(:,i)
   enddo
  else
   do j = 1,num 
    states(:,j) = states(:,j) + vec(i,j)*c1
   enddo
   call multiply_dvector_by_dH(c1,ctmp,N,Hi,Hj,Hel,Hlen,maxHlen,basis,hs_size)
   W = ctmp - beta(i)*c0
   alpha(i) = dot_product(c1,W)
   W = W - alpha(i)*c1
   if(i.ne.nstates)then
    beta(i+1) = sqrt(dot_product(W,W))
    C2 = W/beta(i+1)
   endif
   C0 = C1
   c1 = c2/sqrt(dot_product(c2,c2))
  endif
 enddo

 deallocate(work) !free the memory
 deallocate(vec)

 if(ORTHOGONALIZE)then
  deallocate(l_vecs_c)
 endif

 return
 end subroutine get_ground_state 

 subroutine get_frequ_vec_lanczos(cgs,N,hs_size,basis,Hi,Hj,Hval,Hlen,maxHlen,&
                                  w,Eg,vec_out,nistates,ORTHOGONALIZE)
 use mod_parameters, only: igam
 use mod_hamiltonian, only: multiply_cvector_by_cH
 use mod_random
 implicit none
 logical loop
 logical ORTHOGONALIZE
 integer(KIND=8) cnt, NLanczos
 integer(KIND=8) info      !output status for lapack calls
 integer(KIND=8) i,j,k     !dummy counters 
 integer(KIND=8) N, hs_size, Hlen, nistates, maxHlen
 integer(KIND=8), dimension(1:hs_size) :: basis
 integer(KIND=8), dimension(1:maxHLen) :: Hi, Hj
 double complex, dimension(1:maxHLen) :: Hval
 double complex, dimension(1:N) :: cgs, C0, C1, C2
 double precision, allocatable, dimension(:,:) :: l_vecs_c, tmpmat, tmp2d
 double precision, allocatable, dimension(:) :: work
 complex(kind=8), allocatable, dimension(:) :: phivec
 complex(kind=8), dimension(1:N) :: vec_out
 double complex, dimension(1:N) :: fvecs
 double precision, dimension(1:nimax) :: enew, eold
 double precision w, eg, norm, rtmp
 double precision, allocatable, dimension(:) :: alpha, beta, tmp1d
 double precision, allocatable, dimension(:) :: alpha_temp, beta_temp
 norm = sqrt(dot_product(cgs,cgs))
 C0 = cgs/norm

 Nlanczos = min(nistates,nimax)
 allocate(alpha(1:nlanczos))
 allocate(beta(1:nlanczos))

 if(ORTHOGONALIZE)then
  allocate(l_vecs_c(1:N,1:nlanczos))
  l_vecs_c(:,1) = c0(:)
  cnt = 0
 endif

 Eold = 1e4
 Enew = 0.0d0
 call multiply_cvector_by_cH(c0,c1,N,Hi,Hj,Hval,Hlen,maxHlen,basis,hs_size)

 alpha(1) = dot_product(c0,c1)
 beta(1) = 0.0d0
 c1 = c1 - alpha(1)*C0
 c1 = c1/sqrt(dot_product(c1,c1))
 i = 1

 loop = .true.
 do while(loop)
  i = i + 1
  call multiply_cvector_by_cH(c1,c2,N,Hi,Hj,Hval,Hlen,maxHlen,basis,hs_size)
  alpha(i) = dot_product(c1,c2)
  beta(i) = dot_product(c0,c2)

  c2 = c2 - alpha(i)*c1 - beta(i)*C0
  C0 = c1
  if(ORTHOGONALIZE)then
   cnt = cnt + 1
   l_vecs_c(:,i) = c1
   rtmp = abs(dot_product(c1,cgs))
   if(cnt.ge.northo.or.rtmp.gt.orth_limit)then
    do j = 1,i
     rtmp = dot_product(l_vecs_c(:,j),C2)
     c2 = c2 - rtmp*l_vecs_c(:,j)
    enddo
   endif
  endif
  c1 = c2/sqrt(dot_product(c2,c2))
  !now check the convergence of the lowest states. We do this when i = NLanczos
  if(i.eq.nlanczos)then
   !diagonalize the problem
   !first allocate space for the problem
   allocate(work(1:2*Nlanczos-2))
   allocate(tmpmat(1:nlanczos,1:nlanczos))
   allocate(beta_temp(1:nlanczos))
   allocate(alpha_temp(1:nlanczos))
   alpha_temp = alpha
   beta_temp = beta
   !lapack call
   Eold(1:nimax) = Enew(1:nimax)
   call dstev('N',nlanczos,alpha_temp,beta_temp(2:),tmpmat,nlanczos,work,info)
   Enew(1:nlanczos) = alpha_temp(1:nlanczos)   !store new eigenvalues for the ground states

   deallocate(work)       !get rid of memory we need to reuse
   deallocate(tmpmat)
   deallocate(beta_temp)
   deallocate(alpha_temp)

   if(maxval(abs(Enew(1:nconvergei)-Eold(1:nconvergei))).le.1e-8)then
    loop = .false.
   else
    !increment nlanczos
    k = nlanczos
    nlanczos = min(nlanczos+ndiag,nimax)

    !reallocate memory...
    !first we redimension alpha
    allocate(tmp1d(1:k))
    tmp1d(1:k) = alpha(:)
    deallocate(alpha)
    allocate(alpha(1:nlanczos))
    alpha(1:k) = tmp1d;
    alpha(k+1:) = 0.0d0

    !redimension beta
    tmp1d(1:k) = beta(:)
    deallocate(beta)
    allocate(beta(1:nlanczos))
    beta(1:k) = tmp1d;
    beta(k+1:) = 0.0d0

    deallocate(tmp1d)
    !redimension l_vecs_c
    if ( orthogonalize ) then
     allocate(tmp2d(1:N,1:k))
     do j = 1,k
      tmp2d(:,j) = l_vecs_c(:,j)
     enddo
     deallocate(l_vecs_c)
     allocate(l_vecs_c(1:N,1:nlanczos))
     l_vecs_c = 0.d0
     do j = 1,k
      l_vecs_c(:,j) = tmp2d(:,j)
     enddo
     deallocate(tmp2d)
    endif
   endif
  endif
  if(i.eq.nimax) loop = .false.
 enddo

 !now solve the problem
 allocate(work(1:2*Nlanczos-2))
 allocate(tmpmat(1:nlanczos,1:nlanczos))
 allocate(phivec(1:nlanczos))
 allocate(alpha_temp(lbound(alpha,1):ubound(alpha,1)))
 allocate(beta_temp(lbound(beta,1):ubound(beta,1)))
 alpha_temp = alpha
 beta_temp = beta
 call dstev('V',nlanczos,alpha_temp,beta_temp(2:),tmpmat,nlanczos,work,info)
 Enew = 0.d0
 Enew(1:nlanczos) = alpha_temp
 deallocate(alpha_temp)
 deallocate(beta_temp)

 vec_out = dcmplx(0.0d0,0.0d0)
 if ( orthogonalize ) then
  do j = 1,i
   !put the norm back in and divide out the energy denominator
   phivec(:) = norm*dcmplx(tmpmat(1,:),0.d0)/(Eg - Enew(1:nlanczos) + w + igam)
   phivec = matmul(tmpmat,phivec)
   vec_out(:) = vec_out(:) + phivec(j)*l_vecs_c(:,j)
  enddo
 else
  c1 = cgs/norm
  do j = 1,i
   !put the norm back in and divide out the energy denominator
   phivec(:) = norm*dcmplx(tmpmat(1,:),0.d0)/(Eg - Enew(1:nlanczos) + w + igam)
   phivec = matmul(tmpmat,phivec)
   vec_out(:) = vec_out(:) + phivec(j)*c1
   call multiply_cvector_by_cH(c1,c2,N,Hi,Hj,Hval,Hlen,maxHlen,basis,hs_size)
   c2 = c2 - alpha(j)*c1 - beta(j)*c0
   c0 = c1
   c1 = c2/sqrt(dot_product(c2,c2))
  enddo
 endif

 deallocate(work) !free the memory
 deallocate(tmpmat)
 deallocate(phivec)
 deallocate(alpha)
 deallocate(beta)
 !free all memory no longer needed
 if(ORTHOGONALIZE)then
  deallocate(l_vecs_c)
 endif

 return
 end subroutine get_frequ_vec_lanczos


end module mod_lanczos
