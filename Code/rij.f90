Rij = 0.0d0

do  j =  1,ncu
 do orb = 0,2
  i = j + orb*Ncu
  ibl = 3*Ncu + nox + j 
  ibr = 3*Ncu + nox + j + 1
  itl = 3*Ncu + nox + j + nox/2
  itr = 3*Ncu + nox + j + nox/2 + 1
  Rij( i,itr,:) = [ alat, blat]/2;              Rij(itr, i,:) =-[ alat, blat]/2
  Rij( i,ibr,:) = [ alat,-blat]/2;              Rij(ibr, i,:) =-[ alat,-blat]/2
  Rij( i,ibl,:) = [-alat,-blat]/2;              Rij(ibl, i,:) =-[-alat,-blat]/2
  Rij( i,itl,:) = [-alat, blat]/2;              Rij(itl, i,:) =-[-alat, blat]/2

  ibl = 3*Ncu + j
  ibr = 3*Ncu + j + 1
  itl = 3*Ncu + j + nox/2
  itr = 3*Ncu + j + nox/2 + 1
  Rij( i,itr,:) = [ alat, blat]/2;              Rij(itr, i,:) =-[ alat, blat]/2
  Rij( i,ibr,:) = [ alat,-blat]/2;              Rij(ibr, i,:) =-[ alat,-blat]/2
  Rij( i,ibl,:) = [-alat,-blat]/2;              Rij(ibl, i,:) =-[-alat,-blat]/2
  Rij( i,itl,:) = [-alat, blat]/2;              Rij(itl, i,:) =-[-alat, blat]/2
 enddo
enddo

!============================================================
!Oxygen- Oxygen
do i = 3*Ncu+nox+1,3*Ncu+nox+nox/2-1
  j = i + 1
  Rij(i,j,:) = [alat,0.0d0];                 
  Rij(j,i,:) =-[alat,0.0d0] 
enddo
do i = 3*Ncu+nox+nox/2+1,3*Ncu+2*nox-1
  j = i + 1
  Rij(i,j,:) = [alat,0.0d0];                 
  Rij(j,i,:) =-[alat,0.0d0] 
enddo

do i = 3*Ncu+1,3*Ncu+nox/2
  j = i + nox/2
  Rij(i,j,:) = [0.0d0,blat];                 
  Rij(j,i,:) =-[0.0d0,blat] 
enddo

do i = 3*Ncu+1,3*Ncu+nox/2-1
  j = i + 1
  Rij(i,j,:) = [alat,0.0d0];                 
  Rij(j,i,:) =-[alat,0.0d0]
enddo
do i = 3*Ncu+nox/2+1,3*Ncu+nox-1
  j = i + 1
  Rij(i,j,:) = [alat,0.0d0];                 
  Rij(j,i,:) =-[alat,0.0d0]
enddo

do i = 3*Ncu+nox+1,3*Ncu+Nox+nox/2
  j = i + nox/2
  Rij(i,j,:) = [0.0d0,blat];                 
  Rij(j,i,:) =-[0.0d0,blat] 
enddo
