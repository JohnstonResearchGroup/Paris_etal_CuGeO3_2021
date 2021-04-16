hopping = 0.0d0

!3d_{xy} -> 2px
do  i =  1,ncu
 ibl = 3*Ncu + nox + i
 ibr = 3*Ncu + nox + i + 1
 itl = 3*Ncu + nox + i + nox/2 
 itr = 3*Ncu + nox + i + nox/2 + 1
 hopping( i,itr) = tpd1x;              hopping(itr, i) = tpd1x
 hopping( i,ibr) =-tpd1x;              hopping(ibr, i) =-tpd1x
 hopping( i,ibl) =-tpd1x;              hopping(ibl, i) =-tpd1x
 hopping( i,itl) = tpd1x;              hopping(itl, i) = tpd1x

 ibl = 3*Ncu + i
 ibr = 3*Ncu + i + 1
 itl = 3*Ncu + i + nox/2 
 itr = 3*Ncu + i + nox/2 + 1
 hopping( i,itr) = tpd1y;              hopping(itr, i) = tpd1y
 hopping( i,ibr) = tpd1y;              hopping(ibr, i) = tpd1y
 hopping( i,ibl) =-tpd1y;              hopping(ibl, i) =-tpd1y
 hopping( i,itl) =-tpd1y;              hopping(itl, i) =-tpd1y
enddo

do  i =  1,ncu
 ibl = 3*Ncu + nox + i
 ibr = 3*Ncu + nox + i + 1
 itl = 3*Ncu + nox + i + nox/2 
 itr = 3*Ncu + nox + i + nox/2 + 1
 hopping( i+ncu,itr) =-tpd2x;              hopping(itr, i+ncu) =-tpd2x
 hopping( i+ncu,ibr) =-tpd2x;              hopping(ibr, i+ncu) =-tpd2x
 hopping( i+ncu,ibl) = tpd2x;              hopping(ibl, i+ncu) = tpd2x
 hopping( i+ncu,itl) = tpd2x;              hopping(itl, i+ncu) = tpd2x

 ibl = 3*Ncu + i
 ibr = 3*Ncu + i + 1
 itl = 3*Ncu + i + nox/2 
 itr = 3*Ncu + i + nox/2 + 1
 hopping( i+ncu,itr) = tpd2y;              hopping(itr, i+ncu) = tpd2y
 hopping( i+ncu,ibr) =-tpd2y;              hopping(ibr, i+ncu) =-tpd2y
 hopping( i+ncu,ibl) =-tpd2y;              hopping(ibl, i+ncu) =-tpd2y
 hopping( i+ncu,itl) = tpd2y;              hopping(itl, i+ncu) = tpd2y
enddo

do  i =  1,ncu
 ibl = 3*Ncu + nox + i
 ibr = 3*Ncu + nox + i + 1
 itl = 3*Ncu + nox + i + nox/2 
 itr = 3*Ncu + nox + i + nox/2 + 1
 hopping( i+2*ncu,itr) =-tpd3x;              hopping(itr, i+2*ncu) =-tpd3x
 hopping( i+2*ncu,ibr) =-tpd3x;              hopping(ibr, i+2*ncu) =-tpd3x
 hopping( i+2*ncu,ibl) = tpd3x;              hopping(ibl, i+2*ncu) = tpd3x
 hopping( i+2*ncu,itl) = tpd3x;              hopping(itl, i+2*ncu) = tpd3x

 ibl = 3*Ncu + i
 ibr = 3*Ncu + i + 1
 itl = 3*Ncu + i + nox/2 
 itr = 3*Ncu + i + nox/2 + 1
 hopping( i+2*ncu,itr) =-tpd3y;              hopping(itr, i+2*ncu) =-tpd3y
 hopping( i+2*ncu,ibr) = tpd3y;              hopping(ibr, i+2*ncu) = tpd3y
 hopping( i+2*ncu,ibl) = tpd3y;              hopping(ibl, i+2*ncu) = tpd3y
 hopping( i+2*ncu,itl) =-tpd3y;              hopping(itl, i+2*ncu) =-tpd3y
enddo

!============================================================
!Oxygen- Oxygen
do i = 3*Ncu+nox+1,3*Ncu+nox+nox/2-1
  j = i + 1
  hopping(i,j) = -tppx
  hopping(j,i) = -tppx
enddo
do i = 3*Ncu+nox+nox/2+1,3*Ncu+2*nox-1
  j = i + 1
  hopping(i,j) = -tppx
  hopping(j,i) = -tppx
enddo

do i = 3*Ncu+1,3*Ncu+nox/2
  j = i + nox/2
  hopping(i,j) = -tppy
  hopping(j,i) = -tppy
enddo

do i = 3*Ncu+1,3*Ncu+nox/2-1
  j = i + 1
  hopping(i,j) = tpppy
  hopping(j,i) = tpppy
enddo
do i = 3*Ncu+nox/2+1,3*Ncu+nox-1
  j = i + 1
  hopping(i,j) = tpppy
  hopping(j,i) = tpppy
enddo

do i = 3*Ncu+nox+1,3*Ncu+Nox+nox/2
  j = i + nox/2
  hopping(i,j) = tpppx
  hopping(j,i) = tpppx
enddo
