Uprime = 0.0d0
jd = 0.0d0

Uprime(1,2) = racahA+4.0d0*racahB+racahC
Uprime(2,1) = racahA+4.0d0*racahB+racahC

Uprime(1,3) = racahA-4.0d0*racahB+racahC
Uprime(3,1) = racahA-4.0d0*racahB+racahC

Uprime(2,3) = racahA-4.0d0*racahB+racahC 
Uprime(3,2) = racahA-4.0d0*racahB+racahC 

jd(1,2) = racahC
jd(2,1) = racahC

jd(1,3) = 4.0d0*racahB + racahC
jd(3,1) = 4.0d0*racahB + racahC

jd(2,3) = 4.0d0*racahB + racahC
jd(3,2) = 4.0d0*racahB + racahC
