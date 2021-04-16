! mod_parameters is a module that holds all parameters for the simulation. 
module mod_parameters
 integer(KIND=8), parameter :: create=1
 integer(KIND=8), parameter :: annhil=0
 integer(KIND=8), parameter :: spinup = 1
 integer(KIND=8), parameter :: spindn = 0
 integer(KIND=8), parameter :: nquanta = 0   !# phonon quanta
 integer(KIND=8), parameter :: nup = 2       !# up electrons
 integer(KIND=8), parameter :: ndn = 1       !# down electrons
 integer(KIND=8), parameter :: ncu = 3       !# Cu atoms
 integer(KIND=8), parameter :: nox = 2*ncu+2     !# oxygen atoms
 integer(KIND=8), parameter :: N = 3*Ncu+2*Nox   !# orbitals in the cluster
 double precision, parameter :: alat = 1.0d0     !lattice constants
 double precision, parameter :: blat = 1.0d0
 double precision, parameter :: racahA = 6.45d0
 double precision, parameter :: racahB = 0.25d0
 double precision, parameter :: racahC = 0.35d0
 double precision, parameter :: edz2_r2 = 1.20d0
 double precision, parameter :: edx2_y2 = 1.60d0
 double precision, parameter :: edxy = 0.0d0
 double precision, parameter :: deltaE_Cu = 0.50d0
 double precision, parameter :: deltaE_O = 0.98d0
 double precision, parameter :: X0 = 0.0d0
 double precision, parameter :: XAmp = 0.07d0
 double precision, parameter :: hbar = 0.658211951d0 !eV.fs
 double precision, parameter :: omegamode = 0.15d0*0.0136d0/hbar 
 double precision, parameter :: tau = 500d0
 double precision, parameter :: decu = 4.0d0
 double precision, parameter :: Jdimer = 1.00d0
 double precision, parameter :: epx = 5.60d0!-0.7d0           
 double precision, parameter :: epy = 5.25d0!-0.7d0        
 double precision, parameter :: Up = 4.1d0
 double precision, parameter :: Ud = 10.0d0
 double precision, parameter :: Jp = 0.6d0
 double precision, parameter :: Jpd = 0.11d0
 double precision, parameter :: Upd = 1.2d0
 double precision, parameter :: Vdd = 0.5d0
 double precision, parameter :: tpd1x = 1.1152d0    !These were computed using (pds) =1.6 and (pdp) = -0.72
 double precision, parameter :: tpd1y = 0.8011d0    !0.6382d0
 double precision, parameter :: tpd2x = 0.2944d0    !0.3175d0
 double precision, parameter :: tpd2y = 0.6832d0    !0.4405d0
 double precision, parameter :: tpd3x = 0.6088d0    !0.4132d0
 double precision, parameter :: tpd3y = 0.5191d0    !0.3853d0
 double precision, parameter :: tppx = 0.96d0  !0.840d0
 double precision, parameter :: tppy = 1.0d0 !0.960d0
 double precision, parameter :: tpppx = 0.250d0
 double precision, parameter :: tpppy = 0.240d0
 double precision, parameter :: wph = 0.074d0
 double precision, parameter :: Uq = 5.0d0
 double precision, parameter :: gz = 0.00d0
 double precision, parameter :: sig = 0.26d0
 double precision, parameter :: Amp = 5.0d0
 double precision, parameter :: eV_to_fs = 4.13567d0
 double precision, parameter :: pulsewidth = 50.0d0/hbar
 double precision, parameter :: pump_freq = 4.65d0                 !hbar*w_pump in eV
 double precision, parameter :: eps_pumpAx = 0.0d0    !polarization of the pump
 double precision, parameter :: eps_pumpAy = 1.0d0
 double complex, parameter :: igam = (0.0d0,0.200d0)
contains
end module mod_parameters
