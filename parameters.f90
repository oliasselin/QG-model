MODULE parameters

   IMPLICIT NONE

    integer, parameter :: n1=128, n2=128, n3=128
    integer, parameter :: npe=8

    integer, parameter :: n1d=n1+2, n2d=n2, n3d=n3
    integer, parameter :: n3h0=n3/npe, n3h1=n3/npe+2, n3h2=n3/npe+4
    integer, parameter :: n3d1=n3d+2*npe                !for transposed f, all z-lev + 1-lev halos stacked

    integer, parameter :: izbot1=2,izbot2=3
    integer, parameter :: iztop1=n3h1-1,iztop2=n3h2-2
    

    integer, parameter :: ktx = n1/2,  kty =n2/2
    integer, parameter :: iktx= ktx+1, ikty=n2, iktyp=n2/npe

    double complex :: i = (0.,1.)
    double precision, parameter :: twopi=4.D0*asin(1.D0)

    double precision, parameter :: L1=twopi, L2=twopi, L3=twopi               !Domain size
    double precision, parameter :: dx=L1/n1,dy=L2/n2,dz=L3/n3                 !Cell dimensions  

    real, parameter :: ktrunc_x = twopi/L1 * float(n1)/3.           ! dimensional truncation wavenumber (x)
    real, parameter :: ktrunc_z = twopi/L3 * float(n3)/3.           ! dimensional truncation wavenumber (x)


    !Tags to specify run!
    !-------------------!
    
    integer, parameter :: linear=0                      !1: set the nonlinear terms (advection) to 0. 
    integer, parameter :: inviscid=0                    !1: No dissipation, otherwise: dissipation
    integer, parameter :: init_wageo=0                  !1: Initialize wk with Ro*wak

    integer :: dealiasing=1         ! 1: Dealias, 0: don't. May not work though...

    !Should eventually plot both energies
!    integer, parameter :: plot_energy=1      !Use 1: energy_linear (equivalent to boussinesq including variable density, 2: energy_lipps)


    !Initial structure!
    !-----------------!


    integer, parameter :: generic=1 
    integer, parameter :: init_vertical_structure=1
    integer, parameter :: linear_vert_structure=0          !1: LINEAR psi as in SB2013      2: kz=1 for all k's       OTHERWISE: kz consistent with QG.  
    integer, dimension(1) :: seed = 2     ! Seed for the random number generator                                                                                                    
    integer, parameter :: initial_k = 5                    !Peak of k-spec
    integer, parameter :: enveloppe = 0                    !1: Enveloppe allowing b=0 at the boundaries
    double precision, parameter :: z_env   = twopi/8!twopi/3.!twopi/8         !Center of the tanh enveloppe
    double precision, parameter :: sig_env = twopi/24!twopi/6.!twopi/24      !Width  of the tanh enveloppe
    double precision, parameter :: z0  = L3/2                   !Position of the tropopause (between 0 and L3)

    
    !Normalization at the tropopause!
    !-------------------------------!
    
    integer, parameter :: norm_trop = 1                !Normalize (1) (or not: 0) with the RMS value of U at the tropopause (overrides normalize=1, normalization from total energy)                                                        
    integer, parameter :: trop_height = n3/2         !Position (in stag) where we want U_RMS to be computed for normalization                                                      
    double precision, parameter :: URMS = 1.D0

    !Or normalization from total energy!
    !----------------------------------!

    integer, parameter :: norm_energy=1      !Use 1: energy_linear (equivalent to boussinesq including variable density, 2: energy_lipps to normalize fields.)
    integer, parameter :: normalize=0
    real, parameter :: e_init=0.175!0.007            !0.175 for  generic psi and 0.007 for TS-type of init give u'~1.
    real, parameter :: k_init=(2./3.)*e_init            
    real, parameter :: p_init=e_init-k_init    


    !Base-state!
    !----------!

    integer, parameter :: boussinesq=1
    integer, parameter :: constant_N=0
    integer, parameter :: base_state=1
    integer, parameter :: fraction=128                   !If h#=150m, then fraction=133.333333~128
    double precision :: H_N = L3/fraction                          !Caracteristic length scale of N^2 for the TANH prof. (1/alpha...)
    double precision, parameter :: N_2_trop = 0.0001 !(2.*grav/10000.)*(1-t0_bot)/(1+t0_bot)             !BV frequency at the tropopause + in the tropsphere
    double precision, parameter :: N_2_stra = 0.0004                    !BV frequency in the stratosphere
    double precision, parameter :: gamma_N1=(sqrt(N_2_stra)-sqrt(N_2_trop))/(sqrt(N_2_stra)+sqrt(N_2_trop))       !This is alpha for N~1+alpha tanh(z/h)


   ! USEFUL INDEX !                                                                                                                          
   ! ------------ !                                                                                                                         

    integer :: ikx,iky,ikyp,izh0,izh1,izh2,izth   
    integer :: ix,iy,iz,izs
    integer :: kx,ky,kh2
    integer :: jj
    double precision :: x,y,z,zs

    ! USEFUL ARRAYS !                                                                                                                                                 
    ! ------------- !                                                                                                   
                                        
    integer, save :: kxa(iktx),kya(ikty)
    integer, save :: L(iktx,ikty)
    integer, save :: zath(n3)

    double precision, save :: xa(n1),ya(n2)
    double precision, save :: za(n3) ,zah2(n3h2) ,zah1(n3h1) ,zah0(n3h0)
    double precision, save :: zas(n3),zash2(n3h2),zash1(n3h1),zash0(n3h0)  !staggered version zasX(iz)=zaX(iz) + dz/2

    double precision, save :: r_1(n3h2),r_2(n3h2),r_3(n3h2)  !z-dependent r coefficients (r_1,2 unstag, while r_3 is stag)
    double precision, save :: r_1st(n3),r_2st(n3)            !Staggered and transposed versions of r_1 and r_2 for the elliptic equation.   (This could be a single value, not an array...) 
    double precision, save :: r_3u(n3h2),r_5u(n3h2)          !Unstaggered versions of r_3 and r_5 for the omega-equation verification.   
    double precision, save :: r_1ut(n3),r_2ut(n3)            !Unstaggered and transposed versions of r_1 and r_2 for the omega equation.  
    double precision, save :: r_3t(n3)                       !Transposed verion of r_3   (contains all n3 z-levels) for the pressure solver (still stag) 
    double precision, save :: r_3ut(n3)                      !Unstag version of r_3t
    double precision, save :: r_5ut(n3)                      
    double precision, save :: rho_st(n3)                     !Transposed verion of rho_s (contains all n3 z-levels) for diagnostics         (still stag) 
    double precision, save :: rho_ut(n3)                     !Transposed verion of rho_u (contains all n3 z-levels) for omega eqn (unstag)
    double precision, save :: a_ell_t(n3),b_ell_t(n3)        !coefficients of the elliptic equation for psi (LHQG)  --- transposed (for elliptic.f90)
    double precision, save :: a_ell_ut(n3)                   !coefficient of omega eqn --- transposed and UNstag - only computed for smooth_N for now
    double precision, save :: a_helm(n3),b_helm(n3)          !coefficients of the elliptic equation for phi ( QG )  --- transposed (for elliptic.f90)
    double precision, save :: a_ell(n3h2),b_ell(n3h2)        !coefficients of the elliptic equation for psi  --- not transposed (for setting the initial q in QG, and recover q from psi in LH)
    double precision, save :: a_ell_u(n3h2),b_ell_u(n3h2)    !coefficients of the elliptic equation for psi  --- unstaggered
    double precision, save :: rho_s(n3h2),rho_u(n3h2)        !Staggered and unstaggered versions of rho_0 (BS density)
    double precision, save :: r_1s(n3h2),r_2s(n3h2)          !Staggered version of r_1,2 (necessary to plot PV and initialize q) +  also conditions of integrability...
    double precision, save :: pi_0(n3h2)                     !For computing E_lh, Exner function's base-state (unstaggered) 
    double precision, save :: pi_0s(n3h2)                    !Staggered version of pi_0 (useful in two_exp BS
    double precision, save :: pi_0st(n3)                     !Staggered version of transposed pi_0 (useful in two_exp BS



    !I choose:
    !x N=0.03 s^-1 (to get N_troposphere ~ 0.01 right)
    !f=0.0001 s^-1 (real value of earth)
    !H=20 000m / 2pi
    !L=H/Ar 
    !U=1m/s 
    ! => Fr = pi/3
    ! => Ro = pi Ar
    !Keep in mind that: N/f = 1/Ar * (Ro/Fr).


    !Primary parameters!
    !------------------!

    double precision, parameter :: H_scale=20000./L3                 !Actual H in m ( z_real = H z' where z' in [0:L3]  is the nondim z.)
    double precision, parameter :: cor=0.0001!0.00000000001!0.0005 !0.0001                           !Actual f = 0.0001 s^-1 (real value of planet Earth)
    double precision, parameter :: L_scale=5.*(0.5*(sqrt(N_2_trop)+sqrt(N_2_stra))*H_scale/cor)  !200*H_scale                 !Actual L in m ( x_real = L x' where x' in [0:2pi] is the nondim x.)
    double precision, parameter :: U_scale=1.                        !Actual U in m/s (u_real = U u' where u' is the nondim velocity ur implemented in the code)

   
    double precision, parameter :: Ar2= (H_scale/L_scale)**2!(1./64.)**2!(1./10.)**2 !0.01     !Aspect ratio squared = (H/L)^2     
    double precision, parameter :: Ro = U_scale/(cor*L_scale)                                  !Rossby number  U/fL
    double precision, parameter :: Fr = U_scale/(0.5*(sqrt(N_2_trop)+sqrt(N_2_stra))*H_scale)  !Froude number  U/N(z0)H

    double precision, parameter :: big_F = Fr*Fr/(Ro*Ro)                    ! (Fr/Ro)^2



    !Timestepping!
    !------------!

    real :: time=0.
    integer :: iter
    integer :: itermax=1000000000
    real :: maxtime=40                      
    double precision, parameter :: delt=0.05*U_scale*dz    !0.0005*U_scale*dz                ! T_visc = 0.25D0*dz*dz/nu
    double precision, parameter :: gamma=1e-2!4e-3!1e-2!7.e-3            !Robert filter parameter

    !Other successful viscosity: 5e-2 * (10./ktrunc_x ) **2. 
    !PERFECT VISCOSITY: 0.01 * (64./(1.*n1)) **(4./3.)
    !In reality, nuh is 1/Re and nuz is 1/(Ar2*Re) with 1/Re = UL/nu

    double precision, parameter :: coeff =0.4!0.4!0.1!0.075
    double precision, parameter :: coeffz=0.!coeff!/10.!/1000!/10.

    integer, parameter :: ilap = 8                   !horizontal viscosity = nuh nabla^(2*ilap). So ilap =1 is regular viscosity. ilap>1 is hyperviscosity

    !General dissipation! (test for hyperviscosity: see Oct 10 2014 toread)
    double precision, parameter :: nuh  =  coeff * (64./(1.*n1)) **(4./3.) * (3./n1)**(2*(ilap-1))             !6e-2 * (10./ktrunc_x ) **2. ! horizontal visc coeff (regular viscosity)
    double precision, parameter :: nuz  = (coeffz* (64./(1.*n1)) **(4./3.) )                                      ! horizontal visc coeff (regular viscosity)
    double precision, parameter :: nuth =  coeff * (64./(1.*n1)) **(4./3.) * (3./n1)**(2*(ilap-1))             ! horizontal visc coeff (regular viscosity)
    double precision, parameter :: nutz = (coeffz* (64./(1.*n1)) **(4./3.) )                                      ! horizontal visc coeff (regular viscosity)

    !"Exact" dissipation:!
!    double precision, parameter :: nuh  =  0.015 * sqrt(Ar2) * (64./(1.*n1)) **(4./3.)              !6e-2 * (10./ktrunc_x ) **2. ! horizontal visc coeff (regular viscosity)
!    double precision, parameter :: nuz  = (0.015 * sqrt(Ar2) * (64./(1.*n1)) **(4./3.) )/Ar2       ! horizontal visc coeff (regular viscosity)
!    double precision, parameter :: nuth =  0.015 * sqrt(Ar2) * (64./(1.*n1)) **(4./3.)              ! horizontal visc coeff (regular viscosity)
!    double precision, parameter :: nutz = (0.015 * sqrt(Ar2) * (64./(1.*n1)) **(4./3.) )/Ar2       ! horizontal visc coeff (regular viscosity)



    !Output!
    !------!

    integer, parameter :: out_etot   = 1, freq_etot   = INT(0.01/delt)!50!346!n3/64!n3!64!n3!50*n3/64      !Total energy                                                    
    integer, parameter :: out_hspec  = 1, freq_hspec  = 5*freq_etot!n3/64!n3!freq_etot*10     !Horizontal energy spectrum at various heights 
    integer, parameter :: out_hg     = 0                 !Output geostrophic horizontal spectrum as well?
    integer, parameter :: out_vspec  = 0, freq_vspec =  freq_hspec
    integer, parameter :: out_vbuoy  = 1, freq_vbuoy =  freq_hspec
    integer, parameter :: out_vbuoyr = 0, freq_vbuoyr=  freq_etot
    integer, parameter :: out_ens    = 0, freq_ens   =  3*n3!freq_etot*10
    integer, parameter :: out_pv     = 0, freq_pv    =  3*n3!freq_etot*10

    integer, parameter :: out_ez     = 1, freq_ez    =  freq_etot        !E(z) (freq has to be a multiple of that of etot) 
    integer, parameter :: out_rotz   = 0, freq_rotz  =  freq_etot 
    integer, parameter :: out_ensz   = 0, freq_ensz  =  3*n3!freq_ens
    integer, parameter :: out_pvz    = 0, freq_pvz   =  freq_pv
    integer, parameter :: out_cond   = 0, freq_cond  =  5*freq_etot!*10        !Plot the conditions of integrability of the balance equations.
    integer, parameter :: out_grow   = 1, freq_grow  =  5*freq_etot!*10        !Plot the conditions of integrability of the balance equations.
    integer, parameter :: out_omega  = 1, freq_omega =  5*freq_etot!*10        !Compute the QG ageotrophic vertical velocity wak and war
    integer, parameter :: out_condwz = 0, freq_condwz=  freq_omega!*10        !Plot the w_z condition (requires out_omega = 1)
    integer, parameter :: out_cont   = 0, freq_cont  =  freq_etot!*10        !Plot the anelastic divergence (should be 0 because of the proj method)



    

    !For conditions:
    double precision :: jump_region_width = 5.

    !For vspec
    integer, parameter :: num_couples=n1 !Should be well enough to have a good estimate of the vertical spectrum.
    integer, save :: x0(num_couples) 
    integer, save :: y0(num_couples) 
    integer, parameter :: variance_spectrum=1               !1: Plots variance spectrum (no rho(z) or other z-factors, just fields squared). 0: regular spectra

    integer, parameter :: parabolic_L_peak = 1              !1: Refines L_peak by using npt (odd number>=3) to fit a parabola.  0: Use the peak simply (causes discontinuous L_peak)
    integer, parameter :: npt = 5

    !For slices                                                                                                                     
    integer, parameter :: stag=1,unstag=2
    integer, parameter :: where_bz=unstag
    integer, parameter :: num_spec = 10

    integer, parameter :: height(num_spec)=[n3/4, n3/2-25*n3/fraction, n3/2-10*n3/fraction, n3/2-5*n3/fraction, n3/2-2*n3/fraction, n3/2-1,  n3/2+2*n3/fraction, n3/2+5*n3/fraction, n3/2+10*n3/fraction, 3*n3/4]
!HERE MAX IS H5!                              0              1                  2                    3                     4           5             6                   7                   8               9


    !Slices
    integer, parameter :: max_slices = 99     
    integer, parameter :: nfields  = 7         !Don't forget to change tag_slice_xz(nfields) accordingly in "mpi.f90"
    integer, parameter :: nfields2 = 5         !Don't forget to change tag_slice_xz(nfields) accordingly in "mpi.f90"
    integer :: count_slice(nfields) = 0       !number of slices
    integer :: count_slice2(nfields2) = 0       !number of slices
    integer :: zval=n3/2                      !z-level at which we wish to plo a slice                                                                                                                               
    integer :: yval=n2/2
    integer :: xval=n1/2
    integer :: hlvl(nfields)=[2,1,2,1,2,1,1]                                   
    integer :: hlvl2(nfields2)=[1,1,1,1,0]                                   

    integer, parameter :: bot_height = n3/2+5*n3/fraction
    integer, parameter :: mid_height = n3/2+n3/fraction
    integer, parameter :: top_height = n3/2

    integer, parameter :: out_slab = 0, freq_slab = 1
    integer, parameter :: slab_mype   = npe/2-1 
    integer :: count_eta = 0
    
                                              !halo levels (u=2,zz=1...)                                                                                                                                                     
    integer :: id_field                       !dummy index to differenciate fields plotted  

    integer, parameter :: out_slice   = 0, freq_slice =  5* freq_etot
    integer, parameter :: out_eta     = 0, freq_eta   =  freq_hspec
    integer, parameter :: out_tspec   = 0

END MODULE parameters
