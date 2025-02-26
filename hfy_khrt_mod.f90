
!-------------------------------------------------------------------------------
! KHRT model: a hybrid secondary breakup model, including two submodel:
! i : KH, i.e. Kelvin-Helmholtz model, the mechanism of that is the same with
!     WAVE model which models process of first atomization.
! ii: RT, i.e. Rayleigh-Taylor model, modeling the process of light and heavy
!     fluids merge into each other

! author    date            .org
! snx       2021.05.10      simpop.cn
!-------------------------------------------------------------------------------

module hfy_khrt_mod
  implicit none
  private
  intrinsic :: acos, selected_real_kind
  !-----------------------------------------------------------------------------
  ! Precision Choice, by default
  ! (1) p=10 => RKD = 8, double precision for real/floating variables
  integer, parameter :: RKD = selected_real_kind( p=10 )
  !-----------------------------------------------------------------------------
  ! Named-Constants only used within this module
  ! Real
  real(RKD), parameter :: QUARTER = 0.25_RKD
  real(RKD), parameter :: THIRD   = 1.0_RKD/3.0_RKD
  real(RKD), parameter :: HALF    = 0.5_RKD

  real(RKD), parameter :: ZERO  =  0.0_RKD
  real(RKD), parameter :: ONE   =  1.0_RKD
  real(RKD), parameter :: TWO   =  2.0_RKD
  real(RKD), parameter :: THREE =  3.0_RKD
  real(RKD), parameter :: FOUR  =  4.0_RKD
  real(RKD), parameter :: FIVE  =  5.0_RKD
  real(RKD), parameter :: SIX   =  6.0_RKD
  real(RKD), parameter :: SEVEN =  7.0_RKD
  real(RKD), parameter :: EIGHT =  8.0_RKD
  real(RKD), parameter :: NINE  =  9.0_RKD
  real(RKD), parameter :: TEN   = 10.0_RKD

  ! other parameters

  real(RKD), parameter :: PI = acos( -ONE )

  ! surface tension between gaseous phase and liquid parcel
  real(RKD), parameter :: SIGMA_GP = 25.0e-3_RKD

  ! dynamic viscosity of liquid parcels
  real(RKD), parameter :: MU_PARC = 0.75e-3_RKD

  ! density of liquid droplets
  real(RKD), parameter :: DENSLIQ = 790.0_RKD
  !-----------------------------------------------------------------------------

  ! For Interacting and info transfer between subroutines in this Module.
  ! Specificly,
  ! These 4 global variables are computed in hfy_khrt_choose_submodel, and then
  !  used in the subsequent hfy_khrt_kh & hfy_khrt_rt submodel.
  real(RKD) :: Lambda_RT = ZERO  ! 3.2.2.4-27
  real(RKD) :: radi_s_KH = ZERO  ! 3.2.2.4-2,  rs
  real(RKD) :: radi_b_KH = ZERO  ! 3.2.2.4-12, rb
  real(RKD) :: Lambda_kh = ZERO  ! 3.2.2.4-4
  real(RKD) :: Omega_kh  = ZERO  ! 3.2.2.4-5

  !-----------------------------------------------------------------------------
  ! Number of newly-generated small son parcels in this rank.
  ! After each KH breakup, numb_parc_local_new += 1
  integer :: numb_newly_gen_parc = 0

  ! Local number of parcels already in this rank
  integer :: numb_local_old_parc = 0
  !-----------------------------------------------------------------------------

  ! Store properties/attributes, such as velocity, mass, diameter etc., from
  ! eptp/pepa in CS platform. And also the two flags for KHRT model, which have
  ! already been added to pepa array in platform
  real(RKD), dimension(:,:), allocatable :: prop_old_parc
  real(RKD), dimension(:,:), allocatable :: prop_new_parc

  ! auxiliary array index. The number must be <= NIND, for both arrays:
  ! prop_old_parc, and prop_new_parc

  integer, parameter :: IKDP = 1         ! diameter
  integer, parameter :: IKNP = IKDP + 1  ! number-attribute
  integer, parameter :: IKMP = IKNP + 1  ! mass

  integer, parameter :: IKUP = IKMP + 1  ! u of parcel
  integer, parameter :: IKVP = IKUP + 1  ! v of parcel
  integer, parameter :: IKWP = IKVP + 1  ! w of parcel

  integer, parameter :: IKUF = IKWP + 1  ! u of fluid seen
  integer, parameter :: IKVF = IKUF + 1  ! v of fluid seen
  integer, parameter :: IKWF = IKVF + 1  ! w of fluid seen

  ! density of gas/continous phase just there at the parcel's position
  integer, parameter :: IKRF = IKWF + 1
  ! molecular dynamical viscosity locally seen by parcel
  integer, parameter :: IMUF = IKRF + 1
  ! molecular dynamical viscosity of the liquid parcel
  integer, parameter :: IMUP = IMUF + 1
  ! surface tension: \sigma
  integer, parameter :: ISGM = IMUP + 1

  ! flags of KH/RT submodels
  integer, parameter :: IFKH = ISGM + 1 ! relate to radius of big son parcel
  integer, parameter :: IFRT = IFKH + 1 ! relate to cumulative time

  ! maximum number of user-defined array indexes
  integer(RKD), parameter :: NIND = 15
  !-----------------------------------------------------------------------------

  ! initialized from Code_Saturne Lagrangian parcel
  ! Lagrangian time step
  real(RKD) :: dt_lagr = ZERO
  !-----------------------------------------------------------------------------

  ! public procedures, called out of module
  public :: hfy_khrt_read,            &
            hfy_khrt_write,           &
            hfy_khrt_choose_submodel, &
            hfy_khrt_finalize

contains
  subroutine hfy_khrt_read()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! read info from cs platfrom, before applied KHRT model to the parcels
    !---------------------------------------------------------------------------
    use lagran, only: dtp, nbpart,                                &
                      lagr_resize_particle_set,                   &
                      ipepa, jisor,                               &
                      pepa, jrpoi, jrkh, jrrt,                    &
                      eptp, jdp, jmp, jup, jvp, jwp, juf, jvf, jwf
    use numvar, only: icrom, iprpfl, iviscl
    use field,  only: field_get_val_s
    !---------------------------------------------------------------------------
    intrinsic :: size, allocated

    integer, parameter :: NOVERSIZE = 5000
    integer :: ndim, ndim1, ndim2, newdim, olddim, nsize
    integer :: ip, iel
    logical :: yes
    real(RKD), dimension(:), pointer :: crom
    real(RKD), dimension(:), pointer :: viscl
    !---------------------------------------------------------------------------
    ! initialization of KHRT model

    ! Lagrange time step
    dt_lagr = dtp

    ! The number of parcels already exist in domain
    numb_local_old_parc = nbpart

    numb_newly_gen_parc = 0
    !---------------------------------------------------------------------------

    ! request memory, for both old and new parcels
    yes = allocated( prop_old_parc )

    if( .not. yes ) then
      ndim = nbpart + NOVERSIZE

      allocate( prop_old_parc(1:NIND, 1:ndim), source=ZERO )
      allocate( prop_new_parc(1:NIND, 1:ndim), source=ZERO )
    else
      olddim = size( prop_old_parc, dim=2 )

      if( olddim < numb_local_old_parc ) then
        deallocate( prop_old_parc )
        deallocate( prop_new_parc )

        newdim = nbpart + NOVERSIZE

        allocate( prop_old_parc(1:NIND, 1:newdim), source=ZERO )
        allocate( prop_new_parc(1:NIND, 1:newdim), source=ZERO )
      else
        prop_old_parc(:, :) = ZERO
        prop_new_parc(:, :) = ZERO
      end if
    end if

    !---------------------------------------------------------------------------
    ! read parcels properties

    ! density field of continuous phase
    call field_get_val_s( icrom, crom )
    ! Molecular dynamic viscosity, ref user/cs_user_physical_parameters.f90
    call field_get_val_s( iprpfl(iviscl), viscl )

    do ip = 1, numb_local_old_parc
      ! dimameter, number attribute, mass
      prop_old_parc(IKDP, ip) = eptp(jdp,   ip)
      prop_old_parc(IKNP, ip) = pepa(jrpoi, ip)
      prop_old_parc(IKMP, ip) = eptp(jmp,   ip) * pepa(jrpoi, ip)

      ! velocity of this parcel
      prop_old_parc(IKUP, ip) = eptp(jup, ip)
      prop_old_parc(IKVP, ip) = eptp(jvp, ip)
      prop_old_parc(IKWP, ip) = eptp(jwp, ip)

      ! fluid velocity seen by parcel
      prop_old_parc(IKUF, ip) = eptp(juf, ip)
      prop_old_parc(IKVF, ip) = eptp(jvf, ip)
      prop_old_parc(IKWF, ip) = eptp(jwf, ip)

      ! gas density, dynamic viscosity of gas and parcel, surface tension
      ! locally seen by parcel
      iel = ipepa(jisor, ip)

      prop_old_parc(IKRF, ip) = crom(iel)   ! gas density
      prop_old_parc(IMUF, ip) = viscl(iel)  ! gas dyn viscosity
      prop_old_parc(IMUP, ip) = MU_PARC     ! parcel dyn viscosity
      prop_old_parc(ISGM, ip) = SIGMA_GP    ! surface tension

      ! flags of KHRT model
      prop_old_parc(IFKH, ip) = pepa(jrkh, ip)  ! radius of big son parcel
      prop_old_parc(IFRT, ip) = pepa(jrrt, ip)  ! cumulative time in RT submodel
    end do

    ! check for particle count limit, assuming that all parcels will break up
    ! by KH submodel
    nsize = numb_local_old_parc * 2
    if( lagr_resize_particle_set(nsize) < 0 ) then
      ! btfprint instead
      print *, "Reach the Limit of Number of Parcels Set, KHRT Model."
      return
    end if

  end subroutine hfy_khrt_read

  subroutine hfy_khrt_write()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! write the new parcels, if any, to cs platfrom, after applied KHRT model to
    ! the parcels
    !---------------------------------------------------------------------------
    use lagran, only: nbpart, nbptot, dnbpar,                                   &
                      ipepa, jisor,                                             &
                      pepa,  jrpoi, jrkh, jrrt, jrtsp,                          &
                      eptp,  jdp,   jmp,  jup,  jvp,  jwp, juf, jvf, jwf

    use parall, only: parcpt, irangp
    !---------------------------------------------------------------------------
    integer :: ip, nip, iel, nbp_gen_sum
    real(RKD) :: dnbp_gen
    !---------------------------------------------------------------------------

    ! for old parcels
    ! if the parent parcel breaks up with KH submodel, now these arrays contain
    !   properties of its big son parcel
    ! if the parent parcel breaks up with RT submodel, now these arrays contain
    !   the properties of its only son parcel
    do ip = 1, numb_local_old_parc
      ! dimameter, number attribute, mass
      eptp(jdp,   ip) = prop_old_parc(IKDP, ip)
      pepa(jrpoi, ip) = prop_old_parc(IKNP, ip)
      eptp(jmp,   ip) = prop_old_parc(IKMP, ip) / prop_old_parc(IKNP, ip)

      ! velocity of this parcel
      eptp(jup, ip) = prop_old_parc(IKUP, ip)
      eptp(jvp, ip) = prop_old_parc(IKVP, ip)
      eptp(jwp, ip) = prop_old_parc(IKWP, ip)

      ! fluid velocity seen by parcel
      !eptp(juf, ip) = prop_old_parc(IKUF, ip)
      !eptp(jvf, ip) = prop_old_parc(IKVF, ip)
      !eptp(jwf, ip) = prop_old_parc(IKWF, ip)

      ! gas density and viscosity seen by parcel, and parcel's viscosity, and
      ! surface tension between gas and parcel, unchage.

      ! flags of KHRT model
      pepa(jrkh, ip) = prop_old_parc(IFKH, ip)  ! radius of big son parcel
      pepa(jrrt, ip) = prop_old_parc(IFRT, ip)  ! cumulative of time
    end do

    ! for newly generated parcels
    do ip = 1, numb_newly_gen_parc
      nip = ip + numb_local_old_parc

      ! you should make sure that all properties other than those relating KHRT
      ! have been written in eptp/pepa

      ! dimameter, number attribute, mass
      eptp(jdp,   nip) = prop_new_parc(IKDP, ip)
      pepa(jrpoi, nip) = prop_new_parc(IKNP, ip)
      eptp(jmp,   nip) = prop_new_parc(IKMP, ip) / prop_new_parc(IKNP, ip)

      ! velocity of this parcel
      eptp(jup, nip) = prop_new_parc(IKUP, ip)
      eptp(jvp, nip) = prop_new_parc(IKVP, ip)
      eptp(jwp, nip) = prop_new_parc(IKWP, ip)

      ! density of gas seen by parcel, viscosities of gas and liquid parcel,
      ! and surface tension between parcel and gas there, not change

      ! flags of new parcels
      pepa(jrkh, nip) = prop_new_parc(IFKH, ip) ! not ZERO
      pepa(jrrt, nip) = prop_new_parc(IFRT, ip) ! or ZERO

      pepa(jrtsp, nip) = ZERO
    end do

    ! update nbpart
    nbpart = nbpart + numb_newly_gen_parc

    ! update nbptot
    nbp_gen_sum = numb_newly_gen_parc

    if( irangp .ge. 0 ) then
      ! MPI_Reduce and then broadcast if parallel computing
      call parcpt( nbp_gen_sum )
    end if

    nbptot = nbptot + nbp_gen_sum

    ! update dnbpar
    dnbp_gen = ZERO
    do ip = 1, numb_newly_gen_parc
      dnbp_gen = dnbp_gen + prop_new_parc(IKNP, ip)
    end do

    dnbpar = dnbpar + dnbp_gen

  end subroutine hfy_khrt_write

  subroutine hfy_khrt_choose_submodel()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! According to the flags computed in this subroutine, check which submodel
    ! to work or neither work.
    !---------------------------------------------------------------------------
    use lagran, only: ipepa, pepa, eptp
    use cstphy, only: gravx => gx, gravy => gy, gravz => gz
    !---------------------------------------------------------------------------
    intrinsic :: min, exp, dot_product, abs, sqrt

    ! local name-constants
    real(RKD), parameter :: B0   = 0.61_RKD
    real(RKD), parameter :: B1   = 1.73_RKD ! KH, maybe [1.0, 60.0]
    real(RKD), parameter :: CRT  = 0.13_RKD  !0.1_RKD  ! RT
    real(RKD), parameter :: CTAU = HALF     ! RT
    real(RKD), parameter :: WECR = 12.0_RKD
    ! if detached-mass / parent parcel mass >= YCR, trigger KH submodel
    real(RKD), parameter :: YCR = 0.05_RKD

    ! local variables
    integer :: ip, nip, iel
    real(RKD) :: Ta, Oh, Weg, Wel, Rel
    real(RKD) :: tau_KH, Omega_RT, tau_RT
    real(RKD) :: diam_p, radi_p, diam_b, rho_l, magn, gt, krt, mass_p, numb_p
    real(RKD) :: ratio_m, rtemp, t1, t2, t3, t4, min_l, min_r, rho_g, mu_g
    real(RKD) :: rep, rep1, rep2, sigma_local, mu_this_parc
    real(RKD), dimension(1:3) :: u_relat
    !---------------------------------------------------------------------------
! *** debug sta ***
integer, save :: ipass = 0
ipass = ipass + 1
! *** debug end ***
    ! loop over each parcel

    do ip = 1, numb_local_old_parc

      ! some local temporary variable

      ! surface tension for this parcel
      sigma_local  = prop_old_parc(ISGM, ip)
      mu_this_parc = prop_old_parc(IMUP, ip)
      mu_g         = prop_old_parc(IMUF, ip)
      rho_g        = prop_old_parc(IKRF, ip)

      diam_p = prop_old_parc(IKDP, ip)
      radi_p = prop_old_parc(IKDP, ip) * HALF

      mass_p = prop_old_parc(IKMP, ip)
      numb_p = prop_old_parc(IKNP, ip)
      !-------------------------------------------------------------------------
      ! *** check if update flags of RT submodel ***

      ! RT branch first

      ! vel of parcel relative to fluid vel seen by this parcel
      u_relat = [prop_old_parc(IKUP, ip) - prop_old_parc(IKUF, ip), &
                 prop_old_parc(IKVP, ip) - prop_old_parc(IKVF, ip), &
                 prop_old_parc(IKWP, ip) - prop_old_parc(IKWF, ip)]

      magn = magn_vec_( u_relat )

      rep1 = ONE
      rep2 = 1000.0_RKD
      rep  = rho_g * magn * diam_p / mu_g

      ! density of this parcel

      rho_l = DENSLIQ

      ! drag, ref 3.2.2.4-25, we have rearanged the expression, no explicit
      ! expression of CD

      if( rep < rep1 ) then
        gt = 18_RKD * magn * mu_g / (rho_l * diam_p*diam_p)
      else if( rep1<=rep .and. rep<rep2 ) then
        gt = 18_RKD * magn * mu_g / (rho_l * diam_p*diam_p) &
           * (ONE + 0.15_RKD*rep**0.687_RKD)
      else
        ! rep >= rep2
        gt = 0.33_RKD * rho_g * magn*magn / (rho_l * diam_p)
      end if

      !< if gravity acceleration added

      !< use cstphy in base/
      !< Generally, the gravity should read from Code_Saturne. 
      !< Note: Drag force is along the NEGATIVE position of \vec u_relat
      !< if no gravity, comment on the line below
      gt = gt + dot_product( [gravx, gravy, gravz], u_relat ) / magn

      !< the wave number, 3.2.2.4-25
      t1 = abs( gt )
      t2 = rho_l - rho_g
      t3 = THREE * sigma_local
      krt = sqrt( t1 * t2 / t3 )

      !< wave length, 3.2.2.4-27
      Lambda_RT = TWO * PI * CRT / krt

      if( Lambda_RT < diam_p ) then
        prop_old_parc(IFRT, ip) = prop_old_parc(IFRT, ip) + dt_lagr
      end if

      ! *** add code to compute tau_RT

      Omega_RT                                        &
      = sqrt(2*(abs(gt*(rho_l-rho_g)))**1.5_RKD       &
      / (3*(3*sigma_local)**0.5_RKD * (rho_g + rho_l)))

      tau_RT = CTAU / Omega_RT
      !-------------------------------------------------------------------------

      ! *** check if update flag of KH submodel ***

      Weg = rho_g * (magn**TWO) * radi_p / sigma_local  !< 3.2.2.4-6
      Wel = rho_l * (magn**TWO) * radi_p / sigma_local  !< 3.2.2.4-7
      Rel = rho_l * magn * radi_p / mu_this_parc        !< 3.2.2.4-9

      Oh = sqrt( Wel ) / Rel  !< 3.2.2.4-8
      Ta = Oh * sqrt( Weg )   !< 3.2.2.4-10

      !< 3.2.2.4-4
      t1 = 9.02_RKD * radi_p
      t2 = ONE + 0.45_RKD *Oh**HALF
      t3 = ONE + 0.4_RKD * Ta**0.7_RKD
      t4 = (ONE + 0.865_RKD * Weg**1.67_RKD)**0.6_RKD

      Lambda_KH = t1 * t2 * t3 / t4

      !< 3.2.2.4-5
      t1 = 0.34_RKD + 0.385_RKD * Weg**1.5_RKD
      t2 = ONE + Oh
      t3 = ONE + 1.4_RKD*Ta**0.6_RKD
      t4 = sqrt( sigma_local / (rho_l*radi_p**THREE) )

      Omega_KH = t1 / (t2*t3) * t4

! *** debug sta ***
!print *, 'l450, ipass/ip = ', ipass, ip
!print *, 't1  = ', t1,  '\',&
!         'weg = ', weg, &
!         't2  = ', t2,  &
!         'Oh  = ', Oh,  &
!         'Ta  = ', Ta,  &
!         'sgml = ', sigma_local,  &
!         'rho_l = ', rho_l, &
!         'radi_p = ', radi_p
! *** debug end ***
      ! 3.2.2.4-2

      rtemp = B0 * Lambda_KH

      if( rtemp <= diam_p ) then
        radi_s_KH = rtemp
      else
        !< second part in 3.2.2.4-2

        t1 = THREE * PI * radi_p*radi_p * magn
        t2 = TWO * Omega_KH
        min_l = (t1 / t2)**THIRD

        t1 = THREE * radi_p*radi_p * Lambda_KH
        min_r = (t1 * QUARTER)**THIRD

        !< if B0 * Lambda > radi_p

        radi_s_KH = min( min_l, min_r )
      end if

      ! 3.2.2.4-3
      tau_KH = 3.726_RKD * B1 * radi_p / (Lambda_KH * Omega_KH)

      if( Weg > WECR ) then

        ! 3.2.2.4-12
        prop_old_parc(IFKH, ip)                   &
          = radi_s_KH                             &
          + (prop_old_parc(IFKH, ip) - radi_s_KH) &
          * exp( -dt_lagr/tau_KH )
      end if

      ! radius of big son parcel
      radi_b_KH = prop_old_parc(IFKH, ip)
      !-------------------------------------------------------------------------

      ! *** if choose RT submodel ***

      if( prop_old_parc(IFRT, ip) > tau_RT ) then
        call hfy_khrt_rt_( ip )

        cycle
      end if
      !-------------------------------------------------------------------------

      ! *** if choose KH submodel ***

      ratio_m = ZERO

      diam_b  = TWO * prop_old_parc(IFKH, ip)
      ratio_m = ONE - (diam_b / diam_p)**THREE

      if( ratio_m > YCR ) then
        numb_newly_gen_parc = numb_newly_gen_parc + 1

        ! first copy all infos to this new parcels
        nip = numb_newly_gen_parc + numb_local_old_parc
        pepa(:, nip) = pepa(:, ip)
        eptp(:, nip) = pepa(:, ip)
        ipepa(:, nip) = ipepa(:, ip)

        call hfy_khrt_kh_( ip )

      end if
      !-------------------------------------------------------------------------

      ! *** no choose any of the breaking submodel ***

    end do

  end subroutine hfy_khrt_choose_submodel

  subroutine hfy_khrt_kh_( ip )
  !=============================================================================
    !---------------------------------------------------------------------------
    ! Kelvin-Holmholtz submodel:
    ! ONE parent parcel -> ONE big son parcel + ONE small son parcel

    ! both old and new parcels have some new properties
    !---------------------------------------------------------------------------
    ! dummy arguments
    integer, intent(in) :: ip

    ! local variables
    real(RKD) :: mass_p, numb_p, diam_p, rho_p, dm, radi_p, radi_b, mass_s
    real(RKD) :: diam_b, rho_b, mass_b, diam_s, rho_s
    real(RKD), dimension(1:3) :: vel_p, vel_s
    integer :: nip
    !---------------------------------------------------------------------------

    ! temperarily store the mass, number and diam
    mass_p = prop_old_parc(IKMP, ip)
    numb_p = prop_old_parc(IKNP, ip)
    diam_p = prop_old_parc(IKDP, ip)

    radi_p = diam_p * HALF

    ! density parent parcel
    rho_p = DENSLIQ

    nip = numb_newly_gen_parc
    !---------------------------------------------------------------------------

    ! diameter
    prop_old_parc(IKDP, ip ) = radi_b_KH * TWO
    prop_new_parc(IKDP, nip) = radi_s_KH * TWO

    ! mass
    radi_b = radi_b_KH

    prop_new_parc(IKMP, nip) = mass_p * (ONE - (radi_b/radi_p)**THREE)
    prop_old_parc(IKMP, ip ) = mass_p - prop_new_parc(IKMP, nip)

    ! number-attribute
    mass_s = prop_new_parc(IKMP, nip)
    rho_s  = rho_p
    diam_s = prop_new_parc(IKDP, nip)
    prop_new_parc(IKNP, nip) = SIX * mass_s / (rho_s * PI * diam_s**THREE)

    mass_b = prop_old_parc(IKMP, ip)
    diam_b = prop_old_parc(IKDP, ip)
    rho_b  = rho_p

    prop_old_parc(IKNP, ip) = SIX * mass_b / (rho_b * PI * diam_b**THREE)
    !---------------------------------------------------------------------------

    ! velocities of fluid seen by parcel, these kind of vels of the two son
    ! parcels are both inherited from their parent parcel
    !prop_new_parc(IKUF, nip) = prop_old_parc(IKUF, ip)
    !prop_new_parc(IKVF, nip) = prop_old_parc(IKVF, ip)
    !prop_new_parc(IKWF, nip) = prop_old_parc(IKWF, ip)
    !---------------------------------------------------------------------------

    ! velocities of the two son parcels, 3.2.2.4-17 -- 3.2.2.4-23
    mass_s = prop_new_parc(IKMP, nip)
    vel_p = [prop_old_parc(IKUP, ip),  &
             prop_old_parc(IKVP, ip),  &
             prop_old_parc(IKWP, ip)]

    ! after this subroutine, vel_p will changed to be the value of big son parc
    call vel_kh_( mass_p, mass_s, vel_p, vel_s )

    ! small son parcel
    prop_new_parc(IKUP, nip) = vel_s(1)
    prop_new_parc(IKVP, nip) = vel_s(2)
    prop_new_parc(IKWP, nip) = vel_s(3)

    ! big son parcel
    prop_old_parc(IKUP, ip) = vel_p(1)
    prop_old_parc(IKVP, ip) = vel_p(2)
    prop_old_parc(IKWP, ip) = vel_p(3)
    
    ! flags
    ! new parcel
    prop_new_parc(IFRT, nip) = ZERO
    prop_new_parc(IFKH, nip) = HALF * prop_new_parc(IKDP, nip) ! not ZERO

    ! back to zero for the two flags of old parcel
    prop_old_parc(IFRT, ip) = ZERO
    prop_old_parc(IFKH, ip) = HALF * prop_old_parc(IKDP, ip)

  end subroutine hfy_khrt_kh_

  subroutine hfy_khrt_rt_( ip )
  !=============================================================================
    !---------------------------------------------------------------------------
    ! Rayleigh-Taylor summodel: ONE parent parcel -> ONE son parcel
    ! The parcel would only change its properties below:
    ! (1) diameter
    ! (2) number attribute
    !---------------------------------------------------------------------------

    ! dummy arguments
    integer, intent(in) :: ip
    ! local variables
    real(RKD) :: rho, diam_p, mass_p, numb_p, mass, numb, diam
    !---------------------------------------------------------------------------

    mass_p = prop_old_parc(IKMP, ip)
    numb_p = prop_old_parc(IKNP, ip)
    diam_p = prop_old_parc(IKDP, ip)

    ! density
    rho = DENSLIQ
    !---------------------------------------------------------------------------

    ! new diameter
    prop_old_parc(IKDP, ip) = Lambda_RT

    ! new number attribute of parcel, or 'weight of particle' in Code_Saturne
    mass = mass_p
    diam = prop_old_parc(IKDP, ip)
    prop_old_parc(IKNP, ip) = SIX * mass / (rho * PI * diam**THREE)

    ! set the two flags
    prop_old_parc(IFKH, ip) = HALF * prop_old_parc(IKDP, ip)  !not ZERO
    prop_old_parc(IFRT, ip) = ZERO

  end subroutine hfy_khrt_rt_

  subroutine hfy_khrt_finalize ()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! free memories
    !---------------------------------------------------------------------------

    if( allocated(prop_old_parc) ) deallocate( prop_old_parc )
    if( allocated(prop_new_parc) ) deallocate( prop_new_parc )

  end subroutine hfy_khrt_finalize

  !*****************************************************************************
  ! internal procedures
  function magn_vec_( vec_1d )
  !=============================================================================
    !---------------------------------------------------------------------------
    ! Computing the magnitude of a 3D vector

    ! var         intent  attr         rank       size
    ! vec_1d      in      real(RKD)    (:)        1:3
    ! vec_magn_   out     real(RKD)    -          -
    !---------------------------------------------------------------------------
    intrinsic :: sqrt
    ! dummy argument  and function result return
    real(RKD), dimension(1:3), intent(in) :: vec_1d
    real(RKD) :: magn_vec_
    ! local variables
    integer :: i

    magn_vec_ = ZERO
    do i = 1, 3
      magn_vec_ = magn_vec_ + vec_1d(i) * vec_1d(i)
    end do

    magn_vec_ = sqrt ( magn_vec_ )
  end function magn_vec_

  subroutine vel_kh_( mp, ms, vel_p, vel_s )
  !=============================================================================
    !---------------------------------------------------------------------------
    ! internal private procedure, used only in this module
    ! for KH submodel, to computing the velocity vectors of the two son parcels

    ! var     intent    attr        rank      size    denote
    ! mp      in        real(RKD)   -         -       Mass of Parent parcel
    ! ms      in        real(RKD)   -         -       Mass of Small son parcel
    ! vel_p   inout     real(RKD)   (:)       1:3     VELocity of Parent parcel
    ! vel_s   out       real(RKD)   (:)       1:3     VELocity of Small son parcel

    ! initially vel_p has values of parcel parcel vector, but would change to
    !  velocity vector of big son parcel
    !---------------------------------------------------------------------------
    !< intrinsic proceudre
    intrinsic :: cos, sin

    ! arguments
    real(RKD), intent(in) :: mp, ms
    real(RKD), dimension(1:3), intent(inout) :: vel_p
    real(RKD), dimension(1:3), intent(out) :: vel_s

    ! local parameter
    real(RKD), parameter :: C1 = ZERO !0.188_RKD !10.0_RKD !0.01_RKD !< KH, 3.2.2.4-17
    ! local variables
    real(RKD) :: magn, magnsq, xi, dun, t1, t2
    real(RKD), dimension(1:3) :: l1, l2, ln

    ! first direction vector orthogonal with vel_p, 3.2.2.4-18
    magn = magn_vec_( vel_p ); magnsq = magn * magn

    l1(1) = ONE - vel_p(1) * vel_p(1) / magnsq
    l1(2) =     - vel_p(1) * vel_p(2) / magnsq
    l1(3) =     - vel_p(1) * vel_p(3) / magnsq

    ! normalize \vec l1
    l1(1:3) = l1(1:3) / magn_vec_( l1 )

    ! the second direction vector, orthogonal with \vec vec_p & \vec l1
    ! 3.2.2.4-19
    magn = magn_vec_( vel_p )
    call vec_product_( vel_p/magn, l1, l2 )

    ! a random number, meeting uniform distribution
    call random_number( xi )

    ! 3.2.2.4-20
    t1 = TWO * PI * xi
    ln = cos( t1 )*l1 + sin( t1 )*l2

    ! 3.2.2.4-17
    dun = C1 * Lambda_KH * Omega_KH

    !< 3.2.2.4-21
    magn = magn_vec_( vel_p )
    t1 = dun / magn
    t2 = ONE - t1*t1
    ! velocity of small son parcel, 3.2.2.4-22
    vel_s = sqrt( t2 )*vel_p + dun * ln

    ! velocity of big son parcel, 3.2.2.4-23
    vel_p = (mp*vel_p - ms*vel_s) / (mp-ms)

  end subroutine vel_kh_

  subroutine vec_product_( vec1, vec2, vec_product )
  !=============================================================================
    !---------------------------------------------------------------------------
    ! internal procedure
    ! vector product of \vec vec1 and \vec vec2, the restult is also a vector
    !   \vec vec_product
    !---------------------------------------------------------------------------

    !< dummy argument
    real(RKD), dimension(1:3), intent(in) :: vec1, vec2
    real(RKD), dimension(1:3), intent(out) :: vec_product

    vec_product(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
    vec_product(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
    vec_product(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)
  end subroutine vec_product_

end module hfy_khrt_mod
