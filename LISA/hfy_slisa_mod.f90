
module hfy_slisa_mod
  implicit none
  !-----------------------------------------------------------------------------
  ! author    date        aff
  ! snx       21.06.30    simpop.cn

  ! instruction to certain suffix

  ! IN: used only in this module
  ! G : global for this module
  !-----------------------------------------------------------------------------
  private
  save
  intrinsic :: acos, selected_real_kind
  !-----------------------------------------------------------------------------

  integer, parameter :: RKD = selected_real_kind( p=10 )
  !-----------------------------------------------------------------------------

  ! Named-Constants only used within this module
  ! Real

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

  real(RKD), parameter :: QUART = 0.25_RKD
  real(RKD), parameter :: THIRD = 1.0_RKD / 3.0_RKD
  real(RKD), parameter :: HALF  = 0.5_RKD

  real(RKD), parameter :: PI = acos( -ONE )
  !-----------------------------------------------------------------------------

  ! Abstract Data Type definition

  type, private :: Parcel_t

    ! general attriutes
    real(RKD), dimension(1:3) :: pos_ = [ZERO, ZERO, ZERO]
    real(RKD), dimension(1:3) :: vel_ = [ZERO, ZERO, ZERO]

    real(RKD) :: mass_ = ZERO
    real(RKD) :: diam_ = ZERO
    real(RKD) :: numb_ = ZERO

    !---------------------------------------------------------------------------

    ! These Attributes Below Would Be Added Into eptp/pepa/ipepa

    ! liquidCore = 1 by default
    integer :: iLiqC_ = 1

    ! random number \in [0, 1], attached to parcel, used when breaking with RR
    ! distribution
    real(RKD) :: rand_ = ZERO

    ! breaking length
    !real(RKD) :: bLen_ = ZERO

    ! initial info
    real(RKD) :: dNozIni_ = ZERO  ! diam of noz from which the parc be injected
    real(RKD) :: uConIni_ = ZERO  ! speed of parc along the spray cone
    real(RKD) :: thetIni_ = ZERO  ! the initial half-cone angle of this parc
  end type Parcel_t

  type, private :: Rands_t
    integer :: irand_ = 0 ! index of rands_
    integer :: nsize_ = 0 ! size of rands_

    real(RKD), dimension(:), pointer :: rands_ => NULL()
  end type Rands_t

  type, public :: NozzleSLS_t
    !---------------------------------------------------------------------------
    ! INPUT
    ! from user

    ! time for nozzle beginning and ending to work
    real(RKD) :: staTime_ = ZERO
    real(RKD) :: endTime_ = ZERO

    ! density of liquid droplets
    real(RKD) :: dens_ = ZERO
    ! surface tension
    real(RKD) :: tens_ = ZERO
    ! viscosity of liquid
    real(RKD) :: visc_ = ZERO

    ! time rate of mass of injecting
    real(RKD) :: massRate_ = ZERO
    ! pressure difference
    real(RKD) :: presDiff_ = ZERO

    ! diameter of this nozzle
    real(RKD) :: diam_ = ZERO

    ! position of nozzle face
    real(RKD), dimension(1:3) :: posi_ = [ZERO, ZERO, ZERO]
    ! orientation vector normal to nozzle face: direction
    real(RKD), dimension(1:3) :: dire_ = [ZERO, ZERO, ZERO]

    ! HalF-COne ANGLe
    real(RKD) :: hfCoAngl_ = ZERO
    ! dispersion radius, for each nozzle
    real(RKD) :: dispAngl_ = ZERO

    ! the number of parcels scheduled to be injected
    ! 1, The final value may be marginally altered
    integer :: numbTotParcs_ = 0
    !---------------------------------------------------------------------------

    ! OUTPUT
    ! these are also the properties of this nozzle, but they would be abtained
    ! from properties above input by user

    ! magnitude of injecting velocity
    real(RKD) :: magnVel_ = ZERO

    ! parcels to be injected this time step and its number
    type(Parcel_t), dimension(:), pointer :: parcsInj_ => NULL()
    integer :: numbInjParcs_ = 0

    ! nozzle-attached random number series, and the number of series 
    type(Rands_t), dimension(:), pointer :: rands_ => NULL()
    ! number of random number series needed
    integer :: numbRands_ = 3

  end type nozzleSLS_t
  !-----------------------------------------------------------------------------

  ! number of nozzles, evaluated in uslag1.f90 

  integer, public :: numbNozz_SLS = 0

  ! the nozzles in sLISA, open outside module

  type(NozzleSLS_t), dimension(:), pointer :: nozzles_SLS => NULL()
  public :: nozzles_SLS

  ! Lagrangian time step

  real(RKD) :: dtLagr_ING = ZERO

  ! density field of gaseous phase

  real(RKD), dimension(:), pointer :: densg_ING => NULL()
  !-----------------------------------------------------------------------------

  ! Subroutines are by default private

  ! These Subroutines Are Made Public
  public ::         &
    HFY_SLS_Ini,    &
    HFY_SLS_Fin,    &
    HFY_SLS_Inject, &
    HFY_SLS_Atomize

  !private ::
  ! CheckIn_,
  ! CheckIrand_,
  ! CalcIniFilmThickness_
  ! CalcInjSpeed_
  ! CalcMagnVec_,
  ! CalcVecProduct_,
  ! GenRands_,
  ! Resize_,
  ! UpdateNoz_,
contains
  subroutine HFY_SLS_Ini()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! Allocate nozzles of type NozzleSLS_t: nozzles_SLS
    !---------------------------------------------------------------------------

    use lagran, only: dtp
    use field,  only: field_get_val_s
    use numvar, only: icrom
    !---------------------------------------------------------------------------

    intrinsic :: allocated

    integer, parameter :: NOVERSIZE = 100
    integer, save :: ipass = 0
    integer :: inoz, nnozs
    logical :: yes
    type(NozzleSLS_t), pointer :: thisNoz => NULL()
    !---------------------------------------------------------------------------

    ipass = ipass + 1
    !---------------------------------------------------------------------------

    ! Lagragian time step

    dtLagr_ING = dtp

    ! density field of gaseous phase

    call field_get_val_s(icrom, densg_ING)
    !---------------------------------------------------------------------------

    nnozs = numbNozz_SLS

    ! Init nozzles

    if ( ipass == 1 ) then 

      allocate (nozzles_SLS(1:nnozs))

      do inoz = 1, nnozs

        thisNoz => nozzles_SLS(inoz)

        allocate( thisNoz%parcsInj_(NOVERSIZE) )

        thisNoz%numbRands_ = 3

        allocate( thisNoz%rands_(3) )

      end do
    end if

    ! Read parameters

    if( ipass == 1) then

      call hfy_slisa_in()

      ! Check and modified if need
      ! noz%dire_
      ! noz%numbTotParcs_ ! TODO:

      call CheckIn_()

      ! each nozzles' injecting velocity

! *** DEBUG STA ***
      !call CalcInjSpeed_()
      thisNoz%magnVel_ = 79.556_RKD
! *** DEBUG END ***
      ! Generate a data pool of random number, for sLISA model

      call GenRands_()
    end if

  end subroutine HFY_SLS_Ini

  subroutine HFY_SLS_Inject()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! Inject sheet parcels
    !---------------------------------------------------------------------------

    intrinsic :: abs, atan, cos, exp, int, log, sin, sqrt, tan

    real(RKD), parameter :: EPS = 1.0e-2_RKD

    integer :: inoz, nparcs, nnozs, iptr, iparc

    real(RKD) :: dens, rand, modl
    real(RKD) :: v1, v2, v3, r1, r2, alph, radi, thet, v_swirl, v_cone
    real(RKD) :: pis6, diam, dtLagr, thet1, thet2

    real(RKD), dimension(1:3) :: t1_vec, t2_vec, nr_vec, n_cone_vec, nt_vec

    type(NozzleSLS_t), pointer :: thisNoz => NULL()
    type(Parcel_t),    pointer :: theParc => NULL()
    type(Rands_t),     pointer :: myRands => NULL() 
    !---------------------------------------------------------------------------

    ! Update:
    ! 1) noz%numbInjParcs_
    ! 2) call Resize_(): allocate or resize noz%parcsInj_

    call UpdateNoz_()
    !---------------------------------------------------------------------------

    ! temporarily keep
    ! 1) number of nozzles
    ! 2) Lagrange time step

    nnozs = numbNozz_SLS

    dtLagr = dtLagr_ING
    !---------------------------------------------------------------------------

    ! the mass

    do inoz = 1, nnozs
      thisNoz => nozzles_SLS(inoz)

      nparcs = thisNoz%numbInjParcs_

      thisNoz%parcsInj_(:)%mass_ = ZERO

      ! for each each parcel

      thisNoz%parcsInj_(1:nparcs)%mass_ = thisNoz%massRate_ * dtLagr / nparcs
    end do
    !---------------------------------------------------------------------------

    ! thickness of sheet serves as initial diameter of sheet parcel
    ! nozzles(:)%parcsInj_(:)%diam_

! *** DEBUG STA ***
    !call CalcIniFilmThickness_()
    thisNoz%parcsInj_(1:nparcs)%diam_ = thisNoz%diam_ * 0.98_RKD
! *** DEBUG END ***
    !---------------------------------------------------------------------------

    ! parcel-attached random number

    do inoz = 1, nnozs
      thisNoz => nozzles_SLS(inoz)

      nparcs = thisNoz%numbInjParcs_

      do iparc = 1, nparcs
        theParc => thisNoz%parcsInj_(iparc)

        ! first sequence of random number for diameter
        myRands => thisNoz%rands_(1)

        myRands%irand_ = myRands%irand_ + 1
        iptr = myRands%irand_

        theParc%rand_ = myRands%rands_(iptr)
      end do
    end do
    !---------------------------------------------------------------------------

    call CheckIrand_()
    !---------------------------------------------------------------------------

    ! number attribute, or weight

    pis6 = PI / SIX

    do inoz = 1, nnozs
      thisNoz => nozzles_SLS(inoz)

      dens = thisNoz%dens_

      do iparc = 1, nparcs
        theParc => thisNoz%parcsInj_(iparc)

        diam = theParc%diam_

        theParc%numb_ = theParc%mass_ / (dens * pis6 * diam**THREE)
      end do
    end do
    !---------------------------------------------------------------------------

    ! position and velocity of parcels to be injected

    do inoz = 1, nnozs

      thisNoz => nozzles_SLS(inoz)

      ! to define t1_vec
      ! temperarily store the nozzle face vector: \vec n_ref, see fig 3.6
      v1 = thisNoz%dire_(1)
      v2 = thisNoz%dire_(2)
      v3 = thisNoz%dire_(3)

      ! fix if .%dire_ = (1, 0, 0)
      modl = CalcMagnVec_( [v1-ONE, v2, v3] )

      if( modl < EPS ) then
        t1_vec = [ZERO, ONE, ZERO]
      else
        ! comput \vec t1: see eqn 3.2.2.1-1
        t1_vec = [ONE-v1*v1, -v1*v2, -v1*v3]

        ! normalization
        t1_vec(1:3) = t1_vec(1:3) / CalcMagnVec_( t1_vec(1:3) )
      end if

      ! \vec t2 = \vec n_ref \times \vec t1, see eqn 3.2.2.1-4
      call CalcVecProduct_( thisNoz%dire_, t1_vec, t2_vec )

      ! work log, a little difference with 3.2.2.1-40 & 41

      ! radius of nozzle acts as outer radius r2 of the annulus, inner radius r1
      ! can thus be computed by aditional angles of half-cone and dispersion

      ! r2: outer radius
      r2 = thisNoz%diam_ * HALF

      thet1 = thisNoz%hfCoAngl_ - HALF * thisNoz%dispAngl_
      thet2 = thisNoz%hfCoAngl_ + HALF * thisNoz%dispAngl_

      ! r1: inner radius
      r1 = r2 * tan( thet1 ) / tan( thet2 )

      nparcs = thisNoz%numbInjParcs_

      do iparc = 1, nparcs

        theParc => thisNoz%parcsInj_(iparc)
        !-------------------------------------------------------------------------

        ! position in polar coordinates, see 3.2.2.1-38 & 39

        ! for random circumferential angle
        myRands => thisNoz%rands_(2)

        myRands%irand_ = myRands%irand_ + 1

        iptr = myRands%irand_
        rand = myRands%rands_(iptr)

        alph = TWO * PI * rand

        ! for random length of radius
        myRands => thisNoz%rands_(3)

        myRands%irand_ = myRands%irand_ + 1

        iptr = myRands%irand_
        rand = myRands%rands_(iptr)

        radi = r1 + rand * (r2 - r1)
        !-------------------------------------------------------------------------

        ! 3.2.2.1-5, or 3.2.2.1-43
        nr_vec(1:3) =  cos( alph ) * t1_vec(1:3) + sin( alph ) * t2_vec(1:3)

        ! position of injecting points, see eqn 3.2.2.1-42
        theParc%pos_(1:3) = thisNoz%posi_(1:3) + radi * nr_vec(1:3)
        !-------------------------------------------------------------------------

        ! velocity of parcel

        ! ref 1: worklog 2021.04.02
        ! ref 2: 3.2.2.1-44

        thet = atan( radi / (thisNoz%diam_ * HALF) * tan(thet2) )
        theParc%thetIni_ = thet

        ! 3.2.2.1-45
        n_cone_vec(1:3)                     &
        = cos( thet ) * thisNoz%dire_(1:3)  &
        + sin( thet ) * nr_vec(1:3)

        ! 3.2.2.2.1-32
        v_swirl = thisNoz%magnVel_ * sin( thet )

        ! 3.2.2.2.1-33
        v_cone = thisNoz%magnVel_

        ! eqn 3.2.2.1-6
        nt_vec(1:3) = -sin( alph ) * t1_vec(1:3) + cos( alph ) * t2_vec(1:3)

        ! 3.2.2.2.1-39
        theParc%vel_(1:3) = v_cone * n_cone_vec(1:3) + v_swirl * nt_vec(1:3)
      end do
    end do

    call CheckIrand_()
    !---------------------------------------------------------------------------

    ! Other attributes of parcels to be injected 

    do inoz = 1, nnozs
      thisNoz => nozzles_SLS(inoz)

      nparcs = thisNoz%numbInjParcs_

      ! liquidCore = 1

      thisNoz%parcsInj_(1:nparcs)%iLiqC_ = 1

      ! initial diameter

      thisNoz%parcsInj_(1:nparcs)%dNozIni_ = thisNoz%diam_

      ! initial velocity along spay cone

      thisNoz%parcsInj_(1:nparcs)%uConIni_ = thisNoz%magnVel_
    end do

  end subroutine HFY_SLS_Inject

  subroutine HFY_SLS_Atomize()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! 1, Get the position of each sheet parcel, update its diameter and number

    ! 2, Compute Weg of each sheet parcel, confirm the atomization branch
    !   2.1, if Weg < WECR, long-wave branch, compute
    !     * breaking length, if break, compute
    !       # dL
    !       # dD
    !       # call Rosin-Rammler function
    !   2.2, if Weg > WECR, short-wave branch, compute
    !     * breaking length, if break, compute
    !       # dL
    !       # dD
    !       # call Rosin-Rammler function
    !---------------------------------------------------------------------------

    use lagran, only: nbpart,                                                 &
                      eptp,  jdp,    jmp,    jxp,    jyp,    jzp,    jup,     &
                             jvp,    jwp,    juf,    jvf,    jwf,             &
                      pepa,  jrpoi,  jrxp0,  jryp0,  jrzp0,  jrdnoz, jrsht0,  &
                             jrtens, jrvisc, jrucon, jrthet, jrrand, jrdens,  &
                      ipepa, jiliqc, jisor
    !---------------------------------------------------------------------------

    intrinsic :: exp, log, sin, sqrt

    real(RKD), parameter :: WECR = 1.6875_RKD !27.0_RKD/16.0_RKD
    real(RKD), parameter :: CTAU = 12.0_RKD

    integer :: iparc, nparcs, iliqc, iel, nexpn
    real(RKD) :: dist, t1, t2, mass, diam, sht0, diam0, numb0, coef, dens
    real(RKD) :: weg, densg, spdr, tens, dnoz, uConIni, thetIni, xj, xq
    real(RKD) :: densl, tau, expn, bLen, ks, omgs, dl, viscl, oh, dd
    real(RKD) :: rand, dmin, dmax, dbar, k_temp
    real(RKD), dimension(1:3) :: pos0, pos1, vgas, vliq
    !---------------------------------------------------------------------------

    nparcs = nbpart

    do iparc = 1, nparcs

      ! We only Manipulate the Sheet Parcels here

      iliqc = ipepa(jiliqc, iparc)

      if( iliqc == 0 ) cycle
      !-------------------------------------------------------------------------

      ! some common variables, used in this subroutine

      dnoz = pepa(jrdnoz, iparc)  ! diameter of nozzle
      sht0 = pepa(jrsht0, iparc)  ! sheet thickness, i.e. initial diam

      uConIni = pepa(jrucon, iparc)
      thetIni = pepa(jrthet, iparc) ! initial injecting angle

      ! surface tension, sheet parcels' viscocity
      tens  = pepa(jrtens, iparc)
      viscl = pepa(jrvisc, iparc)

      ! the density of gas there: densg
      iel = ipepa(jisor, iparc)
      densg = densg_ING(iel)

      !density of liquid (sheet) parcel: densl
      densl = pepa(jrdens, iparc)

      ! initial position
      pos0(1) = pepa(jrxp0, iparc)
      pos0(2) = pepa(jryp0, iparc)
      pos0(3) = pepa(jrzp0, iparc)

      ! present position
      pos1(1) = eptp(jxp, iparc)
      pos1(2) = eptp(jyp, iparc)
      pos1(3) = eptp(jzp, iparc)

      ! velocity of gas there
      vgas(1) = eptp(juf, iparc)
      vgas(2) = eptp(jvf, iparc)
      vgas(3) = eptp(jwf, iparc)

      ! velocity of the parcel
      vliq(1) = eptp(jup, iparc)
      vliq(2) = eptp(jvp, iparc)
      vliq(3) = eptp(jwp, iparc)

      ! speed of relative velocity
      spdr = CalcMagnVec_( vliq - vgas )

      ! old diameter of sheet parcel

      diam0 = eptp(jdp,   iparc)
      numb0 = pepa(jrpoi, iparc)

      ! mass, note the diffs between this module & CS platform

      mass = eptp(jmp, iparc) * pepa(jrpoi, iparc)

      ! random number attached to this sheet parcel

      rand = pepa(jrrand, iparc)
      !-------------------------------------------------------------------------

      ! Update the Thickness of Liquid Film

      ! the distance this parcel has flowed away

      dist = CalcMagnVec_( pos1 - pos0 )
      !-------------------------------------------------------------------------

      ! update diameter of this sheet parcel, by stretching
      ! 3.2.2.2.3-13 or 3.2.2.2.3-22

      t1 = QUART * (dnoz - sht0) * sht0
      t2 = HALF  * (dnoz - sht0) + dist * sin( thetIni )

      eptp(jdp, iparc) = TWO * t1 / (t2 + 1.0E-18)

      ! new number attribute
      ! density does not change, so n * d^3 = Const, 

      diam = eptp(jdp, iparc)

      pepa(jrpoi, iparc) = numb0 * (diam0/diam)**THREE
      !-------------------------------------------------------------------------

      ! According the local gas Weber number, confirm the breaking branch:
      ! Long wave or otherwise Short wave ?

      ! local gas Weber number, in Fluent, depends on
      ! 1, the liquid velocity,
      ! 2, gas density,
      ! 3, sheet half-thickness

      ! but in StarCCM+, it depends instead
      ! 1, the relative velocity of liquid sheet parcel

      ! the latter is exploited here

      diam = eptp(jdp, iparc)

      weg  = densg * spdr**TWO * (HALF*diam) / tens

      ! which breaking-mode ?

      if( weg < WECR ) then
        ! Long Wave Mode: L, dL, dD

        ! Breaking Time and Length

        ! worklog, 2021.04.09; in 3.2.2.2.3-9

        t1 = (dnoz - sht0) * sht0
        t2 = FOUR * uConIni * sin( thetIni )
        xj = t1 / t2

        xq = densg / densl

        ! breaking time

        expn = THIRD * TWO

        t1 = xj * tens
        t2 = xq**TWO * spdr**FOUR * densl

        ! 3.2.2.2.3-9
        tau = (THREE*CTAU)**expn * (t1/t2)**THIRD

        ! breaking distance this time step, 3.2.2.2.3-11
        bLen = uConIni * tau

        if( dist < blen ) cycle

        ! 3.2.2.2.3-6
        ks = densg * spdr**TWO / (TWO*tens)

        ! diameter of liquid ligment: 3.2.2.2.3-12
        diam = eptp(jdp, iparc)
        dl = sqrt( FOUR * diam / ks )

        ! Ohnesorge number, 3.2.2.2.3-24
        oh = viscl / sqrt( densl * tens * dl )

        ! 3.2.2.2.3-23
        expn = ONE / SIX; coef = 1.88_RKD
        dd = coef * dl * (ONE + THREE*oh)**expn

        ! Rosin-Rammler distribution

        nexpn = 3
        dbar  = dd
        dmin  = dd / TEN
        dmax  = dd

        k_temp = ONE - exp( (dmin/dbar)**nexpn - (dmax/dbar)**nexpn )

        eptp(jdp, iparc)                                                    &
        = dbar * ((dmin / dbar)**nexpn - log(ONE - k_temp*rand))**(ONE/nexpn)

        ipepa(jiliqc, iparc) = 0

      else
        ! Short Wave Mode: L, dL, dD

        xq = densg / densl

        call SolveOmegaKs_(ks, omgs)

        ! Breaking Time and Length

        ! breaking time, 3.2.2.2.3-19
        tau = CTAU / omgs

        ! breaking length, 3.2.2.2.3-20
        bLen = uConIni * tau

        if( dist < bLen ) cycle

        ! present diameter of sheet parc
        diam = eptp(jdp, iparc)

        ! diameter of liquid ligament

        dl = sqrt( EIGHT * diam / ks )

        ! Ohnesorge number, 3.2.2.2.3-24

        oh = viscl / sqrt( densl * tens * dl )

        ! 3.2.2.2.3-23

        expn = ONE / SIX
        coef = 1.88_RKD

        dd = coef * dl * (ONE + THREE*oh)**expn

        ! Rosin-Rammler distribution

        nexpn = 3
        dbar  = dd
        dmin  = dd / TEN
        dmax  = dd

        k_temp = ONE - exp( (dmin/dbar)**nexpn - (dmax/dbar)**nexpn )

        t1 = (dmin / dbar)**nexpn - log(ONE - k_temp*rand)

        eptp(jdp, iparc) = dbar * t1**(ONE/nexpn)

        ipepa(jiliqc, iparc) = 0
      end if

      ! other attributes: such as number attribute

      mass = eptp(jmp,    iparc) * pepa(jrpoi, iparc)
      diam = eptp(jdp,    iparc)
      dens = pepa(jrdens, iparc)

      pepa(jrpoi, iparc) = SIX * mass / (PI * dens * diam**THREE)
      eptp(jmp,   iparc) = mass / pepa(jrpoi, iparc)

    end do

  contains
    subroutine SolveOmegaKs_( oKs, oOmgs )
    !===========================================================================
      !-------------------------------------------------------------------------
      ! Calculate Ks and Omegas in 3.2.2.2.3-14
      ! Ref 1: worklog 2021.04.19
      ! Ref 2: T. Sauer, Numerical Analysis, E3, p28
      !-------------------------------------------------------------------------
      real(RKD), intent(out) :: oKs, oOmgs

      real(RKD), parameter :: TOL = 1.e-4_RKD

      real(RKD) :: ksta, kend, xk
      !-------------------------------------------------------------------------

      ksta = ZERO
      kend = xq * spdr**TWO * densl / tens

      do while ( HALF*(kend-ksta) > TOL )
        xk = HALF * (ksta + kend)

        if( dfunc_(xk) < TOL ) exit

        if( dfunc_(ksta) * dfunc_(xk) < ZERO ) then
          kend = xk
        else
          ksta = xk
        end if
      end do

      oKs = xk
      oOmgs = func_(xk)
    end subroutine SolveOmegaKs_

    function dfunc_(x) result(y)
    !===========================================================================
      !-------------------------------------------------------------------------
      ! Abstract function of dOmegas: 3.2.2.2.3-15
      !-------------------------------------------------------------------------
      real(RKD), intent(in) :: x
      real(RKD) :: y

      real(RKD) :: c1, c2, c3, c4, c5, c6, c7, tt1, tt2, tt3
      !-------------------------------------------------------------------------

      c1 =  EIGHT * viscl**TWO
      c2 =  xq    * spdr**TWO
      c3 = -THREE * tens / (TWO * densl)
      c4 =  FOUR  * viscl**TWO
      c5 =  xq    * spdr**TWO
      c6 = -tens  / densl
      c7 = -FOUR  * viscl

      tt1 = c1 * x**THREE + c2 * x + c3 * x**TWO
      tt2 = sqrt( c4 * x**FOUR + c5 * x**TWO + c6 * x**THREE )
      tt3 = c7 * x

      y = tt1 / tt2 + tt3
    end function dfunc_

    function func_(x) result(y)
    !===========================================================================
      !-------------------------------------------------------------------------
      ! 3.2.2.2.3-14
      !-------------------------------------------------------------------------
      intrinsic :: sqrt

      real(RKD), intent(in) :: x
      real(RKD) :: y

      real(RKD) :: c1, c2, c3, c4, tt1, tt2, tt3, tt4
      !-------------------------------------------------------------------------

      c1 = -TWO  * viscl
      c2 =  FOUR * viscl**TWO
      c3 =  xq   * spdr**TWO
      c4 = -tens / densl

      tt1 = c1 * x**TWO
      tt2 = c2 * x**FOUR
      tt3 = c3 * x**TWO
      tt4 = c4 * x**THREE

      y = tt1 + sqrt( tt2 + tt3 + tt4 )

    end function func_
  end subroutine HFY_SLS_Atomize

  subroutine HFY_SLS_Fin()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! de-allocate dynamical arrays which only depend on number of nozzles
    !---------------------------------------------------------------------------

    integer :: inoz, iuse, nnozs

    type(NozzleSLS_t), pointer :: thisNoz => NULL()
    !---------------------------------------------------------------------------

    nnozs = numbNozz_SLS

    do inoz = 1, nnozs
      thisNoz => nozzles_SLS(inoz)

      do iuse = 1, 3
        deallocate( thisNoz%rands_(iuse)%rands_ )
      end do

      deallocate( thisNoz%rands_ )

      deallocate( thisNoz%parcsInj_ )
    end do

    nullify(thisNoz)

    deallocate( nozzles_SLS )

  end subroutine HFY_SLS_Fin

  !*****************************************************************************
  ! INTERNAL PRIVATE PROCEDURES
  !*****************************************************************************

  subroutine UpdateNoz_()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! Update
    ! 1, noz%numbInjParcs_
    ! 2, call Resize_()

    !---------------------------------------------------------------------------

    intrinsic :: cos, max, sqrt

    integer :: inoz, nnozs
    real(RKD) :: meff, dens, diam, thet, dpres, kv, dtLagr
    real(RKD) :: t1, t2, t3, maxl, maxr, worktime, xnumb

    type(NozzleSLS_t), pointer :: thisNoz
    !---------------------------------------------------------------------------

    nnozs = numbNozz_SLS

    ! number of parcels to be injected this time step

    do inoz = 1, nnozs
      thisNoz => nozzles_SLS(inoz)

      dtLagr = dtLagr_ING

      ! 3.2.2.1-9
      worktime = thisNoz%endTime_ - thisNoz%staTime_
      xnumb = thisNoz%numbTotParcs_ / worktime * dtLagr

      ! Guarantee that at least one parcel would be injected one time
      if( xnumb < ONE ) then
        thisNoz%numbInjParcs_ = 1
      else
        thisNoz%numbInjParcs_ = int( xnumb )
      end if
    end do

    ! Allocate or resize noz%parcsInj_. This array depends on the number of
    ! parcels injected this time

    call Resize_()

  end subroutine UpdateNoz_

  subroutine Resize_()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! Allocate or resize nozzles%parcelInj_ which depends on the number of
    ! injected parcels this time

    ! After the injected number of parcels known,
    ! 1, if not allocated, allocate it
    ! 2, if allocated, check the size of the arrays,
    !   2.1, if parc number > size(.%parcelsInj_), resize it
    !---------------------------------------------------------------------------
    intrinsic :: associated, size

    integer, parameter :: NOVERSIZE = 1000

    integer :: nsize, nparcs, nnozs, inoz
    logical :: yes

    type(NozzleSLS_t), pointer :: thisNoz => NULL()
    !---------------------------------------------------------------------------

    ! temporarily store the number of nozzles
    nnozs  = numbNozz_SLS

    do inoz = 1, nnozs
      thisNoz => nozzles_SLS(inoz)

      nparcs = thisNoz%numbInjParcs_

      yes = associated( thisNoz%parcsInj_ )

      if( .not. yes ) then
        allocate( thisNoz%parcsInj_(nparcs+NOVERSIZE) )
      else

        nsize = size( thisNoz%parcsInj_, dim=1 )

        if( nsize < nparcs ) then
          deallocate( thisNoz%parcsInj_ )

          allocate( thisNoz%parcsInj_(1:nparcs+NOVERSIZE) )
        end if
      end if
    end do

  end subroutine Resize_

  subroutine GenRands_()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! Allocate and generate random number sequences. They are attached to each
    ! nozzle
    !---------------------------------------------------------------------------

    intrinsic :: mod, associated

    integer, parameter :: NOVERSIZE = 1000000

    integer :: iuse, inoz, nuses, nnozs, nparcs!, ndigits, i, j
    !integer, dimension(:), allocatable :: seed

    type(NozzleSLS_t), pointer :: thisNoz => NULL()
    type(Rands_t),     pointer :: myRands => NULL()
    !---------------------------------------------------------------------------

    nnozs = numbNozz_SLS

    do inoz = 1, nnozs
      thisNoz => nozzles_SLS(inoz)

      nuses = thisNoz%numbRands_

      do iuse = 1, nuses
        myRands => thisNoz%rands_(iuse)

        !call random_seed( size=ndigits )

        !allocate( seed(1:ndigits) )

        ! the seed is computed by repeated the values in intrinsic subroutine
        !   'date_and_time(values=xxx)'
        !call date_and_time( values=seed(1:8) )

        !do i = 9, ndigits
        !  j = mod(i, 8); if( j==0 ) j = 8
        !  seed(i) = seed(j)
        !end do

        !call random_seed( put=seed )
        call random_seed()

        nparcs = thisNoz%numbTotParcs_

        allocate( myRands%rands_(1:nparcs+NOVERSIZE) )

        ! the length or size of %rands_(:), used to check if %irand_ beyond
        myRands%nsize_ = nparcs + NOVERSIZE

        call random_number( myRands%rands_(:) )

        ! initialize the index indicator
        myRands%irand_ = 0

        !deallocate( seed )
      end do
    end do

  end subroutine GenRands_

  subroutine CheckIn_()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! Check two parameters in input
    ! 1, noz%direc_
    ! 2, noz%numbTotParcs_  ! TODO:

    !---------------------------------------------------------------------------

    real(RKD), parameter :: EPS = 1.0e-2_RKD

    integer :: inoz, nnozs
    real(RKD) :: magn

    type(NozzleSLS_t), pointer :: thisNoz => NULL()
    !---------------------------------------------------------------------------

    nnozs = numbNozz_SLS

    do inoz = 1, nnozs
      thisNoz => nozzles_SLS(inoz)

      magn = CalcMagnVec_( thisNoz%dire_ )

      if( magn < EPS ) then
        ! FIXME: btf_print()
        print *, "WARNING: nozzle ", inoz, " n_ref 's Magnitude Is Too Small"
      end if

      thisNoz%dire_(:) = thisNoz%dire_(:) / magn
    end do

  end subroutine CheckIn_

  function CalcMagnVec_( aVec ) result( magn )
  !=============================================================================
    !---------------------------------------------------------------------------
    ! Compute the magnitude of a vector
    !---------------------------------------------------------------------------

    intrinsic :: sqrt

    real(RKD), dimension(1:3), intent(in) :: aVec
    real(RKD) :: magn

    integer :: i
    !---------------------------------------------------------------------------

    magn = ZERO

    do i = 1, 3
      magn = magn + aVec(i) * aVec(i)
    end do

    magn = sqrt( magn )

  end function CalcMagnVec_

  subroutine CalcVecProduct_( vec1, vec2, vec_product )
  !=============================================================================
    !---------------------------------------------------------------------------
    ! internal procedure

    ! vector product of \vec vec1 and \vec vec2, the restult is also a vector
    !   \vec vec_product
    !---------------------------------------------------------------------------

    real(RKD), dimension(1:3), intent(in) :: vec1, vec2
    real(RKD), dimension(1:3), intent(out) :: vec_product
    !---------------------------------------------------------------------------

    vec_product(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
    vec_product(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
    vec_product(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)

  end subroutine CalcVecProduct_

  subroutine CheckIrand_()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! Check whether the parcel's random number index indicator exeeds its bound
    !---------------------------------------------------------------------------

    integer :: nsize, inoz, nnoz, iptr, iuse
    type(NozzleSLS_t), pointer :: thisNoz => NULL()
    type(Rands_t),     pointer :: myRands => NULL()
    !---------------------------------------------------------------------------

    nnoz = numbNozz_SLS

    do inoz = 1, nnoz
      thisNoz => nozzles_SLS(inoz)

      do iuse = 1, 3
        myRands => thisNoz%rands_(iuse)

        iptr  = myRands%irand_
        nsize = myRands%nsize_

        if( iptr > nsize ) then
          ! FIXME: btf_print()
          print *, "Parcel Index of rands Beyond! Nozzle/iuse: ", inoz, iuse
        end if
      end do
    end do

    nullify(thisNoz)
    nullify(myRands)

  end subroutine CheckIrand_

  subroutine CalcIniFilmThickness_()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! Calculating thickness of the liquid film, as initial diam of sheet parcels
    !---------------------------------------------------------------------------
    intrinsic :: cos, sqrt

    real(RKD), parameter :: EPS = 1.e-18_RKD

    integer :: inoz, nnozs, nparcs
    real(RKD) :: meff, dens, diam, thet, t1, su

    type(NozzleSLS_t), pointer :: thisNoz => NULL()
    !---------------------------------------------------------------------------

    nnozs = numbNozz_SLS

    do inoz = 1, nnozs

      thisNoz => nozzles_SLS(inoz)

      meff = thisNoz%massRate_
      dens = thisNoz%dens_
      diam = thisNoz%diam_
      thet = thisNoz%hfCoAngl_

      ! 3.2.2.2.1-3
      ! Small u: inj vel along nozzle axis

      su = thisNoz%magnVel_ * cos( thet )

      ! 3.2.2.2.1-2

      t1 = diam * diam - FOUR * meff / (PI*dens*su)

      ! Fixed by shuning@simpop.cn
      ! Avoiding the minus 0, say, -5.12e-23
      t1 = t1 + EPS

      if( t1 < 0 ) then
        ! FIXME: btf_print()
        print *, "WARNING: 3.2.2.2.1-2, sqrt < 0 "
      end if

      nparcs = thisNoz%numbInjParcs_
      thisNoz%parcsInj_(1:nparcs)%diam_ = HALF * ( diam - sqrt(t1) )
    end do

    thisNoz => NULL()

  end subroutine CalcIniFilmThickness_

  subroutine CalcInjSpeed_()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! Calculate the injecting velocity of each nozzle:
    ! noz%magnVel_
    !---------------------------------------------------------------------------
    integer :: inoz, nnozs
    real(RKD) :: meff, dens, diam, thet, dpre, ul, ur

    type(NozzleSLS_t), pointer :: thisNoz => NULL()
    !---------------------------------------------------------------------------

    nnozs = numbNozz_SLS

    do inoz = 1, nnozs
      ! 1st: Magnitude of injecting velocity

      thisNoz => nozzles_SLS(inoz)

      meff = thisNoz%massRate_
      dens = thisNoz%dens_
      diam = thisNoz%diam_
      thet = thisNoz%hfCoAngl_

      dpre = thisNoz%presDiff_

      ! 3.2.2.2.1-4 & 3.2.2.2.1-5

      ul = 0.7_RKD * sqrt( TWO * dpre / dens )

      ur = FOUR * meff / (PI * diam**TWO * dens * cos(thet))

      thisNoz%magnVel_ = max( ul, ur )
    end do

  end subroutine CalcInjSpeed_

  ! ...
end module hfy_slisa_mod
