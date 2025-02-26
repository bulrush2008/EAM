
module hfy_wave_mod
  implicit none
  !-----------------------------------------------------------------------------
  ! author    date        aff
  ! snx       21.06.18    simpop.cn
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

    real(RKD) :: dens_ = ZERO
    !---------------------------------------------------------------------------

    ! These Attributes Below Would Be Added Into eptp/pepa/ipepa

    ! liquidCore = 1 by default
    integer :: iBlob_ = 1

    ! random number \in [0, 1], attached to parcel, needed by new velocity after
    ! breaking up
    real(RKD) :: rand_ = ZERO

    ! breaking length
    real(RKD) :: bLen_ = ZERO
  end type Parcel_t

  type, private :: Rands_t
    integer :: irand_ = 0 ! index of rands_
    integer :: nsize_ = 0 ! size of rands_

    real(RKD), dimension(:), pointer :: rands_ => NULL()
  end type Rands_t

  type, public :: NozzleWAV_t
    !---------------------------------------------------------------------------
    ! INPUT
    ! from user

    ! time for nozzle beginning and ending to work
    real(RKD) :: staTime_ = ZERO
    real(RKD) :: endTime_ = ONE

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

    ! position of nozzle face
    real(RKD), dimension(1:3) :: posi_ = [ZERO, ZERO, ZERO]
    ! orientation vector normal to nozzle face: direction
    real(RKD), dimension(1:3) :: dire_ = [ZERO, ZERO, ZERO]

    ! HalF-COne ANGLe
    real(RKD) :: hfCoAngl_ = ZERO

    ! the number of parcels scheduled to be injected
    ! 1, The final value may be marginally altered
    integer :: numbTotParcs_ = 0

    ! nozzle diameter, used to calculate Cd: 3.2.2.3-9
    real(RKD) :: diam_ = ZERO
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
    ! 3 random number series needed:
    ! 2 for velocity direction
    ! 1 for changing velocity after breaking
    integer :: numbRands_ = 3

  end type nozzleWAV_t
  !-----------------------------------------------------------------------------

  ! number of nozzles, defined in uslag1.f90 
  integer, public :: numbNozz_WAV = 0

  ! nozzles in WAVE model
  type(NozzleWAV_t), dimension(:), pointer :: nozzles_WAV => NULL()
  public :: nozzles_WAV

  ! Lagrangian time step
  real(RKD) :: dtLagr = ZERO

  ! density of liquid blob parcel
  real(RKD) :: densLiq = ZERO

  ! breaking length
  real(RKD), dimension(:), allocatable :: bLenNoz

  ! density field of gas field
  real(RKD), dimension(:), pointer :: cRomCS => NULL()
  !-----------------------------------------------------------------------------

  ! Subroutines are by default private

  ! These Subroutines Are Made Public
  public ::           &
    HFY_WAV_Ini,      &
    HFY_WAV_Inject,   &
    HFY_WAV_Atomize,  &
    HFY_WAV_Fin

  !private ::
  ! CalcBreakingLength_,
  ! CalcMagnVec_,
  ! CalcVecProduct_,
  ! CheckIrand_,
  ! CheckNoz_,
  ! GenRands_,
  ! Resize_,
  ! UpdateNoz_
contains
  subroutine HFY_WAV_Ini()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! 1, Allocate nozzles of type NozzleWAV_t: nozzles_WAV if ipass=1
    ! 2, allocate bLenNoz if ipass=1
    ! 3, read user-setting parameters and then call checkNoz_, if ipass = 1
    ! 4, generated nozzle-attached random sequences, if ipass = 1
    ! 5, Calculate breaking length each time step
    !   5.1 bLenNoz(:)
    ! 6, init some private but global variables, each time step
    !   6.1 dtLagr
    !   6.2 densLiq
    ! 7, update the density field
    !   7.1 cRomCS(:)
    !---------------------------------------------------------------------------
    use lagran, only: dtp
    use field,  only: field_get_val_s
    use numvar, only: icrom
    !---------------------------------------------------------------------------
    intrinsic :: allocated, associated

    integer :: nnozs
    integer, save :: ipass = 0
    !---------------------------------------------------------------------------

    ipass = ipass + 1
    !---------------------------------------------------------------------------

    ! numbNozz_WAV has been initialized in uslag1.f90

    nnozs = numbNozz_WAV

    ! request memory for nozzles

    if( .not. associated(nozzles_WAV) ) then 
      allocate( nozzles_WAV(1:nnozs) )
    end if
    !---------------------------------------------------------------------------
    ! read user-setting parameters

    if( ipass == 1 ) then

      call hfy_wave_in()

      ! nomalize the nozzle direction vector, if need
      call CheckNoz_()
    end if
    !---------------------------------------------------------------------------

    ! Generate a data pool of random number, for sLISA model

    if( ipass == 1) call GenRands_()
    !---------------------------------------------------------------------------

    ! keeping the density field of gas phase this time step
    ! NOTE: cRomCS is global arrays in module

    !cRomCS => NULL()
    call field_get_val_s( icrom, cRomCS )
    !---------------------------------------------------------------------------

    ! for breaking length array
    if( .not. allocated(bLenNoz) ) then
      allocate( bLenNoz(1:nnozs) )
    end if

    ! breaking length may be different bwtween nozzles
    call CalcBreakingLength_()

    !---------------------------------------------------------------------------
    ! the common variables alive in WAVE model
    ! NOTE: dtLagr is global variables in module

    dtLagr = dtp

    ! density of liquid blob. It is temporarily assumed that density of injected
    ! blob parcels for all nozzles is same.

    densLiq = nozzles_WAV(1)%dens_

  end subroutine HFY_WAV_Ini

  subroutine HFY_WAV_Inject()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! Inject blob parcels
    !---------------------------------------------------------------------------

    intrinsic :: cos, sin

    ! contraction coeffcient, 3.2.2.3-7
    real(RKD), parameter :: CA  = 0.5_RKD
    real(RKD), parameter :: EPS = 1.E-2_RKD

    integer :: inoz, nnozs, iparc, nparcs, iptr
    real(RKD) :: masr, mass, diam, dens, modl, v1, v2, v3, thet, alph
    real(RKD) :: rand
    real(RKD), dimension(1:3) :: t1Vec, t2Vec, nRefVec, nRadVec

    type(NozzleWAV_t), pointer :: thisNoz => NULL()
    type(Parcel_t),    pointer :: theParc => NULL()
    type(Rands_t),     pointer :: myRands => NULL()
    !---------------------------------------------------------------------------

    ! Update nozzles,
    ! 1, injecting speed
    ! 2, number blob parcels to be injected
    ! 3, resize noz%parcsInj_(:)

    call UpdateNoz_()
    !---------------------------------------------------------------------------

    ! for attributes: mass, diameter, number, and also rand, iblob & bLen

    nnozs = numbNozz_WAV

    do inoz = 1, nnozs
      thisNoz => nozzles_WAV(inoz)

      nparcs = thisNoz%numbInjParcs_

      do iparc = 1, nparcs
        theParc => thisNoz%parcsInj_(iparc)
        !-----------------------------------------------------------------------

        ! mass attribute

        masr = thisNoz%massRate_

        theParc%mass_ = masr * dtLagr / nparcs  ! 3.2.2.3-6
        !-----------------------------------------------------------------------

        ! diameter attribute, 3.2.2.3-7

        theParc%diam_ = sqrt( CA ) * thisNoz%diam_
        !-----------------------------------------------------------------------

        ! parcel's density

        theParc%dens_ = thisNoz%dens_
        !-----------------------------------------------------------------------

        ! number attribute

        mass = theParc%mass_
        dens = theParc%dens_
        diam = theParc%diam_

        ! FIXME: avoiding big error when diameter is two small !
        theParc%numb_ = SIX * mass / (dens * PI * diam * diam * diam)
        !-----------------------------------------------------------------------

        ! attached random number, to disturb the velocity after breaking

        myRands => thisNoz%rands_(3)
        myRands%irand_ = myRands%irand_ + 1

        iptr = myRands%irand_

        theParc%rand_ = myRands%rands_(iptr)
        !-----------------------------------------------------------------------

        ! denoting the blob parcel

        theParc%iBlob_ = 1
        !-----------------------------------------------------------------------

        ! breaking length, bLenNoz(1:nnoz) was computed when initializing

        theParc%bLen_ = bLenNoz(inoz)
        !-----------------------------------------------------------------------

        theParc => NULL()
        myRands => NULL()
      end do

      thisNoz => NULL()
    end do

    call CheckIrand_()
    !---------------------------------------------------------------------------

    ! position of parcels

    nnozs = numbNozz_WAV

    do inoz = 1, nnozs
      thisNoz => nozzles_WAV(inoz)

      nparcs = thisNoz%numbInjParcs_

      do iparc = 1, nparcs
        theParc => thisNoz%parcsInj_(iparc)

        theParc%pos_(1:3) = thisNoz%posi_(1:3)

        theParc => NULL()
      end do

      thisNoz => NULL()
    end do
    !---------------------------------------------------------------------------

    ! velocity

    nnozs = numbNozz_WAV

    do inoz = 1, nnozs
      thisNoz => nozzles_WAV(inoz)

      ! temperarily store the nozzle face vector: \vec n_ref, see fig 3.6
      nRefVec(1:3) = thisNoz%dire_(1:3)

      modl = CalcMagnVec_( nRefVec - [ONE, ZERO, ZERO] )

      if( modl < EPS ) then
        t1Vec = [ZERO, ONE, ZERO]
      else
        v1 = nRefVec(1)
        v2 = nRefVec(2)
        v3 = nRefVec(3)

        ! comput \vec t1: see eqn 3.2.2.1-1
        t1Vec = [ONE-v1*v1, -v1*v2, -v1*v3]

        ! normalization
        t1Vec(1:3) = t1Vec(1:3) / CalcMagnVec_( t1Vec )
      end if

      ! \vec t2 = \vec n_ref \times \vec t1, see eqn 3.2.2.1-4
      call CalcVecProduct_( nRefVec, t1Vec, t2Vec )

      nparcs = thisNoz%numbInjParcs_

      do iparc = 1, nparcs
        theParc => thisNoz%parcsInj_(iparc)

        ! radial angle

        myRands => thisNoz%rands_(1)
        myRands%irand_ = myRands%irand_ + 1
        iptr = myRands%irand_

        rand = myRands%rands_(iptr)
        thet = thisNoz%hfCoAngl_ * rand

        ! circumferential angle

        myRands => thisNoz%rands_(2)
        myRands%irand_ = myRands%irand_ + 1
        iptr = myRands%irand_

        rand = myRands%rands_(iptr)
        alph = TWO * PI * rand

        ! vector along radial direction, 3.2.2.3-12
        nRadVec = cos( alph ) * t1Vec + sin( alph ) * t2Vec

        ! final velocity vector

        theParc%vel_(:)                                     &
        = thisNoz%magnVel_                                  &
        * ( cos(thet)*nRefVec(1:3) + sin(thet)*nRadVec(1:3) )

        myRands => NULL()
        theParc => NULL()
      end do

      thisNoz => NULL()
    end do

    call CheckIrand_()

  end subroutine HFY_WAV_Inject

  subroutine HFY_WAV_Atomize()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! Check each parcel in computing domain,
    ! - if iBlob == 1, i.e. it is a blob parcel, operate it
    ! - if iBlob == 0, i.e. it is a real parcel, next one

    ! NOTE:
    ! Logically speaking, each blob parcel suffers either no breaking or at most
    ! breaking once
    ! 1, if no breaking within breaking length, the parcel would change directly
    ! to a real parcel by changing iBlob = 0
    ! 2, if a blob parcel breaks in breaking length, the two son parcels would
    ! both be real ones
    !---------------------------------------------------------------------------

    use lagran, only: nbpart, dnbpar, nbptot, lagr_resize_particle_set,         &
                      ipepa,  jiblob, jisor,                                    &
                      pepa,   jrxp0,  jryp0, jrzp0, jrblen, jrtens, jrvisc,     &
                              jrfrb,  jrdens,jrrand,jrtsp,                      &
                      eptp,   jxp,    jyp,   jzp,   jdp,    jup,    jvp,    jwp,&
                              juf,    jvf,   jwf,   jmp,    jrpoi,              &
                      eptpa
    use parall, only: parcpt, irangp
    !---------------------------------------------------------------------------

    intrinsic :: cos, exp, min, sin, sqrt

    real(RKD), parameter :: B0   = 0.61_RKD ! in 3.2.2.3-18
    real(RKD), parameter :: B1   = 1.73_RKD ! in 3.2.2.3-19
    real(RKD), parameter :: CC1  = ZERO     ! in 3.2.2.3-33, disturb vel
    real(RKD), parameter :: WECR = 12.0_RKD ! same with KH submodel in KHRT
    real(RKD), parameter :: YCR  = 0.03_RKD ! critical detached mass, or 0.05

    integer :: iparc, nParcsIn, iblob, iel, nParcsNew, nip, nSize, nbp_gen_sum
    real(RKD) :: diamp, rand, angl, vvert, vpara, massb, masss, ratio
    real(RKD) :: dist, blen, densg, tensp, weg, spdr, diams, diamb
    real(RKD) :: wel, densl, rel, radip, viscl, oh, ta, lambda_kh, omega_kh
    real(RKD) :: c1, c2, c3, c4, c5, c6, c7, blbd, radis, minl, minr, tau_kh
    real(RKD) :: radib, massp, dmass, upp, vpp, wpp, spdp, spdpsq
    real(RKD) :: dnbp
    real(RKD), dimension(1:3) :: pos, pos0, velg, vell, l1v, l2v, lvert
    real(RKD), dimension(1:3) :: vels, velb
    !---------------------------------------------------------------------------

    ! number of newly-generated parcels this time step
    nParcsNew = 0

    ! number of current parcels in domain
    nParcsIn = nbpart

    ! check particle count limit, assuming that all parcels would break
    nsize = nParcsIn * 2
    if( lagr_resize_particle_set(nsize) < 0 ) then
      ! btfprint instead
      print *, "WARNING:"
      print *, "Parcel Set Reachs the Number Limit, in WAVE Model."
      return
    end if

    do iparc = 1, nParcsIn
      !-------------------------------------------------------------------------
      ! temporily keeping the variables, used in this subroutine

      ! vector attributes

      ! gas velocity there
      velg(1) = eptp(juf, iparc)
      velg(2) = eptp(jvf, iparc)
      velg(3) = eptp(jwf, iparc)

      ! parcel velocity
      vell(1) = eptp(jup, iparc)
      vell(2) = eptp(jvp, iparc)
      vell(3) = eptp(jwp, iparc)

      ! position when injecting
      pos0(1) = pepa(jrxp0, iparc)
      pos0(2) = pepa(jryp0, iparc)
      pos0(3) = pepa(jrzp0, iparc)

      ! current position
      pos(1) = eptp(jxp, iparc)
      pos(2) = eptp(jyp, iparc)
      pos(3) = eptp(jzp, iparc)

      ! breaking length attached to each blob parcel

      blen = pepa(jrblen, iparc)

      ! density of gas
      iel = ipepa(jisor, iparc)
      densg = cRomCS(iel)

      ! diameter & radius
      diamp = eptp(jdp, iparc)
      radip = eptp(jdp, iparc) * HALF

      ! mass
      massp = eptp(jmp, iparc) * pepa(jrpoi, iparc)

      ! surface tension
      tensp = pepa(jrtens, iparc)

      ! density of liquid blob or parcel
      densl = pepa(jrdens, iparc)

      ! viscosity of liquid blob or parcel
      viscl = pepa(jrvisc, iparc)
      !-------------------------------------------------------------------------

      ! *** 1st: All parcels exceeding breaking length, must be real parcels ***

      ! 1st judgement from breaking length

      dist = CalcMagnVec_( pos - pos0 )

      ! all parcels beyond the Breaking Length must be real parcels

      if( dist > bLen ) then
        ipepa(jiblob, iparc) = 0

        cycle
      end if
      !-------------------------------------------------------------------------

      ! *** 2nd: further, parcels with iBlob = 0 are all real parcels ***

      iblob = ipepa(jiblob, iparc)

      if( iblob == 0) then
        cycle
      end if
      !-------------------------------------------------------------------------

      ! *** 3rd: all other parcels are blobs. ***

      ! 3.1 Gas Weber: weg

      ! relative speed

      spdr = CalcMagnVec_( vell - velg )

      ! gas Weber number, 3.2.2.3-22

      weg = densg * spdr**TWO * radip / tensp

      ! do nothing if gas weber number is low

      if( weg < WECR ) cycle

      ! cumulate the mass for KH-mode break up

      ! Weber number of liquid blob parcel

      ! 3.2.2.3-23, liquid Weber number
      wel = densl * spdr**TWO * radip / tensp

      ! 3.2.2.3-25, liquid Reynolds number
      rel = densl * spdr * radip / viscl

      ! Ohnesorge number, 3.2.2.3-24
      oh = sqrt( wel ) / rel

      ! Taylor number, 3.2.2.3-26
      ta = oh * sqrt( weg )

      ! length of wave of least stability, in 3.2.2.3-20

      c1 = 9.02_RKD
      c2 = 0.45_RKD
      c3 = 0.4_RKD
      c4 = 0.7_RKD
      c5 = 0.865_RKD
      c6 = 1.67_RKD
      c7 = 0.6_RKD

      ! 3.2.2.3-20
      lambda_kh                   &
      = c1 * radip                &
      * (ONE + c2 * oh **HALF)    &
      * (ONE + c3 * ta **c4  )    &
      / (ONE + c5 * weg**c6  )**c7

      ! rate of increasing of wave of least stability, in 3.2.2.3-21

      c1 = 0.34_RKD
      c2 = 0.385_RKD
      c3 = 1.5_RKD
      c4 = 1.4_RKD
      c5 = 0.6_RKD

      ! 3.2.2.3-21
      omega_kh                              &
      = (c1  + c2 * weg**c3)                &
      / (ONE + oh)                          &
      / (ONE + c4 * ta **c5)                &
      * sqrt( tensp / densl / radip**THREE )

      ! radius of small son parcel, whether the blob parcel breaks this time

      ! blbd: B0 * Lambda_KH, in 3.2.2.3-18

      blbd = B0 * lambda_kh

      if( blbd < radip ) then
        radis = blbd
        diams = radis * TWO
      else
        minl = (THREE * PI * radip**TWO * spdr / omega_kh * HALF)**THIRD
        minr = (THREE * radip**TWO * lambda_kh * QUART)**THIRD

        radis = min( minl, minr )
        diams = radis * TWO
      end if

      ! characteristic time for break up

      c1     = 3.726_RKD
      tau_kh = c1 * B1 * radip / (lambda_kh * omega_kh) ! 3.2.2.3-19

      ! update the flag of diameter of big son parcel

      ! fdb: flags of diameter of big son parcel
      radib = pepa(jrfrb, iparc) 

      radib = radis + (radib - radis) * exp( -dtLagr / tau_kh )
      diamb = radib * TWO

      ! cumulate the detached mass for small son parcel

      ratio = ONE - (diamb / diamp)**THREE
      dmass = massp * ratio

      if( ratio < YCR ) then
        ! update jrfrb
        pepa(jrfrb, iparc) = radib

        cycle
      end if

      ! perform KH break mode

      ! generate a new parcel

      nParcsNew = nParcsNew + 1

      ! the index the the new son parcel
      nip = nParcsIn + nParcsNew

      eptp(:, nip) = eptp(:, iparc)
      pepa(:, nip) = pepa(:, iparc)

      ipepa(:, nip) = ipepa(:, iparc)

      eptpa(:, nip) = eptpa(:, iparc)

      ! define properties of the two son parcels

      ! diameter
      eptp(jdp, iparc) = diamb
      eptp(jdp, nip  ) = diams

      ! number attributes

      massb = massp - dmass
      masss = dmass

      pepa(jrpoi, iparc) = SIX * massb / ( densl * PI * diamb**THREE)
      pepa(jrpoi, nip  ) = SIX * masss / ( densl * PI * diams**THREE)

      ! mass

      eptp(jmp, iparc) = massb / pepa(jrpoi, iparc)
      eptp(jmp, nip  ) = masss / pepa(jrpoi, nip  )

      ! position
      !eptp(jxp, iparc) = pos(1)
      !eptp(jyp, iparc) = pos(2)
      !eptp(jzp, iparc) = pos(3)

      eptp(jxp, nip  ) = pos(1)
      eptp(jyp, nip  ) = pos(2)
      eptp(jzp, nip  ) = pos(3)

      ! speed of the parcel

      spdp   = CalcMagnVec_( vell )
      spdpsq = spdp * spdp

      upp = vell(1)
      vpp = vell(2)
      wpp = vell(3)

      ! 3.2.2.3-34
      l1v = [ONE - upp * upp / spdpsq,  &
                 - upp * vpp / spdpsq,  &
                 - upp * wpp / spdpsq]

      l1v = l1v / CalcMagnVec_( l1v )

      ! 3.2.2.3-35
      call CalcVecProduct_( vell, l1v, l2v )

      rand = pepa(jrrand, iparc)
      angl = TWO * PI * rand

      ! 3.2.2.3-36
      lvert = cos( angl ) * l1v + sin( angl ) * l2v

      ! 3.2.2.3-33
      vvert = CC1 * lambda_kh * omega_kh

      ! 3.2.2.3-37
      vpara = sqrt( ONE - (vvert/spdp)**TWO )

      ! for small son parcel, 3.2.2.3-38
      vels(1:3) = vpara * vell(1:3) + vvert * lvert(1:3)

      ! for small son parcel

      eptp(jup, nip) = vels(1)
      eptp(jvp, nip) = vels(2)
      eptp(jwp, nip) = vels(3)

      ! for big son parcel

      massb = eptp(jmp, iparc)
      masss = eptp(jmp, nip  )

      ! momentum conservation, 3.2.2.3-39
      velb(1:3) = massp * vell(1:3) - masss * vels(1:3)
      velb(1:3) = velb(1:3) / massb

      eptp(jup, iparc) = velb(1)
      eptp(jvp, iparc) = velb(2)
      eptp(jwp, iparc) = velb(3)
      !-------------------------------------------------------------------------

      ! other attributes of the two son parcels

      ! the two son parcels are now real rather blob parcels

      ipepa(jiblob, iparc) = 0
      ipepa(jiblob, nip  ) = 0

      ! jrfrb for two new son parcels
      pepa(jrfrb, iparc) = eptp(jdp, iparc) * HALF
      pepa(jrfrb, nip  ) = eptp(jdp, nip  ) * HALF

      ! tension
      pepa(jrtens, iparc) = tensp
      pepa(jrtens, nip  ) = tensp

      ! viscosity
      pepa(jrvisc, iparc) = viscl
      pepa(jrvisc, nip  ) = viscl

      ! density
      pepa(jrdens, iparc) = densl
      pepa(jrdens, nip  ) = densl

      ! life
      pepa(jrtsp, nip) = ZERO

    end do

    ! update some other variables: nbpart, dnbpar, nbpartall

    ! nbpart
    nbpart = nbpart + nParcsNew

    ! update nbptot
    nbp_gen_sum = nParcsNew

    if( irangp .ge. 0 ) then
      ! MPI_Reduce and then broadcast if parallel computing
      call parcpt( nbp_gen_sum )
    end if

    nbptot = nbptot + nbp_gen_sum

    ! update dnbpar

    dnbp = ZERO

    ! loop on all parcels, including newly generated ones

    do iparc = 1, nbpart 
      dnbp = dnbp + pepa(jrpoi, iparc)
    end do

    dnbpar = dnbp

  end subroutine HFY_WAV_Atomize

  subroutine HFY_WAV_Fin()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! de-allocate dynamical arrays which only depend on number of nozzles
    !---------------------------------------------------------------------------

    integer :: inoz, iuse, nnozs
    type(NozzleWAV_t), pointer :: thisNoz => NULL()
    !---------------------------------------------------------------------------

    nnozs = numbNozz_WAV

    do inoz = 1, nnozs
      thisNoz => nozzles_WAV(inoz)

      do iuse = 1, 3
        deallocate( thisNoz%rands_(iuse)%rands_ )
      end do

      deallocate( thisNoz%rands_ )

      deallocate( thisNoz%parcsInj_ )
    end do

    deallocate( nozzles_WAV )

    deallocate( bLenNoz )

  end subroutine HFY_WAV_Fin

  !*****************************************************************************
  ! INTERNAL PRIVATE PROCEDURES
  !*****************************************************************************

  subroutine CalcBreakingLength_()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! Calculating breaking length: 3.2.2.3-8 in Proposal

    ! out: bLenNoz(:)
    !---------------------------------------------------------------------------

    use mesh, only: ncelet, ncel, xyzcen

    intrinsic :: sqrt

    ! Levich Constant, in 3.2.2.3-8, 7 to 16, or even 30, see RC Liu's thesis
    real(RKD), parameter :: CL = 6.0_RKD  

    integer :: inoz, nnozs, iel, iproc
    real(RKD) :: diam, densl, densg, xnoz, ynoz, znoz

    type(NozzleWAV_t), pointer :: thisNoz => NULL()

    interface    
      subroutine findpt( ncelet, ncel, xyzcen, xx, yy, zz, node, ndrang )
        integer :: ncelet, ncel, node, ndrang
        double precision :: xyzcen(3, ncelet)
        double precision :: xx, yy, zz
      end subroutine findpt
    end interface
    !---------------------------------------------------------------------------

    nnozs = numbNozz_WAV

    do inoz = 1, nnozs
      thisNoz => nozzles_WAV(inoz)

      diam  = thisNoz%diam_
      densl = thisNoz%dens_

      xnoz = thisNoz%posi_(1)
      ynoz = thisNoz%posi_(2)
      znoz = thisNoz%posi_(3)

      iel = 0
      iproc = 0
      call findpt(ncelet, ncel, xyzcen, xnoz, ynoz, znoz, iel, iproc )

      densg = cRomCS(iel)
! *** DEBUG STA ***
!print *, densl, densg
! *** DEBUG END ***
      bLenNoz(inoz) = CL * diam * sqrt( densl / densg )

      thisNoz => NULL()
    end do

  end subroutine CalcBreakingLength_

  subroutine CheckNoz_()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! 1, If the vec's magnitude /= unit, normalize the nozzle face vector
    ! 2, If the vec's magnitude \approx 0, print warning message
    !---------------------------------------------------------------------------

    real(RKD), parameter :: EPS = 1.0e-2_RKD

    integer :: inoz, nnozs
    real(RKD) :: magn
    type(NozzleWAV_t), pointer :: thisNoz => NULL()
    !---------------------------------------------------------------------------

    nnozs = numbNozz_WAV

    do inoz = 1, nnozs
      thisNoz => nozzles_WAV(inoz)

      magn = CalcMagnVec_( thisNoz%dire_ )

      if( magn < EPS ) then
        print *, "WARNING: nozzle ", inoz, " n_ref 's Magnitude Is Too Small"
      end if

      thisNoz%dire_(:) = thisNoz%dire_(:) / magn
    end do

  end subroutine CheckNoz_

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

  subroutine GenRands_()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! Allocate and generate random number sequences. They are attached to each
    ! nozzle
    !---------------------------------------------------------------------------

    intrinsic :: mod, associated

    integer, parameter :: NOVERSIZE = 10000

    integer :: iuse, inoz, nuses, nnozs, nparcs!, ndigits, i, j
    !integer, dimension(:), allocatable :: seed

    type(NozzleWAV_t), pointer :: thisNoz => NULL()
    type(Rands_t),     pointer :: myRands => NULL()
    !---------------------------------------------------------------------------

    nnozs = numbNozz_WAV

    do inoz = 1, nnozs
      thisNoz => nozzles_WAV(inoz)

      thisNoz%numbRands_ = 3

      nuses = thisNoz%numbRands_

      if( .not. associated(thisNoz%rands_) ) then
        allocate( thisNoz%rands_(1:nuses) )
      end if

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

  subroutine UpdateNoz_()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! Update nozzles
    ! 1, .%magnVel_
    ! 2, .%numbParcsInj_
    ! 3, call Resize_()
    !---------------------------------------------------------------------------

    intrinsic :: int, sqrt

    integer :: inoz, nnozs
    real(RKD) :: Cd, masr, diam, dens, dprs, area, wtim, xnmb

    type(NozzleWAV_t), pointer :: thisNoz => NULL()
    !---------------------------------------------------------------------------

    nnozs = numbNozz_WAV

    do inoz = 1, nnozs
      thisNoz => nozzles_WAV(inoz)

      !-------------------------------------------------------------------------
      ! the injecting speed

      masr = thisNoz%massRate_
      diam = thisNoz%diam_
      dens = thisNoz%dens_
      dprs = thisNoz%presDiff_

      ! area of this nozzle
      area = PI * diam * diam * QUART

      ! 3.2.2.3-9
      Cd = masr / (area * sqrt( TWO * dens * dprs ))

      ! 3.2.2.3-10
      thisNoz%magnVel_ = Cd * sqrt( TWO * dprs / dens )
      !-------------------------------------------------------------------------
! *** DEBUG STA ***
print *, 'Inj Speed, Cd, dprs, dens: ', thisNoz%magnVel_, Cd, dprs, dens
! *** DEBUG END ***
      ! the number of blob parcels to be injected this time step

      wtim = thisNoz%endTime_ - thisNoz%staTime_

      xnmb = thisNoz%numbTotParcs_ / wtim * dtLagr

      ! Guarantee that at least one parcel would be injected one time
      if( xnmb < ONE ) then
        thisNoz%numbInjParcs_ = 1
      else
        thisNoz%numbInjParcs_ = int( xnmb )
      end if

      nullify(thisNoz)
    end do

    ! Allocate or resize noz%parcsInj_. This array depends on the number of
    ! parcels injected this time

    call Resize_()

  end subroutine UpdateNoz_

  subroutine Resize_()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! Allocate or resize nozzles%parcelInj_(:) which depends on the number of
    ! injected parcels this time

    ! After the injected number of parcels known,
    ! 1, if not allocated, allocate it
    ! 2, if allocated, check the size of the arrays,
    !   2.1, if parc number > size(noz%parcelsInj_), resize it
    !---------------------------------------------------------------------------
    intrinsic :: associated, size

    integer, parameter :: NOVERSIZE = 1000

    integer :: nsize, nparcs, nnozs, inoz
    logical :: yes

    type(NozzleWAV_t), pointer :: thisNoz => NULL()
    !---------------------------------------------------------------------------

    ! temporarily store the number of nozzles

    nnozs  = numbNozz_WAV

    do inoz = 1, nnozs
      thisNoz => nozzles_WAV(inoz)

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

      thisNoz => NULL()
    end do

  end subroutine Resize_

  subroutine CheckIrand_()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! Check whether the parcel's random number index indicator exeeds its bound
    !---------------------------------------------------------------------------

    integer :: nsize, inoz, nnoz, iptr, iuse
    type(NozzleWAV_t), pointer :: thisNoz => NULL()
    type(Rands_t),     pointer :: myRands => NULL()
    !---------------------------------------------------------------------------

    nnoz = numbNozz_WAV

    do inoz = 1, nnoz
      thisNoz => nozzles_WAV(inoz)

      do iuse = 1, 3
        myRands => thisNoz%rands_(iuse)

        iptr  = myRands%irand_
        nsize = myRands%nsize_

        if( iptr > nsize ) then
          ! FIXME: btf_print()
          print *, "Parcel Index of rands Beyond! Nozzle/iuse: ", inoz, iuse
        end if

        myRands => NULL()
      end do

      thisNoz => NULL()
    end do

  end subroutine CheckIrand_

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

end module hfy_wave_mod
