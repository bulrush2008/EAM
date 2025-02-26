
module hfy_srr_mod
  implicit none
  !-----------------------------------------------------------------------------
  ! author    date        aff
  ! snx       21.05.25    simpop.cn
  !-----------------------------------------------------------------------------
  private
  save
  intrinsic :: selected_real_kind, acos
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
  real(RKD), parameter :: THIRD = 1.0_RKD/3.0_RKD
  real(RKD), parameter :: HALF  = 0.5_RKD

  real(RKD), parameter :: PI = acos(-1.0_RKD)
  !-----------------------------------------------------------------------------

  ! Abstract Data Type definition

  type, private :: parcel_t
    !private
    real(RKD), dimension(1:3) :: pos_ = [ZERO, ZERO, ZERO]
    real(RKD), dimension(1:3) :: vel_ = [ZERO, ZERO, ZERO]

    real(RKD) :: mass_ = ZERO
    real(RKD) :: diam_ = ZERO
    real(RKD) :: numb_ = ZERO
  end type parcel_t

  type, private :: rands_t
    integer :: irand_ = 0 ! index of rands_
    integer :: nSize_ = 0 ! size of rands_

    real(RKD), dimension(:), pointer :: rands_ => NULL()
  end type rands_t

  type, public :: nozzle_t
    !private
    ! for INPUT:

    ! time for nozzle beginning and ending to work
    real(RKD) :: staTime_ = ZERO
    real(RKD) :: endTime_ = ONE

    ! density of liquid droplets
    real(RKD) :: dens_ = ZERO

    ! time rate of mass of injecting
    real(RKD) :: massRate_ = ZERO

    ! position of nozzle face: Cartesian coordinates
    real(RKD), dimension(1:3) :: posi_ = [ZERO, ZERO, ZERO]
    ! vector normal to nozzle face, better nomalized
    real(RKD), dimension(1:3) :: dire_ = [ZERO, ZERO, ZERO]

    ! HalF-COne ANGLe
    real(RKD) :: hfCoAngl_ = ZERO
    ! dispersion radius, for each nozzle
    real(RKD) :: dispAngl_ = ZERO
    ! radius of nozzle
    real(RKD) :: dispRadi_ = ZERO

    !---------------------------------------------------------------------------
    ! RR model: exponent, mean/min/max diameters
    integer :: expnRR_ = 0

    real(RKD) :: meanDiamRR_ = ZERO
    real(RKD) :: miniDiamRR_ = ZERO
    real(RKD) :: maxiDiamRR_ = ZERO
    !---------------------------------------------------------------------------

    ! swirling fraction
    real(RKD) :: fracSwl_ = ZERO
    ! magnitude of injecting velocity
    real(RKD) :: magnVel_ = ZERO

    ! the number of parcels scheduled to be injected
    ! 1, The final value may be marginally altered
    integer :: numbTotParcs_ = 0
    !---------------------------------------------------------------------------
    ! for OUTPUT:

    ! the parcels injected each time step
    integer :: numbInjParcs_ = 0
    ! properties of parcel
    ! rank 1: diameter, number, mass, position and velocity
    ! rank 1: each parcel
    ! rank 2: each nozzle
    type(parcel_t), dimension(:), pointer :: parcsInj_ => NULL()

    ! nozzle-attached the random number series
    type(rands_t), dimension(:), pointer :: rands_ => NULL()
  end type nozzle_t
  !-----------------------------------------------------------------------------

  ! number of nozzles 
  integer, public :: numb_nozz_srr = 0

  ! the nozzles in RR model
  type(nozzle_t), dimension(:), pointer :: nozzles_srr => NULL()
  public :: nozzles_srr
  !-----------------------------------------------------------------------------

  ! subroutines are by default private

  ! public subroutines
  public ::         &
    hfy_srr_ini,    &
    hfy_srr_fin,    &
    hfy_srr_inject 

  !private ::
  ! resize_,
  ! gen_rands_,
  ! check_nozzle_,
  ! chekc_iptr_rands_,
  ! vec_magn_,
  ! sign_,
  ! vec_product_,
contains
  subroutine hfy_srr_ini()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! allocate nozzles_srr 
    !---------------------------------------------------------------------------
    intrinsic :: allocated

    integer :: nnozs

    integer, save :: ipass = 0
    !---------------------------------------------------------------------------

    ipass = ipass + 1

    ! 初始化，仅调用一次
    if( ipass >= 2 ) return
    !---------------------------------------------------------------------------

    ! the dim of the arrays below only dependent on nozzle number

    nnozs = numb_nozz_srr

    if( .not. associated(nozzles_srr) )  &
      allocate( nozzles_srr(1:nnozs) )

    ! read in parameters set by user
    call hfy_srr_in()

    ! check and nomalize \vec n_ref, i.e. dire_nozz_srr if needed
    call checkNoz_()
    !---------------------------------------------------------------------------

    ! the random number pool, size of the pool
    !TODO: 在需要的时候生成，十最好的方式。但如何解决这样的生成，是否仍满足需要的分布，
    ! 仍未解决
    call GenRands_()

  end subroutine hfy_srr_ini

  subroutine hfy_srr_fin()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! de-allocate dynamical arrays which only depend on number of nozzles
    !---------------------------------------------------------------------------

    integer :: inoz, iuse, nnozs
    type(nozzle_t), pointer :: thisNoz => NULL()
    !---------------------------------------------------------------------------

    nnozs = numb_nozz_srr

    do inoz = 1, nnozs
      thisNoz => nozzles_srr(inoz)

      do iuse = 1, 3
        deallocate( thisNoz%rands_(iuse)%rands_ )
      end do

      deallocate( thisNoz%rands_ )

      deallocate( thisNoz%parcsInj_ )
    end do

    deallocate( nozzles_srr )

  end subroutine hfy_srr_fin

  subroutine hfy_srr_inject()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! inject parcels with Rosin-Rammler distribution
    !---------------------------------------------------------------------------
    use lagran, only: dtp
    !---------------------------------------------------------------------------

    intrinsic :: int, exp, log, abs, tan, atan, sin, cos, sqrt

    real(RKD), parameter :: EPS_RR = 1.0e-4_RKD

    integer :: inoz, nparcs, nnozs, iptr, iparc
    real(RKD) :: worktime, xnumb, tdmin, tdmax, tdmean, tq, dens, k_temp, rand
    real(RKD) :: v1, v2, v3, r1, r2, alpha_rr, r_rr, theta_rr, v_swirl, v_cone
    real(RKD) :: pi6, diam, frac, dt_lagr
    real(RKD), dimension(1:3) :: t1_vec, t2_vec, nr_vec, n_cone_vec, nt_vec
    type(nozzle_t), pointer :: thisNoz  => NULL()
    type(parcel_t), pointer :: theParc  => NULL()
    type(rands_t),  pointer :: rands => NULL() 
    !---------------------------------------------------------------------------
    ! number of parcels to be injected this time

    nnozs = numb_nozz_srr

    do inoz = 1, nnozs

      thisNoz => nozzles_srr(inoz)

      dt_lagr = dtp

      ! 3.2.2.1-9
      worktime = thisNoz%endTime_ - thisNoz%staTime_
      xnumb = thisNoz%numbTotParcs_ / worktime * dt_lagr

      if( xnumb < ONE ) then
        thisNoz%numbInjParcs_ = 1
      else
        thisNoz%numbInjParcs_ = int( xnumb )
      end if
    end do

    ! allocate or resize noz%parcsInj_
    call resize_()
    !---------------------------------------------------------------------------

    ! the diameters of parcels to be injected

    do inoz = 1, nnozs

      thisNoz => nozzles_srr(inoz)

      dt_lagr = dtp
      nparcs = thisNoz%numbInjParcs_

      ! mass of each parcel
      thisNoz%parcsInj_(:)%mass_ = thisNoz%massRate_ * dt_lagr / nparcs

      ! the diameters
      tdmin  = thisNoz%miniDiamRR_
      tdmax  = thisNoz%maxiDiamRR_
      tdmean = thisNoz%meanDiamRR_
      tq     = thisNoz%expnRR_

      k_temp = ONE - exp( (tdmin/tdmean)**tq - (tdmax/tdmean)**tq )

      do iparc = 1, nparcs

        theParc => thisNoz%parcsInj_(iparc)

        rands => thisNoz%rands_(1)
        rands%irand_ = rands%irand_ + 1

        iptr = rands%irand_
        rand = rands%rands_(iptr)

        theParc%diam_ &
        = tdmean * ((tdmin / tdmean)**tq - log(ONE - k_temp*rand))**(ONE/tq)

        diam = theParc%diam_; dens = thisNoz%dens_; pi6 = PI/SIX

        theParc%numb_ = theParc%mass_ / (dens * pi6 * diam**THREE)
      end do
    end do

    call CheckIrands_()
    !---------------------------------------------------------------------------

    ! position and velocity of parcels to be injected

    do inoz = 1, nnozs

      thisNoz => nozzles_srr(inoz)

      ! to define t1_vec
      ! temperarily store the nozzle face vector: \vec n_ref, see fig 3.6
      v1 = thisNoz%dire_(1)
      v2 = thisNoz%dire_(2)
      v3 = thisNoz%dire_(3)

      ! comput \vec t1: see eqn 3.2.2.1-1
      t1_vec(1) = ONE - v1 * v1
      t1_vec(2) =     - v1 * v2
      t1_vec(3) =     - v1 * v3

      ! normalization
      t1_vec(1:3) = t1_vec(1:3) / vec_magn_( t1_vec(1:3) )

      ! need being fixed if %dire_ = (1, 0, 0)
      if(     abs(v1 - ONE) < EPS_RR  &
        .and. abs(v2      ) < EPS_RR  &
        .and. abs(v3      ) < EPS_RR ) then

        t1_vec = [ZERO, ONE, ZERO]
      end if

      ! \vec t2 = \vec n_ref \times \vec t1, see eqn 3.2.2.1-4
      call vec_product_( thisNoz%dire_, t1_vec, t2_vec )

      ! 3.2.2.1-40 & 41
      r1                                                    &
      = thisNoz%dispRadi_                                   &
      * tan( thisNoz%hfCoAngl_ - HALF * thisNoz%dispAngl_ ) &
      / tan( thisNoz%hfCoAngl_ )

      r2                                                    &
      = thisNoz%dispRadi_                                   &
      * tan( thisNoz%hfCoAngl_ + HALF * thisNoz%dispAngl_ ) &
      / tan( thisNoz%hfCoAngl_ )

      nparcs = thisNoz%numbInjParcs_

      do iparc = 1, nparcs

        theParc => thisNoz%parcsInj_(iparc)
        !-------------------------------------------------------------------------

        ! position in polar coordinates, see eqn 3.2.2.1-38 & 39

        rands => thisNoz%rands_(2)

        rands%irand_ = rands%irand_ + 1

        iptr = rands%irand_
        rand = rands%rands_(iptr)

        alpha_rr = TWO * PI * rand

        rands => thisNoz%rands_(3)

        rands%irand_ = rands%irand_ + 1

        iptr = rands%irand_
        rand = rands%rands_(iptr)

        r_rr = r1 + rand*(r2 - r1)
        !-------------------------------------------------------------------------

        ! 3.2.2.1-5, or 3.2.2.1-43
        nr_vec(1:3) =  cos(alpha_rr)*t1_vec(1:3) + sin(alpha_rr)*t2_vec(1:3)

        ! position of injecting points, see eqn 3.2.2.1-42
        theParc%pos_(1:3) = thisNoz%posi_(1:3) + r_rr*nr_vec(1:3)
        !-------------------------------------------------------------------------

        ! velocity of parcel

        ! eqn 3.2.2.1-44
        theta_rr = atan( r_rr / thisNoz%dispRadi_ * tan(thisNoz%hfCoAngl_) )

        ! eqn 3.2.2.1-45
        n_cone_vec(1:3)                       &
        = cos(theta_rr) * thisNoz%dire_(1:3)  &
        + sin(theta_rr) * nr_vec(1:3)

        ! eqn 3.2.2.1 - 35
        frac = thisNoz%fracSwl_

        v_swirl                                          &
        = thisNoz%magnVel_                               &
        / sqrt( ONE + ((ONE - abs(frac))/abs(frac)**TWO) )

        ! eqn 3.2.2.1 - 36
        v_cone                                         &
        = thisNoz%magnVel_                             &
        / sqrt( ONE + (abs(frac)/(ONE-abs(frac)))**TWO )

        ! eqn 3.2.2.1-6
        nt_vec(1:3) = -sin(alpha_rr)*t1_vec(1:3) + cos(alpha_rr)*t2_vec(1:3)

        ! eqn 3.2.2.1 - 37
        theParc%vel_(1:3)                                        &
        = v_cone*n_cone_vec(1:3) + sign_(frac)*v_swirl*nt_vec(1:3)
      end do
    end do

    call CheckIrands_()

  end subroutine hfy_srr_inject

  !*****************************************************************************
  ! INTERNAL PRIVATE PROCEDURES
  !*****************************************************************************

  subroutine resize_()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! allocate or resize nozzles%parcelInj

    ! after the injected number of parcels,
    ! 1, if not allocated, allocate it
    ! 2, if allocated, check the size of the arrays,
    !   2.1, if parc number > size of array, resize it
    !---------------------------------------------------------------------------
    intrinsic :: size

    integer, parameter :: NOVERSIZE = 1000

    integer :: nsize, nparcs, nnozs, inoz
    logical :: yes
    type(nozzle_t), pointer :: thisNoz => NULL()
    !---------------------------------------------------------------------------

    nnozs  = numb_nozz_srr

    do inoz = 1, nnozs
      thisNoz => nozzles_srr(inoz)

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

  end subroutine resize_

  subroutine GenRands_()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! allocate and generate random number sequences
    !---------------------------------------------------------------------------

    intrinsic :: mod, associated

    integer, parameter :: NOVERSIZE = 10000

    integer :: i, iuse, inoz, nnozs, nparcs, ndigits, nuses
    integer, dimension(:), allocatable :: seed
    type(nozzle_t), pointer :: thisNoz => NULL()
    type(rands_t),  pointer :: rands => NULL()
    !---------------------------------------------------------------------------

    nnozs = numb_nozz_srr

    nuses = 3

    do inoz = 1, nnozs
      thisNoz => nozzles_srr(inoz)

      if( .not. associated(thisNoz%rands_) )  &
        allocate( thisNoz%rands_(1:nuses) )

      do iuse = 1, nuses
        rands => thisNoz%rands_(iuse)

        call random_seed( size=ndigits )

        allocate( seed(1:ndigits) )

        do i = 1, ndigits
          seed(i) = i + (iuse+(inoz-1)*nuses) * 10000
        end do

        call random_seed( put=seed )
        !call random_seed()

        nparcs = thisNoz%numbTotParcs_

        allocate( rands%rands_(1:nparcs+NOVERSIZE) )

        ! the length or size of %rands_(:), used to check if %irand_ beyond
        rands%nSize_ = nparcs + NOVERSIZE

        call random_number( rands%rands_(:) )

        ! initialize the index indicator
        rands%irand_ = 0

        deallocate( seed )
      end do
    end do

  end subroutine GenRands_

  subroutine checkNoz_()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! 1, If the vec's magnitude /= unit, normalize the nozzle face vector
    ! 2, If the vec's magnitude \approx 0, print warning message
    !---------------------------------------------------------------------------

    real(RKD), parameter :: EPS = 1.0e-2_RKD

    integer :: inoz, nnozs
    real(RKD) :: magn
    type(nozzle_t), pointer :: thisNoz => NULL()
    !---------------------------------------------------------------------------

    nnozs = numb_nozz_srr

    do inoz = 1, numb_nozz_srr
      thisNoz => nozzles_srr(inoz)

      magn = vec_magn_( thisNoz%dire_ )

      if( magn < EPS ) then
        print *, "WARNING: nozzle ", inoz, " n_ref 's Magnitude Is Too Small"
      end if

      thisNoz%dire_(:) = thisNoz%dire_(:) / magn
    end do

  end subroutine checkNoz_

  function vec_magn_( vec_1d_ )
  !=============================================================================
    !---------------------------------------------------------------------------
    ! compute the magnitude of a 3d vector

    ! var         intent  attr          rank
    ! vec1d       in      real(RKD)     (:)
    ! vec_magn_   out     real(RKD)     0
    !---------------------------------------------------------------------------
    intrinsic :: sqrt

    real(RKD), dimension(1:3), intent(in) :: vec_1d_
    real(RKD) :: vec_magn_
    integer :: i
    !---------------------------------------------------------------------------

    vec_magn_ = ZERO

    do i = 1, 3
      vec_magn_ = vec_magn_ + vec_1d_(i) * vec_1d_(i)
    end do

    vec_magn_ = sqrt(vec_magn_)

  end function vec_magn_

  function sign_( x_ )
  !=============================================================================
    !---------------------------------------------------------------------------
    ! internal function
    ! return the sign of a variable of double precision

    ! var     intent  attr
    ! sign_   out     real(RKD)
    ! x_      in      real(RKD)
    !---------------------------------------------------------------------------

    real(RKD), intent(in) :: x_
    real(RKD) :: sign_
    !---------------------------------------------------------------------------

    if( x_ > ZERO ) then
      sign_ = ONE
    else
      sign_ = -ONE
    end if

  end function sign_

  subroutine vec_product_( vec1, vec2, vec_product )
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

  end subroutine vec_product_

  subroutine CheckIrands_()
  !=============================================================================
    !---------------------------------------------------------------------------
    ! check if the parcel index indicator exeeds its bound
    !---------------------------------------------------------------------------
    intrinsic :: size

    integer :: nsize, inoz, nnoz, iptr, iuse
    type(nozzle_t), pointer :: thisNoz => NULL()
    type(rands_t),  pointer :: rands => NULL()
    !---------------------------------------------------------------------------

    nnoz = numb_nozz_srr

    do inoz = 1, nnoz
      thisNoz => nozzles_srr(inoz)

      do iuse = 1, 3
        rands => thisNoz%rands_(iuse)

        iptr  = rands%irand_
        nsize = rands%nSize_

        if( iptr > nsize ) then
          print *, "Parcel Index of rands_srr Beyond! Nozzle: ", inoz
        end if
      end do
    end do

  end subroutine CheckIrands_

  ! ...
end module hfy_srr_mod
