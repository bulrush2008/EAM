
subroutine hfy_srr_in()
!===============================================================================
  !-----------------------------------------------------------------------------
  ! author    date        aff
  ! snx       21.05.24    simpop.cn

  ! this subroutine is exposed to user
  !-----------------------------------------------------------------------------
  use hfy_srr_mod
  !-----------------------------------------------------------------------------
  implicit none

  integer, parameter :: RKD = selected_real_kind( p=10 )

  integer :: inoz, nnozs
  type(nozzle_t), pointer :: thisNoz => NULL()
  !-----------------------------------------------------------------------------

  nnozs = numb_nozz_srr

  do inoz = 1, nnozs

    thisNoz => nozzles_srr(inoz)

    thisNoz%dens_ = 790.0_RKD

    ! mass flow rate for each nozzle
    thisNoz%massRate_ = 0.1_RKD

    ! ratio of swirling flow
    thisNoz%fracSwl_ = 0.1_RKD

    ! magnitude of vel at injecting points
    thisNoz%magnVel_ = 15.1_RKD
    !---------------------------------------------------------------------------

    ! positions of nozzles
    thisNoz%posi_ = [0.0_RKD, 0.03_RKD, 0.0_RKD]

    ! normal vectors of faces, should be unitary
    thisNoz%dire_ = [0.0_RKD, -1.0_RKD, 0.0_RKD]

    ! dispersion radius for each nozzle
    thisNoz%dispRadi_ = 0.0003_RKD

    ! half-cone angle in radian
    thisNoz%hfCoAngl_ = 0.2_RKD

    ! dispersion angle
    thisNoz%dispAngl_ = 0.003_RKD
    !---------------------------------------------------------------------------

    ! in RR model
    thisNoz%expnRR_ = 3

    ! representative diameter
    thisNoz%meanDiamRR_ = 0.04_RKD

    ! lower and upper limit of diameter of parcels in RR
    thisNoz%miniDiamRR_ = 1.e-4_RKD
    thisNoz%maxiDiamRR_ = 0.1_RKD
    !---------------------------------------------------------------------------

    ! number of parcels per inje point over the whole injecting period
    thisNoz%numbTotParcs_ = 10000
  end do

end subroutine hfy_srr_in
