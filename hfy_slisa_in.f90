
subroutine hfy_slisa_in()
!===============================================================================
  !-----------------------------------------------------------------------------
  ! author    date        aff
  ! snx       21.06.03    simpop.cn

  ! This file is used by the user to input the setting parameters
  !-----------------------------------------------------------------------------

  use hfy_slisa_mod
  !-----------------------------------------------------------------------------
  implicit none

  integer, parameter :: RKD = selected_real_kind( p=10 )

  integer :: inoz, nnozs
  type(NozzleSLS_t), pointer :: thisNoz => NULL()
  !-----------------------------------------------------------------------------

  nnozs = numbNozz_SLS

  do inoz = 1, nnozs

    thisNoz => nozzles_SLS(inoz)

    ! working time
    thisNoz%staTime_ = 0.1_RKD
    thisNoz%endTime_ = 0.2_RKD

    ! density, surface tension and visocity of fuel oil
    thisNoz%dens_ = 790.0_RKD
    thisNoz%tens_ = 25.0e-3_RKD
    thisNoz%visc_ = 1.8e-5_RKD

    ! mass flow rate
    thisNoz%massRate_ = 0.1_RKD

    ! pressure difference
    thisNoz%presDiff_ = 0.1e-5_RKD
    !---------------------------------------------------------------------------

    ! diameter of this nozzle, serving as out radius of spray dispersion
    thisNoz%diam_ = 0.0006_RKD

    ! positions of this nozzle
    thisNoz%posi_ = [0.0_RKD, 0.03_RKD, 0.0_RKD]

    ! vector denoting direction of nozzle face
    thisNoz%dire_ = [0.0_RKD, -1.0_RKD, 0.0_RKD]

    ! half-cone angle in radian
    thisNoz%hfCoAngl_ = 0.2_RKD

    ! dispersion angle in radian
    thisNoz%dispAngl_ = 0.003_RKD
    !---------------------------------------------------------------------------

    ! number of parcels per inje point over the whole injecting period
    thisNoz%numbTotParcs_ = 10000
  end do

end subroutine hfy_slisa_in
