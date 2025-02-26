
subroutine hfy_wave_in()
!===============================================================================
  !-----------------------------------------------------------------------------
  ! author    date        aff
  ! snx       21.06.16    simpop.cn

  ! This file is used by the user to input the setting parameters
  !-----------------------------------------------------------------------------

  use hfy_wave_mod
  !-----------------------------------------------------------------------------
  implicit none

  integer, parameter :: RKD = selected_real_kind( p=10 )

  integer :: nnozs, inoz

  type(nozzleWAV_t), pointer :: thisNoz => NULL()
  !-----------------------------------------------------------------------------

  nnozs = numbNozz_WAV

  do inoz = 1, nnozs
    thisNoz => nozzles_WAV(inoz)

    if( inoz == 1 ) then
      ! 1st nozzle

      ! working time
      thisNoz%staTime_ = 0.0_RKD
      thisNoz%endTime_ = 1.0_RKD

      ! position and orientation
      thisNoz%posi_ = [0.0_RKD, 0.1_RKD, 0.0_RKD]
      thisNoz%dire_ = [0.0_RKD, 1.0_RKD, 0.0_RKD]

      ! density, viscosity and surface tension
      thisNoz%dens_ = 790.0_RKD
      thisNoz%tens_ = 25.0e-3_RKD
      thisNoz%visc_ = 1.8e-5_RKD

      ! presure difference
      thisNoz%presDiff_ = 0.4e5_RKD

      ! massrate
      thisNoz%massRate_ = 6.19e-3_RKD

      ! half-cone angle
      thisNoz%hfCoAngl_ = 0.15_RKD

      ! noz diameter
      thisNoz%diam_ = 0.6e-3_RKD

      ! total number of blob parcels to be inject
      thisNoz%numbTotParcs_ = 10000
    else
      ! other nozzles
      ! ...
    end if
  end do

  nullify(thisNoz)

end subroutine hfy_wave_in