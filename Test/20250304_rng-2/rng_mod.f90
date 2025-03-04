
! File: independent_rng.f90
module rng_mod
  implicit none
  private

  type, public :: RNG
    private
    integer(8)        :: state
    real              :: spare = 0.0  ! 高斯分布需要，均匀分布不需要
    logical           :: has_spare = .false.
  contains
    procedure :: init => seed_rng
    procedure :: normal => rng_normal
    procedure :: uniform => rng_uniform
  end type RNG

  integer(8), parameter :: LCG_MULT = 6364136223846793005_8
  integer(8), parameter :: LCG_INCR = 1442695040888963407_8
  integer(8), parameter :: MODULUS  = 9223372036854775807_8  ! 2^63-1

contains

  subroutine seed_rng(this, seed_val)
    class(RNG), intent(inout) :: this
    integer(8), intent(in)    :: seed_val
    this%state = mod(abs(seed_val), MODULUS)
    this%has_spare = .false.
  end subroutine seed_rng

  function rng_uniform(this) result(u)
    class(RNG), intent(inout) :: this
    real :: u
    integer(8) :: next_state
    
    ! 更新状态并生成均匀分布随机数
    ! 注意，mod(x,y) 中 x 可能由于过大，被机器自动识别为负数，因此需要 abs()
    next_state = mod(LCG_MULT * this%state + LCG_INCR, MODULUS)
    u = real(next_state) / real(MODULUS)  ! [0,1) 范围
    u = abs(u)
    this%state = next_state
  end function rng_uniform

  function rng_normal(this) result(n)
    class(RNG), intent(inout) :: this
    real :: n
    real :: u, v, s

    if (this%has_spare) then
      this%has_spare = .false.
      n = this%spare
    else
      do  ! Marsaglia polar method
        ! Generate in [-1, 1) range
        u = 2.0 * real(mod(LCG_MULT * this%state + LCG_INCR, MODULUS)) / real(MODULUS) - 1.0
        this%state = mod(LCG_MULT * this%state + LCG_INCR, MODULUS)
        v = 2.0 * real(mod(LCG_MULT * this%state + LCG_INCR, MODULUS)) / real(MODULUS) - 1.0
        this%state = mod(LCG_MULT * this%state + LCG_INCR, MODULUS)
        
        s = u*u + v*v
        if (s > 0.0 .and. s < 1.0) exit
      end do

      s = sqrt(-2.0 * log(s) / s)
      n = u * s     ! First normal value
      this%spare = v * s  ! Store second value
      this%has_spare = .true.
    end if
  end function rng_normal

end module rng_mod